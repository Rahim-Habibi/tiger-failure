#include "libmesh/libmesh_config.h"
#include "TigerTerminator.h"
#include "RankTwoScalarTools.h"
#include "MooseEnum.h"
#include "Executioner.h"
#include "FEProblem.h"
#include <cfloat>
#include "Function.h"
#define PI 3.141592653589793238462643383279502884

registerMooseObject("TigerApp", TigerTerminator);

InputParameters
TigerTerminator::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Requests termination of the current solve based on the evaluation of"
                             " a parsed logical expression of the Postprocessor value(s).");

  MooseEnum failModeOption("HARD SOFT", "HARD");
  params.addParam<MooseEnum>(
      "fail_mode",
      failModeOption,
      "Abort entire simulation (HARD) or just the current time step (SOFT).");
  params.addParam<std::string>("message",
                               "An optional message to be output instead of the default message "
                               "when the termination condition is triggered");

  MooseEnum errorLevel("INFO WARNING ERROR NONE", "INFO");
  params.addParam<MooseEnum>(
      "error_level",
      errorLevel,
      "The error level for the message. A level of ERROR will always lead to a hard "
      "termination of the entire simulation.");
  params.addParam<Real>("cohesion", 0.0, "bulck cohesion should be given by user");
  params.addParam<Real>("phi", 35.0, "internal friction coefficient in degree");
  params.addParam<Real>("mu_s", 0.6, "static friction coefficient");
  params.addParam<Real>(
      "dt_min", 43200, "min dt below which the code is assumed to have converged to the threshold");
  MooseEnum criterion_type_options("Mohr_Coulomb Slip_Tendency Stress_Excess", "Slip_Tendency", "Fracture_Potential");
  params.addParam<MooseEnum>(
      "criterion_type", criterion_type_options, "Criterion to use for the threshold");
  return params;
}

TigerTerminator::TigerTerminator(const InputParameters & parameters)
  : AuxKernel(parameters),
    _fail_mode(getParam<MooseEnum>("fail_mode").getEnum<FailMode>()),
    _error_level(getParam<MooseEnum>("error_level").getEnum<ErrorLevel>()),
    _cohesion(getParam<Real>("cohesion")),
    _phi(getParam<Real>("phi")),
    _mu_s(getParam<Real>("mu_s")),
    _dt_min(getParam<Real>("dt_min")),
    _problem(_c_fe_problem),
    _dt(_problem.dt()),
    _criterion_type(getParam<MooseEnum>("criterion_type").getEnum<CriterionType>()),
    _sig_n(getMaterialProperty<Real>("NORMAL_STRESS")),
    _sig_t(getMaterialProperty<Real>("SHEAR_STRESS")),
    _tensor(nullptr)
{
  if (hasMaterialProperty<RankTwoTensor>("rank_two_tensor"))
    _tensor = &getMaterialProperty<RankTwoTensor>("rank_two_tensor");
  else
    mooseError("Error in RankTwoBasedFailureCriteriaNOSPD! Required rank two tensor is not "
               "available for current peridynamics model!");

  if (isNodal())
    paramError("variable", "This AuxKernel only supports Elemental fields");
}

void
TigerTerminator::handleMessage()
{
  std::string message;
  if (!isParamValid("message"))
  {
    message = "Terminator '" + name() + "' is causing ";
    if (_fail_mode == FailMode::HARD)
      message += "the execution to terminate.\n";
    else
      message += "a time step cutback by marking the current step as failed.\n";
  }
  else
    message = getParam<std::string>("message");

  switch (_error_level)
  {
    case ErrorLevel::INFO:
      mooseInfo(message);
      break;

    case ErrorLevel::WARNING:
      mooseWarning(message);
      break;

    case ErrorLevel::ERROR:
      mooseError(message);
      break;

    default:
      break;
  }
}

Real
TigerTerminator::computeValue()
{
  Real ratio = 1.0;
  switch (_criterion_type)
  {
    case CriterionType::Mohr_Coulomb:
    {
      Real _tau_strength = _cohesion - (std::max(0.0, -_sig_n[_qp]) * tan(_phi * PI / 180));
      ratio = _sig_t[_qp] / _tau_strength;
      if (abs(ratio) >= 0.48)
      {
        handleMessage();
        if (_fail_mode == FailMode::HARD)
          _c_fe_problem.terminateSolve();
        else
        {
          if (_dt >= _dt_min)
            getMooseApp().getExecutioner()->fixedPointSolve().failStep();
          // Real old_ratio = ratio;
          else
            _c_fe_problem.terminateSolve();
        } // end of else
      }
      // return tau_over_sigma;
      break;
    }
    case CriterionType::Slip_Tendency:
    {
      ratio = _sig_t[_qp] / (std::max(0.0, -_sig_n[_qp]) * _mu_s);
      if (abs(ratio) >= 0.99)
      {
        handleMessage();
        if (_fail_mode == FailMode::HARD)
          _c_fe_problem.terminateSolve();
        else
        {
          if (_dt >= _dt_min)
            getMooseApp().getExecutioner()->fixedPointSolve().failStep();
          else
            _c_fe_problem.terminateSolve();
        }
      }
      // return tau_over_normal;
      break;
    }
    case CriterionType::Stress_Excess:
    {
      ratio = _sig_t[_qp] - _mu_s * std::max(0.0, -_sig_n[_qp]);

      if (ratio >= 0.02e6)
      {
        handleMessage();
        if (_fail_mode == FailMode::HARD)
          _c_fe_problem.terminateSolve();
        else
        {
          if (_dt >= _dt_min)
            getMooseApp().getExecutioner()->fixedPointSolve().failStep();
          else
            _c_fe_problem.terminateSolve();
        }
      }
      // return tau_over_normal;
      break;
    }
    case CriterionType::Fracture_Potential:
    {
      RankTwoTensor avg_tensor = 0.5 * ((*_tensor)[0] + (*_tensor)[1]);
      Real _max_shear = 2.0 * RankTwoScalarTools::maxShear(avg_tensor);         // calculates 'sigma 1 - sigma 3' in tems of mohr circle
      Real _distance = _sig_n[_qp] * (tan(_phi * PI / 180) * cos(_phi * PI /180));
      ratio = _max_shear / (_max_shear + 2 * _distance);
        if (abs(ratio) >= 0.99)
      {
        handleMessage();
        if (_fail_mode == FailMode::HARD)
          _c_fe_problem.terminateSolve();
        else
        {
          if (_dt >= _dt_min)
            getMooseApp().getExecutioner()->fixedPointSolve().failStep();
          else
            _c_fe_problem.terminateSolve();
        }
      }
      break;
    }
  }
  return ratio;
}
