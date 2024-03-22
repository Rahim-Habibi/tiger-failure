
#pragma once

#include "libmesh/libmesh_config.h"
#include "AuxKernel.h"
#include "libmesh/fparser.hh"
#include "RankTwoTensor.h"

/**
 * This AuxKernel requests termination of the current solve based on
 * the values of failure criterion defined as a ratio
 *
 *                     <((((((\\\
 *                     /      . }\
 *                     ;--..--._|}
 *  (\                 '--/\--'  )
 *   \\                | '-'  :'|
 *    \\               . -==- .-|
 *     \\               \.__.'   \--._
 *     [\\          __.--|       //  _/'--.
 *     \ \\       .'-._ ('-----'/ __/      \
 *      \ \\     /   __>|      | '--.       |
 *       \ \\   |   \   |     /    /       /
 *        \ '\ /     \  |     |  _/       /
 *         \  \       \ |     | /        /
 *          \  \      \        /
 */
class TigerTerminator : public AuxKernel
{
public:
  static InputParameters validParams();

  TigerTerminator(const InputParameters & parameters);

protected:
  /// handle output of the optional message
  void handleMessage();
  virtual Real computeValue() override;
  const enum class FailMode { HARD, SOFT } _fail_mode;
  const enum class ErrorLevel { INFO, WARNING, ERROR, NONE } _error_level;

  /// cohesion of the material in Pa
  Real _cohesion;
  /// friction coefficient of material in degree
  Real _phi;
  Real _mu_s;
  /// to call main object of the problem
  FEProblemBase & _problem;
  /// a criterion in form of time step
  Real & _dt;
  ///  min dt below which the code is assumed to have converged to the threshold
  const Real & _dt_min;
  /// normal stress
  const MaterialProperty<Real> & _sig_n;
  /// shear stress
  const MaterialProperty<Real> & _sig_t;
  /// stress tensor
  const MaterialProperty<RankTwoTensor> * _tensor;

private:
  // enum to select failure criterion
  const enum class CriterionType { Mohr_Coulomb, Slip_Tendency, Stress_Excess, Fracture_Potential } _criterion_type;
};
