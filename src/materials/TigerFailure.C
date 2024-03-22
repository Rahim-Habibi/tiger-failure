

#include "TigerFailure.h"

registerMooseObject("TigerApp", TigerFailure);

InputParameters
TigerFailure::validParams()
{
  InputParameters params = Material::validParams();

      params.addParam<std::string>("base_name", "the identical base name provided "
            "in TensorMechanics Action");
      params.set<bool>("use_displaced_mesh") = true;
  params.addClassDescription(
      "Material to compute Normal and shear stress on the 2D element.  This is mostly designed for 2D "
      "elements living in 3D space, however, the 1D-element and 3D-element cases are set  "
      "zero intentionally.  The Variable for this Material must be an elemental Variable");

  return params;
}

TigerFailure::TigerFailure(const InputParameters & parameters)
  : Material(parameters),
    _use_displaced_mesh(getParam<bool>("use_displaced_mesh")),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _sig_n(declareProperty<Real>("NORMAL_STRESS")),
    _sig_t(declareProperty<Real>("SHEAR_STRESS"))
{

}

void
TigerFailure::computeQpProperties()
{
	_ratio = _sig_t [_qp] - 0.6 * (- _sig_n [_qp]);
}
