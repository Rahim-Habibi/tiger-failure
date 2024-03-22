#pragma once
#include "Material.h"
#include "libmesh/quadrature_gauss.h"
#include "RankTwoTensor.h"

class TigerFailure : public Material
{
public:
  static InputParameters validParams();

TigerFailure (const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

const bool _use_displaced_mesh;
 std::string _base_name;
/// normal stress perpendicular on the element
const MaterialProperty <Real> & _sig_n;
/// shear stress parallel to the element
const MaterialProperty<Real>  &  _sig_t;
};
