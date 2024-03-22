//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TigerInvCriticalNucleationAreaAux.h"

registerMooseObject("MooseApp", TigerInvCriticalNucleationAreaAux);

InputParameters
TigerInvCriticalNucleationAreaAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("compute critical nucleation area distribution");
  params.addRequiredCoupledVar("normal_stress", "the on-fault normal stress");
  params.addParam<Real>("mu_s", 0.6, "static friction coefficient");
  params.addParam<Real>("mu_d", 0.2, "dynamic friction coefficient");
  params.addParam<Real>(
      "CG", 9.2e+09, "C(nu)*G in Uenishi 2009, see Galis 2015, eq 19, with G rigidity");
  params.addParam<Real>("d_c", 0.1, "linear slip weakening distance");

  return params;
}

TigerInvCriticalNucleationAreaAux::TigerInvCriticalNucleationAreaAux(
    const InputParameters & parameters)
  : AuxKernel(parameters),
    _normal_stress(coupledValue("normal_stress")),
    _mu_s(getParam<Real>("mu_s")),
    _mu_d(getParam<Real>("mu_d")),
    _CG(getParam<Real>("CG")),
    _d_c(getParam<Real>("d_c"))
{
}

Real
TigerInvCriticalNucleationAreaAux::computeValue()
{
  Real PI = std::atan(1.0) * 4;
  Real W = (_mu_s - _mu_d) * std::max(0.0, -_normal_stress[_qp]) / _d_c;
  Real L = 0.624 * _CG / W;
  return 1. / (PI * L * L);
}
