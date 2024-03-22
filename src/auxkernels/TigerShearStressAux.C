//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "TigerShearStressAux.h"

registerMooseObject("MooseApp", TigerShearStressAux);

InputParameters
TigerShearStressAux::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("compute stress excess from stress tensor and friction parameters");
  params.addRequiredParam<MaterialPropertyName>("total_stress",
                                                "The rank two material tensor name");

  return params;
}

TigerShearStressAux::TigerShearStressAux(const InputParameters & parameters)
  : AuxKernel(parameters), _TenMech_total_stress(getMaterialProperty<RankTwoTensor>("total_stress"))
{
}

RealVectorValue
TigerShearStressAux::compute_normal_vector()
{
  if (_current_elem->dim() != 2)
    mooseError("trying to normal vector, but dim !=2");
  RealVectorValue v1 = _current_elem->point(1) - _current_elem->point(0);
  RealVectorValue v2 = _current_elem->point(2) - _current_elem->point(1);
  return (v1.cross(v2)).unit();
}

Real
TigerShearStressAux::computeValue()
{
  RealVectorValue normal_vector = compute_normal_vector();
  RealVectorValue traction = _TenMech_total_stress[_qp] * normal_vector;
  // note that sig_n is negative (in compression)
  Real sig_n = traction * normal_vector;

   RealVectorValue shearTraction = traction - sig_n * normal_vector;
  Real sig_t = sqrt(shearTraction * shearTraction);
  return sig_t;
}
