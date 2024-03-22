//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"

/**
 * This auxiliary kernel computes its value by dividing "numerator" by
 * "denominator.  For efficiency, it doesn't check the denominator for
 * zero before dividing.  Perhaps a derived class, CheckedTigerNormalStressAux,
 * could be added for people who want this feature.
 */
class TigerNormalStressAux : public AuxKernel
{
public:
  static InputParameters validParams();

  TigerNormalStressAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  const MaterialProperty<RankTwoTensor> & _TenMech_total_stress;
  /// a function to calculate the components of the normal vector
  RealVectorValue compute_normal_vector();
};
