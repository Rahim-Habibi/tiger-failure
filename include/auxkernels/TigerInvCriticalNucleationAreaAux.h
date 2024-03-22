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
 * zero before dividing.  Perhaps a derived class, CheckedTigerInvCriticalNucleationAreaAux,
 * could be added for people who want this feature.
 */
class TigerInvCriticalNucleationAreaAux : public AuxKernel
{
public:
  static InputParameters validParams();

  TigerInvCriticalNucleationAreaAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;
  const VariableValue & _normal_stress;
  /// static friction
  Real _mu_s;
  // dynamic friction
  Real _mu_d;
  // C(nu)*G in Uenishi 2009, see Galis 2015, eq 19, with G rigidity (in Pa)
  Real _CG;
  // linear slip weakening distance (in m)
  Real _d_c;
};
