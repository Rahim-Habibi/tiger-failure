/**************************************************************************/
/*  TIGER - THMC sImulator for GEoscience Research                        */
/*                                                                        */
/*  Copyright (C) 2017 by Maziar Gholami Korzani, Robert Egert            */
/*  Karlsruhe Institute of Technology, Institute of Applied Geosciences   */
/*  Division of Geothermal Research                                       */
/*                                                                        */
/*  This file is part of TIGER App                                        */
/*                                                                        */
/*  This program is free software: you can redistribute it and/or modify  */
/*  it under the terms of the GNU General Public License as published by  */
/*  the Free Software Foundation, either version 3 of the License, or     */
/*  (at your option) any later version.                                   */
/*                                                                        */
/*  This program is distributed in the hope that it will be useful,       */
/*  but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          */
/*  GNU General Public License for more details.                          */
/*                                                                        */
/*  You should have received a copy of the GNU General Public License     */
/*  along with this program.  If not, see <http://www.gnu.org/licenses/>  */
/**************************************************************************/

#pragma once

#include "GeneralUserObject.h"
#include "RankTwoTensor.h"

class TigerFrictionUO : public GeneralUserObject
{
public:
  static InputParameters validParams();

  TigerFrictionUO(const InputParameters & parameters);

  virtual void execute();
  virtual void initialize();
  virtual void finalize();

  /// friction matrix; to be called in Material
  virtual RankTwoTensor Friction(const int & dim, const Real & porosity, const Real & scale_factor, const std::vector<Real> kmat) const = 0;

protected:
  // Creates the friction tensor as function of input and dimension
  virtual RankTwoTensor FrictionTensorCalculator(const int & dim, const std::vector<Real> & mu, const MooseEnum & _friction_type) const;

  enum FT {isotropic, orthotropic, anisotropic};
};