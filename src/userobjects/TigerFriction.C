/**************************************************************************/
/*  TIGER - THMC sImulator for GEoscience Research                        */
/*                                                                        */
/*  Copyright (C) 2017 by Maziar Gholami Korzani                          */
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

#include "TigerFrictionUO.h"

InputParameters
TigerFrictionUO::validParams()
{
  InputParameters params = GeneralUserObject::validParams();

  return params;
}

TigerFrictionUO::TigerFrictionUO(const InputParameters & parameters)
  : GeneralUserObject(parameters)
{
}

void
TigerFrictionUO::execute(){}

void
TigerFrictionUO::initialize(){}

void
TigerFrictionUO::finalize(){}

RankTwoTensor
TigerFrictionUO::FrictionTensorCalculator(const int & dim, const std::vector<Real> & mu, const MooseEnum & _friction_type) const
{
  RealVectorValue mu_x;
  RealVectorValue mu_y;
  RealVectorValue mu_z;

  if (dim == 1)
  {
    switch (_friction_type)
    {
      case FT::isotropic:
        if (mu.size() != 1)
          mooseError(name(),": One input value is needed for isotropic distribution of friction! You provided ", mu.size(), " values.\n");
        mu_x = RealVectorValue(mu[0], 0.0, 0.0);
        mu_y = RealVectorValue(0.0  , 0.0, 0.0);
        mu_z = RealVectorValue(0.0  , 0.0, 0.0);
        break;
      case FT::orthotropic:
      case FT::anisotropic:
        mooseError(name(),": One dimensional elements cannot have non-isotropic friction values.\n");
        break;
    }
  }
  else if (dim == 2)
  {
    switch (_friction_type)
    {
      case FT::isotropic:
        if (mu.size() != 1)
          mooseError(name(),": One input value is needed for isotropic distribution of friction! You provided ", mu.size(), " values.\n");
        mu_x = RealVectorValue(mu[0], 0.0  , 0.0);
        mu_y = RealVectorValue(0.0  , mu[0], 0.0);
        mu_z = RealVectorValue(0.0  , 0.0  , 0.0);
        break;

      case FT::orthotropic:
        if (mu.size() != 2)
          mooseError(name(),": Two input values are needed for orthotropic distribution of friction in two dimensional elements! You provided ", mu.size(), " values.\n");
        mu_x = RealVectorValue(mu[0], 0.0  , 0.0);
        mu_y = RealVectorValue(0.0  , mu[1], 0.0);
        mu_z = RealVectorValue(0.0  , 0.0  , 0.0);
        break;

      case FT::anisotropic:
        if (mu.size() != 4)
          mooseError(name(),": Four input values are needed for anisotropic distribution of friction in two dimensional elements! You provided ", mu.size(), " values.\n");
        mu_x = RealVectorValue(mu[0], mu[1], 0.0);
        mu_y = RealVectorValue(mu[2], mu[3], 0.0);
        mu_z = RealVectorValue(0.0  , 0.0  , 0.0);
        break;
    }
  }
  else if (dim == 3)
  {
    switch (_friction_type)
    {
      case FT::isotropic:
        if (mu.size() != 1)
          mooseError(name(),": One input value is needed for isotropic distribution of friction! You provided ", mu.size(), " values.\n");
        mu_x = RealVectorValue(mu[0], 0.0, 0.0);
        mu_y = RealVectorValue(0.0, mu[0], 0.0);
        mu_z = RealVectorValue(0.0, 0.0, mu[0]);
        break;

      case FT::orthotropic:
        if (mu.size() != 3)
          mooseError(name(),": Three input values are needed for orthotropic distribution of friction! You provided ", mu.size(), " values.\n");
        mu_x = RealVectorValue(mu[0], 0.0, 0.0);
        mu_y = RealVectorValue(0.0, mu[1], 0.0);
        mu_z = RealVectorValue(0.0, 0.0, mu[2]);
        break;

      case FT::anisotropic:
        if (mu.size() != 9)
          mooseError(name(),": Nine input values are needed for anisotropic distribution of friction! You provided ", mu.size(), " values.\n");
        mu_x = RealVectorValue(mu[0], mu[1], mu[2]);
        mu_x = RealVectorValue(mu[3], mu[4], mu[5]);
        mu_x = RealVectorValue(mu[6], mu[7], mu[8]);
        break;
    }
  }
  return RankTwoTensor::initializeFromRows(mu_x, mu_y, mu_z);

}