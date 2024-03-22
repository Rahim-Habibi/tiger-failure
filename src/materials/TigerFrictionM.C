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

#include "TigerFrictionM.h"
#include "MooseMesh.h"
#include "Function.h"

registerMooseObject("TigerApp", TigerFrictionM);

InputParameters
TigerFrictionM::validParams()
{
  InputParameters params = Material::validParams();

  params.addParam<std::string>("base_name", "the identical base name provided "
        "in TensorMechanics Action");
  params.addParam<std::vector<FunctionName>>("extra_friction_vector",
        "Vector of values defining the extra stress "
        "to add, in order 11, 22, 33. Functions can be provided as well.");
  params.addClassDescription("friction coefficient for rupture");

  return params;
}

TigerFrictionM::TigerFrictionM(const InputParameters & parameters)
  : Material(parameters),
    _base_name(isParamValid("base_name") ? getParam<std::string>("base_name") + "_" : ""),
    _TenMech_extra_friction(declareProperty<RankTwoTensor>(_base_name + "extra_friction"))
{
}

void
TigerFrictionM::computeQpProperties()
{

// Extra friction can be added and included in TensorMechanics Action
  const std::vector<FunctionName> & friction_fct(
      getParam<std::vector<FunctionName>>("extra_friction_vector"));
  const unsigned num = friction_fct.size();
  if (!(num == 0 || num == 3 || num == 9))
    mooseError("Please supply either zero or 3 or 9 extra stresses. "
               "You supplied ", num,".\n");

  _extra_friction.resize(num);
  for (unsigned i = 0; i < num; ++i)
    _extra_friction[i] = &getFunctionByName(friction_fct[i]);

  if (_extra_friction.size() == 3)
    for (unsigned i = 0; i < num; ++i)
      _TenMech_extra_friction[_qp](i, i) = _extra_friction[i]->value(_t, _q_point[_qp]);

if (_extra_friction.size() == 9) {
    for (unsigned i = 0; i < 3; ++i) {
        for (unsigned j = 0; j < 3; ++j) {
            _TenMech_extra_friction[_qp](i, j) = _extra_friction[i * 3 + j]->value(_t, _q_point[_qp]);
        }
    }
}
}
