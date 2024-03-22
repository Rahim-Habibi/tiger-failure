[Mesh]
  [./file]
    type = FileMeshGenerator
    file = box.msh
  [../]
[]
[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  pressure = pressure
  #has_supg = true
  temperature = temperature
  stress_free_temperature = 'T0'
  thermal_expansion_coeff = 1e-5 # 8e-6  engineering toolbox # linear
[]
# properties of fluid is defined here, and later will be call by material block
[Modules]
  [./TensorMechanics]
    [./Master]
      [./all]
      add_variables = true
      strain = SMALL
      incremental = true
      #temperature = temperature
      #use_displaced_mesh = true
      eigenstrain_names = 'reduced_eigenstrain'
      volumetric_locking_correction = true
      #additional_generate_output = 'stress_xx stress_xy stress_xz stress_yy stress_yz stress_zz'
      additional_generate_output = 'stress_xx stress_xy stress_xz stress_yx stress_yy stress_yz stress_zx stress_zy stress_zz'
      [../]
    [../]
  [../]
 []
  [FluidProperties]
   [./water_uo]
  type = TigerIdealWater
  cp = 3800
  thermal_conductivity = 0.6
  bulk_modulus = 2.16356e+09
  [../]
[]


# Permeability Assigning
[UserObjects]
  [./rock_uo0]
    type = TigerPermeabilityConst
    permeability_type = isotropic
     k0 = '3.0e-17'
  [../]
  [./rock_uo1]
    type = TigerPermeabilityConst
    permeability_type = isotropic
    k0 = '1.5e-14'
  [../]
  [./rock_uo2]
    type = TigerPermeabilityConst
    permeability_type = isotropic
    k0 = '5.0e-5'
  [../]
  [./rock_uo3]
    type = TigerPermeabilityConst
    permeability_type = isotropic
    k0 = '4.0e-10'
  [../]
  [./supg_w]
    type = TigerSUPG
    effective_length = min
    supg_coeficient = transient_tezduyar
  [../]
  [./supg_f]
    type = TigerSUPG
    effective_length = average
    supg_coeficient = transient_tezduyar
  [../]
  [./supg_m]
    type = TigerSUPG
    effective_length = directional_average
    supg_coeficient = transient_tezduyar
  [../]
[]
[AuxVariables]
  [./p0]
    family = LAGRANGE
    order = FIRST
  [../]
  [./dP]
    family = LAGRANGE
    order = FIRST
  [../]
  [./T0]
    family = LAGRANGE
    order = FIRST
  [../]
  [./dT]
    family = LAGRANGE
    order = FIRST
  [../]
  [./stress_excess]
    family = MONOMIAL
    order = CONSTANT
    #initial_condition = '-1e10'
    block = 'middle right1 right2 left1 left2'
  [../]
  [./invCritNucArea]
    family = MONOMIAL
    order = CONSTANT
    #initial_condition = '-1e10'
    block = 'middle right1 right2 left1 left2'
  [../]
  [./normal_stress]
    family = MONOMIAL
    order = CONSTANT
    #initial_condition = '-1e10'
    block = 'middle right1 right2 left1 left2'
  [../]
  [./shear_stress]
    family = MONOMIAL
    order = CONSTANT
    #initial_condition = '-1e10'
    block = 'middle right1 right2 left1 left2'
  [../]
  [./aux_equals_1]
 initial_condition = 1e7
[../]
[]
[AuxKernels]
  [./p0_ker]
    type = FunctionAux
    variable = 'p0'
    function = '(1000*9.81*(500-(z+200)))'
    execute_on = 'initial'
  [../]
  [./stress_excess_ker]
    type = ParsedAux
    variable = 'stress_excess'
    coupled_variables = 'shear_stress normal_stress pressure'
    expression = 'shear_stress - 0.6 * max(0.0, -normal_stress - pressure)'
    execute_on = 'TIMESTEP_END'
    block = 'middle right1 right2 left1 left2'
  [../]

  [./dP_ker]
    type = ParsedAux
    variable = 'dP'
    coupled_variables = 'pressure p0'
    expression = 'pressure-p0'
    execute_on = 'TIMESTEP_END'
  [../]
  [./T0_ker]
    type = FunctionAux
    variable = 'T0'
    function = '((0.03*(500-z))+283.15)'
    execute_on = 'initial'
  [../]

  [./dT_ker]
    type = ParsedAux
    variable = 'dT'
    coupled_variables = 'temperature T0'
    expression = 'temperature-T0'
    execute_on = 'TIMESTEP_END'
  [../]


  [./normalstress_ker]
    type = TigerNormalStressAux
    variable = 'normal_stress'
    total_stress = stress
    execute_on = 'TIMESTEP_END'
    block = 'middle right1 right2 left1 left2'
  [../]
  [./shearstress_ker]
    type = TigerShearStressAux
    variable = 'shear_stress'
    total_stress = stress
    execute_on = 'TIMESTEP_END'
    block = 'middle right1 right2 left1 left2'
  [../]
[]
# Materials Assigning
[Materials]
  [./Elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    youngs_modulus = 15e9
    poissons_ratio = 0.2
    block = ' box middle right1 right2 left1 left2'
  [../]
  [./stress]
    type = ComputeFiniteStrainElasticStress
    block = ' box middle right1 right2 left1 left2'
  [../]
  [./thermal_expansion]
    type = ComputeThermalExpansionEigenstrain
    eigenstrain_name = 'thermal_eigenstrain'
    block = ' box middle right1 right2 left1 left2'
  [../]
  [./reduced_order_eigenstrain]
    type = ComputeReducedOrderEigenstrain
    input_eigenstrain_names = 'thermal_eigenstrain'
    eigenstrain_name = 'reduced_eigenstrain'
    block = ' box middle right1 right2 left1 left2'
  [../]
  [./rock_m]
    type = TigerMechanicsMaterialM
    disps = 'disp_x disp_y disp_z'
    incremental = true
    biot_coefficient = 1
    solid_bulk_modulus = 8e9
    extra_stress_vector = '0 0 0'
    # output_properties = 'extra_stress'
    # outputs = 'exodus'
    block = ' box middle right1 right2 left1 left2'
  [../]
#   [./strain]
#   type = ComputeSmallStrain
#   eigenstrain_names = ini_stress
# [../]
# [./ini_stress]
#   type = ComputeEigenstrainFromInitialStress
#   initial_stress = 'ini_xx 0 0  0 0 0  0 0 0'
#   eigenstrain_name = ini_stress
# [../]
  [./rock_g0]
    type = TigerGeometryMaterial
     gravity = '0 0 -9.8'
     scale_factor = 1 # 100
     block = ' box'
  [../]
  [./rock_fault]
    type = TigerGeometryMaterial
     gravity = '0 0 -9.8'
     scale_factor = 0.01 # Fault Apperture
     block = 'middle right1 right2 left1 left2'
  [../]
  [./rock_p0]
    type = TigerPorosityMaterial
    porosity = 0.15
    specific_density = 2500
    block = ' box'
  [../]
  [./rock_p1]
    type = TigerPorosityMaterial
    porosity = 0.15
    specific_density = 2500
    block = ' middle right1 right2 left1 left2'
  [../]
  [./rock_f]
    type= TigerFluidMaterial
    fp_uo = water_uo
    block = '  box middle right1 right2 left1 left2'
  [../]
  [./rock_h0]
    type = TigerHydraulicMaterialH
    pressure = pressure
    compressibility = 5.0e-10
    kf_uo = rock_uo0
    block = ' box'
  [../]
  [./FaultActivation]
    type = TigerHydraulicMaterialH
    pressure = pressure
    compressibility = 5.0e-10
    kf_uo = rock_uo3
    block = 'middle right1 right2 left1 left2'
  [../]
  [./Thermalsediment]
    type = TigerThermalMaterialT
    conductivity_type = isotropic
    mean_calculation_type = geometric
    lambda = 3
    specific_heat = 800
    has_supg = true
    supg_uo = supg_f
    advection_type = darcy_velocity
    block = ' box'
  [../]
  [./ThermalFaults]
    type = TigerThermalMaterialT
    conductivity_type = isotropic
    mean_calculation_type = geometric
    lambda = 3
    specific_heat = 800
    has_supg = true
    advection_type = darcy_velocity
    supg_uo = supg_f
    block = 'middle right1 right2 left1 left2'
  [../]
[]
# Boundary Conditions
[BCs]
  [./no_x]
    type = DirichletBC
    variable = disp_x
    boundary = 'right left'
    value = 0.0
  [../]
  [./no_y]
    type = DirichletBC
    variable = disp_y
    boundary = 'back front'
    value = 0.0
  [../]
  [./no_z]
    type = DirichletBC
    variable = disp_z
    boundary = 'bottom'
    value = 0.0
  [../]
  # [./overburdenz]
  #   type = Pressure
  #   boundary = 'top'
  #   variable = disp_z
  #   function = ini_zz_force
  # [../]
  # [./overburdenx]
  #   type = Pressure
  #   boundary = 'right'
  #   variable = disp_x
  #   function = ini_xx_force
  # [../]
  # [./overburdeny]
  #   type = Pressure
  #   boundary = 'back'
  #   variable = disp_y
  #   function = ini_xx_force
  # [../]
  [./Top]
    type = FunctionDirichletBC
    variable = pressure
    boundary = 'top bottom'
    function = hydrostatic
  [../]
  [./TempBottom]
    type = FunctionDirichletBC
    variable = temperature
    boundary = 'top bottom'
    function = Thermal
  [../]
  [./t_inject]
    type = DirichletBC
    variable = temperature
    boundary = 'inj'
    value = 313.15 # 40 C
  [../]
[]
# [Controls]
#   [./flowrate_well_injection_year1]
#     type = TimePeriod
#       enable_objects = 'BCs::t_inject'
#       start_time = '35033400'
#       end_time = 77583600 # 492 days
#   [../]
#   []
[Functions]
  [./ini_xx]
    type = ParsedFunction
    expression = '-((-z) * 9.81 * 2500) '
  [../]
  [./ini_yy]
    type = ParsedFunction
    expression = '-((-z) * 9.81 * 2500) '
  [../]
  [./ini_zz]
    type = ParsedFunction
    expression = '10e6' # total stress is negative, pore pressure is positive. so the effective can be read in this way
    # expression = 'if(z >= 0, 0,
    #           if(z < 0, (-((-z) * 9.81 * 2500) + (1000*9.81*(-z))),0))'
  [../]
  [./ini_xx_force]
    type = ParsedFunction
    value = '((-z) * 9.81 * 2500)'
  [../]
  [./ini_yy_force]
    type = ParsedFunction
    value = '((-z) * 9.81 * 2500)'
  [../]
  [./ini_zz_force]
    type = ParsedFunction
    value = '((-z) * 9.81 * 2500)-9000000'
  [../]
  [./hydrostatic]
    type = ParsedFunction
    expression = '(1000*9.81*(-z))'
  [../]
  [./Thermal]
    type = ParsedFunction
    expression = '((0.03*(-z))+283.15)'
  [../]
[]
[ICs]
  [./hydrostatic_ic]
    type = FunctionIC
    variable = pressure
    function = hydrostatic
  [../]
  [./temperature_ic]
    type = FunctionIC
    variable = temperature
    function = Thermal
  [../]
[]
# Variable Definitions
[Variables]
  [./pressure]
  [../]
  [./temperature]
  [../]
  [./disp_x]
    initial_condition = 0
  [../]
  [./disp_y]
    initial_condition = 0
  [../]
  [./disp_z]
    initial_condition = 0
  [../]
[]
# DirecKernels
# [DiracKernels]
#   [./pump_in]
#     type = TigerHydraulicPointSourceH
#     point = '5000 5000 -5000'
#     mass_flux = 100.0 # Mass Flow (kg/s)
#     variable = pressure
#   [../]
#   # [./pump_out]
#   #   type = TigerHydraulicPointSourceH
#   #   point = '5289.5 5000 -2052.38'
#   #   mass_flux = -100.0 # Mass Flow (kg/s)
#   #   variable = pressure
#   # [../]
# []

# Kernels
[Kernels]
    #   [./gravity]
    # type = Gravity
    # use_displaced_mesh = false
    # variable = disp_z
    # value = -10e-6 # remember this is in MPa
    # [../]
  [./gravity_z]
    type = TigerMechanicsGravityM
    variable = 'disp_z'
    component = 2
    use_displaced_mesh = false
  [../]
  [./hdiff]
    type = TigerHydraulicKernelH
    variable = pressure
  [../]
  [./htime]
    type = TigerHydraulicTimeKernelH
    variable = pressure
  [../]
  [./T_diff]
    type = TigerThermalDiffusionKernelT
    variable = temperature
  [../]
  [./T_advect]
    type = TigerThermalAdvectionKernelT
    variable = temperature
    pressure = pressure
  [../]
  [./T_time]
    type = TigerThermalTimeKernelT
    variable = temperature
  [../]
  # [./hm]
  #   type = TigerHydroMechanicsKernelHM
  #   variable = pressure
  #   displacements = 'disp_x disp_y disp_z'
  # [../]
  [./poro_x]
    type = PoroMechanicsCoupling
    variable = disp_x
    porepressure = pressure
    component = 0
    use_displaced_mesh = false
  [../]
  [./poro_y]
    type = PoroMechanicsCoupling
    variable = disp_y
    porepressure = pressure
    component = 1
    use_displaced_mesh = false
  [../]
  [./poro_z]
    type = PoroMechanicsCoupling
    variable = disp_z
    porepressure = pressure
    component = 2
    use_displaced_mesh = false
  [../]
[]
[Preconditioning]
  active = 'p1'
  [./p1]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -pc_hypre_type -snes_type -snes_linesearch_type -sub_pc_factor_shift_type'
    petsc_options_value = 'hypre boomeramg newtonls basic NONZERO'
  [../]
  [./p2]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -sub_pc_type -snes_type -snes_linesearch_type -sub_pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'asm lu newtonls basic NONZERO 51'
  [../]
  [./p3]
    type = SMP
    full = true
    #petsc_options = '-snes_ksp_ew'
    petsc_options_iname = '-pc_type -ksp_type -sub_pc_type -snes_type -snes_linesearch_type -pc_asm_overlap -sub_pc_factor_shift_type -ksp_gmres_restart'
    petsc_options_value = 'asm gmres lu newtonls basic 2 NONZERO 51'
  [../]
  [./p4]
    type = FSP
    full = true
    topsplit = pT
    [./pT]
      splitting = 'p T'
      splitting_type = multiplicative
      petsc_options_iname = '-ksp_type -pc_type -snes_type -snes_linesearch_type'
      petsc_options_value = 'fgmres lu newtonls basic'
    [../]
    [./p]
      vars = 'pressure'
      petsc_options_iname = '-ksp_type -pc_type -sub_pc_type'
      petsc_options_value = 'fgmres asm ilu'
    [../]
    [./T]
      vars = 'temperature'
      petsc_options_iname = '-ksp_type -pc_type -pc_hypre_type'
      petsc_options_value = 'preonly hypre boomeramg'
    [../]
  [../]
[]
# Executioners
[Executioner]
  type = Transient
 l_tol = 1e-11
 nl_rel_tol = 1e-6
 nl_abs_tol = 1e-10
 l_max_its = 20
 nl_max_its = 20
 [./TimeStepper]
   type = IterationAdaptiveDT
   dt = 86400
   growth_factor = 10
 [../]
 dtmax = 10512000
 end_time = 315360000
 #end_time = 946080000
  automatic_scaling = true
  #auto_preconditioning = true
  solve_type = 'NEWTON'
[]
#Steady_state detection
# [Executioner]
#   type = Transient
#   dt = 86400
#   #num_steps = 30
#   #steady_state_start_time = -1000
#   steady_state_detection = true
#   steady_state_tolerance = 1e-6
#   check_aux = true
# []
[Postprocessors]
  [./dt]
    type = TimestepSize
  [../]
  [./max_dP]
    type = ElementExtremeValue
    variable = dP
    execute_on = 'initial timestep_end'
    block = 'middle right1 right2 left1 left2'
  [../]
  [./max_stress_excess]
    type = ElementExtremeValue
    variable = stress_excess
    execute_on = 'initial timestep_end'
    block = 'middle right1 right2 left1 left2'
  [../]
[]

[UserObjects]
  # [./terminator]
  #   type = Terminator
  #   expression = 'max_stress_excess > 0.9'
  #   fail_mode = SOFT
  #   execute_on = ' timestep_end'
  #   enable = True
  # []
  # [./terminator2]
  #   type = Terminator
  #   expression = 'dt < 3600'
  #   fail_mode = HARD
  #   execute_on = 'timestep_end'
  #   enable = True
  # []
[]


# [Postprocessors]
#   [./production_temp]
#     type = PointValue
#     variable = temperature
#     point = '4856.916 5000 -1988.58'
#   [../]
#   [./Injection_temp]
#     type = PointValue
#     variable = temperature
#     point = '4019.506 5000 -1961.17'
#   [../]
#   [./production_pressure]
#     type = PointValue
#     variable = pressure
#     point = '4856.916 5000 -1988.58'
#   [../]
#   [./Injection_pressure]
#     type = PointValue
#     variable = pressure
#     point = '4019.506 5000 -1961.17'
#   [../]
# []
# Outputs
[Outputs]
  [out]
    type = Exodus #= true
  elemental_as_nodal = true
  []
  csv = true
  checkpoint = true
  print_linear_residuals = false
  print_nonlinear_converged_reason = true
  print_linear_converged_reason = false
[]
