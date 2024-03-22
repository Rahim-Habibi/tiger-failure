[Mesh]
  [./file]
    type = FileMeshGenerator
    file = box.msh
  [../]
  #allow_renumbering = false
[]
[GlobalParams]
  displacements = 'disp_x disp_y disp_z'
  #thermal_expansion_coeff = 1e-5 # 8e-6  engineering toolbox # linear
  PorousFlowDictator = dictator
  gravity = '0 0 -9.8'
  biot_coefficient = 1.0
[]
[FluidProperties]
  [./the_simple_fluid]
    type = SimpleFluidProperties
    thermal_expansion = 0.0
    bulk_modulus = 1e9
    viscosity = 1.0e-3
    density0 = 1000.0
    cv = 4000.0   #Constant specific heat capacity at constant volume (J/kg/K)
    cp = 1000.0   #Constant specific heat capacity at constant pressure (J/kg/K)
    porepressure_coefficient = 1.0  #The enthalpy is internal_energy + P / density * porepressure_coefficient.
  [../]
 []
# Permeability Assigning
[UserObjects]
  [./dictator]
  type = PorousFlowDictator
  porous_flow_vars = 'pressure  temperature disp_x disp_y disp_z'
    number_fluid_phases = 1
    number_fluid_components = 1
  [../]
  [pc]
    type = PorousFlowCapillaryPressureConst
    pc = 10000000
  []
[]
[AuxVariables]
  [./stress_xx]
order = CONSTANT
family = MONOMIAL
[../]
[./stress_xy]
order = CONSTANT
family = MONOMIAL
[../]
[./stress_xz]
order = CONSTANT
family = MONOMIAL
[../]
[./stress_yx]
order = CONSTANT
family = MONOMIAL
[../]
[./stress_yy]
order = CONSTANT
family = MONOMIAL
[../]
[./stress_yz]
order = CONSTANT
family = MONOMIAL
[../]
[./stress_zx]
order = CONSTANT
family = MONOMIAL
[../]
[./stress_zy]
order = CONSTANT
family = MONOMIAL
[../]
[./stress_zz]
order = CONSTANT
family = MONOMIAL
[../]
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
  [./fcn]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./effective_fluid_pressure]
    family = MONOMIAL
    order = CONSTANT
  [../]
  [./porosity]
  family = MONOMIAL
  order = CONSTANT
  [../]
  [mass_frac_phase0_species0]
    initial_condition = 1 # all water in phase=0
  []
  # [mass_frac_phase1_species0]
  #   initial_condition = 1 # no water in phase=1
  # []
  [sgas]
    family = MONOMIAL
    order = CONSTANT
  []
  [swater]
    family = MONOMIAL
    order = CONSTANT
  []
[]
[AuxKernels]
  [./stress_xx]
type = RankTwoAux
rank_two_tensor = stress
variable = stress_xx
index_i = 0
index_j = 0
[../]
[./stress_xy]
type = RankTwoAux
rank_two_tensor = stress
variable = stress_xy
index_i = 0
index_j = 1
[../]
[./stress_xz]
type = RankTwoAux
rank_two_tensor = stress
variable = stress_xz
index_i = 0
index_j = 2
[../]
[./stress_yx]
type = RankTwoAux
rank_two_tensor = stress
variable = stress_yx
index_i = 1
index_j = 0
[../]
[./stress_yy]
type = RankTwoAux
rank_two_tensor = stress
variable = stress_yy
index_i = 1
index_j = 1
[../]
[./stress_yz]
type = RankTwoAux
rank_two_tensor = stress
variable = stress_yz
index_i = 1
index_j = 2
[../]
[./stress_zx]
type = RankTwoAux
rank_two_tensor = stress
variable = stress_zx
index_i = 2
index_j = 0
[../]
[./stress_zy]
type = RankTwoAux
rank_two_tensor = stress
variable = stress_zy
index_i = 2
index_j = 1
[../]
[./stress_zz]
type = RankTwoAux
rank_two_tensor = stress
variable = stress_zz
index_i = 2
index_j = 2
[../]
  [./p0_ker]
    type = FunctionAux
    variable = 'p0'
    function = '(1000*9.81*(500-(z+200)))'
    execute_on = 'initial'
  [../]
  # [./dP_ker]
  #   type = ParsedAux
  #   variable = 'dP'
  #   args = 'pressure p0'
  #   function = 'pressure-p0'
  #   execute_on = 'TIMESTEP_END'
  # [../]

  [./T0_ker]
    type = FunctionAux
    variable = 'T0'
    function = '((0.03*(500-z))+283.15)'
    execute_on = 'initial'
  [../]
  # [./dT_ker]
  #   type = ParsedAux
  #   variable = 'dT'
  #   args = 'temperature T0'
  #   function = 'temperature-T0'
  #   execute_on = 'TIMESTEP_END'
  # [../]
  # [./terminator]
  #   type = TigerTerminator
  #   total_stress = stress
  #   criterion_type = Slip_Tedency
  #   cohesion = 0
  #   phi = 35
  #   variable = fcn
  #   fail_mode = SOFT
  #   use_displaced_mesh = false
  #   execute_on = ' TIMESTEP_END'
  #   block = '11 12'
  # [../]
  [./effective_fluid_pressure]
  type = ParsedAux
  coupled_variables = 'pressure  swater '
    expression = 'pressure * swater '
    variable = effective_fluid_pressure
  [../]
  # [swater]
  #   type = PorousFlowPropertyAux
  #   variable = swater
  #   property = saturation
  #   phase = 0
  #   execute_on = timestep_end
  # []
  # [sgas]
  #   type = PorousFlowPropertyAux
  #   variable = sgas
  #   property = saturation
  #   phase = 1
  #   execute_on = timestep_end
  # []
  # [./porosity]
  # type = PorousFlowPropertyAux
  # variable = porosity
  # property = porosity
  # execute_on = timestep_end
  # [../]
[]
# Materials Assigning
[Materials]
  [./temperature]
    type = PorousFlowTemperature  #to provide temperature at the quadpoints or nodes and derivatives of it with respect to the PorousFlow variables
    #temperature = temperature
    PorousFlowDictator = dictator
    #at_nodes = True
  [../]
  [massfrac]
    type = PorousFlowMassFraction
  []

 [relperm_water]
   type = PorousFlowRelativePermeabilityConst
   phase = 0
 []

  [./water]
    type = PorousFlowSingleComponentFluid  #This Material calculates fluid properties at the quadpoints or nodes for a single component fluid
    fp = the_simple_fluid
    phase = 0
    #at_nodes = true
  [../]
  # [./normalshear]
  #   type = TigerNormalShearStressM
  #   total_stress = stress
  #   DoUWantShear = true
  #   block = '11 12'
  #   #outputs = 'exodus'
  #   []
  [./Elasticity_tensor]
    type = ComputeIsotropicElasticityTensor  #Compute a constant isotropic elasticity tensor.
    youngs_modulus = 15e10
    poissons_ratio = 0.0
    block = 'box middle right1 right2 left1 left2'
  [../]
  [./strain]
    type = ComputeSmallStrain
    eigenstrain_names = ' initial_stress'  #List of eigenstrains to be applied in this strain calculation
  [../]
  [./stress]
    type = ComputeLinearElasticStress # Compute stress using elasticity for small strains
    block = 'box middle right1 right2 left1 left2'
  [../]
  [thermal_expansion]
  type = PorousFlowConstantThermalExpansionCoefficient  #Computes the effective thermal expansion coefficient
  fluid_coefficient = 5e-6
  drained_coefficient = 2e-4
   []
  [./initial_strain]
    type = ComputeEigenstrainFromInitialStress
    initial_stress = 'ini_xx 0 0  0 ini_yy 0  0 0 ini_zz'
    eigenstrain_name = initial_stress  # it is called by "strain" to compute eigenstrain resulted from initial stress
    # output_properties = 'initial_stress'
    # outputs = exodus
  [../]
  [saturation_calculator]
    type = PorousFlow1PhaseFullySaturated #This Material is used for the fully saturated single-phase situation where porepressure is the primary variable
    porepressure = pressure
  []
  # [ppss]
  #   type = PorousFlow1PhaseP   #used for the fully saturated single-phase situation where porepressure is the primary variable, it can be replaced by "PorousFlow1PhaseFullySaturated" when you do not know the capillar pressure
  #   porepressure = pressure
  #   capillary_pressure = pc
  # []
  [./effective_fluid_pressure_mat]
      type = PorousFlowEffectiveFluidPressure  #calculates an effective fluid pressure: effective_stress = total_stress + biot_coeff*effective_fluid_pressure. The effective_fluid_pressure = sum_{phases} and provides effective pressure for PoroMechanicsCoupling
      at_nodes = true
      block = 'box middle right1 right2 left1 left2'
  [../]
  [./effective_fluid_pressure_qp]
      type = PorousFlowEffectiveFluidPressure  #calculates an effective fluid pressure: effective_stress = total_stress + biot_coeff*effective_fluid_pressure. The effective_fluid_pressure = sum_{phases} and provides effective pressure for PoroMechanicsCoupling
      block = 'box middle right1 right2 left1 left2'
  [../]
  [./volumetric_strain]
      type = PorousFlowVolumetricStrain  # computes volumetric strain and volumetric strain rate linked to TensorMechanics strain calculator
      block = 'box middle right1 right2 left1 left2'
  [../]
  #Scaling Geometry
  [./rock_g0]
    type = TigerGeometryMaterial
     scale_factor = 1 # 100
     block = 'box'
  [../]
  [./rock_fault]
    type = TigerGeometryMaterial
     scale_factor = 0.01 # Fault Apperture
     block = ' middle right1 right2 left1 left2'
  [../]
  #Prorosity
  [./rock_p0]
    type = PorousFlowPorosityConst
    #fluid = true       # If true, porosity will be a function of effective porepressure, volumetric_strain and temperature
    porosity = 0.01
    #mechanical = true
    #thermal = true
    # reference_temperature = 'T0'
    # reference_porepressure = 20e6
    block = 'box'
  [../]
  [./rock_p1]
    type = PorousFlowPorosityConst
    porosity = 0.15
    # fluid = true
    # mechanical = true
    # thermal = true
    # reference_temperature = 'T0'
    # reference_porepressure = 20e6
    block = ' middle right1 right2 left1 left2'
  [../]
  #Permeability     the value and type of permeability has huge effect on residual
  [./rock_h0]
    type = PorousFlowPermeabilityConst
    permeability = '0 0 0  0 1e-15 0 0 0 0'
    block = 'box'
  [../]
  [./FaultActivation]
    type = PorousFlowPermeabilityConst
       permeability = '0 0 0  0 1e-10 0  0 0 0'
    block = ' middle right1 right2 left1 left2'
  [../]
  #ThermalConductivity
  [./ThermalConductivity]
    type = PorousFlowThermalConductivityIdeal  #calculates rock-fluid combined thermal conductivity
    dry_thermal_conductivity = '1 0 0  0 1 0  0 0 1'
    block = ' box'
  [../]
  [./ThermalFaults]
    type = PorousFlowThermalConductivityIdeal
    dry_thermal_conductivity = '1 0 0  0 1 0  0 0 1'
    block = ' middle right1 right2 left1 left2'
  [../]
  [./rock_internal_energy_sediment]
   type = PorousFlowMatrixInternalEnergy   #calculates the internal energy of solid rock grains, which is specific_heat_capacity * density * temperature
    specific_heat_capacity = 1000
    density = 2300
    block = 'box'
  [../]
  [./rock_internal_energy_base]
   type = PorousFlowMatrixInternalEnergy
    specific_heat_capacity = 1050
    density = 2400
    block = ' middle right1 right2 left1 left2'
  [../]
  [biot_modulus]
  type = PorousFlowConstantBiotModulus  #Computes the Biot Modulus
  fluid_bulk_modulus = 1e9
  solid_bulk_compliance = 1e-8
[]
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
  [./overburden]
    type = Pressure
    boundary = 'top'
    variable = disp_z
    component = 2
    function = ini_zz_force
  [../]
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
  [./t_inject2]
    type = DirichletBC
    variable = temperature
    boundary = 'inj'
    value = 283.15 # 60 C
  [../]
  # [./production]
  #   type = PorousFlowSink
  #   variable = pressure
  #   flux_function = 1e-2
  #   boundary = 'Right' # could be an area or line
  # []
  # [./injection]
  #   type = PorousFlowSink
  #   variable = pressure
  #   fluid_phases = 0
  #   flux_function = -1e-2
  #   boundary = 'inj' # could be an area or line
  # []
  []
[Functions]
  [./constrain_effective_fluid_pressure]
    type = ParsedFunction
    symbol_names = effective_fluid_pressure_at_wellbore
    symbol_values = effective_fluid_pressure_at_wellbore
    expression = 'max(effective_fluid_pressure_at_wellbore, 20e6)'
  [../]
  [./ini_xx]
    type = ParsedFunction
    value = '-((-z + 500) * 9.81 * 2500) *  0.65'
  [../]
  [./ini_yy]
    type = ParsedFunction
    value = '-((-z + 500) * 9.81 * 2500) * 1.25'
  [../]
  [./ini_zz_force]
    type = ParsedFunction
    value = '((-z + 500) * 9.81 * 2500)'
  [../]
  [./ini_zz]
    type = ParsedFunction
    value = '-((-z + 500) * 9.81 * 2500)'
  [../]
  [./hydrostatic]
    type = ParsedFunction
    value = '(1000*9.81*(500-(z+200)))'
    #function = '(1000*9.81*(500-(z+200)))'
  [../]
  [./Thermal]
    type = ParsedFunction
    value = '((0.03*(500-z))+283.15)'
  [../]
[]
[ICs]
  [./hydrostatic_ic]
    type = FunctionIC
    variable = p0
    function = hydrostatic
    block = 'box middle right1 right2 left1 left2'
  [../]
  [./temperature_ic]
    type = FunctionIC
    variable = temperature
    function = Thermal
    block = 'box middle right1 right2 left1 left2'
  [../]
[]
# Variable Definitions
[Variables]
  [./pressure]
    scaling = 1e-6
  [../]
  [./temperature]
    #scaling = 1e-7
  [../]
  [./disp_x]
    scaling = 1e-15
  [../]
  [./disp_y]
    scaling = 1e-16
  [../]
  [./disp_z]
    scaling = 1e-20
  [../]
[]

# Kernels
[Kernels]
  # [gravity]   # needs density block
  #   type = Gravity
  #   use_displaced_mesh = false
  #   variable = disp_z
  #   value = -10e6 # remember this is in MPa
  # []
  [./weight]
    type = BodyForce
    variable = disp_z
    value = 25000 # this is density*gravity
  [../]
  # [mass_water_dot]
  #  type = PorousFlowFullySaturatedMassTimeDerivative #Derivative of fluid-component mass with respect to time, it works on qp rather than node
  #  variable = pressure
  #  #fluid_component = 0
  #  multiply_by_density = false
  #  coupling_type = ThermoHydroMechanical
  # []
 [flux_water]
   type = PorousFlowAdvectiveFlux #advective flux of the component given by fluid_component
   use_displaced_mesh = false
   fluid_component = 0
   variable = pressure
 []
 # [vol_strain_rate_water]
 #   type = PorousFlowMassVolumetricExpansion #Energy-density*rate_of_solid_volumetric_expansion (node)
 #   variable = pressure
 #   fluid_component = 0
 # []

 # [energy_dot]
 #   type = PorousFlowEnergyTimeDerivative #Derivative of heat-energy-density wrt time
 #   variable = temperature
 #  []
 # [advection]
 #   type = PorousFlowHeatAdvection # heat flux advected by the fluid
 #   use_displaced_mesh = false
 #   variable = temperature
 #  []
 [conduction]
   type = PorousFlowHeatConduction #Heat conduction
   use_displaced_mesh = false
   variable = temperature
   []
 # [vol_strain_rate_heat]
 #   type = PorousFlowHeatVolumetricExpansion #Energy-density*rate_of_solid_volumetric_expansion (node)
 #   variable = temperature
 #   []
 [grad_stress_x]
   type = StressDivergenceTensors
   temperature = temperature
   variable = disp_x
   eigenstrain_names = strain
   use_displaced_mesh = false
   component = 0
   []
 [poro_x]
   type = PorousFlowEffectiveStressCoupling #Implements the weak form of the expression biot_coefficient * grad(effective fluid pressure)
   variable = disp_x
   use_displaced_mesh = false
   component = 0
   []
 [grad_stress_y]
   type = StressDivergenceTensors
   temperature = temperature
   variable = disp_y
   eigenstrain_names = strain
   use_displaced_mesh = false
   component = 1
   []
 [poro_y]
   type = PorousFlowEffectiveStressCoupling
   variable = disp_y
   use_displaced_mesh = false
   component = 1
   []
 [grad_stress_z]
   type = StressDivergenceTensors
   temperature = temperature
   variable = disp_z
   eigenstrain_names = strain
   use_displaced_mesh = false
   component = 2
   []
 [poro_z]
   type = PorousFlowEffectiveStressCoupling
   variable = disp_z
   use_displaced_mesh = false
   component = 2
 []
[]
[Debug]
  show_var_residual_norms = true
[]
[Preconditioning]
  active = 'p1'
    [basic]
      type = SMP
      full = true
      petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
      petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
      petsc_options_value = ' asm      lu           NONZERO                   2'
    []
    [preferred_but_might_not_be_installed]
      type = SMP
      full = true
      petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
      petsc_options_value = ' lu       mumps'
    []
    [smp]   # this is Jacobian debugger uses internal PETSc functionality to create a finite differenced Jacobian matrix from the residuals and compares
    type = SMP  # it to the implemented Jacobian (usually found in computeQpJacobian and computeQpOffDiagJacobian)
    full = true
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = 'gmres      asm      lu           NONZERO                   2             '
  []
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
  [typically_efficient]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_hypre_type'
    petsc_options_value = ' hypre    boomeramg'
  []
  [strong]
    type = SMP
    full = true
    petsc_options = '-ksp_diagonal_scale -ksp_diagonal_scale_fix'
    petsc_options_iname = '-pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap'
    petsc_options_value = ' asm      ilu           NONZERO                   2'
  []
  [probably_too_strong]
    type = SMP
    full = true
    petsc_options_iname = '-pc_type -pc_factor_mat_solver_package'
    petsc_options_value = ' lu       mumps'
  []
  [ali]
    type = SMP
    full = true
  #petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol -snes_max_it'
    petsc_options_value = 'gmres      asm      lu           NONZERO                   2               1E2       1E-5        500'
  []
[]
# Executioners Transient
# [Executioner]
#   type = Transient
#  l_tol = 1e-11
#  nl_rel_tol = 1e-6
#  nl_abs_tol = 1e-10
#  l_max_its = 20
#  nl_max_its = 20
#  start_time = -10
#  [./TimeStepper]
#    type = IterationAdaptiveDT
#    dt = 10
#    growth_factor = 10
#  [../]
#  dtmax = 10512000
#  end_time = 315360000
#  #end_time = 946080000
#   automatic_scaling = true
#   #auto_preconditioning = true
#   solve_type = 'NEWTON'
# []
# Executioners Steady
[Executioner]
  type = Steady

  solve_type = 'NEWTON' # default = PJFNK | NEWTON

  l_max_its  = 20
  l_tol      = 1e-4
  nl_max_its = 500
  nl_rel_tol = 1e-12
  nl_abs_tol = 1e-6
[]
[Postprocessors]
  [./effective_fluid_pressure_at_wellbore]
  type = PointValue
  variable = effective_fluid_pressure
  point = '3639.94 5000 -1953.78'
  use_displaced_mesh = false
 [../]
 [./constrained_effective_fluid_pressure_at_wellbore]
  type = FunctionValuePostprocessor
  function = constrain_effective_fluid_pressure
[../]
[]

# Outputs
[Outputs]
  exodus = true
  print_linear_residuals = true
  print_nonlinear_converged_reason = true
  print_linear_converged_reason = true
[]
