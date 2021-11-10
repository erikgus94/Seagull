# to run:     ../seagull-opt -i boussinesq5_thin_detailed_no_slip.i
# peacock:    ~/projects/moose/python/peacock/peacock -r boussinesq5_thin_detailed_no_slip_out.e


[Mesh]
  [cmg]
    type = CartesianMeshGenerator
    dim = 2
    dx = '0.0375 0.425 0.0375'
    dy = '0.0375 0.425 0.0375'
    ix = '8 90 8'
    iy = '8 90 8'
    subdomain_id = '0 0 0
                    0 1 0
                    0 0 0'
  []
  [delete_center]
    type = BlockDeletionGenerator
    input = 'cmg'
    block = 1
  []
  [add_inlet1]
    type = ParsedGenerateSideset
    input = 'delete_center'
    combinatorial_geometry = 'x<0.0375 & y>0.49999'
    new_sideset_name = 'inlet1'
  []
  [add_inlet2]
    type = ParsedGenerateSideset
    input = 'add_inlet1'
    combinatorial_geometry = 'x>0.4625 & y>0.49999'
    new_sideset_name = 'inlet2'
  []
  [add_outlet]
    type = ParsedGenerateSideset
    input = 'add_inlet2'
    combinatorial_geometry = '(x>0.4625 | x<0.0375) & y<0.00001'
    new_sideset_name = 'outlet'
  []
  [add_wall_top]
    type = ParsedGenerateSideset
    input = 'add_outlet'
    combinatorial_geometry = 'x>0.0375 & x<0.4625 & (y>0.49999 | (y<0.46250001 & y>0.3))'
    new_sideset_name = 'walls_top'
  []
  [add_wall_bot]
    type = ParsedGenerateSideset
    input = 'add_wall_top'
    combinatorial_geometry = 'x>0.0375 & x<0.4625 & (y<0.00001 | (y>0.03749999 & y<0.1))'
    new_sideset_name = 'walls_bottom'
  []
  [add_wall_left]
    type = ParsedGenerateSideset
    input = 'add_wall_bot'
    combinatorial_geometry = '(x<0.00001 | (y>0.0374999 & y<0.46249999 & x>0.03749999 & x<0.2))'
    new_sideset_name = 'walls_left'
  []
  [add_wall_right]
    type = ParsedGenerateSideset
    input = 'add_wall_left'
    combinatorial_geometry = '(x>0.49999 | (y>0.0374999 & y<0.46249999 & x>0.3 & x<0.46250001))'
    new_sideset_name = 'walls_right'
    show_info = true
  []
[]

mu = 0.001241                   # The viscosity [Pa*s]
rho = 1611                      # Density [kg/m³]
k = 1.10                        # Thermal conductivity [W/m*K]
cp = 2097.8                     # Specific heat capacity [J/kg*K]
alpha = 0.0003035               # Thermal expansion
temp_ref = 873.15               # Reference temperature for the thermal expansion law
#von_karman_const = 0.4

vel = 'velocity'
velocity_interp_method = 'rc'
advected_interp_method = 'upwind'
cold_temp = 660
hot_temp = 740

# initial conditions
IC_u = 6.6e-3
IC_v = 5.8e-3

IC_pressure = 4188

[Variables]
  [u]
    type = INSFVVelocityVariable
  []
  [v]
    type = INSFVVelocityVariable
  []
  [pressure]
    type = INSFVPressureVariable
  []
  [T]
    type = INSFVEnergyVariable
  []
  [lambda]
    family = SCALAR
    order = SECOND
    initial_condition = 4.299839e-10
  []
[]

[AuxVariables]
  [U]
    order = CONSTANT
    family = MONOMIAL
    fv = true
  []
  [vel_x]
    order = SECOND
    family = MONOMIAL
  []
  [vel_y]
    order = SECOND
    family = MONOMIAL
  []
  [viz_T]
    order = SECOND
    family = MONOMIAL
  []
  [mixing_len]
    type = MooseVariableFVReal
  []
[]

[ICs]
  [temperature_gradient]
    type = FunctionIC
    variable = T
    function = 'gradient_temp'
  []
  [v_gradient]
    type = FunctionIC
    variable = v
    function = 'gradient_v'
  []
  [u_gradient]
    type = FunctionIC
    variable = u
    function = 'gradient_u'
  []
  [pressure_gradient]
    type = FunctionIC
    variable = pressure
    function = 'gradient_pressure'
  []
[]

[Functions]
  # linear temperature profile for initialization
  [gradient_temp]
    type = ParsedFunction
    value = 'if(y<0.0375, hot_temp, if(y>0.4625, cold_temp, hot_temp - (hot_temp-cold_temp) * y / 0.4625))'
    vars = 'cold_temp hot_temp'
    vals = '${cold_temp} ${hot_temp}'
  []
  [gradient_v]
    type = ParsedFunction       #0.0046875
    value = 'if(x>0.0046875 & x<0.0375-0.0046875  & y>0.0375 & y<0.4625, IC_v,
    if(x>0.4625+0.0046875   & x<0.5-0.0046875     & y>0.0375 & y<0.4625, -IC_v,
    if(x>0.0046875          & x<0.0375-0.0046875  & y>0.009  & y<0.0375, IC_v * y / 0.0375,
    if(x>0.4625+0.0046875   & x<0.5-0.0046875     & y>0.009  & y<0.0375, -IC_v * y / 0.0375,
    if(x>0.0046875          & x<0.0375-0.0046875  & y>0.4625 & y<0.49, IC_v - IC_v * (y - 0.4625) / 0.0375,
    if(x>0.4625+0.0046875   & x<0.5-0.0046875     & y>0.4625 & y<0.49, -IC_v + IC_v * (y - 0.4625) / 0.0375,
    0))))))'
    vars = 'IC_v'
    vals = '${IC_v}'
  []
  [gradient_u]
    type = ParsedFunction
    value = 'if(y>0.0046875 & y<0.0375-0.0046875  & x>0.0375  & x<0.4625, -IC_u * (1 - (y/0.0375)/3),
    if(y>0.4625+0.0046875   & y<0.5-0.0046875     & x>0.0375  & x<0.4625, IC_u * (2/3 - ((0.4625 - y)/0.0375)/3),
    if(y>0.0046875          & y<0.0375-0.0046875  & x>0.009   & x<0.0375, -IC_u * x / 0.0375,
    if(y>0.4625+0.0046875   & y<0.5-0.0046875     & x>0.009   & x<0.0375, IC_u * x / 0.0375,
    if(y>0.0046875          & y<0.0375-0.0046875  & x>0.4625  & x<0.49, -IC_u + IC_u * (x - 0.4625) / 0.0375,
    if(y>0.4625+0.0046875   & y<0.5-0.0046875     & x>0.4625 & x<0.49, IC_u - IC_u * (x - 0.4625) / 0.0375,
    0))))))'
    vars = 'IC_u'
    vals = '${IC_u}'
  []
  [gradient_pressure]
    type = ParsedFunction
    value = 'IC_pressure - 2 * IC_pressure * y / 0.5'
    vars = 'IC_pressure'
    vals = '${IC_pressure}'
  []
[]

[AuxKernels]
  [mag]
    type = VectorMagnitudeAux
    variable = U
    x = u
    y = v
    execute_on = 'initial timestep_end'
  []
  [vel_x]
    type = ParsedAux
    variable = vel_x
    function = 'u'
    execute_on = 'initial timestep_end'
    args = 'u'
  []
  [vel_y]
    type = ParsedAux
    variable = vel_y
    function = 'v'
    execute_on = 'initial timestep_end'
    args = 'v'
  []
  [viz_T]
    type = ParsedAux
    variable = viz_T
    function = 'T'
    execute_on = 'initial timestep_end'
    args = 'T'
  []
  inactive = 'mixing_len'
  [mixing_len]
    type = WallDistanceMixingLengthAux
    walls = 'walls_top walls_right walls_bottom walls_left inlet1 inlet2 outlet'
    variable = mixing_len
    execute_on = 'initial'
    delta = 0.5
  []
[]

[FVKernels]
  [mass]
    type = INSFVMassAdvection
    variable = pressure
    vel = ${vel}
    advected_interp_method = ${advected_interp_method}
    velocity_interp_method = ${velocity_interp_method}
    u = u
    v = v
    pressure = pressure
    mu = 'mu'
    rho = ${rho}
  []
  [mean_zero_pressure]
    type = FVScalarLagrangeMultiplier
    variable = pressure
    lambda = lambda
  []

  [u_time]
    type = INSFVMomentumTimeDerivative
    variable = 'u'
    rho = ${rho}
  []
  [u_advection]
    type = INSFVMomentumAdvection
    variable = u
    advected_quantity = 'rhou'
    vel = ${vel}
    velocity_interp_method = ${velocity_interp_method}
    advected_interp_method = ${advected_interp_method}
    pressure = pressure
    u = u
    v = v
    mu = 'mu'
    rho = ${rho}
  []
  [u_viscosity]
    type = FVDiffusion
    variable = u
    coeff = ${mu}
  []
  [u_pressure]
    type = INSFVMomentumPressure
    variable = u
    momentum_component = 'x'
    pressure = pressure
  []
  [u_turb]
    type = INSFVMixingLengthReynoldsStress
    variable = u
    rho = ${rho}
    mixing_length = mixing_len
    momentum_component = 'x'
    u = u
    v = v
  []

  [v_time]
    type = INSFVMomentumTimeDerivative
    variable = v
    rho = ${rho}
  []
  [v_advection]
    type = INSFVMomentumAdvection
    variable = v
    advected_quantity = 'rhov'
    vel = ${vel}
    velocity_interp_method = ${velocity_interp_method}
    advected_interp_method = ${advected_interp_method}
    pressure = pressure
    u = u
    v = v
    mu = 'mu'
    rho = ${rho}
  []
  [v_viscosity]
    type = FVDiffusion
    variable = v
    coeff = ${mu}
  []
  [v_pressure]
    type = INSFVMomentumPressure
    variable = v
    momentum_component = 'y'
    pressure = pressure
  []
  [v_turb]
    type = INSFVMixingLengthReynoldsStress
    variable = v
    rho = ${rho}
    mixing_length = mixing_len
    momentum_component = 'y'
    u = u
    v = v
  []
  [v_buoyancy]
    type = INSFVMomentumBoussinesq
    variable = v
    T_fluid = T
    gravity = '0 -9.82 0'
    rho = ${rho}
    ref_temperature = ${temp_ref}
    momentum_component = 'y'
  []
  [v_gravity]
    type = INSFVMomentumGravity
    variable = v
    gravity = '0 -9.82 0'
    rho = ${rho}
    momentum_component = 'y'
  []

  [temp_time]
    type = INSFVEnergyTimeDerivative
    variable = T
    rho = ${rho}
    cp_name = 'cp'
  []
  [temp_conduction]
    type = FVDiffusion
    coeff = 'k'
    variable = T
  []
  [temp_advection]
    type = INSFVEnergyAdvection
    variable = T
    vel = ${vel}
    velocity_interp_method = ${velocity_interp_method}
    advected_interp_method = ${advected_interp_method}
    pressure = pressure
    u = u
    v = v
    mu = 'mu'
    rho = ${rho}
  []
  [temp_turb]
    type = WCNSFVMixingLengthEnergyDiffusion
    variable = T
    rho = ${rho}
    cp = 'cp'
    mixing_length = mixing_len
    schmidt_number = 0.5
    u = u
    v = v
  []
[]

[ScalarKernels]
  [mean_zero_pressure_lm]
    type = AverageValueConstraint
    variable = lambda
    pp_name = pressure_integral
    value = 0
  []
[]

[FVBCs]
  [no_slip_x]
    type = INSFVNoSlipWallBC
    variable = u
    boundary = 'walls_top walls_right walls_bottom walls_left inlet1 inlet2 outlet'
    function = 0
  []
  [no_slip_y]
    type = INSFVNoSlipWallBC
    variable = v
    boundary = 'walls_top walls_right walls_bottom walls_left inlet1 inlet2 outlet'
    function = 0
  []

  [T_hot]
    type = FVDirichletBC
    variable = 'T'
    boundary = 'walls_bottom'
    value = ${hot_temp}
  []

  [T_cold]
    type = FVDirichletBC
    variable = 'T'
    boundary = 'walls_top'
    value = ${cold_temp}
  []
[]

[Functions]
  [mu_rampdown]
    type = PiecewiseLinear
    x = '0.1 1 2 3'
    y = '${fparse 1e3*mu} ${fparse 1e2*mu} ${fparse 1e1*mu} ${mu}'
  []
[]

[Postprocessors]
  # To watch the rampdown
  [mu]
    type = FunctionValuePostprocessor
    function = mu_rampdown
  []
  # To evaluate convergence (visual inspection still necessary!)
  [average_temperature]
    type = ElementAverageValue
    variable = T
  []
  [average_velocity_v]
    type = ElementAverageValue
    variable = v
  []
  [average_velocity_u]
    type = ElementAverageValue
    variable = u
  []
  [average_pressure]
    type = ElementAverageValue
    variable = pressure
  []
  [pressure_integral]
    type = ElementIntegralVariablePostprocessor
    variable = pressure
    execute_on = linear
  []
[]

[Materials]
  [const_functor]
    type = ADGenericConstantFunctorMaterial
    prop_names = 'cp k alpha_b'
    prop_values = '${cp} ${k} ${alpha}'
  []
  [mu_functor]
    type = ADGenericFunctionFunctorMaterial
    prop_names = 'mu'
    prop_values = 'mu_rampdown'
  []
  [ins_fv]
    type = INSFVMaterial
    u = 'u'
    v = 'v'
    pressure = 'pressure'
    temperature = 'T'
    rho = ${rho}
  []
[]

[Executioner]
  type = Transient
  dt = 0.01
  dtmax = 1
  # 1s after the end of the mu_rampdown
  # going further could be required to relax the system to steady state further
  start_time = 0
  end_time = 6
  timestep_tolerance = 1e-5
  solve_type = 'NEWTON'
  petsc_options_iname = '-pc_type -pc_asm_overlap -ksp_gmres_restart -sub_pc_type -sub_pc_factor_shift_type'
  petsc_options_value = 'asm              2              1200                lu           NONZERO'
  nl_rel_tol = 1e-5
  nl_abs_tol = 1e-5
  line_search = 'none'
  automatic_scaling = true
  off_diagonals_in_auto_scaling = true
[]
[Debug]
  show_var_residual_norms = true
[]
[Outputs]
  #file_base = 'reduced_new'
  exodus = true
[]
