%% Koopman MPC controller setup
function ocp = kmpc_setup_acados(mpc_setup)
prediction_model = kmpc_prediction_model(mpc_setup);

%% OCP solver model
ocp_model = acados_ocp_model();
ocp_model.set('name', 'kmpc_acados');
N = mpc_setup.N;
Ts = mpc_setup.Ts;
ocp_model.set('T', Ts);  % [s] prediction horizon (just one step in this formulation)

% symbolics
ocp_model.set('sym_x', prediction_model.sym_x);
ocp_model.set('sym_u', prediction_model.sym_u);
ocp_model.set('sym_p',prediction_model.sym_p);

% cost
cost_type = 'auto';  % auto, linear_ls, nonlinear_ls, ext_cost
ocp_model.set('cost_type_0', cost_type);
ocp_model.set('cost_type', cost_type);
ocp_model.set('cost_type_e', cost_type);
ocp_model.set('cost_expr_ext_cost_0', prediction_model.cost_expr_ext_cost_0);
ocp_model.set('cost_expr_ext_cost', prediction_model.cost_expr_ext_cost_0);
ocp_model.set('cost_expr_ext_cost_e', prediction_model.cost_expr_ext_cost_e);

% dynamics
ocp_model.set('dyn_type', 'discrete');
ocp_model.set('dyn_expr_phi', prediction_model.dyn_expr_phi);

% constraints
ocp_model.set('constr_x0', zeros(length(prediction_model.sym_x),1));  % initialize with zeros, change later

ocp_model.set('constr_expr_h_0', prediction_model.constr_expr_h_0);
ocp_model.set('constr_lh_0', prediction_model.constr_lh_0);
ocp_model.set('constr_uh_0', prediction_model.constr_uh_0);

% Jbu = zeros(nu); for ii=1:nu Jbu(ii,ii)=1.0; end
% lbu = 0*ones(nu, 1);
% ubu = T_max*ones(nu, 1);
% ocp_model.set('constr_Jbu', Jbu);
% ocp_model.set('constr_lbu', lbu);
% ocp_model.set('constr_ubu', ubu);

%Jbx = zeros(nx); for ii=1:nx Jbx(ii,ii)=1.0; end
%lbx = -ones(nx, 1);
%ubx = ones(nx, 1);
%	ocp_model.set('constr_Jbx', Jbx);
%	ocp_model.set('constr_lbx', lbx);
%	ocp_model.set('constr_ubx', ubx);

%% OCP solver settings
ocp_opts = acados_ocp_opts();
ocp_opts.set('param_scheme_N', 1);
ocp_opts.set('parameter_values', ones(size(prediction_model.sym_p)));  % initialize with zero, change later

% QP solver
qp_solver = 'partial_condensing_osqp';
ocp_opts.set('qp_solver', qp_solver);
% full_condensing_hpipm
% partial_condensing_hpipm
% full_condensing_qpoases
% partial_condensing_osqp

if (contains(qp_solver, 'hpipm'))
	ocp_opts.set('qp_solver_cond_N', 2);  % ?
	ocp_opts.set('qp_solver_cond_ric_alg', 1);  % 0: dont factorize hessian in the condensing; 1: factorize
	ocp_opts.set('qp_solver_ric_alg', 1);
    ocp_opts.set('qp_solver_iter_max', 50);
end

if (contains(qp_solver, 'osqp'))
    ocp_opts.set('qp_solver_warm_start', 2);
    ocp_opts.set('nlp_solver_warm_start_first_qp',1);
    ocp_opts.set('qp_solver_iter_max', 4000);
end

ocp_opts.set('sim_method', 'discrete');

%% Simulink block settings
simulink_opts = get_acados_simulink_opts;

% inputs
simulink_opts.inputs.y_ref_0 = 0;
simulink_opts.inputs.y_ref = 0;
simulink_opts.inputs.y_ref_e = 0;
simulink_opts.inputs.lbx = 0;
simulink_opts.inputs.ubx = 0;
simulink_opts.inputs.lbx_e = 0;
simulink_opts.inputs.ubx_e = 0;
simulink_opts.inputs.lbu = 0;
simulink_opts.inputs.ubu = 0;
simulink_opts.inputs.lg = 0;
simulink_opts.inputs.ug = 0;
simulink_opts.inputs.lh = 0;
simulink_opts.inputs.uh = 0;
simulink_opts.inputs.lh_e = 0;
simulink_opts.inputs.uh_e = 0;
% simulink_opts.inputs.x_init = 1;  ?
% simulink_opts.inputs.reset_solver = 1;  ?

% outputs
% simulink_opts.outputs.utraj = 1;
% simulink_opts.outputs.xtraj = 1;
simulink_opts.outputs.cost_value = 1;
simulink_opts.outputs.KKT_residual = 0;
% simulink_opts.outputs.KKT_residuals = 1;  ?
simulink_opts.outputs.x1 = 0;  % ?
simulink_opts.outputs.sqp_iter = 0;

simulink_opts.samplingtime = 't0';
    % 't0' (default) - use time step between shooting node 0 and 1
    % '-1' - inherit sampling time from other parts of simulink model

%% Create the OCP solver
ocp = acados_ocp(ocp_model, ocp_opts, simulink_opts);
disp(['Created the OCP solver.',newline,newline])

%% Compile for Simulink if needed
% if mpc_setup.compile_for_simulink
%     cd setup/c_generated_code
%     make_sfun;
% end