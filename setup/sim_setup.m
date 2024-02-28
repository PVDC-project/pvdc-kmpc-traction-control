%% Simulation model setup
function simulation = sim_setup(h,mu_x)
sim_model = acados_sim_model();
model = simulation_model();

% symbolics
sim_model.set('sym_x', model.sym_x);
sim_model.set('sym_u', model.sym_u);
sim_model.set('sym_xdot', model.sym_xdot);
sim_model.set('sym_p', model.sym_p);
sim_model.set('dyn_type', 'implicit');
sim_model.set('dyn_expr_f', model.expr_f_impl);

h_sim = h;  % simulation step size
sim_model.set('T', h_sim);

% acados sim opts
sim_opts = acados_sim_opts();
sim_opts.set('num_steps', 4);
sim_opts.set('parameter_values', mu_x);  % initial mu_x

% create sim
simulation = acados_sim(sim_model, sim_opts);

disp(['Created the simulation model.',newline,newline])
end

%{
% Simulink block settings
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
%}