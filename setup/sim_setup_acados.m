%% Simulation model setup
function simulation = sim_setup_acados(Ts,mu_x)
sim_model = acados_sim_model();
model = simulation_model_acados();

% symbolics
sim_model.set('sym_x', model.sym_x);
sim_model.set('sym_u', model.sym_u);
sim_model.set('sym_xdot', model.sym_xdot);
sim_model.set('sym_p', model.sym_p);
sim_model.set('dyn_type', 'implicit');
sim_model.set('dyn_expr_f', model.expr_f_impl);

sim_model.set('T', Ts);

% acados sim opts
sim_opts = acados_sim_opts();
sim_opts.set('num_steps', 4);
sim_opts.set('parameter_values', mu_x);  % initial mu_x

% create sim
simulation = acados_sim(sim_model, sim_opts);

% add parameters
simulation.model_struct.R = model.wheel_radius;   % [m] wheel radius
simulation.model_struct.T_max = model.max_torque;  % [Nm] maximum wheel torque

disp(['Created the simulation model.',newline])
end