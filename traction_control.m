%% Traction control using NMPC and KMPC
clear all;close all;clc;

%% Environment setup
controller_type = 1;  % 0 - off, 1 - NMPC, 2 - KMPC, 3 - PID
compile_for_simulink = 0;
model_name = 'tc';

addpath('./models')             % prediction and simulation models
addpath('./setup')              % controller and simulation setup
addpath('./functions')          % utility functions
addpath('./postprocessing')    % plotting

%% Simulation setup
Ts = 2e-3;  % [s] sampling time
kappa_ref = 0.1;  % slip reference
v0 = 5;  % [km/Ts] initial car speed
v0 = v0 / 3.6;  % convert to m/s
sim_model = simulation_model();
R = sim_model.wheel_radius;   % [m] wheel radius
w0 = v0/R;  % [rad/s] initial wheel speed
T_max = sim_model.max_torque;  % [Nm] maximum wheel torque

%% Create the controller
N = 4;      % prediction horizon length
T = N*Ts;   % [s] horizon length time

if controller_type == 1
    nmpc_setup(N,Ts,R,kappa_ref,compile_for_simulink);
elseif controller_type == 2 
    kmpc_setup(N,Ts,R,kappa_ref,compile_for_simulink);
    load setup/kmpc_data.mat PU  % for input scaling
elseif controller_type == 3
    Kp = 7500;  % Proportional gain
    Ki = 1000;  % Integral gain
    Kd = 0;     % Derivative gain    
end

% state-space dimensions
ocp_nx = 3;
ocp_nu = 1;

%% Create the simulation model
% must be after FORCES controller setup, otherwise MATLAB crashes (?)
simulation = sim_setup(Ts);

%% Closed loop simulation
sim_x0 = [v0; w0];  % initial state
ocp_x0 = [w0*R-v0; w0; 0];  % wheel slip velocity, wheel speed, integral state
T_sim = 3;  % [s] simulation time
N_sim = T_sim/Ts;  % number of simulation steps
print_step = ceil(N_sim/10);  % print stats 10 times during the simulation

% state and input logging
sim_nx = sim_model.nx;
x_sim = nan(sim_nx, N_sim+1);  % simulated state log
x_sim(:,1) = sim_x0;
x_ocp = nan(ocp_nx, N_sim+1);  % predicted state log
x_ocp(:,1) = ocp_x0;
u_sim = nan(ocp_nu, N_sim);  % control input log

% NMPC performance logging
solve_time_log = nan(1,N_sim);
status_log = nan(1,N_sim);

T_ref = torque_ramp(0.1,0.6,0,T_max,Ts,T_sim);  % motor torque reference
e_int = 0;  % for the integral state calculation
distance = 0;  % maximize final distance

for ii=1:N_sim
    % piecewise varying friction coefficient
    if distance >= 0; mu_x = 0.3; end
    if distance > 4; mu_x = 0.15; end
    if distance > 7.5; mu_x = 0.6; end

    if ~controller_type
        u_sim(:,ii) = T_ref(ii);
    else
        if controller_type == 1  % NMPC
            problem.xinit = x_ocp(:,ii);
            problem.x0 = repmat(ones(4,1), N, 1);  % solver initial guess
            warning('Fix solver initial guess')
            problem.hu = repmat(T_ref(ii), N, 1);
            problem.hl = zeros(N, 1);
            problem.all_parameters = repmat([T_ref(ii);mu_x], N, 1);
            [output, exitflag, info] = nmpc(problem);
            if ~isfield(output,'u0') && isfield(output,'zopt')
                output.u0 = output.zopt(1);
            end
        end
        
        if controller_type == 2  % KMPC
            problem.minus_x0 = -lifting_function(x_ocp(:,ii));  % (negative) initial state
            problem.T_ref = T_ref(ii)*ones(N,1);  % set for all stages
            f = [-2*w_u*T_ref(ii); w_x1; -kappa_ref*w_x1*R; 0;  % linear cost term
                 zeros(size(problem.minus_x0,1)-3,1)];  % extend for the lifted state
            problem.linear_cost = repmat(f,N,1);  % set for all stages
            [output, exitflag, info] = kmpc(problem);
            info.fevalstime = 0;  % add to struct for statistics
            output.u0 = mapminmax('reverse',output.u0,PU);  % unscale the input
        end
        
        if controller_type == 3  % PID
            v = x_sim(1,ii);
            w = x_sim(2,ii);
            kappa = 1 - v / (w*R);
            output.u0 = pid_controller(kappa_ref, kappa, Kp, Ki, Kd, T_ref(ii));
            info.it = 0;
            info.solvetime = 0;
            info.fevalstime = 0;
            exitflag = 1;
        end

        % check for failure and display some statistics
        status = exitflag;  % https://forces.embotech.com/Documentation/exitflags/index.html#tab-exitflag-values
        status_log(ii) = status;
        solve_time_log(ii) = info.solvetime + info.fevalstime;
        if (status==1 && ~mod(ii,print_step))
            num_iter = info.it;
            time_tot = info.solvetime + info.fevalstime;
            fprintf('\nstatus = %d, num_iter = %d, time_tot = %f [ms]\n',status, num_iter, time_tot*1e3);
        end
        if status~=1
            warning(['solver failed with status ',num2str(status)]);
        end

        % get solution for sim
        u_sim(:,ii) = output.u0;        
    end

	% set initial state of sim
	simulation.set('x', x_sim(:,ii));
	% set input in sim
	simulation.set('u', u_sim(:,ii));
    % set parameter
    simulation.set('p', mu_x);

	% simulate state
	simulation.solve();

	% get new sim state
	x_sim(:,ii+1) = simulation.get('xn');

    % get new ocp state
    v = x_sim(1,ii+1);
    w = x_sim(2,ii+1);
    s = w*R-v;
    e0 = 0.1;  % for slip modification
    kappa = (w*R-v)*w*R / ((w*R)^2 + e0);
    if T_ref(ii)>0  && ...              % integrate only after the torque ramp starts
       abs(u_sim(:,ii)-T_ref(ii))>1     % anti-windup
        e_int = e_int + (kappa_ref-kappa) * Ts;
    end
    x_ocp(:,ii+1) = [s; w; e_int];
    
    % update covered distance
    distance = distance + v*Ts;
end

disp([newline,newline,'Simulation done.'])
disp(['Covered distance: ',num2str(distance),' m'])
disp(['Average solve time: ',num2str(1e3*mean(solve_time_log)),' ms'])
disp(['Maximum solve time: ',num2str(1e3*max(solve_time_log)),' ms'])
disp(['Solver error rate: ',num2str(round(100*sum(status_log~=1)/N_sim)),'%'])
disp(['Real time feasible: ',num2str(round(100*sum(solve_time_log<Ts)/N_sim)),'%'])

%% plot the results
plot_results;

%% Literature
% Explicit Nonlinear Model Predictive Control for Electric Vehicle Traction Control
% https://ieeexplore.ieee.org/document/8390933
% Adaptive model predictive traction control for electric vehicles
% https://ieeexplore.ieee.org/document/8795687
% Model-based current limiting for traction control of an electric four-wheel drive race car
% https://ieeexplore.ieee.org/document/6862532
% An MPC/hybrid system approach to traction control
% https://ieeexplore.ieee.org/document/1624479
% Comparison of Centralized and Multi-Layer Architectures for Nonlinear Model Predictive Torque-Vectoring and Traction Control
% https://link.springer.com/article/10.1007/s12239-023-0090-x
% Optimal slip based traction control for electric vehicles using feedback linearization
% https://ieeexplore.ieee.org/document/7231734
% A Novel Design of Traction Control Based on a Piecewise-Linear Parameter-Varying Technique for Electric Vehicles With In-Wheel Motors
% https://ieeexplore.ieee.org/document/8424922
% Integrated Braking and Traction Torque Vectoring Control Based on Vehicle Yaw Rate for Stability Improvement of All-Wheel-Drive Electric Vehicles
% https://ieeexplore.ieee.org/document/10114899
% Integrated stability and traction control for electric vehicles using model predictive control
% https://sci-hub.st/https://doi.org/10.1016/j.conengprac.2016.06.005
% Advanced slip ratio for ensuring numerical stability of low-speed driving
% simulation. Part I: Longitudinal slip ratio
% https://sci-hub.st/https://doi.org/10.1177/0954407018759738
% Model-based Traction Control for Electric Vehicles
% https://link.springer.com/article/10.1365/s38314-014-0239-5
% C-GLISp
% https://ieeexplore.ieee.org/document/9667199
% https://github.com/bemporad/GLIS_MATLAB/tree/main
% A Survey of Traction Control and Antilock Braking Systems
% https://ieeexplore.ieee.org/document/6928502
% MPC tuning guidelines
% https://pubs.acs.org/doi/10.1021/acs.iecr.9b05931