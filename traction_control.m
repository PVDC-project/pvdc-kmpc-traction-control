%% Traction control using NMPC and KMPC
clear all;close all;clc;

%% Environment setup
addpath('./models')             % prediction and simulation models
addpath('./setup')              % controller and simulation setup
addpath('./functions')          % utility functions
addpath('./postprocessing')     % plotting

%% Simulation setup
Ts = 2e-3;  % [s] sampling time

% inputs: [v;w], T, mu_x
% output: [v;w] after Ts
simulation = @(x,u,p) simulation_model(x,u,p,Ts);
sim_nx = 2;

VEHICLE = vehicle_parameters();
R = VEHICLE.WHEEL_RADIUS;
T_max = VEHICLE.MAX_MOTOR_TORQUE;

kappa_ref = 0.1;        % slip reference
mu_x = 0.3;             % [-] road-tire friction coefficient
change_friction = 1;    % test the controller on different surfaces

%% Controller setup
% 0 - off
% 1 - NMPC FORCES high-level interface
% 2 - KMPC FORCES low-level interface
% 3 - PID
% 4 - KMPC YALMIP
% 5 - KMPC FORCES Y2F interface
% 6 - adaptive KMPC FORCES Y2F (prediction model changes with vehicle speed)
% 7 - adaptive KMPC YALMIP (quadprog, DAQP or OSQP)
% 8 - KMPC acados
controller_type = 8;
N = 5;                      % prediction horizon length
compile_for_simulink = 0;   % create the S-function block?
use_yalmip = controller_type == 4 || controller_type == 7;
is_adaptive = controller_type == 6 || 7;
mpc_setup = struct('N',N,'Ts',Ts,'R',R,'kappa_ref',kappa_ref',...
                   'compile_for_simulink',compile_for_simulink,...
                   'use_yalmip',use_yalmip);

switch controller_type
    case 1
        nmpc_setup(mpc_setup);
        warning('Fix solver initial guess'); input('Continue?')
    case 2
        kmpc_setup_low_level_interface(mpc_setup);
        load models/dense_Fwu.mat F w_u  % for low-level interface formulation
        load models/kmpc_data.mat PX PU  % for state and input scaling
    case 3
        Kp = 750;  % proportional gain
        Ki = 100;  % integral gain
        Kd = 0;    % derivative gain
    case 4
        controller = kmpc_setup_y2f(mpc_setup);
        load models/kmpc_data.mat PX PU  % for state and input scaling
    case 5
        kmpc_setup_y2f(mpc_setup);
        load models/kmpc_data.mat PX PU  % for state and input scaling
    case {6,7}
        controllers = kmpc_setup_y2f_adaptive(mpc_setup);
        load models/kmpc_data_adaptive.mat vx kmpc_datas
        speed_limits = vx;  % to avoid name clashes later
    case 8
        controller = kmpc_setup_acados(mpc_setup);
        load models/kmpc_data.mat PX PU  % for state and input scaling
    otherwise
        disp('controller type not recognized, running without control')
        controller_type = 0;
end

% state-space dimensions
ocp_nx = 3;
ocp_nu = 1;

%% Closed loop simulation
v0 = 2;         % [km/h] initial car speed
v0 = v0/3.6;    % convert to m/s
w0 = v0/R;      % [rad/s] initial wheel speed (no slip)

sim_x0 = [v0; w0];              % initial state
ocp_x0 = [w0*R-v0; w0; 0];      % wheel slip velocity, wheel speed, integral state
T_sim = 3;                      % [s] simulation time
N_sim = T_sim/Ts;               % number of simulation steps
print_step = ceil(N_sim/10);    % print stats 10 times during the simulation

% state and input logging
x_sim = nan(sim_nx, N_sim+1);  % simulated state log
x_sim(:,1) = sim_x0;
x_ocp = nan(ocp_nx, N_sim+1);  % controller state log
x_ocp(:,1) = ocp_x0;
u_sim = nan(ocp_nu, N_sim);    % control input log

% MPC performance logging
solve_time_log = nan(1,N_sim);
status_log = nan(1,N_sim);
solver_time_log = nan(1,N_sim);

% motor torque reference
T_ref = torque_ramp(0.1,0.6,0,T_max,Ts,T_sim);
% T_ref = T_max/2 * custom_prbs([1 T_sim/Ts], 0.5);
% T_ref = repelem(T_ref,50);  % hold the input for a few samples (more realistic)
% T_ref = T_ref(1:T_sim/Ts);  % trim to match the desired length

e_int = 0;      % for the integral state calculation
distance = 0;   % compare final distances
mu_x_log = nan(1,N_sim);

vid = zeros(1,N_sim);
disp([newline,'Starting the closed-loop simulation...'])
for ii=1:N_sim
    if change_friction  % varying friction coefficient
        if distance >= 0; mu_x = 0.3; end
        if distance > 4; mu_x = 0.15; end
        if distance > 7.5; mu_x = 0.6; end
    end
    mu_x_log(ii) = mu_x;

    % get the control input
    switch controller_type
        case 0  % no control
            u_sim(:,ii) = T_ref(ii);
        case 1  % NMPC FORCES high-level interface
            tic
            problem.xinit = x_ocp(:,ii);
            problem.x0 = repmat(ones(4,1), N, 1);  % TODO: fix solver initial guess
            problem.hu = repmat(T_ref(ii), N, 1);
            problem.hl = zeros(N, 1);
            problem.all_parameters = repmat([T_ref(ii);0.3], N, 1);  % TODO: use mu_x info?
            [output, exitflag, info] = nmpc(problem);
            info.totaltime = toc;
            if exitflag~=1
                warning(['solver failed with status ',num2str(exitflag)]);
                output.u0 = T_ref(ii);  % don't modify the torque if failed
            else
                output.u0 = output.zopt(1);
            end
        case 2  % KMPC FORCES low-level interface
            tic
            state_to_lift = mapstd_custom('apply',x_ocp(1:2,ii),PX);  % don't scale e_int
            z0 = [lifting_function(state_to_lift); x_ocp(3,ii); kappa_ref];
            problem.F_times_z0 = F*z0;
            T = mapstd_custom('apply',T_ref(ii),PU);
            problem.T_ref = T*ones(N,1);
            f = [-2*w_u*T; 0; 0; 0];
            problem.linear_cost = repmat(f,N,1);

            [output, exitflag, info] = kmpc(problem);
            info.totaltime = toc;
            info.fevalstime = 0;
            if exitflag~=1
                warning(['solver failed with status ',num2str(status)]);
                output.u0 = T_ref(ii);  % don't modify the torque if failed
            else
                output.u0 = mapstd_custom('reverse',output.u0,PU);  % unscale the input
            end
        case 3  % PID
            tic
            v = x_sim(1,ii);
            w = x_sim(2,ii);
            e0 = 0.1;  % for slip modification
            kappa = (w*R-v)*w*R / ((w*R)^2 + e0);
            output.u0 = pid_controller(kappa_ref, kappa, Kp, Ki, Kd, T_ref(ii));
            info.totaltime = toc;
            info.it = 0;
            info.solvetime = 0;
            info.fevalstime = 0;
            exitflag = 1;
        case 4  % KMPC YALMIP
            tic
            state_to_lift = mapstd_custom('apply',x_ocp(1:2,ii),PX);  % don't scale e_int
            problem{1} = [lifting_function(state_to_lift); x_ocp(3,ii); kappa_ref];
            problem{2} = mapstd_custom('apply',T_ref(ii),PU);

            [yalmip_output, exitflag, a, b, c, info] = controller(problem);
            info.totaltime = toc;
            if isfield(info.solveroutput,'output')
                info.it = info.solveroutput.output.iterations; % quadprog
            elseif isfield(info.solveroutput,'iter')
                info.it = info.solveroutput.iter;              % DAQP
            elseif isfield(info.solveroutput,'info')
                info.it = info.solveroutput.info.iter;         % OSQP
            else
                info.it = 1;                                   % other solvers
            end
            info.solvetime = info.solvertime;   % FORCES stats
            info.fevalstime = 0;                % FORCES stats
            solver_time_log(ii) = info.solvertime;

            if exitflag  % with YALMIP 0 means OK
                warning(['solver failed with flag ',num2str(exitflag),' - ', a{1}])
                output.u0 = T_ref(ii);
            else
                output.u0 = mapstd_custom('reverse',yalmip_output{1}(:,1),PU);
            end
        case 5  % KMPC FORCES Y2F interface
            tic
            state_to_lift = mapstd_custom('apply',x_ocp(1:2,ii),PX);  % don't scale e_int
            problem{1} = [lifting_function(state_to_lift); x_ocp(3,ii); kappa_ref];
            problem{2} = mapstd_custom('apply',T_ref(ii),PU);

            [forces_output, exitflag, info] = kmpc(problem);
            info.totaltime = toc;
            info.fevalstime = 0;  % add to struct for statistics
            if exitflag~=1
                warning(['solver failed with status ',num2str(status)]);
                output.u0 = T_ref(ii);  % don't modify the torque if failed
            else
                output.u0 = mapstd_custom('reverse',forces_output{1}(:,1),PU);  % unscale the input
            end
        case 6  % adaptive KMPC Y2F
            error('not implemented correctly')
            tic
            state_to_lift = mapstd_custom('apply',x_ocp(1:2,ii),PX);  % don't scale e_int
            problem{1} = [lifting_function(state_to_lift); x_ocp(3,ii); kappa_ref];
            problem{2} = mapstd_custom('apply',T_ref(ii),PU);
            if v < 10; [forces_output, exitflag, info] = kmpc10(problem); end
            if v >= 10 && v < 20; [forces_output, exitflag, info] = kmpc20(problem); end
            if v >= 20 && v < 30; [forces_output, exitflag, info] = kmpc30(problem); end
            if v >= 30; [forces_output, exitflag, info] = kmpc40(problem); end
            info.totaltime = toc;
            info.fevalstime = 0;  % add to struct for statistics
            if exitflag~=1
                warning(['solver failed with status ',num2str(status)]);
                output.u0 = T_ref(ii);  % don't modify the torque if failed
            else
                output.u0 = mapstd_custom('reverse',forces_output{1}(:,1),PU);  % unscale the input
            end
        case 7  % adaptive KMPC YALMIP
            v = x_sim(1,ii)*3.6;                    % [km/h]
            ctrl_id = find(speed_limits > v, 1);    % find which controller to use
            if isempty(ctrl_id); ctrl_id = length(speed_limits); end
            tic
            state_to_lift = mapstd_custom('apply',x_ocp(1:2,ii),kmpc_datas{ctrl_id}.PX);  % don't scale e_int
            problem{1} = [lifting_function(state_to_lift,kmpc_datas{ctrl_id}); x_ocp(3,ii); kappa_ref];
            problem{2} = mapstd_custom('apply',T_ref(ii),kmpc_datas{ctrl_id}.PU);
            controller = controllers{ctrl_id};
            [yalmip_output, exitflag, a, b, c, info] = controller(problem);
            info.totaltime = toc;
            if isfield(info.solveroutput,'output')
                info.it = info.solveroutput.output.iterations; % quadprog
            elseif isfield(info.solveroutput,'iter')
                info.it = info.solveroutput.iter;              % DAQP
            elseif isfield(info.solveroutput,'info')
                info.it = info.solveroutput.info.iter;         % OSQP
            else
                info.it = 1;                                   % other solvers
            end
            info.solvetime = info.solvertime;   % FORCES stats
            info.fevalstime = 0;                % FORCES stats

            if exitflag  % with YALMIP 0 means OK
                warning(['solver failed with flag ',num2str(exitflag),' - ', a{1}])
                output.u0 = T_ref(ii);
            else
                output.u0 = mapstd_custom('reverse',yalmip_output{1}(:,1),kmpc_datas{ctrl_id}.PU);
            end
        case 8  % KMPC acados
            tic
            state_to_lift = mapstd_custom('apply',x_ocp(1:2,ii),PX);  % don't scale e_int
            z0 = [lifting_function(state_to_lift); x_ocp(3,ii); kappa_ref];
            T = mapstd_custom('apply',T_ref(ii),PU);
            controller.set('p',[z0;T]);
            controller.solve();
            info.totaltime = toc;
            
            exitflag = controller.get('status');
            if exitflag ~= 0
                warning(['acados solver failed with status ',num2str(exitflag)]);
                output.u0 = T_ref(ii);
            else
                Uopt = controller.get('u', 0);
                output.u0 = mapstd_custom('reverse',Uopt(1),PU);  % unscale the input
            end
            % TODO: fix
            info.it = 0;
            info.solvetime = 0;
            info.fevalstime = 0;
    end

    % check for failure and display some statistics
    status_log(ii) = exitflag;
    solve_time_log(ii) = info.totaltime;  % info.solvetime + info.fevalstime;
    if ~mod(ii,print_step)
        disp(['Simulated ',num2str(10*(ii/print_step)),'%'])
        num_iter = info.it;
        time_tot = info.totaltime;  % info.solvetime + info.fevalstime;
        fprintf('status = %d, num_iter = %d, time_tot = %f [ms]\n',exitflag, num_iter, time_tot*1e3);
    end

    % simulate one closed-loop timestep
    u_sim(:,ii) = output.u0;
    x_sim(:,ii+1) = simulation(x_sim(:,ii),u_sim(:,ii),mu_x);

    % get new controller state
    v = x_sim(1,ii+1);
    w = x_sim(2,ii+1);
    s = w*R-v;
    e0 = 0.1;  % for slip modification
    kappa = s*w*R / ((w*R)^2 + e0);
    if T_ref(ii)>0  && ...              % integrate only after the torque ramp starts
       abs(u_sim(:,ii)-T_ref(ii))>1     % anti-windup
        e_int = e_int + (kappa_ref-kappa) * Ts;
    end
    x_ocp(:,ii+1) = [s; w; e_int];
    
    % update covered distance
    distance = distance + v*Ts;
end

disp([newline,'Simulation done.'])

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