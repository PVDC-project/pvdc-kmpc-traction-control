%% ------ Traction control using MPC and Carmaker ------
%% Environment setup
cd C:\Projects_Josip\tc\carmaker\src_cm4sl
clear;clc;close all;
addpath('../../models')             % prediction and simulation models
addpath('../../setup')              % controller and simulation setup
addpath('../../functions')          % utility functions
addpath('../../postprocessing/')    % plotting

% if not empty, simout will be saved in data/<save_filename>.mat
save_filename = 'nmpc_2';

%% Simulation setup
Ts = 2e-3;          % [s] sampling time
VEHICLE = vehicle_parameters();
R = VEHICLE.WHEEL_RADIUS;
Tmax = VEHICLE.MAX_MOTOR_TORQUE;
kappa_ref = 0.1;    % slip reference
v0 = 2;             % [km/h] initial car speed
w0 = v0/3.6/R;      % [rad/s] initial wheel speed

%% Controller setup
% 0 - off
% 1 - NMPC FORCES high-level interface
% 2 - KMPC FORCES low-level interface
% 3 - PID
% 4 - PID + random (data collection)
% 5 - KMPC YALMIP
% 6 - KMPC YALMIP adaptive
% 7 - open-loop random inputs (data collection)
% 8 - KMPC FORCES Y2F interface
controller_type = 1;
N = 5;                      % prediction horizon length
compile_for_simulink = 1;   % create the S-function block?
use_yalmip = controller_type == 5 || controller_type == 6;

mpc_setup = struct('N',N,'Ts',Ts,'R',R,'kappa_ref',kappa_ref',...
                   'compile_for_simulink',compile_for_simulink,...
                   'use_yalmip',use_yalmip);

% cost function weights
k = 0.001;
mpc_setup.w_p = k*750;    % slip tracking error weight
mpc_setup.w_i = k*500000; % integral state weight
mpc_setup.w_u = 1e-6;     % torque reduction weight

%w_p = 1e4;     % slip tracking error weight
%w_i = 1e3;      % integral state weight
%w_u = 1e-7;     % torque reduction weight
switch controller_type
    case 1
        nmpc_setup(mpc_setup);
    case 2
        kmpc_setup_low_level_interface(mpc_setup);
        load ../../models/dense_Fwu.mat F w_u   % for the problem formulation
        load ../../models/kmpc_data.mat PX PU   % for state and input scaling
    case {3,4}
        % discrete PI(D) parameters
        P = 750;
        I = 50000;
    case 5
%       controller = kmpc_setup_y2f(mpc_setup);  
%       loading/saving with optimizer doesn't always work
%       https://groups.google.com/g/yalmip/c/qK_uFt942Yo/m/IoD234uC22oJ
        save ../../data/kmpc_yalmip.mat mpc_setup   % for KMPC setup within the MATLAB System block
        load ../../models/kmpc_data.mat PX PU Alift % for state and input scaling
    case 6
        save ../../data/kmpc_yalmip.mat mpc_setup
    case 7
        Ts_rand = 50*Ts;    % how often to change the random input torque
        rng(42)             % fix the seed for reproducibility
        Tsim = 100;         % sim time max. 100 s
        Nsim = Tsim/Ts_rand;% number of random inputs
        T_max = 250;        % [Nm]
        T_rand = 0.5 * T_max * rand(1,Nsim);  % uniform distribution
        time_vector = 0:Ts_rand:(Nsim-1)*Ts_rand;
        sim_input = [time_vector' T_rand'];
    case 8
        kmpc_setup_y2f(mpc_setup);
        load ../../models/kmpc_data.mat PX PU Alift % for state and input scaling
    otherwise
        warning('controller type not recognized, running without control')
        controller_type = 0;
end

%% Run the simulation
disp('Starting simulation...')
sim('generic');
disp('Simulation done.')
%% Postprocessing
% TODO: plot the results
% get_mpc_inputs(sigsOut);

if controller_type == 4
    save_collected_data(sigsOut);
end

if ~isempty(save_filename)
    save_filename = ['../../data/',save_filename,'.mat'];
    save(save_filename,'sigsOut')
    disp(['Data saved in ',save_filename])
end
