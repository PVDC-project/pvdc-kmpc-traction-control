%% ------ Traction control using MPC and Carmaker ------
%% Environment setup
clear;clc;close all;
addpath('../../models')             % prediction and simulation models
addpath('../../setup')              % controller and simulation setup
addpath('../../functions')          % utility functions
addpath('../../postprocessing/')    % plotting

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
% 1 - NMPC
% 2 - KMPC
% 3 - PID
% 4 - PID + random (data collection)
% 5 - KMPC YALMIP
controller_type = 5;
N = 5;                      % prediction horizon length
compile_for_simulink = 1;   % create the S-function block?
use_yalmip = controller_type == 5;

mpc_setup = struct('N',N,'Ts',Ts,'R',R,'kappa_ref',kappa_ref',...
                   'compile_for_simulink',compile_for_simulink,...
                   'use_yalmip',use_yalmip);

% cost function weights
mpc_setup.w_p = 1e4;      % slip tracking error weight
mpc_setup.w_i = 1e3;      % integral state weight
mpc_setup.w_u = 1e-7;     % torque reduction weight

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
        load ../../models/kmpc_data.mat PX PU       % for state and input scaling
    otherwise
        disp('controller type not recognized, running without control')
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
