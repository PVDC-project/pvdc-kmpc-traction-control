%% ------ Traction control using NMPC and dSPACE ------
%% Environment setup
clear all; clc; close all;
addpath('../models')             % prediction and simulation models
addpath('../setup')              % controller and simulation setup
addpath('../functions')          % utility functions

%% Experiment setup
controller_type = 2;  % 0 - off, 1 - NMPC, 2 - KMPC, 3 - PID
Ts = 2e-3;  % [s] sampling time
kappa_ref = 0.1;  % slip reference
Tmax = 250;  % [Nm] maximum engine torque
R = 0.318;   % [m] wheel radius
% w0 = 4.5;  % [rad/s] initial wheel speed

N = 5;                      % prediction horizon length
compile_for_simulink = 1;   % create the S-function block?
use_yalmip = controller_type == 5 || controller_type == 6;
compile_for_dspace = 1;     % generate code for embedded hardware

mpc_setup = struct('N',N,'Ts',Ts,'R',R,'kappa_ref',kappa_ref',...
                   'compile_for_simulink',compile_for_simulink,...
                   'use_yalmip',use_yalmip);

if controller_type == 1
    load ../data/mpc_inputs_nmpc.mat
    x0 = squeeze(x0)';
    nmpc_setup(mpc_setup,compile_for_simulink);
elseif controller_type == 2
    load ../data/mpc_inputs_kmpc.mat
    mpc_setup.w_p = 0.75;   % slip tracking error weight
    mpc_setup.w_i = 500;    % integral state weight
    mpc_setup.w_u = 1e-6;   % torque reduction weight
    mpc_setup.compile_for_dspace = compile_for_dspace;
    kmpc_setup_y2f(mpc_setup);
end

%% Build the code
build_dir = './build';
if ~exist(build_dir,'dir')
    mkdir(build_dir)
end
addpath(build_dir);
cd(build_dir);

if controller_type == 1
    [errorFlag, errorMsg] = rti_build2('nmpc_traction_control_dspace', 'Command', 'CM');
elseif controller_type == 2
    [errorFlag, errorMsg] = rti_build2('kmpc_traction_control_dspace', 'Command', 'CM');
end
cd('../');