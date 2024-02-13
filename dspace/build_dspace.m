%% ------ Traction control using NMPC and dSPACE ------
%% Environment setup
clear all; clc; close all;
addpath('../models')             % prediction and simulation models
addpath('../setup')              % controller and simulation setup
addpath('../functions')          % utility functions

%% Experiment setup
controller_type = 1;  % 0 - off, 1 - NMPC, 2 - KMPC, 3 - PID
Ts = 2e-3;  % [s] sampling time
kappa_ref = 0.1;  % slip reference
Tmax = 250;  % [Nm] maximum engine torque
w0 = 4.5;  % [rad/s] initial wheel speed

if controller_type == 1
    R = 0.318;  % [m] wheel radius
    N = 4;  % prediction horizon length
    compile_for_simulink = 1;  % for use with Carmaker
    compile_for_dspace = 1;  % generate code for embedded hardware
    nmpc_setup(N,Ts,R,kappa_ref,compile_for_simulink,compile_for_dspace);
elseif controller_type == 2
    error('KMPC not implemented yet')
%     kmpc_setup(N,Ts,R,kappa_ref,compile_for_simulink,w_u,w_x1,w_x3);
%     load setup/kmpc_data.mat PU  % for input scaling
elseif controller_type == 3
    % PI(D) parameters
    P = 2500;
    I = 5000;   
end

%% Build the code
build_dir = './build';
if ~exist(build_dir,'dir')
    mkdir(build_dir)
end
addpath(build_dir);
cd(build_dir);

[errorFlag, errorMsg] = rti_build2('nmpc_traction_control_dspace', 'Command', 'CM');
cd('../');