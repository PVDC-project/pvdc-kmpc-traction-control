%% ------ Traction control using NMPC and Carmaker ------
%% Environment setup
clear;clc;close all;
addpath('../../models')             % prediction and simulation models
addpath('../../setup')              % controller and simulation setup
addpath('../../functions')          % utility functions
addpath('../../postprocessing/')    % plotting

%% Simulation setup
controller_type = 1;  % 0 - off, 1 - NMPC, 2 - KMPC, 3 - PID
Ts = 2e-3;  % [s] sampling time
kappa_ref = 0.1;  % slip reference
Tmax = 250;  % [Nm] maximum engine torque
R = 0.318;  % [m] wheel radius
v0 = 2;  % [km/h] initial car speed
w0 = v0/3.6/R;  % [rad/s] initial wheel speed

if controller_type == 1
    N = 4;  % prediction horizon length
    compile_for_simulink = 1;  % for use with Carmaker
    nmpc_setup(N,Ts,R,kappa_ref,compile_for_simulink);
elseif controller_type == 2
    error('KMPC not implemented yet')
%     kmpc_setup(N,Ts,R,kappa_ref,compile_for_simulink,w_u,w_x1,w_x3);
%     load setup/kmpc_data.mat PU  % for input scaling
elseif controller_type == 3
    % PI(D) parameters
    P = 800;
    I = 700;   
end

%% Run the simulation
sim('generic');

%% TODO: plot the results
% get_mpc_inputs(sigsOut);