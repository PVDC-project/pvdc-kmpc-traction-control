clear;clc;close all;

controller_type = 1;  % 0 - off, 1 - NMPC, 2 - KMPC, 3 - PID
Ts = 2e-3;  % [s] sampling time
R = 0.318;  % [m] wheel radius
x0 = [1.4; 4.4];  % [v0;w0]
N = 4;  % prediction horizon length

sim('traction_control_simulink');

% TODO: plot the results
