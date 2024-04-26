%% Koopman operator identification for traction control
clear;clc;close all;
cd C:\Projects_Josip\tc\setup

addpath('../models')
addpath('../functions')
addpath('../data')

data_filename = '../models/kmpc_data.mat';

load ../data/collected_data.mat

%% Simulation setup
Ts = 2e-3;  % [s] sampling time
mu_x = 0.3; % [-] road-tire friction coefficient
Nsim = length(vx);
Ntraj = 1;

VEHICLE = vehicle_parameters();
R = VEHICLE.WHEEL_RADIUS;
% T_max = 0.5 * VEHICLE.MAX_MOTOR_TORQUE;
nx = 2;
nu = 1;

%% Visualize collected data
figure;
tsim = 0:Ts:(Nsim-1)*Ts;  % time vector for plotting
plot_traj_no = 1;  % number of trajectories to plot

Xtraj = [1/3.6*vx';w(:,1)'];
Xs = Xtraj;
Ys = circshift(Xtraj,1,2);
Us = tc_torque(:,1)';
for ii=1:plot_traj_no
    % vehicle speed
    subplot(4,1,1); hold on;
    plot(tsim, Xtraj(1,:) * 3.6);
    ylabel('$v_x$ [km/h]')
    title('Dataset (partial)')
    % wheel speed
    subplot(4,1,2); hold on;
    plot(tsim, Xtraj(2,:) * 30/pi);
    ylabel('$\omega$ [rpm]')
    % wheel torque
    subplot(4,1,3); hold on;
    stairs(tsim, Us(:,(ii-1)*Nsim+1:ii*Nsim));
    ylabel('$T$ [Nm]')
    % slip
    subplot(4,1,4); hold on;
    plot(tsim, (Xtraj(2,:)*R-Xtraj(1,:)) ./ (Xtraj(2,:)*R))
    ylabel('$\kappa$ [-]')
    xlabel('Time [s]')
end

%% Transform states from (v,w) to (s,w); s = w*R-v
% the integral state is handled manually
Xs = [Xs(2,:)*R-Xs(1,:); Xs(2,:)];
Ys = [Ys(2,:)*R-Ys(1,:); Ys(2,:)];

%% Scale data for better accuracy
XY = [Xs, Ys];                  % single scaling for the states
[XY,PX] = mapstd_custom(XY);    % map to zero mean and unit covariance
Xs = XY(:,1:Nsim*Ntraj);        % retrieve the matrices
Ys = XY(:,Nsim*Ntraj+1:end);

[Us,PU] = mapstd_custom(Us);    % scale the input

disp('Data scaled.')

%% Koopman operator approximation (EDMD)
% basis function selection and lifting
nrbf = 20;                  % number of basis functions
cent = rand(nx,nrbf)*2 - 1; % RBF centers, uniform in [-1,1]
rbf_type = 'thinplate';         % thinplate, gauss, invquad, invmultquad, polyharmonic

save(data_filename,'PX','PU','cent','rbf_type');

Xlift = lifting_function(Xs);
Ylift = lifting_function(Ys);
nz = size(Xlift,1);  % size of the lifted state

disp(['Data lifted.',newline])

disp('Starting regression for A,B,C...')
% EDMD
W = [Ylift ; Xs];
V = [Xlift ; Us];
VVt = V*V';
WVt = W*V';
ABC = WVt * pinv(VVt);
Alift = ABC(1:nz,1:nz);
Blift = ABC(1:nz,nz+1:end);
Clift = [eye(nx), zeros(nx,nz-nx)];  % set the output matrix manually
% Clift = ABC(nz+1:end,1:nz);
disp('Regression for A, B, C done.');

regression_residual = norm(Ylift - Alift*Xlift - Blift*Us,'fro') / norm(Ylift,'fro');
disp(['Regression residual: ', num2str(regression_residual),newline]);

%% Test predictor performance 
%{
N_test = 100;  % number of prediction steps for testing

% random control input
U = T_max * custom_prbs([nu N_test], 0.5);
% U = torque_ramp(0,0.8*N_test*Ts,0,T_max,Ts,N_test*Ts);  % ramp

% initial state
% TODO: random?
vx0s = [vx_low; (vx_low + vx_high)/2; vx_high];
figure
plot_ids = reshape(1:12, 3, 4).';  % subplot expects row-first indices

for jj = 1:3
    vx0 = vx0s(jj);
    kappa0 = 0.1;
    x0 = [vx0; vx0/(R*(1-kappa0))];

    % sanity check
    % warning('remove sanity check')
    % U = zeros(1,N_test);
    % x0 = [vx0;vx0/R];

    % logs
    X_true = nan(nx,N_test);  % true system state
    Z = nan(nz,N_test);       % Koopman internal state

    % initialize
    X_true(:,1) = x0;               % [v;w]
    s0 = [x0(2)*R-x0(1); x0(2)];    % [s;w]
    s0 = mapstd_custom('apply',s0,PX);     % scale the initial state
    Z(:,1) = lifting_function(s0);  % lift the initial state

    for ii = 1:N_test-1
        % true dynamics
        simulation.set('x', X_true(:,ii));
        simulation.set('u', U(:,ii));
        if (strcmp(simulation.opts_struct.method, 'irk'))
            simulation.set('xdot', zeros(nx,1));
        elseif (strcmp(simulation.opts_struct.method, 'irk_gnsf'))
            n_out = simulation.model_struct.dim_gnsf_nout;
            simulation.set('phi_guess', zeros(n_out,1));
        end
        simulation.solve();
        X_true(:,ii+1) = simulation.get('xn');

        % Koopman predictor
        u = mapstd_custom('apply',U(:,ii),PU);  % scale the input
        Z(:,ii+1) = Alift * Z(:,ii) + Blift * u;
    end

    % extract original state estimates from the Koopman predictor
    X_koop = Clift * Z;
    % reverse the scaling
    X_koop = mapstd_custom('reverse',X_koop,PX);

    % transform (v,w) into (s,w) for the true system
    X_true(1,:) = X_true(2,:)*R - X_true(1,:);  % s = w*R-v

    % plot the results
    t_test = 1e3 * (0:Ts:(N_test-1)*Ts);  % time vector
    lw_koop = 2;  % line width for plots
    
    subplot(4,3,plot_ids((jj-1)*4+1))
    plot(t_test,X_true(1,:)*3.6,'-b','linewidth', lw_koop); hold on
    plot(t_test,X_koop(1,:)*3.6, '--r','linewidth',lw_koop)
    legend('True','Koopman');
    ylabel('$s$ [km/h]')
    title(['Prediction accuracy with $v_{x0}$ = ',num2str(3.6*vx0),' km/h'])

    subplot(4,3,plot_ids((jj-1)*4+2))
    plot(t_test,X_true(2,:)*30/pi,'-b','linewidth', lw_koop); hold on
    plot(t_test,X_koop(2,:)*30/pi, '--r','linewidth',lw_koop)
    legend('True','Koopman');
    ylabel('$\omega$ [rpm]')

    subplot(4,3,plot_ids((jj-1)*4+3))
    plot(t_test,X_true(1,:)./R./X_true(2,:),'linewidth',lw_koop)
    hold on
    plot(t_test,X_koop(1,:)./R./X_koop(2,:),'linewidth',lw_koop)
    legend('True','Koopman')
    ylabel('$\kappa$ [-]')

    subplot(4,3,plot_ids((jj-1)*4+4))
    stairs(t_test,U,'linewidth',lw_koop);
    ylabel('T [Nm]')
    xlabel('time [ms]')
end

% maximize the created figure (requires MATLAB R2018a)
hFig = gcf;
hFig.WindowState = 'maximized';
%}

%% save the approximated system matrices and lifting data
save(data_filename,'Alift','Blift','Clift','PX','PU','cent','rbf_type');
disp(['Koopman model saved in ',data_filename])