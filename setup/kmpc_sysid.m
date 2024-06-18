%% Koopman operator identification for traction control using EDMD
%% Environment setup
clear all;clc;close all;
cd C:\Projects_Josip\tc\setup

addpath('../models')
addpath('../functions')
addpath('../data')

rng(42)  % fix  the seed for reproducibility

%% SYSID setup
% set the vehicle speed limits for data collection
vx_lows = 0;
vx_highs = 40;
% vx_lows = [0 10 20 30];
% vx_highs = [10 20 30 40];

if length(vx_lows) == 1
    adaptive = 0;    
    data_filename = '../models/kmpc_data.mat';
else
    adaptive = 1;
    data_filename = '../models/kmpc_data_adaptive.mat';
end

Ntraj = 1000;     % total trajectories (per speed interval)
Nsim = 250;       % samples per trajectory

% choose the amplitude and DC offset of random excitation
VEHICLE = vehicle_parameters();
U_amplitude = 0.5 * VEHICLE.MAX_MOTOR_TORQUE;
U_offset = 50;  % [Nm]

% basis function selection
nrbf = 30;                  % number of basis functions
rbf_type = 'polynomial';    % polynomial, thinplate, gauss, invquad, invmultquad, polyharmonic
if strcmp(rbf_type,'polynomial')
    order = 4;  % maximum order for the polynomial
    nx = 2;     % original state size
    [~,nrbf] = create_poly_basis(nx,order,true);  % last input creates an .m file
end

vx = vx_highs;                  % speed limits for model selection
kmpc_datas = cell(size(vx));    % models and scaling data

%% Simulation setup
Ts = 2e-3;  % [s] sampling time
mu_x = 0.3; % [-] road-tire friction coefficient
R = VEHICLE.WHEEL_RADIUS;

% use acados for creating the dataset (much faster than ode45)
simulation = sim_setup_acados(Ts,mu_x);
nx = 2;
nu = 1;
print_step = Ntraj/10;  % print ten times during data collection

%% Data collection (v,w)
for k = 1:length(vx_lows)
    vx_low = vx_lows(k);
    vx_high = vx_highs(k);
    disp(['Collecting data for speeds: ',num2str(vx_low),' to ',num2str(vx_high),' km/h.'])

    % random control input forcing
    U = U_amplitude * rand(1,Ntraj*Nsim);   % uniform distribution
    U = U + U_offset;                       % add the offset
    U = repelem(U,1);                       % hold the input for a few samples (more realistic)
    U = U(1:Ntraj*Nsim);                    % trim to match the desired length
    % U = T_max*custom_prbs([nu Ntraj*Nsim], 0.5);  % PRBS
    % U = torque_ramp(0.1,0.6,0,T_max,h_sim,h_sim*Ntraj*Nsim);  % ramp

    % random initial condition from the given interval
    vx_low = vx_low/3.6;    % [m/s]
    vx_high = vx_high/3.6;  % [m/s]
    kappa_low = 0;          % [-] 
    kappa_high = 0.2;       % [-]
    vx_init = vx_low + (vx_high-vx_low) * rand(1,Ntraj);                % uniform distribution
    kappa_init = kappa_low + (kappa_high-kappa_low) * rand(1,Ntraj);    % uniform distribution
    w_init = vx_init ./ (R*(1-kappa_init));                             % defined by slip and speed
    Xinit = [vx_init; w_init];                                          % simulation model uses [v;w]

    % logs
    Xs = nan(nx,Ntraj*Nsim);
    Ys = nan(nx,Ntraj*Nsim);
    Us = nan(nu,Ntraj*Nsim);

    % simulate Ntraj trajectories of Nsim samples
    for ii=1:Ntraj
        if ~mod(ii,print_step)
            disp(['Collected ',num2str(10*(ii/print_step)),'%'])
        end
        Xcurrent = Xinit(:,ii);
        for jj=1:Nsim
            current_index = (ii-1)*Nsim + jj;  % index in the merged dataset

            % set initial state and input
            simulation.set('x', Xcurrent);
            simulation.set('u', U(:,current_index));

            % initialize implicit integrator
            if (strcmp(simulation.opts_struct.method, 'irk'))
                simulation.set('xdot', zeros(nx,1));
            elseif (strcmp(simulation.opts_struct.method, 'irk_gnsf'))
                n_out = simulation.model_struct.dim_gnsf_nout;
                simulation.set('phi_guess', zeros(n_out,1));
            end

            % simulate one step and get the next state
            simulation.solve();
            Xnext = simulation.get('xn');

            % log states and inputs
            Xs(:,current_index) = Xcurrent;
            Ys(:,current_index) = Xnext;
            Us(:,current_index) = U(:,current_index);

            % advance to the next step
            Xcurrent = Xnext;
        end
    end
    disp(['Data collection done.',newline])

    %% Visualize collected data
    figure;
    tsim = 0:Ts:(Nsim-1)*Ts;  % time vector for plotting
    plot_traj_no = 10;  % number of trajectories to plot

    for ii=1:plot_traj_no
        Xtraj = Xs(:,(ii-1)*Nsim+1:ii*Nsim);
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
        ylabel('$\sigma_x$ [-]')
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
    cent = rand(nx,nrbf)*2 - 1; % RBF centers, uniform in [-1,1]
    kmpc_data = struct('PX',PX,'PU',PU,'cent',cent,'rbf_type',rbf_type);

    Xlift = lifting_function(Xs,kmpc_data);
    Ylift = lifting_function(Ys,kmpc_data);
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
    
    kmpc_data.Alift = Alift;
    kmpc_data.Blift = Blift;
    kmpc_data.Clift = Clift;
    
    kmpc_datas{k} = kmpc_data;
    
    %% Test predictor performance 
    N_test = 100;  % number of prediction steps for testing

    % random control input
    U = U_amplitude * custom_prbs([nu N_test], 0.5);
    U = U + U_offset;
    % U = torque_ramp(0,0.8*N_test*Ts,0,T_max,Ts,N_test*Ts);  % ramp

    % initial state
    % TODO: random?
    vx0s = [vx_low; (vx_low + vx_high)/2; vx_high];
    figure
    plot_ids = reshape(1:12, 3, 4).';  % subplot expects row-first indices

    MSE = nan(1,3);
    for jj = 1:3
        vx0 = vx0s(jj);
        kappa0 = 0.1;
        x0 = [vx0; vx0/(R*(1-kappa0))];

        % logs
        X_true = nan(nx,N_test);  % true system state
        Z = nan(nz,N_test);       % Koopman internal state

        % initialize
        X_true(:,1) = x0;                           % [v;w]
        s0 = [x0(2)*R-x0(1); x0(2)];                % [s;w]
        s0 = mapstd_custom('apply',s0,PX);          % scale the initial state
        Z(:,1) = lifting_function(s0,kmpc_data);    % lift the initial state

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
        
        % compare the true slip ratio and predictions
        MSE(jj) = mse(X_true(1,:)./R./X_true(2,:),X_koop(1,:)./R./X_koop(2,:));
        
        % save the approximation data for later plotting
        if (vx0 == (vx_low + vx_high)/2)
            kappa_true = X_true(1,:)./R./X_true(2,:);
            kappa_koop = X_koop(1,:)./R./X_koop(2,:);
            save ../data/sysid.mat kappa_true kappa_koop vx0
            disp(['Saved the approximation data for vx0 = ',num2str(vx0*3.6),' km/h in data/sysid.mat'])
        end
    end
    
    disp(['MSE: ', num2str(sum(MSE)),newline]);

    % maximize the created figure (requires MATLAB R2018a)
    hFig = gcf;
    hFig.WindowState = 'maximized';
    
    % wait a bit before moving on
    pause(2)
end

%% save the approximated system matrices and lifting data
if ~adaptive
    save(data_filename,'Alift','Blift','Clift','PX','PU','cent','rbf_type');
else
    save(data_filename,'vx','kmpc_datas')
end
disp(['Koopman model saved in ',data_filename])

%% remove unnecessary files
delete *.json
disp('Deleted acados .json file(s).')

%% Literature
% Linear predictors for nonlinear dynamical systems: Koopman operator meets model predictive control
% https://arxiv.org/pdf/1611.03537.pdf
% https://github.com/MilanKorda/KoopmanMPC/
% Generalizing Koopman Theory to Allow for Inputs and Control
% https://epubs.siam.org/doi/pdf/10.1137/16M1062296
% Modern Koopman Theory for Dynamical Systems
% https://arxiv.org/abs/2102.12086
% https://github.com/eurika-kaiser/koopman-mpc-comparison
% Experimental physical parameter estimation of a thyristor driven DC-motor using the HMF-method
% https://www.sciencedirect.com/science/article/pii/S0967066198000367
% Optimal Construction of Koopman Eigenfunctions for Prediction and Control
% https://ieeexplore.ieee.org/abstract/document/9022864
% Model Predictive Control of a Vehicle using Koopman Operator
% https://www.sciencedirect.com/science/article/pii/S2405896320331748
% Predictive Direct Yaw Moment Control Based on the Koopman Operator
% https://ieeexplore.ieee.org/abstract/document/10122214