%% Koopman MPC controller setup
function controllers = kmpc_setup_y2f_adaptive(mpc_setup)
N = mpc_setup.N; Ts = mpc_setup.Ts;

load kmpc_data_adaptive.mat vx kmpc_datas
speeds = vx;
controllers = cell(size(speeds));
for speed_idx = 1:length(speeds)
    speed = speeds(speed_idx);
    kmpc_data = kmpc_datas{speed_idx};
    Alift = kmpc_data.Alift;
    Blift = kmpc_data.Blift;
    PU = kmpc_data.PU;

    % add the integral state (e_int_dot = kappa_ref - kappa) using Euler method
    % also add kappa_ref to states (for easier problem formulation)
    nz = size(Alift,1);
    e_int_row = zeros(1,nz+2);
    e_int_row(3) = -Ts;         % slip is the third state
    e_int_row(end-1) = 1;       % integral state is second-to-last
    e_int_row(end) = Ts;        % slip reference is the last state
    A = [Alift, zeros(nz,2);
        e_int_row;
        zeros(1,nz+1), 1];      % slip reference has no dynamics

    B = [Blift;
        zeros(2,size(Blift,2))];

    % kappa, e_int and kappa_ref are needed to form the cost
    C = zeros(3,nz+2);
    C(1,3) = 1;     % slip
    C(2,end-1) = 1; % integral state
    C(3,end) = 1;   % slip reference

    ny = size(C,1);         % number of outputs
    [nz, nu] = size(B);     % number of (extended) states and inputs

    % get dense form matrices (Nc=Np=N)
    [F,Phi] = dense_prediction_matrices(A,B,C,N);

    % cost matrices
    w_p = 1e-1;     % slip tracking error weight
    w_i = 1e3;      % integral state weight
    w_u = 1e-5;     % torque reduction weight

    if isfield(mpc_setup,'w_p')  % cost weights set outside
        w_p = mpc_setup.w_p;
        w_i = mpc_setup.w_i;
        w_u = mpc_setup.w_u;
    end

    % w_p = 1e-2; w_i = 1e-2; w_u = 1e2;  % test input reference tracking

    % define problem using YALMIP
    u = sdpvar(nu,N,'full');
    z0 = sdpvar(nz,1,'full');
    T_ref = sdpvar;         % for online torque limiting (should be scaled)
    Y = F*z0 + Phi*u(:);    % dense-form prediction

    umin = mapstd_custom('apply',0,PU);  % scale minimum motor torque

    constraints = [];
    objective = 0;

    for k = 1:N
        y = Y((k-1)*ny+1:k*ny);                             % output from k-th prediction step

        objective = objective + w_p * (y(3)-y(1))^2;        % slip tracking cost
        objective = objective + w_i * y(2)^2;               % integral state cost
        objective = objective + w_u * (T_ref-u(k))^2;       % input torque reduction cost
        constraints = [constraints, umin <= u(k) <= T_ref]; % input constraint
    end

    % solver settings
    algorithm = 'PDIP';  % https://forces.embotech.com/Documentation/solver_options/index.html#solve-methods
    codeoptions = ForcesGetDefaultOptions(['kmpc',num2str(speed)],algorithm,'double');
    codeoptions.printlevel = 0; % summary line after each solve
    codeoptions.overwrite = 1;  % always overwrite the solver

    codeoptions.init = 1;  % centered start

    % Numerical instability debugging
    %https://forces.embotech.com/Documentation/high_level_interface/index.html#call-callback-ref
    % codeoptions.MEXinterface.dynamics = 1;
    % codeoptions.MEXinterface.inequalities = 1;
    % codeoptions.MEXinterface.objective = 1;

    % codeoptions.accuracy.ineq = 1e-8;  % infinity norm of residual for inequalities
    % codeoptions.accuracy.eq = 1e-8;    % infinity norm of residual for equalities
    % codeoptions.accuracy.mu = 1e-8;    % absolute duality gap

    % Simulink block options
    if ~mpc_setup.compile_for_simulink
        codeoptions.BuildSimulinkBlock = 0;
    end
    codeoptions.showinfo = 1;  % https://forces.embotech.com/Documentation/solver_options/index.html#solver-info-in-simulink-block

    params = {z0,T_ref};
    outputs = {u,Y};
    parameter_names = {'z0','T_ref'};
    output_names = {'uopt','yopt'};

    % generate code
    solver_dir = './codegen';
    if ~exist(solver_dir,'dir')
        mkdir(solver_dir)
    end
    addpath(solver_dir);
    cd(solver_dir);

    if ~mpc_setup.use_yalmip  % skip codegen if not needed
        optimizerFORCES(constraints, objective, codeoptions, params, outputs, parameter_names, output_names);
    end
    disp(['Generated KMPC',num2str(speed),'.',newline,'------------------------------',newline])

    cd('../');
    
    % YALMIP for comparison
    options = sdpsettings('solver','quadprog','savesolverinput',1,'savesolveroutput',1);
    controllers{speed_idx} = optimizer(constraints, objective, options, params, outputs);
end

disp([num2str(length(speeds)),' controllers created.'])
    %% test call
%     disp([newline,'Testing the generated solver...'])
%     state_to_lift = mapstd_custom('apply',[0;10],PX);
%     problem{1} = [lifting_function(state_to_lift); 0; 0.1];  % lift, e_int, kappa_ref
%     problem{2} = mapstd_custom('apply',250,PU);  % T_ref (scaled)
% 
%     [forces_sol, forces_flag, forces_info] = kmpc10(problem);
%     disp(['FORCES test solution: ',num2str(forces_sol{1})])
%     assert(forces_flag == 1, 'Test call of FORCESPRO solver failed');
%     [forces_sol, forces_flag, forces_info] = kmpc20(problem);
%     disp(['FORCES test solution: ',num2str(forces_sol{1})])
%     assert(forces_flag == 1, 'Test call of FORCESPRO solver failed');
%     [forces_sol, forces_flag, forces_info] = kmpc30(problem);
%     disp(['FORCES test solution: ',num2str(forces_sol{1})])
%     assert(forces_flag == 1, 'Test call of FORCESPRO solver failed');
%     [forces_sol, forces_flag, forces_info] = kmpc40(problem);
%     disp(['FORCES test solution: ',num2str(forces_sol{1})])
%     assert(forces_flag == 1, 'Test call of FORCESPRO solver failed');
% 
%     [yalmip_sol, yalmip_flag, ~, ~, ~, yalmip_info] = controller(problem);
%     disp(['YALMIP test solution: ',num2str(yalmip_sol{1})]);
%     assert(yalmip_flag == 0, 'Test call with YALMIP failed')
%     disp(['Test ok.',newline])
% 
%     assert(max(abs(yalmip_sol{1}-forces_sol{1})) < 1e-6, 'FORCES and YALMIP solutions differ too much.')
end