%% Koopman MPC controller setup
function [] = kmpc_setup_low_level_interface(mpc_setup)
N = mpc_setup.N; Ts = mpc_setup.Ts;
compile_for_simulink = mpc_setup.compile_for_simulink;

% system matrices
load kmpc_data.mat Alift Blift PX PU;

% add the integral state (e_int_dot = kappa_ref - kappa) using Euler method
% also add kappa_ref to states (for easier problem formulation)
nz = size(Alift,1);
e_int_row = zeros(1,nz+2);
e_int_row(3) = -Ts;         % slip is the third state
e_int_row(end-1) = 1;       % integral state is second-to-last
e_int_row(end) = Ts;        % slip reference is the last state
A = [Alift, zeros(nz,2);
    e_int_row;
    zeros(1,nz+1), 1];     % slip reference has no dynamics

B = [Blift;
    zeros(2,size(Blift,2))];

% kappa, e_int and kappa_ref are needed to form the cost
C = zeros(3,nz+2);
C(1,3) = 1;     % slip
C(2,end-1) = 1; % integral state
C(3,end) = 1;   % slip reference

ny = size(C,1);     % number of outputs
nu = size(B,2);     % number of inputs

% get dense form matrices (Nc=Np=N)
[F,Phi] = dense_prediction_matrices(A,B,C,N);

%% Quadratic cost parameters
% z = [u y1 y2 y3]';
% Hessian matrix for the penalized states
w_p = 1e-1;     % slip tracking error weight
w_i = 1e3;      % integral state weight
w_u = 1e-5;     % torque reduction weight
 
H = [w_u 0 0 0;
     0 w_p 0 -w_p;
     0 0 w_i 0;
     0 -w_p 0 w_p];
H_u = H(1);
H_x = H(2:end,2:end);
H = blkdiag(kron(eye(N),H_u), kron(eye(N),H_x));  % Hessian for the stacked optimizer

% Linear term vector for the penalized states
f = [];  % set as a parameter below
% f = [-2*wu*T_ref; 0; 0; 0];

%% FORCESPRO multistage form
% note: in order to be compatible with a condensing-based solver,
%       variable ordering has to be z(i) = [u_i; x_i] for i=1...N         
stages = MultistageProblem(1);  % dense formulation, one stage

% dimension
stages(1).dims.n = N*(nu+ny);   % number of stage variables
stages(1).dims.r = N*ny;        % number of equality constraints ??
stages(1).dims.l = N*nu;        % number of lower bounds (input only)
stages(1).dims.u = N*nu;        % number of upper bounds (input only)

% cost
stages(1).cost.H = H;
stages(1).cost.f = f;

% lower bound
umin = mapstd_custom('apply',0,PU);     % scale minimum motor torque
stages(1).ineq.b.lbidx = 1:N;           % lower bound acts on these indices
stages(1).ineq.b.lb = umin*ones(N,1);   % lower bound for this stage variable

% upper bound (parametric)
stages(1).ineq.b.ubidx = 1:N;   % upper bound acts on these indices
stages(1).ineq.b.ub = [];       % upper bound for this stage variable

% equality constraints (RHS is set online)
stages(1).eq.D = [-Phi, eye(N*ny)];

% parameters
params(1) = newParam('F_times_z0',1,'eq.c');        % RHS of first eq. constr. is a parameter: c1 = -x0
params(2) = newParam('T_ref',[],'ineq.b.ub');       % upper bound on the input
params(3) = newParam('linear_cost', [], 'cost.f');  % linear cost varies with the reference torque
outputs(1) = newOutput('u0',1,1);

%% solver settings
codeoptions = getOptions('kmpc');
codeoptions.printlevel = 0;
% codeoptions.condense = 1; % enable state-elimination
codeoptions.overwrite = 1;

if ~compile_for_simulink
    codeoptions.BuildSimulinkBlock = 0;
end
codeoptions.showinfo = 1;  % https://forces.embotech.com/Documentation/solver_options/index.html#solver-info-in-simulink-block

%% generate code
solver_dir = './codegen';
if ~exist(solver_dir,'dir')
    mkdir(solver_dir)
end
addpath(solver_dir);
cd(solver_dir);

generateCode(stages,params,codeoptions,outputs);

cd('../');

%% test call
state_to_lift = mapstd_custom('apply',[0;10],PX);
z0 = [lifting_function(state_to_lift); 0; 0.1];  % lift, e_int, kappa_ref
problem.F_times_z0 = F*z0;

T_ref = mapstd_custom('apply',250,PU);  % T_ref (scaled)
problem.T_ref = T_ref*ones(N,1);

f = [-2*w_u*T_ref; 0; 0; 0];
problem.linear_cost = repmat(f,N,1);

[~, exitflag, ~] = kmpc(problem);
assert(exitflag == 1, 'Some issue with FORCESPRO solver');
save models/dense_Fwu.mat F w_u
end