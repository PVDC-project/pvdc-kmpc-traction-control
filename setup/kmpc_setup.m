%% Koopman MPC controller setup
function controller = kmpc_setup(N,Ts,Rw,kappa_ref,compile_for_simulink)
%% system dynamics
load kmpc_data.mat Alift Blift Clift PX;

% add the integral state (e_int_dot = s - kappa_ref*R*w) using Euler method
Alift = [Alift(1:2,1:2), zeros(2,1), Alift(1:2,3:end);
     Ts -Ts*kappa_ref*Rw 1 zeros(1,size(Alift,1)-2);
     Alift(3:end,1:2) zeros(size(Alift,1)-2,1) Alift(3:end,3:end)];

Blift = [Blift(1:2,:)
         zeros(1,size(Blift,2));
         Blift(3:end,:)];

Clift = [eye(3), zeros(3,size(Clift,2)-2)];

ny = size(Clift,1);         % number of outputs
[nz, nu] = size(Blift);     % number of states and inputs in the Koopman model

% get dense form matrices (Nc=Np=N)
[F,Phi] = dense_prediction_matrices(Alift,Blift,Clift,N);

%% MPC parameters
% cost matrices
w_x1 = 1e4;           % wheel slip velocity tracking error weight (unscaled)
w_x3 = 0*1e2;      % integral state weight (unscaled)
w_u = 1e-6;         % torque reduction weight (scaled)

% w_x1 = 1e-6; w_x3 = 1e-6; w_u = 1e6;  % test MPC functionality

Q = diag([w_x1 0 w_x3]);    % state cost, no penalty on wheel speed
R = w_u;                    % input cost

%% define problem using YALMIP
u = sdpvar(nu,N);
z0 = sdpvar(nz,1);
T_ref = sdpvar;  % for online torque limiting
Y = F*z0 + Phi*u(:);  % dense-form prediction

constraints = [];
objective = 0;

for k = 1:N
    y = Y((k-1)*ny+1:k*ny);  % output from k-th prediction step
    s = (y(1)-PX.xoffset(1)) / PX.gain(1);  % unscale the first state
    w = (y(2)-PX.xoffset(2)) / PX.gain(2);  % unscale the second state
    es = s - kappa_ref*Rw*w;  % s tracking error (unscaled)
    ey = [es;0;y(3)];  % state error (ignore y(2), drive y(3) to zero)
    objective = objective + ey'*Q*ey;  % tracking cost
    eu = T_ref-u(k);  % input torque reduction
    objective = objective + eu'*R*eu;  % input torque reduction cost
    constraints = [constraints, 0 <= u(k) <= T_ref];  % input constraint
end

%% solver settings
algorithm = 'PDIP';  % https://forces.embotech.com/Documentation/solver_options/index.html#solve-methods
codeoptions = ForcesGetDefaultOptions('kmpc',algorithm,'double');
codeoptions.printlevel = 0;  % summary line after each solve
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
if ~compile_for_simulink
    codeoptions.BuildSimulinkBlock = 0;
end
codeoptions.showinfo = 1;  % https://forces.embotech.com/Documentation/solver_options/index.html#solver-info-in-simulink-block

params = {z0,T_ref};
outputs = {u,Y};
parameter_names = {'z0','T_ref'};
output_names = {'uopt','yopt'};

%% generate code
solver_dir = './codegen';
if ~exist(solver_dir,'dir')
    mkdir(solver_dir)
end
addpath(solver_dir);
cd(solver_dir);

% optimizerFORCES(constraints, objective, codeoptions, params, outputs, parameter_names, output_names);

cd('../');

% YALMIP for comparison
options = sdpsettings('solver','daqp','savesolverinput',1,'savesolveroutput',1);
controller = optimizer(constraints, objective, options, params, outputs);

%% test call
disp([newline,'Testing the generated solver...'])
state_to_lift = mapminmax('apply',[0;10],PX);
problem{1} = [lifting_function(state_to_lift); 0];  % z0 (no integral state)
problem{2} = 1;  % T_ref (scaled)

% [forces_sol, exitflag, info] = kmpc(problem);
% disp(['FORCES test solution: ',num2str(forces_sol{1})])
% assert(exitflag == 1, 'Test call of FORCESPRO solver failed');

[yalmip_sol, yalmip_flag, ~, ~, ~, yalmip_struct] = controller(problem);
disp(['YALMIP test solution: ',num2str(yalmip_sol{1})])
assert(yalmip_flag == 0, 'Test call with YALMIP failed')

% assert(max(abs(yalmip_sol{1}-forces_sol{1})) < 1e-6, 'FORCES and YALMIP solutions differ too much.')
end