%% Nonlinear MPC controller setup
function [] = nmpc_setup(mpc_setup,compile_for_dspace)
N = mpc_setup.N; Ts = mpc_setup.Ts; R = mpc_setup.R; kappa_ref = mpc_setup.kappa_ref;
if nargin < 2
    compile_for_dspace = 0;
end

model = {};
model.N = N;    % prediction horizon
model.nvar = 4; % number of stage variables, nx+nu
model.neq = 3;  % number of equality constraints, nx
model.nh = 1;   % number of nonlinear inequality constraints
model.npar = 2; % number of runtime parameters

% model.objective = @(z,p) eval_obj(z,p,R,kappa_ref);
model.LSobjective = @(z,p) LSobjective(z,p,R,kappa_ref);

model.continuous_dynamics = @(x,u,p) nmpc_prediction_model(x,u,p,kappa_ref);
model.E = [zeros(3,1), eye(3)]; % selection matrix

model.xinitidx = 2:4;  % indices affected by initial condition

model.ineq = @eval_const;  % handle to nonlinear constraints
model.hlidx = 1;  % only the input is constrained
model.huidx = 1;  
model.hu = [];  % set to T_ref at runtime
model.hl = [];  % set to 0 at runtime

algorithm = 'PDIP_NLP';  % https://forces.embotech.com/Documentation/solver_options/index.html#solve-methods
codeoptions = ForcesGetDefaultOptions('nmpc',algorithm,'double');
codeoptions.reinitialize = 0;
codeoptions.printlevel = 0;
codeoptions.nlp.integrator.type = 'ERK4';
codeoptions.nlp.integrator.Ts = Ts;
codeoptions.nlp.integrator.nodes = 4;
codeoptions.nlp.stack_parambounds = 1;
codeoptions.nlp.ad_tool = 'casadi';
codeoptions.nlp.hessian_approximation = 'gauss-newton';
codeoptions.nlp.linear_solver = 'symm_indefinite';
codeoptions.nlp.strictCheckDistinctStages = 1;

% codeoptions.init = 2;  % primal warm start 
% https://forces.embotech.com/Documentation/solver_options/index.html#solver-initialization

% Numerical instability debugging
%https://forces.embotech.com/Documentation/high_level_interface/index.html#call-callback-ref
codeoptions.MEXinterface.dynamics = 1;
codeoptions.MEXinterface.inequalities = 1;
codeoptions.MEXinterface.objective = 1;

% Simulink block options
if ~mpc_setup.compile_for_simulink
    codeoptions.BuildSimulinkBlock = 0;
end
codeoptions.showinfo = 1;  % https://forces.embotech.com/Documentation/solver_options/index.html#solver-info-in-simulink-block

% Codegen options
if compile_for_dspace
    codeoptions.platform = 'dSPACE-MicroLabBox';  % to specify the platform
    codeoptions.printlevel = 0;  % on some platforms printing is not supported
    codeoptions.cleanup = 0;  % to keep necessary files for target compile    
    codeoptions.timing = 1;
    codeoptions.embedded_timing = 1;    
    
    codeoptions.optimize_choleskydivision = 1;
    codeoptions.optimize_registers = 1;
    codeoptions.optimize_uselocalsheavy = 1;
    codeoptions.optimize_operationsrearrange = 1;
    codeoptions.optimize_loopunrolling = 1;
    codeoptions.optimize_enableoffset = 1;
end
codeoptions.overwrite = 1;
codeoptions.nohash = 1;
solver_dir = './codegen';
if ~exist(solver_dir,'dir')
    mkdir(solver_dir)
end
addpath(solver_dir);
cd(solver_dir);

% output = newOutput('u0', 1, 1);
output = newOutput('zopt');
FORCES_NLP(model, codeoptions, output);

cd('../');

%% Test call
problem.xinit = [1; 10; 0];  % initial state
problem.x0 = repmat([250;1;10;0], model.N, 1);  % solver initial guess
problem.hu = repmat(250, model.N, 1);
problem.hl = zeros(model.N, 1);
problem.all_parameters = repmat([250;0.15], model.N, 1);
% problem.z_init_0 = repmat([250;1;1;0], model.N, 1);

[~, exitflag, ~] = nmpc(problem);
assert(exitflag == 1, 'Some issue with FORCESPRO solver');
disp('Solver check succesful.')
end

%% Auxiliary functions
% function f = eval_obj(z,p,R,kappa_ref)
%     w_x1 = 1e3;         % wheel slip velocity tracking error weight
%     w_x3 = 1e1;         % integral state weight
%     w_u = 1e-2;         % torque reduction weight
%     f = w_x1*(z(2)-kappa_ref*z(3)*R) + ...  % tracking error cost
%         w_x3*z(4)^2 + ...                   % integral state cost
%         w_u*(z(1)-p(1))^2;                  % torque reduction cost
% end

function f = LSobjective(z,p,R,kappa_ref)
    w_x1 = 1e3 / 1^2;       % wheel slip velocity tracking error weight
    w_x3 = 1e5 / 1^2;       % integral state weight
    w_u = 1e-2 / (250^2);   % torque reduction weight
    s = z(2); w = z(3);
    e0 = 0.1;  % for slip modification
    kappa = s*w*R / ((w*R)^2 + e0);
%     kappa = s/(w*R);
%     kappa = kappa .* 1./(1+exp(-5*(w-1.5)));
    f = [sqrt(2*w_u)*(z(1)-p(1));                       % torque reduction cost
         sqrt(2*w_x1)*(kappa_ref-kappa);   % tracking error cost  
         sqrt(2*w_x3)*z(4)];                            % integral state cost
% sqrt(2*w_x1)*(z(2)-kappa_ref*z(3)*R);  % tracking error cost
% sqrt(2*w_x1)*(kappa_ref-z(2)/(z(3)*R+1e-3));  % tracking error cost  
end


function h = eval_const(z)
    h = z(1);  % lower and upper bounds on the input
end