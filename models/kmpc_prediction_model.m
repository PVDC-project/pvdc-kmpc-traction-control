%% Koopman MPC prediction model and partial setup (for acados)
function model = kmpc_prediction_model(mpc_setup)
import casadi.*
N = mpc_setup.N;
Ts = mpc_setup.Ts;

%% sytem dynamics
load kmpc_data.mat Alift Blift PU;

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

%% (unnamed) symbolic variables
Y = SX.sym('Y',ny*N,1);     % stacked output vector
U = SX.sym('U',nu*N,1);     % stacked input vector
z0 = SX.sym('z0',nz,1);     % initial (lifted) state
T_ref = SX.sym('T_ref');    % torque reference for online limiting, should be scaled
sym_p = vertcat(z0,T_ref);  % parameter vector

%% discrete system dynamics
dyn_expr_phi = F * z0 + Phi * U;

%% cost definition
w_p = 1e-1;     % slip tracking error weight
w_i = 1e3;      % integral state weight
w_u = 1e-5;     % torque reduction weight

if isfield(mpc_setup,'w_p')  % cost weights set outside
    w_p = mpc_setup.w_p;
    w_i = mpc_setup.w_i;
    w_u = mpc_setup.w_u;
end

cost_expr_ext_cost_0 = 0;  % penalize only the control at the initial node
cost_expr_ext_cost_e = 0;  % penalize only the outputs at the terminal node
for k = 1:N
    y = Y((k-1)*ny+1:k*ny);  % output from the k-th prediction step
    u = U((k-1)*nu+1:k*nu);  % input from the (k-1)-th prediction step
    cost_expr_ext_cost_e = cost_expr_ext_cost_e + w_p * (y(3)-y(1))^2;  % slip tracking cost
    cost_expr_ext_cost_e = cost_expr_ext_cost_e + w_i * y(2)^2;         % integral state cost
    cost_expr_ext_cost_0 = cost_expr_ext_cost_0 + w_u * (T_ref-u)^2;    % input torque reduction cost
end

%% constraints
% torque constraint: 0 (scaled) <= u <= T_ref
umin = mapstd_custom('apply',0,PU);  % scale minimum motor torque (generally not zero)
model.constr_expr_h_0 = vertcat(U,U-T_ref);
model.constr_lh_0 = [umin*ones(N,1); -1e6*ones(N,1)];   % 0 <= u <= inf
model.constr_uh_0 = [1e6*ones(N,1); zeros(N,1)];        % -inf <= u - T_ref <= 0

%% populate structure
model.sym_x = Y;
model.sym_u = U;
model.sym_p = sym_p;
model.dyn_expr_phi = dyn_expr_phi;
model.cost_expr_ext_cost_0 = cost_expr_ext_cost_0;
model.cost_expr_ext_cost_e = cost_expr_ext_cost_e;
% (no intermediate cost)
end
