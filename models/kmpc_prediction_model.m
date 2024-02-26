function model = kmpc_prediction_model()

% not used!

import casadi.*

% system dimensions
nx = 3;  % wheel slip velocity, integral state, wheel speed
nu = 1;  % wheel torque reduction

% system parameters
m = 1600/4;         % [kg] quarter car mass
g = 9.81;           % [m/s^2] gravity constant
R = 0.318;          % [m] wheel radius
Jw = 1.22;          % [kg*m^2] moment of inertia for one wheel and half-axle
T_max = 250;        % [Nm] maximum motor torque
kappa_ref = 0.1;    % [-] longitudinal slip reference
mu_x = 0.15;        % [-] tire-road friction coefficient

% tire model coefficients Fx=Fz*D*sin(C*arctan(B*kappa))
D = 1.2069;
C = 1.35;
K = 30;
B = K/(C*D);

% named symbolic variables
s = SX.sym('s');            % wheel slip velocity (w*R-v) [m/s]
e_int = SX.sym('e_int');    % integral of wheel slip velocity error [m]
w = SX.sym('w');            % rotational speed of the wheel [rad/s]
T_ref = SX.sym('T_ref');    % desired wheel torque [Nm]
dT = SX.sym('dT');          % wheel torque reduction [Nm]

% (unnamed) symbolic variables
sym_x = vertcat(s,e_int,w);
sym_xdot = SX.sym('xdot', nx, 1);
sym_u = dT;
sym_p = T_ref;  % include the reference as a parameter

% dynamics
Fz = m*g;
Fx = mu_x * Fz * D*sin(C*atan(B*s/(w*R)));
sdot = (-R^2/Jw-1/m) * Fx + (T_ref-dT)*R/Jw;
e_int_dot = s - kappa_ref*w*R;
w_dot = 1/Jw * (T_ref-dT-Fx*R);
dyn_expr_f_expl = vertcat(sdot,e_int_dot,w_dot);
dyn_expr_f_impl = dyn_expr_f_expl - sym_xdot;

% cost
w_x1 = 1e3;     % wheel slip velocity tracking error weight
w_x2 = 1e1;     % integral state weight
w_u = 1e-2;     % control input weight
cost_expr_ext_cost_e = w_x1 * (s-kappa_ref*w*R)^2 + ...
                  w_x2 * e_int^2;
cost_expr_ext_cost = cost_expr_ext_cost_e + w_u * dT^2;

% torque constraint: 0 <= dT <= T_ref
model.constr_expr_h = vertcat(dT,dT-T_ref);
model.constr_lh = [0; -1e10];   % 0 <= dT <= inf
model.constr_uh = [1e10; 0];    % -inf <= dT - T_ref <= 0

% populate structure
model.nx = nx;
model.nu = nu;
model.sym_x = sym_x;
model.sym_xdot = sym_xdot;
model.sym_u = sym_u;
model.sym_p = sym_p;
model.dyn_expr_f_expl = dyn_expr_f_expl;
model.dyn_expr_f_impl = dyn_expr_f_impl;
model.cost_expr_ext_cost = cost_expr_ext_cost;
model.cost_expr_ext_cost_e = cost_expr_ext_cost_e;

model.wheel_radius = R;
model.max_torque = T_max;

end
