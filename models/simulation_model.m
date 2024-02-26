function model = simulation_model()
% TODO: find realistic motor torque limit and driving situation for slip
% wet asphalt? (Fig. 4 in https://ieeexplore.ieee.org/document/6043067)
% Burckhardt model: p.24 in https://sci-hub.st/10.1007/978-1-84996-350-3

import casadi.*

% system dimensions
nx = 2;  % vehicle speed, wheel speed
nu = 1;  % wheel torque

% system parameters
m = 1400/4;     % [kg] quarter car mass
g = 9.81;       % [m/s^2] gravity constant
R = 0.318;      % [m] wheel radius
Jw = 1.22;      % [kg*m^2] moment of inertia for one wheel and half-axle
T_max = 500;    % [Nm] maximum motor torque
% mu_x = 0.15;    % [-] tire-road friction coefficient

% tire model coefficients Fx=mu_x*Fz*D*sin(C*arctan(B*kappa))
% Carmaker MF_205_60R15_V91.tir
D = 1.2005;
C = 1.4010;
B = 17.1571;

% Saab?
% D = 1.2069;
% C = 1.35;
% K = 30;
% B = K/(C*D);

% named symbolic variables
v = SX.sym('v');    % longitudinal speed of the vehicle [m/s]
w = SX.sym('w');    % rotational speed of the wheel [rad/s]
T = SX.sym('T');    % motor moment acting on the wheel [Nm]
mu_x = SX.sym('mu_x');  % tire-road friction coefficient [-]

% generic symbolic variables
sym_x = vertcat(v, w);
sym_xdot = SX.sym('xdot', nx, 1);
sym_u = T;
sym_p = mu_x;

% dynamics
Fz = m*g;
e0 = 0.1;  % for slip modification
kappa = (w*R-v)*w*R / ((w*R)^2 + e0);
Fx = mu_x * Fz * D*sin(C*atan(B*kappa));
expr_f_expl = vertcat(Fx/m, ...
                      1/Jw * (T-Fx*R));
expr_f_impl = expr_f_expl - sym_xdot;

% populate structure
model.nx = nx;
model.nu = nu;
model.sym_x = sym_x;
model.sym_xdot = sym_xdot;
model.sym_u = sym_u;
model.sym_p = sym_p;
model.expr_f_expl = expr_f_expl;
model.expr_f_impl = expr_f_impl;

model.wheel_radius = R;
model.max_torque = T_max;

end
