function xdot = simulation_model_ode(x,u,p)
% TODO: find realistic motor torque limit and driving situation for slip
% wet asphalt? (Fig. 4 in https://ieeexplore.ieee.org/document/6043067)

% system dimensions
% nx = 2;  % vehicle speed, wheel speed
% nu = 1;  % wheel torque

% system parameters
m = 1600/4;     % [kg] quarter car mass
g = 9.81;       % [m/s^2] gravity constant
R = 0.318;      % [m] wheel radius
Jw = 1.22;      % [kg*m^2] moment of inertia for one wheel and half-axle
T_max = 500;    % [Nm] maximum motor torque
mu_x = p(1);    % [-] tire-road friction coefficient

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

% state and input variables
v = x(1);           % longitudinal speed of the vehicle [m/s]
w = x(2);           % rotational speed of the wheel [rad/s]
T = min(u,T_max);   % motor torque acting on the wheel [Nm]

% dynamics
Fz = m*g;
kappa = (w*R-v)/(w*R);
Fx = mu_x * Fz * D*sin(C*atan(B*kappa));
vdot = Fx/m;
wdot = 1/Jw * (T-Fx*R);

xdot = [vdot;wdot];
end
