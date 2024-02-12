function xdot = nmpc_prediction_model(x,u,p,kappa_ref)

% system dimensions
% nx = 3;  % wheel slip velocity, wheel speed, integral state
% nu = 1;  % motor torque

% system parameters
m = 1400/4;         % [kg] quarter car mass
g = 9.81;           % [m/s^2] gravity constant
R = 0.318;          % [m] wheel radius
Jw = 2.3;           % [kg*m^2] moment of inertia for one wheel and half-axle
mu_x = p(2);        % [-] tire-road friction coefficient
gear_ratio = 3;     % [-] ratio of wheel torque and engine torque
% kappa_ref = 0.1;    % [-] longitudinal slip reference
% mu_x = 0.15;        
% T_max = 250;        % [Nm] maximum motor torque

% tire model coefficients Fx=Fz*D*sin(C*arctan(B*kappa))
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
s = x(1);           % wheel slip velocity (w*R-v) [m/s]
w = x(2);           % rotational speed of the wheel [rad/s]
% e_int = x(3);     % integral of wheel slip velocity error [m]
T = gear_ratio*u;   % torque acting on the wheel [Nm]
% T_ref = p;        % desired wheel torque [Nm]

% dynamics
Fz = m*g;
Fx = mu_x * Fz * D*sin(C*atan(B*s/(w*R)));
sdot = (-R^2/Jw-1/m) * Fx + T*R/Jw;
w_dot = 1/Jw * (T-Fx*R);
e_int_dot = kappa_ref - s/(w*R);

xdot = [sdot; w_dot; e_int_dot];

end
