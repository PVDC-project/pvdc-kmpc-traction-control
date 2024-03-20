function xdot = nmpc_prediction_model(x,u,p,kappa_ref)

% system dimensions
% nx = 3;  % wheel slip velocity, wheel speed, integral state
% nu = 1;  % motor torque

% load vehicle parameters
VEHICLE = vehicle_parameters();
m = VEHICLE.MASS;
g = VEHICLE.GRAVITY;
R = VEHICLE.WHEEL_RADIUS;
Jw = VEHICLE.WHEEL_INERTIA;
gear_ratio = VEHICLE.GEAR_RATIO;

D = VEHICLE.MF_D;
C = VEHICLE.MF_C;
B = VEHICLE.MF_B;

mu_x = p(2);        % [-] tire-road friction coefficient

% state and input variables
s = x(1);           % wheel slip velocity (w*R-v) [m/s]
w = x(2);           % rotational speed of the wheel [rad/s]
% e_int = x(3);     % integral of wheel slip velocity error [m]
T = gear_ratio*u;   % torque acting on the wheel [Nm]
% T_ref = p;        % desired wheel torque [Nm]

e0 = 0.1;  % for slip modification
kappa = s*w*R / ((w*R)^2 + e0);
% kappa = s/(w*R);
% kappa = kappa .* 1./(1+exp(-5*(w-1.5)));

% dynamics
Fz = m*g;
Fx = mu_x * Fz * D*sin(C*atan(B*kappa));
sdot = (-R^2/Jw-1/m) * Fx + T*R/Jw;
w_dot = 1/Jw * (T-Fx*R);
e_int_dot = kappa_ref - kappa;

xdot = [sdot; w_dot; e_int_dot];

end
