function VEHICLE = vehicle_parameters()
VEHICLE.MASS = 1400/4;          % [kg] quarter car mass
VEHICLE.GRAVITY = 9.81;         % [m/s^2] gravity constant
VEHICLE.WHEEL_RADIUS = 0.318;   % [m] wheel radius
VEHICLE.WHEEL_INERTIA = 1.22;   % [kg*m^2] moment of inertia for one wheel and half-axle
VEHICLE.MAX_MOTOR_TORQUE = 250; % [Nm] maximum motor torque
VEHICLE.GEAR_RATIO = 3;         % [-] ratio of wheel torque and engine torque

% tire model coefficients Fx=mu_x*Fz*D*sin(C*arctan(B*kappa))
% Carmaker MF_205_60R15_V91.tir
VEHICLE.MF_D = 1.2005;
VEHICLE.MF_C = 1.4010;
VEHICLE.MF_B = 17.1571;

% Saab?
% D = 1.2069;
% C = 1.35;
% K = 30;
% B = K/(C*D);
end