%% function for simulating one time step of wheel dynamics
function xnext = simulation_model(x,u,p,Ts)
% define function for ode45, t must be used, dim(dy)=dim(y)
odefun = @(t,y) [simulation_model_ode(y(1:2),y(3),y(4));0;0];
% simulate one timestep
[~,y] = ode45(odefun,[0 Ts],[x;u;p]);
% get the last simulated state
xnext = y(end,1:size(x,1))';
end

function xdot = simulation_model_ode(x,u,p)
% TODO: find realistic motor torque limit and driving situation for slip
% wet asphalt? (Fig. 4 in https://ieeexplore.ieee.org/document/6043067)

% system dimensions
% nx = 2;  % vehicle speed, wheel speed
% nu = 1;  % motor torque

% load vehicle parameters
VEHICLE = vehicle_parameters();
m = VEHICLE.MASS;
g = VEHICLE.GRAVITY;
R = VEHICLE.WHEEL_RADIUS;
Jw = VEHICLE.WHEEL_INERTIA;
T_max = VEHICLE.MAX_MOTOR_TORQUE;
gear_ratio = VEHICLE.GEAR_RATIO;

D = VEHICLE.MF_D;
C = VEHICLE.MF_C;
B = VEHICLE.MF_B;

mu_x = p(1);    % [-] tire-road friction coefficient

% state and input variables
v = x(1);           % longitudinal speed of the vehicle [m/s]
w = x(2);           % rotational speed of the wheel [rad/s]
T = gear_ratio * min(u,T_max);   % torque acting on the wheel [Nm]

% dynamics
Fz = m*g;
e0 = 0.1;  % for slip modification
kappa = (w*R-v)*w*R / ((w*R)^2 + e0);
Fx = mu_x * Fz * D*sin(C*atan(B*kappa));
vdot = Fx/m;
wdot = 1/Jw * (T-Fx*R);

xdot = [vdot;wdot];
end
