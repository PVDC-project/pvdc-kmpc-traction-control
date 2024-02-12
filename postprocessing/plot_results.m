t_sim = 0:Ts:N_sim*Ts;  % time vector for plotting

% controller internal states
figure;
subplot(3,1,1)
plot(t_sim,x_ocp(1,:) * 3.6)
ylabel('$s$ [km/h]')
subplot(3,1,2)
plot(t_sim, x_ocp(2,:) * 30/pi);
ylabel('$\omega$ [rpm]')
subplot(3,1,3)
plot(t_sim,x_ocp(3,:))
ylabel('$e_{int}$ [m]')
xlabel('Time [s]')

% solver performance
figure;
subplot(2,1,1)
plot(t_sim(1:end-1),1e3*solve_time_log,'o')
yline(1e3*Ts,'r--')
ylabel('Solve time [ms]')
subplot(2,1,2)
plot(t_sim(1:end-1),status_log,'o')
ylabel('Solver status')
xlabel('Time [s]')

% closed-loop response
figure;
% vehicle speed
subplot(4,1,1);
plot(t_sim, x_sim(1,:) * 3.6);
ylabel('$v_x$ [km/h]')
% wheel speed
subplot(4,1,2);
plot(t_sim, x_sim(2,:) * 30/pi);
ylabel('$\omega$ [rpm]')
% wheel torque
subplot(4,1,3);
plot(t_sim(1:end-1), u_sim);
hold on
plot(t_sim(1:end-1), T_ref,'r--');
ylabel('$T$ [Nm]')
% slip
subplot(4,1,4);
vx = x_sim(1,:);
omega = x_sim(2,:);
slip = (omega*R-vx) ./ (omega*R);
plot(t_sim,slip)
hold on
yline(kappa_ref,'r--')
ylabel('$\sigma_x$ [-]')
xlabel('Time [s]')