t_sim = 0:Ts:N_sim*Ts;  % time vector for plotting

% controller internal states
% f1 = figure;
% subplot(3,1,1)
% plot(t_sim,x_ocp(1,:) * 3.6)
% ylabel('$s$ [km/h]')
% subplot(3,1,2)
% plot(t_sim, x_ocp(2,:) * 30/pi);
% ylabel('$\omega$ [rpm]')
% subplot(3,1,3)
% plot(t_sim,x_ocp(3,:))
% ylabel('$e_{int}$ [-]')
% xlabel('Time [s]')

% solver performance
f2 = figure;
subplot(2,1,1)
plot(t_sim(1:end-1),1e3*solve_time_log,'o')
yline(1e3*Ts,'r--')
ylabel('Solve time [ms]')
subplot(2,1,2)
plot(t_sim(1:end-1),status_log,'o')
ylabel('Solver status')
xlabel('Time [s]')

% closed-loop response
f3 = figure;
% vehicle speed
subplot(4,1,1);
plot(t_sim, x_sim(1,:) * 3.6);
if is_adaptive
    hold on
    idxs = find(diff(vid))+1;
    for idx = idxs
        xline(t_sim(idx),':',num2str(vid(idx+1)),'LabelVerticalAlignment','bottom','LabelOrientation','horizontal'); 
    end
end
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
v = x_sim(1,:);
w = x_sim(2,:);
e0 = 0.1;  % for slip modification
slip = (w*R-v).*w*R ./ ((w*R).^2 + e0);
plot(t_sim,slip)
hold on
yline(kappa_ref,'r--')
idxs = find(diff(mu_x_log))+1;  % plot jumps in mu_x
for idx = idxs
    xline(t_sim(idx),':',num2str(mu_x_log(idx+1)),'LabelVerticalAlignment','bottom','LabelOrientation','horizontal'); 
end
ylabel('$\kappa$ [-]')
xlabel('Time [s]')

%% Arrange the figures on screen
toolbar_height = 40;    % Windows
header_height = 85;     % MATLAB

monitors = get(0, 'MonitorPositions');
if size(monitors,1) > 1
    monitor = monitors(2, :);
else
    monitor = monitors;
end
x = monitor(1);
y = monitor(2);
w = monitor(3);
h = monitor(4) - toolbar_height;

quarter_x = x + w/2;
quarter_y = y + h/2 + toolbar_height;
quarter_w = w/2;
quarter_h = h/2 - header_height;

set(f2, 'Position', [quarter_x, quarter_y, quarter_w, quarter_h]);
quarter_y = y + toolbar_height;
set(f3, 'Position', [quarter_x, quarter_y, quarter_w, quarter_h]);