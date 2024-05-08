% this script plots data from different simulations for comparison
clear;clc;close all;

lw = 2;  % line width for plots

%% constant mu
nmpc_1 = load('../data/nmpc_1.mat');
kmpc_1 = load('../data/kmpc_1.mat');

kappa_nmpc_1 = get_sim_data(nmpc_1.sigsOut,'kappa');
kappa_kmpc_1 = get_sim_data(kmpc_1.sigsOut,'kappa');
Tm_nmpc_1 = get_sim_data(nmpc_1.sigsOut,'tc_torque');
Tm_kmpc_1 = get_sim_data(kmpc_1.sigsOut,'tc_torque');
Tm_ref = get_sim_data(nmpc_1.sigsOut,'trq_ref');

t_end = 3;  % s
k_end = t_end/2e-3;
k_end_1 = k_end;
t = 0:2e-3:t_end-2e-3;

figure
subplot(2,1,1)
plot(t,kappa_nmpc_1(1:k_end,1),'linewidth', lw)
hold on
plot(t,kappa_kmpc_1(1:k_end,1),'linewidth', lw)
yline(0.1,'k--','linewidth', lw)
ylabel('$\kappa$ [-]')
grid on
legend({'NMPC','KMPC'},'Location','southeast')

subplot(2,1,2)
stairs(t,Tm_nmpc_1(1:k_end,1),'linewidth', lw)
hold on
stairs(t,Tm_kmpc_1(1:k_end,1),'linewidth', lw)
plot(t,Tm_ref(1:k_end),'k--','linewidth', lw)
ylabel('$T_m$ [Nm]')
grid on

xlabel('Time [s]')

%% varying mu
nmpc_2 = load('../data/nmpc_2.mat');
kmpc_2 = load('../data/kmpc_2.mat');

kappa_nmpc_2 = get_sim_data(nmpc_2.sigsOut,'kappa');
kappa_kmpc_2 = get_sim_data(kmpc_2.sigsOut,'kappa');
Tm_nmpc_2 = get_sim_data(nmpc_2.sigsOut,'tc_torque');
Tm_kmpc_2 = get_sim_data(kmpc_2.sigsOut,'tc_torque');
Tm_ref = get_sim_data(nmpc_2.sigsOut,'trq_ref');
mu_x = get_sim_data(nmpc_2.sigsOut,'mu');
mu_x = mu_x(:,1);
idxs = find(diff(mu_x))+1;  % plot jumps in mu_x

t_end = 8;  % s
k_end = t_end/2e-3;
k_end_2 = k_end;
t = 0:2e-3:t_end-2e-3;

figure
subplot(2,1,1)
plot(t,kappa_nmpc_2(1:k_end,1),'linewidth', lw)
hold on
plot(t,kappa_kmpc_2(1:k_end,1),'linewidth', lw)
yline(0.1,'k--','linewidth', lw)
for idx = idxs'
    if idx<length(t)
        xline(t(idx),':',num2str(mu_x(idx+1)),'LabelVerticalAlignment','bottom','LabelOrientation','horizontal','linewidth', lw); 
    end
end
ylabel('$\kappa$ [-]')
limsy=get(gca,'YLim');
set(gca,'Ylim',[0 limsy(2)]);
grid on
legend({'NMPC','KMPC'},'Location','southwest')

subplot(2,1,2)
stairs(t,Tm_nmpc_2(1:k_end,1),'linewidth', lw)
hold on
stairs(t,Tm_kmpc_2(1:k_end,1),'linewidth', lw)
plot(t,Tm_ref(1:k_end),'k--','linewidth', lw)
for idx = idxs'
    if idx<length(t)
        xline(t(idx),':',num2str(mu_x(idx+1)),'LabelVerticalAlignment','bottom','LabelOrientation','horizontal','linewidth', lw); 
    end
end
ylabel('$T_m$ [Nm]')
grid on

xlabel('Time [s]')

%% execution time
nmpc_1_t = get_sim_data(nmpc_1.sigsOut,'mpc_solvetime_ms');
kmpc_1_t = get_sim_data(kmpc_1.sigsOut,'mpc_solvetime_ms');
nmpc_2_t = get_sim_data(nmpc_2.sigsOut,'mpc_solvetime_ms');
kmpc_2_t = get_sim_data(kmpc_2.sigsOut,'mpc_solvetime_ms');

nmpc_1_t = squeeze(nmpc_1_t);
kmpc_1_t = squeeze(kmpc_1_t);
nmpc_2_t = squeeze(nmpc_2_t);
kmpc_2_t = squeeze(kmpc_2_t);

k_start = 10;

nmpc_1_t = nmpc_1_t(k_start:k_end_1);
kmpc_1_t = kmpc_1_t(k_start:k_end_1);
nmpc_2_t = nmpc_2_t(k_start:k_end_2);
kmpc_2_t = kmpc_2_t(k_start:k_end_2);

disp('constant mu_x:')
disp(['NMPC mean: ',num2str(mean(nmpc_1_t))])
disp(['NMPC max: ',num2str(max(nmpc_1_t))])
disp(['KMPC mean: ',num2str(mean(kmpc_1_t))])
disp(['KMPC max: ',num2str(max(kmpc_1_t))])

disp(' ')

disp('varying mu_x:')
disp(['NMPC mean: ',num2str(mean(nmpc_2_t))])
disp(['NMPC max: ',num2str(max(nmpc_2_t))])
disp(['KMPC mean: ',num2str(mean(kmpc_2_t))])
disp(['KMPC max: ',num2str(max(kmpc_2_t))])

%% sysid
sysid = load('../data/sysid.mat');
t = 0:2e-3:(length(sysid.kappa_true)-1)*2e-3;
t = t*1e3;
figure
plot(t,sysid.kappa_true,'linewidth', lw)
hold on
grid on
plot(t,sysid.kappa_koop,'linewidth', lw)
xlabel('Time [ms]')
ylabel('$\kappa$ [-]')
legend({'True','Koopman'})