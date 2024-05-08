%% Save all signals entering the NMPC block for additional evaluation
function get_mpc_inputs_nmpc(sigsOut)
hu = get_sim_data(sigsOut,'mpc_in_hu');
hl = get_sim_data(sigsOut,'mpc_in_hl');
x0 = get_sim_data(sigsOut,'mpc_in_x0');
xinit = get_sim_data(sigsOut,'mpc_in_xinit');
all_parameters = get_sim_data(sigsOut,'mpc_in_all_params');

save('../../data/mpc_inputs_nmpc.mat','hu','hl','x0','xinit','all_parameters')
disp('MPC inputs saved to data/mpc_inputs_nmpc.mat')
end
