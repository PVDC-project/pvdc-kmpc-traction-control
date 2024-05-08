%% Save all signals entering the KMPC block for additional evaluation
function get_mpc_inputs_kmpc(sigsOut)
z0 = get_sim_data(sigsOut,'mpc_in_z0');
T_ref = get_sim_data(sigsOut,'mpc_in_T_ref');

save('../../data/mpc_inputs_kmpc.mat','z0','T_ref')
disp('MPC inputs saved to data/mpc_inputs_kmpc.mat')
end
