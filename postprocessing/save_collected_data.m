%% Save the data for later use in Koopman operator identification
function save_collected_data(sigsOut)
w = get_sim_data(sigsOut,'w');
mu = get_sim_data(sigsOut,'mu');
vx = get_sim_data(sigsOut,'vx');
tc_torque = get_sim_data(sigsOut,'tc_torque');

save('../../data/collected_data.mat','w','mu','vx','tc_torque')
disp('Simulation output saved to data/collected_data.mat')
end