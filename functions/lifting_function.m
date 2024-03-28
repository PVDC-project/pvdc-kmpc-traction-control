function lifted_state = lifting_function(original_state)
% load parameters once for faster evaluation
persistent kmpc_data vehicle_params R
if isempty(kmpc_data)
    kmpc_data = load('kmpc_data.mat');
    vehicle_params = vehicle_parameters();
    R = vehicle_params.WHEEL_RADIUS;
end

% calculate slip
reversed_state = mapstd_custom('reverse',original_state,kmpc_data.PX);
s = reversed_state(1,:);
w = reversed_state(2,:);

e0 = 0.1;  % for slip modification
slip = s.*w*R ./ ((w*R).^2 + e0);

% lifted state = [original state; slip; basis functions]
lifted_state = [original_state;
                slip;
                rbf(original_state,kmpc_data.cent,kmpc_data.rbf_type)];

end