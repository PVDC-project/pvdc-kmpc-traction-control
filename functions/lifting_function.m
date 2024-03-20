function lifted_state = lifting_function(original_state)
persistent kmpc_data
if isempty(kmpc_data)
    kmpc_data = load('kmpc_data.mat');
end

% don't lift the integral state (if present)
if size(original_state,1) == 3
    state_to_lift = original_state(1:2,:);
else
    state_to_lift = original_state;
end

% basis function selection
cent = kmpc_data.cent;
rbf_type = kmpc_data.rbf_type;

% lifted state = [original state; basis functions]
if size(original_state,1) == 3
    lifted_state = [state_to_lift; original_state(3,:); rbf(state_to_lift,cent,rbf_type)];
else
    lifted_state = [state_to_lift; rbf(state_to_lift,cent,rbf_type)];
end

end