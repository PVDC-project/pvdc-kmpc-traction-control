function lifted_state = lifting_function(original_state)
% don't lift the integral state
state_to_lift = original_state(1:2);

% scale the state before lifting
PX = load('kmpc_data.mat','PX').PX;
state_to_lift = mapminmax('apply',state_to_lift,PX);

% basis function selection and lifting (TODO: use monomials?)
cent = load('kmpc_data.mat','cent').cent;
rbf_type = 'thinplate';  % gauss, invquad, invmultquad, polyharmonic

% lifted state = [original state; basis functions]
lifted_state = [state_to_lift; original_state(3); rbf(state_to_lift,cent,rbf_type)];
end