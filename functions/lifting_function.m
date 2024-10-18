function lifted_state = lifting_function(original_state,varargin)
% load parameters once for faster evaluation
persistent kmpc_data vehicle_params R

if isempty(kmpc_data)
    try
        kmpc_data = load('kmpc_data.mat');
    catch
        warning('kmpc_data.mat not found, continuing anyway...')
    end
    vehicle_params = vehicle_parameters();
    R = vehicle_params.WHEEL_RADIUS;
end

% use the specified kmpc_data if present
if nargin > 1
    kmpc_data = varargin{1};
end

% calculate slip
reversed_state = mapstd_custom('reverse',original_state,kmpc_data.PX);
s = reversed_state(1,:);
w = reversed_state(2,:);

e0 = 0.1;  % for slip modification
slip = s.*w*R ./ ((w*R).^2 + e0);

if strcmp(kmpc_data.rbf_type,'polynomial')
    basis = poly_basis(original_state);
%     basis = basis(3:end,:);  % TODO: poly already contains the original state
else
    basis = rbf(original_state,kmpc_data.cent,kmpc_data.rbf_type);
end

lifted_state = [original_state;
                slip;
                basis];
end