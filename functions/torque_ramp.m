function T_ref = torque_ramp(time_start, time_end, value_start, value_end, step_size, time_sim)
% go from VALUE_START to VALUE_END between TIME_START and TIME_END, zero
% before, VALUE_END after
t = (0:step_size:time_sim-step_size);  % time vector
rampOn = @(t) t>=time_start & t<=time_end;  % zero value outside of the specified interval
ramp = @(t) value_start + (t-time_start)/(time_end-time_start) .* (value_end-value_start) .* rampOn(t);
T_ref = ramp(t);
T_ref(time_end/step_size:end) = value_end;  % hold the end value for the rest of the simulation