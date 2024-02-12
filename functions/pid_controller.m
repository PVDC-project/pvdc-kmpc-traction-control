function output = pid_controller(setpoint, process_variable, Kp, Ki, Kd, Tref)
    persistent integral_term prev_error;

    % Initialize variables on the first run
    if isempty(integral_term)
        integral_term = 0;
        prev_error = 0;
    end

    % Calculate error
    error = setpoint - process_variable;

    % Proportional term
    P = Kp * error;

    % Integral term
    integral_term = integral_term + Ki * error;
    
    % Anti-windup (limit the integral term with asymmetric limits)
    if integral_term > Tref
        integral_term = Tref;
    elseif integral_term < 0
        integral_term = 0;
    end

    % Derivative term
    derivative_term = Kd * (error - prev_error);

    % Calculate the control output
    output = P + integral_term + derivative_term;
    
    % Saturate the control output
    if output > Tref
        output = Tref;
    elseif output < 0
        output = 0;
    end

    % Update previous error for the next iteration
    prev_error = error;
end