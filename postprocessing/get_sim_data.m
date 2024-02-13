function Values = get_sim_data(out,str)
    % Extract values from logged simulation data by name
    % Output: rows are samples
    idx = find(strcmp(out.getElementNames,str));
    Values = out{idx}.Values.Data;
end