function [psi,npsi] = create_poly_basis(nx, nd, toFile)
% POLY BASIS
% Function creates polynomial basis up to order nd based on 
% the variable of size nx.
% Function returns the function handle psi and size of the expanded basis npsi

    x = sdpvar(nx,1);       % YALMIP sdp variable
    v = monolist(x,nd);     % create sdpvar with all the monomials
    m = sdisplay(v);        % transforms sdpvar to cell array
    b = m(2:end)';          % exclude constant term
    b_joined = join(b,';'); % joins all elements of b with delimiter
    basis_str = strcat(['[',b_joined{1},';ones(1,size(x,2))]']); % creates basis vector string
    
    pattern = '\((\d+)\)';
    basis_str = regexprep(basis_str, pattern, '($1,:)');
    basis_str = strrep(basis_str, '^', '.^');
    basis_str = strrep(basis_str, '*', '.*');
    
    % create function file
    if(toFile)
        fid = fopen('../functions/poly_basis.m','wt');
        fprintf(fid, 'function z = basis(x)\n\t');
        fprintf(fid, strcat(['z = ', basis_str,';\n']));
        fprintf(fid, 'end');
        fclose(fid);
    end

    func_str = strcat(['@(x)', basis_str]); % create function string
    psi = str2func(func_str);               % create function handle
    
    npsi = size(psi(rand(nx,1)),1);
end