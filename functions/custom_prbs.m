function y = custom_prbs( N, dutyCycle )
% N: number of samples
% dutyCycle: real number betweeb 0 and 1
transpose_output = 0;
if(numel(N) == 1)
    N = [N 1];
    transpose_output = 1;
end
y = rand(N);
y = double(y > 1 - dutyCycle);
if transpose_output
    y = y';
end
end

