clear all; clc; close all;
load ../data/dspace.mat
exitflag = dspace.Y(1).Data;
solve_time_ms = dspace.Y(2).Data;
figure
subplot(2,1,1)
plot(exitflag(1:4085))
ylabel('Solver exitflag')
subplot(2,1,2)
plot(solve_time_ms(1:4085))
ylabel('Solve time [ms]')
xlabel('Sample')