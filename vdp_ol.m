
close all;
clear;
clc;


% dynamics of the Van der Pol Oscillator
x1dot = @(x2) x2;
x2dot = @(x1, x2, u) 4*(1-x1^2)*x2 - x1 + u;


% initial conditions and initialization
t0 = 0;
tend = 100;
tstep = 0.01;
time = t0:tstep:tend;


% initial states
x1(1) = 4;
x2(1) = 4;


%% Run NMPC
% output to console
fprintf('------------------------------------------------\n');
fprintf('-- Open Loop Ietartion Van der Pol Oscillator --\n');
fprintf('------------------------------------------------\n\n');

% closed loop dynamics simulation
for i = 1 : length(time)-1
    x2(i+1) = x2(i) + tstep*x2dot(x1(i), x2(i), 0);
    x1(i+1) = x1(i) + tstep*x1dot(x2(i));
end


%% Postprocessing


figure;
hold on;
title('State Trajectory');
plot(x1, x2);
xlabel('x1');
ylabel('x2');
grid on;
axis equal;
hold off;


figure;
tiledlayout(2, 1);
nexttile;
hold on;
title('x1')
plot(time, x1);
xlabel('Time');
ylabel('x1');
grid on;
hold off;
    nexttile;
    hold on;
    title('x2')
    plot(time, x2);
    xlabel('Time');
    ylabel('x2');
    grid on;
    hold off;
