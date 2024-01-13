close all;
clear;
clc;

%%
% get terminal constraints
loadConstraints = false;
saveConstraints = true;

if not(loadConstraints)
    [K,P,alpha] = getTerminalConstraintsBallBeam(saveConstraints);
else
    load('terminalConstraintsBallBeam.mat');
end



%%
% dynamics of the Van der Pol Oscillator
global m g k
m = 1;
g = 9.81;
k = 5;

x1dot = @(x2) x2;
x2dot = @(x1, x3, x4) k^-1 * (m*x1*x4^2 - m*g*sin(x3));
x3dot = @(x4) x4;
x4dot = @(u) u;


% initial conditions and initialization
t0 = 0;
tend = 30;
tstep = 0.01;
time = t0:tstep:tend;


% initial states
x1(1) = 0.1;
x2(1) = 0.1;
x3(1) = deg2rad(0.5);
x4(1) = 0.01;

x1ol(1) = x1(1);
x2ol(1) = x2(1);
x3ol(1) = x3(1);
x4ol(1) = x4(1);



%% Run NMPC
% output to console
fprintf('---------------------------------------------------\n');
fprintf('-- NMPC Ietartion for the Van der Pol Oscillator --\n');
fprintf('---------------------------------------------------\n\n');

% closed loop dynamics simulation
for i = 1 : length(time)-1
    
    u(i) = -K*[x1(i);x2(i);x3(i);x4(i)];

    % Euler-Cauchy Method for explicit solution of the IVP
    x1(i+1) = x1(i) + tstep*x1dot(x2(i));
    x2(i+1) = x2(i) + tstep*x2dot(x1(i), x3(i), x4(i));
    x3(i+1) = x3(i) + tstep*x3dot(x4(i));
    x4(i+1) = x4(i) + tstep*x4dot(u(i));

end


%% Postprocessing

figure;
hold on;
title('State Trajectory position');
plot(x1, x2);
xlabel('x1');
ylabel('x2');
grid on;
axis equal;
hold off;

figure;
hold on;
title('State Trajectory angle');
plot(x3, x4);
xlabel('x3');
ylabel('x4');
grid on;
axis equal;
hold off;


