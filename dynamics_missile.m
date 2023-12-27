close all;
clear;
clc;

% parameters
global m Iyy;
m = 200;
Iyy = 1000;

% dynamics
u_dot = @(A_xe, q, w, u) A_xe - q*w;
w_dot = @(A_ze, q, u) A_ze + q*u;

drag_x = @(u, w) 0.05*0.7*sqrt(u^2 + w^2)^2*1.225*0.5 * (u/(abs(u)+abs(w)));
drag_z = @(u, w) 0.05*0.7*sqrt(u^2 + w^2)^2*1.225*0.5 * (w/(abs(u)+abs(w)));

% f_A_xe = @(Fx, theta, u) Fx/m - 9.81*sin(theta) - (0.05*0.7*1.225*u^2*0.5)/m;
f_A_xe = @(Fx, theta, u, w) Fx/m - 9.81*sin(theta) - drag_x(u, w)/m;
f_A_ze = @(Fz, theta, u, w) Fz/m + 9.81*cos(theta) - drag_z(u, w)/m;

Xe_dot = @(u, w, theta) u*cos(theta) + w*sin(theta);
Ze_dot = @(u, w, theta) -u*sin(theta) + w*cos(theta);

q_dot = @(My) My/Iyy;
theta_dot = @(q) q;


%% Example maneuvre

t0 = 0;
tf = 60;
tstep = 0.01;
time = t0:tstep:tf;

u0 = 300;
w0 = 0;
theta0 = 0;
q0 = 0;

u  = zeros(length(time), 1);
u(1) = u0;

w  = zeros(length(time), 1);
w(1) = w0;

theta  = zeros(length(time), 1);
theta(1) = theta0;

q  = zeros(length(time), 1);
q(1) = q0;

A_xe = zeros(length(time), 1);
A_xe(1) = 0;

A_ze = zeros(length(time), 1);
A_ze(1) = 0;

Xe = zeros(length(time), 1);
Xe(1) = 0;

Ze = zeros(length(time), 1);
Ze(1) = 0;


for i = 1 : length(time)-1
    % declare control input
    if i > 100 && i < 200
        My = 10;
    elseif i > 200 && i < 400
        My = -10;
    elseif i > 400 && i < 500
        My = 10;
    else
        My = 0;
    end
    % My = 0;
    Fx = 50000;
    Fz = -1.0*(m*9.81);

    % Euler-Cauchy Update for states
    q(i+1) = q(i) + tstep*q_dot(My);
    theta(i+1) = theta(i) + tstep*theta_dot(q(i));
    A_xe(i+1) = f_A_xe(Fx, theta(i), u(i), w(i));
    A_ze(i+1) = f_A_ze(Fz, theta(i), u(i), w(i));
    w(i+1) = w(i) + tstep*w_dot(A_ze(i), q(i), u(i));
    u(i+1) = u(i) + tstep*u_dot(A_xe(i), q(i), w(i), u(i));
    Xe(i+1) = Xe(i) + tstep*Xe_dot(u(i), w(i), theta(i));
    Ze(i+1) = Ze(i) + tstep*Ze_dot(u(i), w(i), theta(i));
end

% Postprocessing
Ze = -Ze; % transform to commen human sense coordinate


%% Analysis

figure;
hold on;
title('xz-plane trajectory');
scatter(Xe, Ze, 10, 'o', 'filled');
xlabel('Xe');
ylabel('Ze');
grid on;
hold off;


figure;
tiledlayout(4, 2);
nexttile;
hold on;
title('u');
scatter(time, u, 3, 'o', 'filled');
xlabel('time');
ylabel('u');
grid on;
hold off;

nexttile;
hold on;
title('w');
scatter(time, w, 3, 'o', 'filled')
xlabel('time');
ylabel('w');
grid on;
hold off;

nexttile;
hold on;
title('q');
scatter(time, q, 3, 'o', 'filled')
xlabel('time');
ylabel('q');
grid on;
hold off;

nexttile;
hold on;
title('theta');
scatter(time, theta, 3, 'o', 'filled')
xlabel('time');
ylabel('theta');
grid on;
hold off;

nexttile;
hold on;
title('Aze');
scatter(time, A_ze, 3, 'o', 'filled')
xlabel('time');
ylabel('Aze');
grid on;
hold off;

nexttile;
hold on;
title('Axe');
scatter(time, A_xe, 3, 'o', 'filled')
xlabel('time');
ylabel('Axe');
grid on;
hold off;

nexttile;
hold on;
title('Ze');
scatter(time, Ze, 3, 'o', 'filled')
xlabel('time');
ylabel('Ze');
grid on;
hold off;

nexttile;
hold on;
title('Xe');
scatter(time, Xe, 3, 'o', 'filled')
xlabel('time');
ylabel('Xe');
grid on;
hold off;


