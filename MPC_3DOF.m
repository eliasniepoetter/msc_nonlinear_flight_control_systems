close all;
clear;
clc;

% parameters
global m Iyy;
m = 200;
Iyy = 1000;

% dynamics of the 3-DOF point mass
u_dot = @(A_xe, q, w) A_xe - q*w;
w_dot = @(A_ze, q, u) A_ze + q*u;

f_A_xe = @(Fx, theta, u) Fx/m - 9.81*sin(theta) - 0.05*0.7*1.225*u^2*0.5/200;
f_A_ze = @(Fz, theta) Fz/m + 9.81*cos(theta);

Xe_dot = @(u, w, theta) u*cos(theta) + w*sin(theta);
Ze_dot = @(u, w, theta) -u*sin(theta) + w*cos(theta);

q_dot = @(My) My/Iyy;
theta_dot = @(q) q;


t0 = 0;
tf = 5;
tstep = 0.1;
time = t0:tstep:tf;


Fx = ones(length(time)-1, 1)*50000;
Fz = zeros(length(time)-1, 1);
My = zeros(length(time)-1, 1);

u0 = 300;
w0 = 0;
theta0 = 0;

theta  = zeros(length(time), 1);
theta(1) = theta0;

u  = zeros(length(time), 1);
u(1) = u0;



q0 = 0;
axe0 = 0;
% axe0 = Fx(1)/m - 9.81*sin(theta(1)) - 0.05*0.7*1.225*u(1)^2*0.5;
aze0 = 0;
xe0 = 0;
ze0 = 0;



w  = zeros(length(time), 1);
w(1) = w0;

q  = zeros(length(time), 1);
q(1) = q0;

A_xe = zeros(length(time), 1);
A_xe(1) = axe0;

A_ze = zeros(length(time), 1);
A_ze(1) = aze0;

Xe = zeros(length(time), 1);
Xe(1) = xe0;

Ze = zeros(length(time), 1);
Ze(1) = ze0;





% output to console
fprintf('-------------------------------------------------\n');
fprintf('-- NMPC Ietartion for a missile-target problem --\n');
fprintf('-------------------------------------------------\n\n');

for i = 1 : length(time)-1
    % calculate control input
    X = [u(i); w(i); A_xe(i); A_ze(i); Xe(i); Ze(i); q(i); theta(i)];
    Tx = 1000;
    Tz = -200;
    % [Fx(i), Fz(i), My(i), flags(i)] = OCP_3DOF(m, Iyy, X, [Tx; Tz]);
    [Fz(i), My(i), flags(i)] = OCP_3DOF(m, Iyy, X, [Tx; Tz]);
    
    % if isnan(Fx(i))
    %     Fx(i) = 0;
    % end

    if isnan(Fz(i))
        Fz(i) = 0;
    end

    if isnan(My(i))
        My(i) = 0;
    end

    % Euler-Cauchy Update for states
    q(i+1) = q(i) + tstep*q_dot(My(i));
    theta(i+1) = theta(i) + tstep*theta_dot(q(i));
    A_xe(i+1) = f_A_xe(Fx(i), theta(i), u(i));
    A_ze(i+1) = f_A_ze(Fz(i), theta(i));
    w(i+1) = w(i) + tstep*w_dot(A_ze(i), q(i), u(i));
    u(i+1) = u(i) + tstep*u_dot(A_xe(i), q(i), w(i));
    Xe(i+1) = Xe(i) + tstep*Xe_dot(u(i), w(i), theta(i));
    Ze(i+1) = Ze(i) + tstep*Ze_dot(u(i), w(i), theta(i));

    if mod(i/(length(time)-1)*100, 10) == 0
        done = i/(length(time)-1)*100;
        disp([num2str(done),'% done']);
    end

    if sqrt((Xe(i+1)-Tx)^2 + (Ze(i+1)-Tz)^2) <= 100
        disp('Hit!');
        break;
    elseif sqrt((Xe(i+1)-Tx)^2 + (Ze(i+1)-Tz)^2) >= 1.5*(sqrt((xe0-Tx)^2 + (ze0-Tz)^2))
        disp('Miss');
        break;
    end

end

% Postprocessing
Ze = -Ze; % transform to commen human sense coordinate


%% Analysis

figure;
hold on;
title('xz-plane trajectory');
scatter(Xe, Ze, 10, 'o', 'filled');
scatter(Tx, -Tz, 30, 'o', 'filled');
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


figure;
tiledlayout(3, 1);
nexttile;
hold on;
title('Fx');
scatter(time(1:end-1), Fx, 3, 'o', 'filled');
xlabel('time');
ylabel('Fx');
grid on;
hold off;

nexttile;
hold on;
title('Fz');
scatter(time(1:end-1), Fz, 3, 'o', 'filled');
xlabel('time');
ylabel('Fz');
grid on;
hold off;

nexttile;
hold on;
title('My');
scatter(time(1:end-1), My, 3, 'o', 'filled');
xlabel('time');
ylabel('My');
grid on;
hold off;









