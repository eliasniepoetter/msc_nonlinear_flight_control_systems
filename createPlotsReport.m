clear;
close all;
clc;


%% Beam & Ball dynamics

global m g k
m = 10;
g = 9.81;
k = 6;

x1dot = @(x2) x2;
% x2dot = @(x1, x3, x4) k^-1 * (m*x1*x4^2 - m*g*sin(x3));
x2dot = @(x1, x3, x4) k^-1 * (m*x1*x4^2 - m*g* ((x3^5)/120 - (x3^3)/6 + x3));
x3dot = @(x4) x4;
x4dot = @(u) u;



%% Generate terminal conditions

[K_pp,P_pp,alpha_pp] = getTerminalConstraintsBallBeam(false);
load('terminalConstraintsBallBeam.mat');


%% LQR vs. Pole plot

x = mpvar('x', [4 1]);

figure;
hold on;

Vx = x(1:2)'*P(1:2,1:2)*x(1:2);
[~,h] = pcontour(Vx, alpha,[-1.5 1.5 -1.5 1.5], '-', 1000);
h.EdgeColor = [0 0 0];
h.LineStyle = "-";
h.LineWidth = 1.125;

Vx_pp = x(1:2)'*P_pp(1:2,1:2)*x(1:2);
[~,h_pp] = pcontour(Vx_pp, alpha_pp,[-1.5 1.5 -1.5 1.5], '--', 1000);
h_pp.EdgeColor = [0 0 0];
h_pp.LineWidth = 1.125;

fs = 20;
xlabel('$x_1$','FontSize', fs, Interpreter='latex');
ylabel('$x_2$','FontSize', fs, Interpreter='latex');
legend('LQR', 'Pole Placement','FontSize', fs, Interpreter='latex');
grid on;
axis equal;
xlim([-1.5 1.5]);
ylim([-1.5 1.5]);
hold off;







%% Quiver plot

x = mpvar('x', [4 1]);
V = x'*P*x;

step_x = 0.01/4;
step_y = 0.03/4;

% Define the grid of points where you want to plot the vectors
[x11, x22] = meshgrid(-0.2:step_x:0.2, -0.6:step_y:0.6);
[x33, x44] = meshgrid(-0.2:step_x:0.2, -0.6:step_y:0.6);


v = zeros(size(x11));
w = zeros(size(x22));
q = zeros(size(x33));
r = zeros(size(x44));
r_closed = zeros(size(x44));
u_open = zeros(size(x11));

for i = 1:numel(x11)
    v(i) = x1dot(x22(i));
    w(i) = x2dot(x11(i), x33(i), x44(i));
    q(i) = x3dot(x44(i));
    r(i) = x4dot(u_open(i));
    r_closed(i) = x4dot(-K*[x11(i); x22(i); x33(i); x44(i)]);
end


% Normalize the vectors for better visualization
magnitude = sqrt(q.^2 + r.^2);
magnitude_closed = sqrt(q.^2 + r_closed.^2);

q_normalized = q ./ magnitude;
r_normalized = r ./ magnitude;
r_closed_normalized = r_closed ./ magnitude_closed;

figure;
hold on;
Vx = x(3:4)'*P(3:4,3:4)*x(3:4);
[C,h] = pcontour(Vx, alpha,[-1.5 1.5 -1.5 1.5], '-', 1000);
h.EdgeColor = [0 0 0];
quiver(x33, x44, q_normalized, r_closed_normalized, 'AutoScale', 'on', 'AutoScaleFactor', 1,Color=[0 0 0]);
xlabel('x3');
ylabel('x4');
grid on;
axis equal;
hold off;



