close all;
clear;
clc;

global g m Iyy ut Fxt Fzt Myt rho S cw


%% Step 1: Calculate linear controller K for linearized system

% parameters
g = 9.91;
m = 200;
Iyy = 5000;
ut = 200;
rho = 1;
S = 0.05;
cw = 0.7;
Fxt = 0.5*rho*ut^2*S*cw;
Fzt = 0;
Myt = 0;

% states
syms u w xe ze q theta
x = [u w xe ze q theta];

% inputs
syms Fx Fz My
u_in = [Fx Fz My];

% dynamics
udot = @(A_xe, q, w) A_xe - q*w;
wdot = @(A_ze, q, u) A_ze + q*u;

drag_x = @(u, w) 0.05*0.7*sqrt(u^2 + w^2)^2*1.225*0.5 * (u/(abs(u)+abs(w)));
drag_z = @(u, w) 0.05*0.7*sqrt(u^2 + w^2)^2*1.225*0.5 * (w/(abs(u)+abs(w)));

f_Axe = @(Fx, theta, u, w) Fx/m - 9.81*sin(theta) - drag_x(u, w)/m;
f_Aze = @(Fz, theta, u, w) Fz/m + 9.81*cos(theta) - drag_z(u, w)/m;

xedot = @(u, w, theta) u*cos(theta) + w*sin(theta);
zedot = @(u, w, theta) -u*sin(theta) + w*cos(theta);

qdot = @(My) My/Iyy;
thetadot = @(q) q;


% calculate jacobians
f_x_udot = jacobian(udot(f_Axe(Fx, theta, u, w), q, w), x);
f_x_wdot = jacobian(wdot(f_Aze(Fz, theta, u, w), q, u), x);
f_x_xedot = jacobian(xedot(u, w, theta), x);
f_x_zedot = jacobian(zedot(u, w, theta), x);
f_x_qdot = jacobian(qdot(My), x);
f_x_thetadot = jacobian(thetadot(q), x);

f_u_udot = jacobian(udot(f_Axe(Fx, theta, u, w), q, w), u_in);
f_u_wdot = jacobian(wdot(f_Aze(Fz, theta, u, w), q, u), u_in);
f_u_xedot = jacobian(xedot(u, w, theta), u_in);
f_u_zedot = jacobian(zedot(u, w, theta), u_in);
f_u_qdot = jacobian(qdot(My), u_in);
f_u_thetadot = jacobian(thetadot(q), u_in);

% evaluate jacobians
eps = 1e-9;
x0 = [ut eps eps eps eps eps];
u0 = [Fxt eps eps];

A = [double(subs(f_x_udot,[x, u_in],[x0, u0]));
    double(subs(f_x_wdot,[x, u_in],[x0, u0]));
    double(subs(f_x_xedot,[x, u_in],[x0, u0])); 
    double(subs(f_x_zedot,[x, u_in],[x0, u0]));
    double(subs(f_x_qdot,[x, u_in],[x0, u0]));
    double(subs(f_x_thetadot,[x, u_in],[x0, u0]));
    ];
A = round(A, 6);

B = [double(subs(f_u_udot,[x, u_in],[x0, u0]));
    double(subs(f_u_wdot,[x, u_in],[x0, u0]));
    double(subs(f_u_xedot,[x, u_in],[x0, u0])); 
    double(subs(f_u_zedot,[x, u_in],[x0, u0]));
    double(subs(f_u_qdot,[x, u_in],[x0, u0]));
    double(subs(f_u_thetadot,[x, u_in],[x0, u0]));
    ];
B = round(B, 6);

C = eye(6);
D = zeros(6, 3);
p = [-200; -200; 0; 0; -200; 0];
K = place(A,B,p);
Acl = A-B*K;

missile_lin = ss(A, B, C, D);
poles_missile_lin = pole(missile_lin);

missile_lin_placed = ss(Acl,B,C,D);
poles_missile_lin_placed = pole(missile_lin_placed);

% plot step responses
% figure;
% step(missile_lin);
% 
% figure;
% step(missile_lin_placed);


%% Step 2: Choose constant kappa

Ak = A+B*K;
Ak_eigen = eig(Ak);
kappa_interval = [0 -max(Ak_eigen)];
kappa = -195;

% MATLAB solver for Lyapunov equation (continous time):
    % AX + XA' + Q = 0
% Allgöwer notation:
    % X = P
    % A = (Ak + kappa*I)'
    % Q = Q_star
    % Q_star = Q + K'RK
    % Q = eye(6) -> weights for states
    % R = eye(3) -> weights for inputs

A_lyap = (Ak + kappa*eye(6))';
Q = eye(6)*1e-5;
R = eye(3)*1e-5;
Q_star = Q + K'*R*K;
P = lyap(A_lyap,Q_star);
chol(P)


% syms Ps;
% eqn = (Ak + kappa*eye(6))'*Ps + Ps*(Ak + kappa*eye(6)) == -(Q + K'*R*K);
% P = solve(eqn,Ps);


%% Step 3: Find largest aplha_1

pvar alpha1;
x = mpvar('x', [6,1]);
% pconst = -0.5*(x'*P*x) + 0.5*alpha1 >= 0;
pconst = -0.5*(x'*P*x) + 100 >= 0;

info = sosopt(pconst, x);
disp(['Feasibility of Optimization Problem: ' num2str(info.feas)]);












%% Step 4: Find largest alpha















