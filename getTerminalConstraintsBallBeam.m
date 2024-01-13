function [K,P,alpha] = getTerminalConstraintsBallBeam(saveConstraints)
%% calculate terminal constraints


%% Step 1: Calculate linear controller K for linearized system

% states
syms x1 x2 x3 x4
x = [x1 x2 x3 x4];

% inputs
syms u

% parameters
global m g k
m = 1;
g = 9.81;
k = 5;

% dynamics
x1dot = @(x2) x2;
x2dot = @(x1, x3, x4) k^-1 * (m*x1*x4^2 - m*g*sin(x3));
x3dot = @(x4) x4;
x4dot = @(u) u;

% calculate jacobians
fx_x1dot = jacobian(x1dot(x2),x);
fx_x2dot = jacobian(x2dot(x1,x3,x4),x);
fx_x3dot = jacobian(x3dot(x4),x);
fx_x4dot = jacobian(x4dot(u),x);

fu_x1dot = jacobian(x1dot(x2),u);
fu_x2dot = jacobian(x2dot(x1,x3,x4),u);
fu_x3dot = jacobian(x3dot(x4),u);
fu_x4dot = jacobian(x4dot(u),u);


% evaluate jacobians
x0 = [0 0 0 0];
u0 = 0;

A = [double(subs(fx_x1dot,[x,u],[x0,u0]));
    double(subs(fx_x2dot,[x,u],[x0,u0]));
    double(subs(fx_x3dot,[x,u],[x0,u0]));
    double(subs(fx_x4dot,[x,u],[x0,u0]));
    ];

B = [double(subs(fu_x1dot,[x,u],[x0,u0]));
    double(subs(fu_x2dot,[x,u],[x0,u0]));
    double(subs(fu_x3dot,[x,u],[x0,u0]));
    double(subs(fu_x4dot,[x,u],[x0,u0]));
    ];

C = eye(4);
D = zeros(4, 1);

% place poles
p = [-5,-10,-6,-12];
K = place(A,B,p);
Acl = A-B*K;

vdp_linearized = ss(A, B, C, D);
vdp_poles = pole(vdp_linearized);

vdp_linearized_placed = ss(Acl,B,C,D);
vdp_poles_placed = pole(vdp_linearized_placed);

% plot step responses for controlled and uncontrolled system
figure;
step(vdp_linearized);

figure;
step(vdp_linearized_placed);


%% Step 2: Choose constant kappa and caclulate matrix P

Ak = A-B*K;
Ak_eigen = eig(Ak);
kappa_interval = [0 -max(Ak_eigen)];
kappa = 4.9;

A_lyap = (Ak + kappa*eye(4))';
Q = 0.25*eye(4);
R = eye(1);
Q_star = Q + K'*R*K;
P = lyap(A_lyap,Q_star);
cholesky = chol(P);



%% Step 3: Find largest alpha1

% lyapunov function candidate
clear x x1 x2 x3 x4;
pvar x1 x2 x3 x4;
x = [x1;x2;x3;x4];
V = x'*P*x;

uLim = 100;
alpha1_min = 0;
alpha1_max = 1000;

% bisection procedure
figure;
hold on;
grid on;
xlabel('x1');
ylabel('x2');
while (alpha1_max-alpha1_min)>0.1
    alpha1 = (alpha1_min + alpha1_max)/2;
    p = x'*P*x - alpha1;
    [xin,xon] = psample(p,x,[0;0;0;0],10000);
    feasible = true;
    for i = 1 : length(xon)
        if abs(K*xon(:,i)) > uLim
            feasible = false;
            break;
        end
    end
    if feasible
        alpha1_min = alpha1;
    else
        alpha1_max = alpha1;
    end
    Vx = x(1:2)'*P(1:2,1:2)*x(1:2);
    pcontour(Vx, alpha1, alpha1*[-1 1 -1 1]);
    pause(0.1)
    alpha1
end


%% Step 4: Find largest aplha

clear pconstr x u

% approx sin through Taylor
% syms xd
% sinTaylor = taylor(sin(xd));


% system dynamics
x = mpvar('x', [4 1]);
f = [x(2);
    % k^-1 * (m*x(1)*x(4)^2 - m*g* (x(3)^5)/120 - (x(3)^3)/6 + x(3));
    k^-1 * (m*x(1)*x(4)^2 - m*g* -(x(3)^3)/6 + x(3));
    x(4);
    -K*x
    ];

V = x'*P*x;

% sos multiplier
z = monomials(x, 0:4);
s = sosdecvar('s',z);

% gamma sublevel set
pvar gama;
alpha_min = 0;
alpha_max = alpha1;
epsilon = 10^-6; % small positive number

% sos constraints
pconstr(1) = V - epsilon*(x'*x) >= 0;
pconstr(2) = s*(V-gama) -jacobian(V,x)*f - epsilon*(x'*x) >= 0;
pconstr(3) = s>= 0;
% pconstr(4) = pi/2 >= x(3);
% pconstr(5) = -pi/2 <= x(3);

% bisection procedure
figure;
hold on;
grid on;
title('Region of Attraction Estimation');
xlabel('x1');
ylabel('x2');
while (alpha_max-alpha_min)>1e-3
    alpha = (alpha_min + alpha_max)/2;
    [info, dopt] = sosopt(subs(pconstr, gama, alpha), x);
        if info.feas == 1
            alpha_min = alpha;
        else
            alpha_max = alpha;
        end
    Vx = x(1:2)'*P(1:2,1:2)*x(1:2);
    pcontour(Vx, alpha);
    pause(0.075)
    disp(alpha);
end


%% save results
if saveConstraints
    save('terminalConstraintsBallBeam', 'K', 'P', 'alpha');
end


end

