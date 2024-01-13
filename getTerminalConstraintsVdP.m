function [K,P,alpha] = getTerminalConstraintsVdP(saveConstraints)
%% calculate terminal constraints


%% Step 1: Calculate linear controller K for linearized system

% states
syms x1 x2
x = [x1 x2];

% inputs
syms u

% dynamics
x1dot = @(x2) x2;
x2dot = @(x1,x2,u) 0.5*(1-x1^2)*x2 - x1 + u;

% calculate jacobians
fx_x1dot = jacobian(x1dot(x2), x);
fx_x2dot = jacobian(x2dot(x1,x2,u), x);
fu_x1dot = jacobian(x1dot(x2), u);
fu_x2dot = jacobian(x2dot(x1,x2,u), u);

% evaluate jacobians
x0 = [0 0];
u0 = 0;

A = [double(subs(fx_x1dot,[x,u],[x0,u0]));
    double(subs(fx_x2dot,[x,u],[x0,u0]));
    ];

B = [double(subs(fu_x1dot,[x,u],[x0,u0]));
    double(subs(fu_x2dot,[x,u],[x0,u0]));
    ];

C = eye(2);
D = zeros(2, 1);

% place poles
p = [-10,-1];
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
kappa = 0.9;

A_lyap = (Ak + kappa*eye(2))';
Q = eye(2);
R = eye(1);
Q_star = Q + K'*R*K;
P = lyap(A_lyap,Q_star);
cholesky = chol(P);



%% Step 3: Find largest alpha1

% lyapunov function candidate
clear x x1 x2;
% x = mpvar('x', [2 1]);
pvar x1 x2;
x = [x1;x2];
V = x'*P*x;

uLim = 100;
alpha1_min = 0;
alpha1_max = 10000;

% bisection procedure
figure;
hold on;
grid on;
xlabel('x1');
ylabel('x2');
while (alpha1_max-alpha1_min)>0.001
    alpha1 = (alpha1_min + alpha1_max)/2;
    p = x'*P*x - alpha1;
    [xin,xon] = psample(p,x,[0;0],10000);
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
    pcontour(V, alpha1, 40*[-1 1 -1 1]);
    pause(0.1)
end


%% Step 4: Find largest aplha

clear pconstr x u

% system dynamics
x = mpvar('x', [2 1]);
f = [x(2); 0.5*(1-x(1)^2)*x(2) - x(1) - K*x];

V = x'*P*x;

% sos multiplier
z = monomials(x, 0:4);
s = sosdecvar('s',z);

% gamma sublevel set
pvar g;
alpha_min = 0;
alpha_max = alpha1;
epsilon = 10^-6; % small positive number

% sos constraints
pconstr(1) = V - epsilon*(x'*x) >= 0;
pconstr(2) = s*(V-g) -jacobian(V,x)*f - epsilon*(x'*x) >= 0;
pconstr(3) = s>= 0;
% pconstr(4) = K*x <= uLim;

% bisection procedure
figure;
hold on;
grid on;
title('Region of Attraction Estimation');
xlabel('x1');
ylabel('x2');
while (alpha_max-alpha_min)>0.001
    alpha = (alpha_min + alpha_max)/2;
    [info, dopt] = sosopt(subs(pconstr, g, alpha), x);
        if info.feas ==1
            alpha_min = alpha;

        else
            alpha_max = alpha;
        end
    pcontour(V, alpha, 30*[-1 1 -1 1]);
    pause(0.075)
end


%% save results
if saveConstraints
    save('terminalConstraintsVdP', 'K', 'P', 'alpha');
end


end

