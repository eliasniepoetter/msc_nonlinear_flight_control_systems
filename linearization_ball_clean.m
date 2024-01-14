% linearization clean

clear

addpath('C:\Users\nicoh\Documents\Stuttgart\3.Semester\systemtheoretsich\Project\msc_nonlinear_flight_control_systems-main\sosopt')
addpath('C:\Users\nicoh\Documents\Stuttgart\3.Semester\systemtheoretsich\Project\msc_nonlinear_flight_control_systems-main\multipoly')


%% K
syms x1 x2 x3 x4 u

m = 10;
g = 9.81;
k = 6;


x = [x1;x2;x3;x4];

f = [x2
    (m*x1*x4^2 - m*g*sin(x3))/(k)
    x4
    u];

A = jacobian(f,x);
B = jacobian(f,u);
C = [1 0 0 0];
D = 0;

Anew = double(subs(A,[x;u],[0;0;0;0;0]));
Bnew = double(subs(B,[x;u],[0;0;0;0;0]));

Q = [2 0 0 0
    0 1 0 0 
    0 0 1 0
    0 0 0 0.1];
R = 1;

[K,~,~] = lqr(Anew,Bnew,Q,R);

%% P

max_eigen = max(eig(Anew - Bnew*K));

kappa = -real(max_eigen) - 0.02;

P = sym('P', [4 4]);

% Lyapunov eq
eq = ((Anew - Bnew*K) + kappa*eye(4))'*P + P * ((Anew - Bnew*K) + kappa*eye(4)) == -(Q + K'*R*K);
p = solve(eq,P);

% Initialize the matrix P
P = zeros(4, 4);

% Populate the matrix P using the struct 'p'
for i = 1:4
    for j = 1:4
        field_name = sprintf('P%d_%d', i, j);
        P(i, j) = p.(field_name);
    end
end

%% set containment x'Px <= alpha contained in K*x < ulimit
% generalised s procedure

x = mpvar('x', [4 1]);

ulimit = 10;

% sos multiplier
z = monomials(x, 0:3);
s = sosdecvar('s',z);

pvar a;
alpha1_min = 0;
alpha1_max = 100000;
epsilon = 10^-6;

pconstr(1) = - s*(-x'*P*x + a) + (-(K*x)^2 + ulimit^2) - epsilon*(x'*x) >= 0;
pconstr(2) = s>=0;

% tic
% [info,dopt,sossol] = gsosopt(pconstr,x,a);
%  toc
tic
while (alpha1_max-alpha1_min)>0.001
    alpha1 = (alpha1_min + alpha1_max)/2;
   
    [info, dopt] = sosopt(subs(pconstr,a, alpha1),x);
    if info.feas == 1
        alpha1_min = alpha1;
    else
        alpha1_max = alpha1;
    end
    % figure(2)
    % hold on
    % pcontour(x'*P*x, alpha1, 100*[-1 1 -1 1])
    % alpha1_max-alpha1_min
end
toc
alpha1_sos = double(alpha1_min);

% %% alpha1
% 
% % system dynamics
% x = mpvar('x', [4 1]);
% 
% ulimit = 10;
% 
% % Number of random points
% num_points = 1000;
% num_alphas = 10000;
% alpha1_max = 100000;
% alpha1_min = 0.01;
% 
% alpha1 = linspace(0.001,alpha1_max,num_alphas);
% 
% while (alpha1_max-alpha1_min)>0.01
%     alpha1 = (alpha1_min + alpha1_max)/2;
%     p = x'*P*x - alpha1;
%     [~,xon] = psample(p,x,[0;0],num_points);
% 
%     for i = 1:num_points
%         flag(i) = abs(K * (xon(:, i))) >= ulimit;
%     end
%     if sum(flag) > 0.5
%         alpha1_max = alpha1;
%     else
%         alpha1_min = alpha1;
%     end
%     % pcontour(x'*P*x, alpha1, 100*[-1 1 -1 1],'r')
% end

%% alpha

clear x
% system dynamics
x = mpvar('x', [4 1]);
f = [x(2)
    (m*x(1)*x(4)^2 - m*g*(x(3)^5/120 - x(3)^3/6 + x(3)))/(k)
    x(4)
    -K*x];
V = x'*P*x;

% sos multiplier
z = monomials(x, 0:3);
s = sosdecvar('s',z);

% gamma sublevel set
pvar alpha_sos;
alpha_min = 0; alpha_max = alpha1_sos;
epsilon = 10^-6; % small positive number

% sos constraints
pconstr(1) = V - epsilon*(x'*x) >= 0;
pconstr(2) = s*(V-alpha_sos) -jacobian(V,x)*f - epsilon*(x'*x) >= 0;
pconstr(3) = s>= 0;
%pconstr(4) = x(3) <= pi/2;
%pconstr(5) = x(3) >= -pi/2;

% bisection procedure
while (alpha_max-alpha_min)>0.0001
    alpha = (alpha_min + alpha_max)/2;
    [info, dopt] = sosopt(subs(pconstr, alpha_sos, alpha), x);
    if info.feas ==1
        alpha_min = alpha;
    else
        alpha_max = alpha;
    end
    % figure(1)
    % pcontour(V, alpha1_sos, 10*[-1 1 -1 1],'r')
    % hold on
    % pcontour(V, alpha, 10*[-1 1 -1 1])
    pcontour(x(1:2)'*P(1:2,1:2)*x(1:2), alpha_min, 2*[-1 1 -1 1])
    alpha_min
end
