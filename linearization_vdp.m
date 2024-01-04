% linearization

clear

syms x1 x2 u

x = [x1;x2];

f = [x2
    0.5*(1-x1^2)*x2 - x1 + u];

A = jacobian(f,x);
B = jacobian(f,u);
C = eye(2);
D = 0;

Anew = double(subs(A,x,[0;0]));
Bnew = double(subs(B,u,0));

Q = [1 0
    0 1];
R = 1;

[K,~,~] = lqr(Anew,Bnew,Q,R);

% sys1 = ss(Anew-Bnew*K,Bnew,C,D);
% step(sys1)

max_eigen = max(eig(Anew - Bnew*K));

kappa = -real(max_eigen) - 0.02;

P = sym('P', [2 2]);

% Lyapunov eq
eq = ((Anew - Bnew*K) + kappa*eye(2))'*P + P * ((Anew - Bnew*K) + kappa*eye(2)) == -(Q + K'*R*K);
p = solve(eq,P);

% Initialize the matrix P
P = zeros(2, 2);

% Populate the matrix P using the struct 'p'
for i = 1:2
    for j = 1:2
        field_name = sprintf('P%d_%d', i, j);
        P(i, j) = p.(field_name);
    end
end


u_limit = 10;

% Cholesky decomposition of P
[V, D] = eig(P);

% Generate 100 values equally distributed along the border
num_points = 100;
theta_values = linspace(0, 2*pi, num_points);
border_points = zeros(2, num_points);

% Objective function for maximization
objective = @(alpha1) -alpha1;

% Initial guess for alpha1
initial_guess = 1.0;

% Constraints
nonlinear_constraint = @(alpha1) deal([], max(abs(K * (sqrt(alpha1) * V * sqrt(D) * [cos(theta_values); sin(theta_values)]))) - u_limit);

% Options for fmincon
options = optimoptions('fmincon', 'Display', 'iter', 'Algorithm', 'interior-point');

% Optimization using fmincon
result = fmincon(objective, initial_guess, [], [], [], [], 0, [], nonlinear_constraint, options);


%% optional plotting
% for i = 1:num_points
%     % Parametric equation of the set
%     x_theta = sqrt(result) * V * sqrt(D) * [cos(theta_values(i)); sin(theta_values(i))];
%     border_points(:, i) = x_theta;
%     %K*x_theta
% end

% % Plot the set
% figure;
% plot(border_points(1, :), border_points(2, :), 'b.');
% title('Set: x''Px < alpha');
% xlabel('x1');
% ylabel('x2');
% axis equal;
% grid on;
%% 

alpha = 2;

clear x

% Dynamics function
fc = @(x) [x(2); 0.5*(1-x(1)^2)*x(2) - x(1) + K*x];
phi = @(x) (fc(x) - (Anew - Bnew * K) * x);

% Objective function for minimization
objective = @(x) -x' * P * phi(x) + kappa * x' * P * x;

% Nonlinear constraint function
nonlinear_constraint = @(x) deal([], x' * P * x - alpha);

% Options for fmincon
options = optimoptions('fmincon', 'Algorithm', 'interior-point');

% Initial guess for x
x0 = ones(2, 1);

% Optimization using fmincon
result2 = fmincon(objective, x0, [], [], [], [], [], [], nonlinear_constraint, options);

% Display the result
disp(['Optimal x: ', num2str(result2')])
disp(['Objective value: ', num2str(-objective(result2))])
