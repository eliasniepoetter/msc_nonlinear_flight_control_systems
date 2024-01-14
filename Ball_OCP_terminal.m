function [u, flag] = Ball_OCP_terminal(x1_current, x2_current, x3_current, x4_current, P,alpha)

global m k g

opti = casadi.Opti();       % Opti Stack object

% parameters
x0 = opti.parameter(4,1);   % initial state values
%N = 500;                     % number of control intervals for shooting
T = 5;                      % Time Horizon
dt = 0.01;                  % Time steps
N = T/dt;

% decision variables
X = opti.variable(4,N+1);   % states
X1 = X(1,:);                % state x1
X2 = X(2,:);                % state x2
X3 = X(3,:);                % state x3
X4 = X(4,:);                % state x4
U = opti.variable(1,N);     % control trajectory
%T = opti.variable();        % time

% objective function
wx1 = 2;                    % weight for x1
wx2 = 1;                    % weight for x2
wx3 = 1;                    % weight for x3
wx4 = 0.1;                    % weight for x4
wu = 1;                     % weight for u

terminalCosts = X(:,N+1)'*P*X(:,N+1);
opti.minimize(terminalCosts + wx1*(X1*X1') + wx2*(X2*X2') + wx4*(X4*X4')+ wx3*(X3*X3') + wu*(U*U'));


% dynamic constraints
%dt = T/N;                   % time step for integration (for optimization)
x1 = casadi.MX.sym('x1');   % helper variables for dynamic constraint
x2 = casadi.MX.sym('x2');
x3 = casadi.MX.sym('x3');
x4 = casadi.MX.sym('x4');
u = casadi.MX.sym('u');


mu = 0.5;
f = casadi.Function('f',{x1, x2, x3, x4, u}, {[x2;(m*x1*x4^2 - m*g*sin(x3))/(k);x4;u]}, {'x1' 'x2' 'x3' 'x4' 'u'}, {'xdot'});
X_next = X(:,1:N) + dt*f(X1(:,1:N),X2(:,1:N),X3(:,1:N),X4(:,1:N),U);

    % k1_X = dt * f(X1(:,1:N),X2(:,1:N),U);
    % 
    % k2_X = dt * f(X1(:,1:N) + 0.5*k1_X(1), X2(:,1:N) + 0.5*k1_X(2), U);
    % 
    % k3_X = dt * f(X1(:,1:N) + 0.5*k2_X(1), X2(:,1:N) + 0.5*k2_X(2), U);
    % 
    % k4_X = dt * f(X1(:,1:N) + k3_X(1), X2(:,1:N) + k3_X(2), U);
    % 
    % X_next = X(:,1:N) + (1/6)*(k1_X + 2*k2_X + 2*k3_X + k4_X);

opti.subject_to(X(:,2:N+1) == X_next);


% additional constraints
% opti.subject_to(X1(N+1)==0);                    % stabilizing terminal constraint 
% opti.subject_to(X2(N+1)==0); % 

opti.subject_to(X(:,N+1)'*P*X(:,N+1) < alpha);  % terminal set              

opti.subject_to(X1(1)==x0(1));                  % initial state
opti.subject_to(X2(1)==x0(2)); 
opti.subject_to(X3(1)==x0(3));
opti.subject_to(X4(1)==x0(4));
opti.subject_to(X1<10);                         % state constraints 
opti.subject_to(X1>-10);                        %
opti.subject_to(X2<20);                         %
opti.subject_to(X2>-20);                        %
opti.subject_to(X3<pi/2);                         %
opti.subject_to(X3>-pi/2);
opti.subject_to(X4<10);                         %
opti.subject_to(X4>-10);
opti.subject_to(U<10);                          % control constraints
opti.subject_to(U>-10);                         %
%opti.subject_to(T>=0);                          % positive Time only


% solve the problem
solver_options = struct;
solver_options.ipopt.print_level = 0;
solver_options.print_time = 0;
solver_options.verbose = 0;
solver_options.ipopt.sb ='yes';
solver_options.ipopt.tol = 1e-6;  % Tolerance for optimality conditions
solver_options.ipopt.constr_viol_tol = 1e-4;  % Tolerance for constraint violations


opti.solver('ipopt', solver_options);           % select solver
opti.set_value(x0,[x1_current;x2_current;x3_current;x4_current]);     % set initial states
try
    sol = opti.solve();
    utrajectory = sol.value(U);
    u = utrajectory(1);
    flag = 0;
    
catch
    u = nan;
    flag = 1;
end

end

