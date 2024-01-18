function [u, flag, utrajectory, xtrajectory, time] = ocp_van_der_pol(x1_current, x2_current, P, alpha)

opti = casadi.Opti();       % Opti Stack object

% parameters
x0 = opti.parameter(2,1);   % initial state values
T = 10;
dt = 0.1;
N = T/dt;                     % number of control intervals for shooting
muvdp = 0.5;

% decision variables
X = opti.variable(2,N+1);   % states
X1 = X(1,:);                % state x1
X2 = X(2,:);                % state x2
U = opti.variable(1,N);     % control trajectory
% T = opti.variable();        % time




% objective function
wx1 = 1;                    % weight for x1
wx2 = 1;                    % weight for x2
wu = 1;                     % weight for u
wT = 1;                     % weight for terminal cost

% terminalCosts = X1(N+1) + X2(N+1);

terminalCosts = [X1(N+1);X2(N+1)]'*P*[X1(N+1);X2(N+1)];

J = wT*terminalCosts + wx1*(X1*X1') + wx2*(X2*X2') + wu*(U*U');
opti.minimize(J);


% dynamic constraints
x1 = casadi.MX.sym('x1');   % helper variables for dynamic constraint
x2 = casadi.MX.sym('x2');
u = casadi.MX.sym('u');

f = casadi.Function('f',{x1, x2, u}, {[x2;muvdp*(1-x1^2)*x2 - x1 + u]}, {'x1' 'x2' 'u'}, {'xdot'});
X_next = X(:,1:N) + dt*f(X1(:,1:N),X2(:,1:N),U);
opti.subject_to(X(:,2:N+1) == X_next);


% additional constraints
% opti.subject_to(X1(N+1)==0);                    % stabilizing terminal constraint 
% opti.subject_to(X2(N+1)==0);                    % 
opti.subject_to(X(:,N+1)'*P*X(:,N+1) < alpha);  % terminal set    
opti.subject_to(X1(1)==x0(1));                  % initial state
opti.subject_to(X2(1)==x0(2));                  % 
opti.subject_to(X1<10);                         % state constraints 
opti.subject_to(X1>-10);                        %
opti.subject_to(X2<10);                         %
opti.subject_to(X2>-10);                        %
opti.subject_to(U<100);                         % control constraints
opti.subject_to(U>-100);
% opti.subject_to(T>=0);                          % positive Time only


% solve the problem
solver_options = struct;
solver_options.ipopt.print_level = 0;
solver_options.print_time = 0;
solver_options.verbose = 0;
solver_options.ipopt.sb ='yes';

opti.solver('ipopt', solver_options);           % select solver
opti.set_value(x0,[x1_current;x2_current]);     % set initial states
try
    sol = opti.solve();
    utrajectory = sol.value(U);
    xtrajectory = sol.value(X);
    time = 0:dt:T;
    u = utrajectory(1);
    flag = 0;
catch
    u = nan;
    utrajectory = zeros(N,1);
    xtrajectory = zeros(N,1);
    time = 0:dt:T;
    flag = 1;
end

end

