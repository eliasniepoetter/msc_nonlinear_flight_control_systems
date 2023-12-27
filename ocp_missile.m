function [control_fz, control_my, flag] = ocp_missile(mass, inertia, X, target_position)

opti = casadi.Opti();       % Opti Stack object

% parameters
N = 10;                     % number of control intervals for shooting
m = mass;
Iyy = inertia;


% decision variables
% states: 
    U = opti.variable(1,N+1);       % x-component of velocity
    W = opti.variable(1,N+1);       % z-component of velocity
    Axe = opti.variable(1,N+1);     % acceleration in x (earth-fixed)
    Aze = opti.variable(1,N+1);     % acceleration in z (earth-fixed)
    Xe = opti.variable(1,N+1);      % earth-fixed x
    Ze = opti.variable(1,N+1);      % earth-fixed z
    Theta = opti.variable(1,N+1);   % pitch angle
    Q = opti.variable(1,N+1);       % pitch rate
    T = opti.variable();            % time
% control input:
    % Fx = opti.variable(1,N);        % force in x
    Fz = opti.variable(1,N);        % force in z
    My = opti.variable(1,N);        % moment around y


% objective function

opti.minimize(T);
% temrinal_dist = sqrt((Xe(N+1)-target_position(1))^2 + (Ze(N+1)-target_position(2))^2);
% opti.minimize(temrinal_dist);


% dynamic constraints
% define helper variables for CasADi Functions
    dt = T/N;
    u = casadi.MX.sym('u');
    w = casadi.MX.sym('w');
    axe = casadi.MX.sym('axe');
    aze = casadi.MX.sym('aze');
    theta = casadi.MX.sym('theta');
    q = casadi.MX.sym('q');

    fx = casadi.MX.sym('fx');
    fz = casadi.MX.sym('fz');
    my = casadi.MX.sym('my');

% definition of state updates (Euler-Cauchy)
    u_dot = casadi.Function('u_dot', {axe, q, w}, {axe - q*w}, {'axe' 'q' 'w'}, {'u_dot'});
    U_next = U(1:N) + dt*u_dot(Axe(1:N),Q(1:N),W(1:N));
    opti.subject_to(U(2:N+1) == U_next);
    
    w_dot = casadi.Function('w_dot', {aze, q, u}, {aze + q*u}, {'aze' 'q' 'u'}, {'w_dot'});
    W_next = W(1:N) + dt*w_dot(Aze(1:N),Q(1:N),U(1:N));
    opti.subject_to(W(2:N+1) == W_next);

    % f_axe = casadi.Function('f_Axe', {fx, theta}, {fx/200 - 9.81*sin(theta)}, {'fx' 'theta'}, {'f_axe'});
    % f_axe = casadi.Function('f_Axe', {theta, u}, {50000/200 - 9.81*sin(theta) - 0.05*0.7*1.225*u^2*0.5/200}, {'theta', 'u'}, {'f_axe'});
    f_axe = casadi.Function('f_Axe', {theta, u, w}, {50000/200 - 9.81*sin(theta) - (0.05*0.7*sqrt((u^2 + w^2) + 1e-9)^2*1.225*0.5 * (u/(abs(u + 1e-9)+abs(w + 1e-9) + 1e-9)))/200}, {'theta', 'u', 'w'}, {'f_axe'});
    Axe_next = f_axe(Theta(1:N), U(1:N),W(1:N));
    opti.subject_to(Axe(2:N+1) == Axe_next);

    % f_aze = casadi.Function('f_Aze', {fz, theta}, {fz/200 + 9.81*cos(theta)}, {'fz' 'theta'}, {'f_aze'});
    f_aze = casadi.Function('f_Aze', {fz, theta, u, w}, {fz/200 + 9.81*cos(theta) - (0.05*0.7*sqrt((u^2 + w^2) + 1e-9)^2*1.225*0.5 * (w/(abs(u + 1e-9)+abs(w + 1e-9) + 1e-9)))/200}, {'fz' 'theta', 'u', 'w'}, {'f_aze'});
    Aze_next = f_aze(Fz,Theta(1:N),U(1:N),W(1:N));
    opti.subject_to(Aze(2:N+1) == Aze_next);

    xe_dot = casadi.Function('xe_dot', {u, w, theta}, {u*cos(theta) + w*sin(theta)}, {'u' 'w' 'theta'}, {'xe_dot'});
    Xe_next = Xe(1:N) + dt*xe_dot(U(1:N),W(1:N),Theta(1:N));
    opti.subject_to(Xe(2:N+1) == Xe_next);

    ze_dot =casadi.Function('ze_dot', {u, w, theta}, {-u*sin(theta) + w*cos(theta)}, {'u' 'w' 'theta'}, {'ze_dot'});
    Ze_next = Ze(1:N) + dt*ze_dot(U(1:N),W(1:N),Theta(1:N));
    opti.subject_to(Ze(2:N+1) == Ze_next);

    q_dot = casadi.Function('q_dot', {my}, {my/1000}, {'my'}, {'q_dot'});
    Q_next = Q(1:N) + dt*q_dot(My);
    opti.subject_to(Q(2:N+1) == Q_next);
    
    theta_dot = casadi.Function('theta_dot', {q}, {q}, {'q'}, {'theta_dot'});
    Theta_next = Theta(1:N) + dt*theta_dot(Q(1:N));
    opti.subject_to(Theta(2:N+1) == Theta_next);


% initial states
opti.subject_to(U(1)==X(1));
opti.subject_to(W(1)==X(2));
opti.subject_to(Axe(1)==X(3));
opti.subject_to(Aze(1)==X(4));
opti.subject_to(Xe(1)==X(5));
opti.subject_to(Ze(1)==X(6));
opti.subject_to(Q(1)==X(7));
opti.subject_to(Theta(1)==X(8));

% input constraints
% opti.subject_to(Fx<100000);
% opti.subject_to(Fx>0);
opti.subject_to(Fz<50000); 
opti.subject_to(Fz>-50000);
opti.subject_to(My<1000); 
opti.subject_to(My>-1000);

% state constraints


% time constraint
opti.subject_to(T>=0);
% opti.subject_to(T<20);

% distance to target within delta
delta = 100;
opti.subject_to(  sqrt((Xe(N+1)-target_position(1))^2 + (Ze(N+1)-target_position(2))^2)  <= delta   );


% solve the problem
solver_options = struct;
solver_options.ipopt.print_level = 0;
solver_options.print_time = 0;
solver_options.verbose = 0;
solver_options.ipopt.sb ='yes';

opti.solver('ipopt', solver_options);           % select solver

try
    sol = opti.solve();
    % Fx_traj = sol.value(Fx);
    Fz_traj = sol.value(Fz);
    My_traj = sol.value(My);
    % control_fx = Fx_traj(1);
    control_fz = Fz_traj(1);
    control_my = My_traj(1);
    flag = 0;
catch
    % control_fx = nan;
    control_fz = nan;
    control_my = nan;
    flag = 1;
end

end

