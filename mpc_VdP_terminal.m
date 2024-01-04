close all;
clear;
clc;

% dynamics of the Van der Pol Oscillator
x1dot = @(x2) x2;
x2dot = @(x1, x2, u) 0.5*(1-x1^2)*x2 - x1 + u;


% initial conditions and initialization
t0 = 0;
tend = 30;
tstep = 0.05;
time = t0:tstep:tend;


% initial states
x1(1) = -3;
x2(1) = 3;

% Feedback gain
K = [0.414213562373093	1.94167511067722];
% level
alpha = 2;
% matrix P
P = [0.906788577923648	0.0848859743757874
    0.0848859743757874	0.731950953837246];


% test robustness agaionst random noise
enableNoise = true;


% output to console
fprintf('---------------------------------------------------\n');
fprintf('-- NMPC Iteration for the Van der Pol Oscillator --\n');
fprintf('---------------------------------------------------\n\n');

controlInput = zeros(length(time)-1,1);
flags = 0;

% closed loop dynamics simulation
for i = 1 : length(time)-1
    if [x1(i);x2(i)]'*P*[x1(i);x2(i)] > alpha
        %disp('>')
        [u, flags(i)] = VdP_OCP_terminal(x1(i), x2(i),P,alpha);
        if isnan(u)
            if i == 1
                % catch if initial problem is infeasible
                controlInput(i) = 0;
            else
                % Handling of infeasible solutions
                controlInput(i) = controlInput(i-1);
                % controlInput(i) = -controlInput(i-1);
                % controlInput(i) = 0;
            end
        else
            % normal MPC iteration
            controlInput(i) = u;
        end
    else
        %disp('<')
        controlInput(i) = -K*[x1(i);x2(i)];
    end

    % Euler-Cauchy Method for explicit solution of the IVP
    x2(i+1) = x2(i) + tstep*x2dot(x1(i), x2(i), controlInput(i));
    x1(i+1) = x1(i) + tstep*x1dot(x2(i));

    if enableNoise
        % add noise
        r = -0.001 + (0.001+0.001).*rand(1,1);
        x2(i+1) = x2(i+1)*(1+r);
        r = -0.001 + (0.001+0.001).*rand(1,1);
        x1(i+1) = x1(i+1)*(1+r);
    end
    if mod(i/(length(time)-1)*100, 10) == 0
        done = i/(length(time)-1)*100;
        disp([num2str(done),'% done']);
    end
end


%% Postprocessing

[V, D] = eig(P);

% Generate 100 values equally distributed along the border
num_points = 100;
theta_values = linspace(0, 2*pi, num_points);
border_points = zeros(2, num_points);
for i = 1:num_points
    x_theta = sqrt(alpha) * V * sqrt(D) * [cos(theta_values(i)); sin(theta_values(i))];
    border_points(:, i) = x_theta;
end

figure;
hold on;
title('State Trajectory');
plot(x1, x2);
xlabel('x1');
ylabel('x2');
grid on;
plot(border_points(1, :), border_points(2, :), 'b.');
hold off;



figure;
tiledlayout(3, 1);
nexttile;
hold on;
title('x1')
plot(time, x1);
xlabel('Time');
ylabel('x1');
grid on;
hold off;
    nexttile;
    hold on;
    title('x2')
    plot(time, x2);
    xlabel('Time');
    ylabel('x2');
    grid on;
    hold off;
        nexttile;
        hold on;
        title('u')
        plot(time(1:end-1), controlInput);
        xlabel('Time');
        ylabel('u');
        grid on;
        hold off;

figure;
hold on;
title('OCP Flags');
plot(flags);
ylabel('Flag');
xlabel('Iteration');
grid on;
hold off;



