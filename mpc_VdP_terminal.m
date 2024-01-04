close all;
clear;
clc;


% dynamics of the Van der Pol Oscillator
x1dot = @(x2) x2;
x2dot = @(x1, x2, u) 0.5*(1-x1^2)*x2 - x1 + u;


% initial conditions and initialization
t0 = 0;
tend = 50;
tstep = 0.05;
time = t0:tstep:tend;


% initial states
x1(1) = -5;
x2(1) = -5;

% Feedback gain
K = [0.414213562373093	2.25454470582718];
% level
alpha = 2000;
% matrix P
P = [261.524853855428	158.946158648386
158.946158648386	181.090662200890];


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
        disp('>')
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
        disp('<')
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

figure;
hold on;
title('State Trajectory');
plot(x1, x2);
xlabel('x1');
ylabel('x2');
grid on;

% Define the quadratic form expression x'Px
% Generate a grid of points
[x1plot, x2plot] = meshgrid(linspace(-5, 5, 100), linspace(-5, 5, 100));

% Evaluate the quadratic form for each point
quad_form = x1plot .* (P(1,1) * x1plot + P(1,2) * x2plot) + x2plot .* (P(2,1) * x1plot + P(2,2) * x2plot);

% Create a contour plot
contour(x1plot, x2plot, quad_form, [alpha, alpha], 'LineWidth', 2);

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


