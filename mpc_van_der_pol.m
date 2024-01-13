close all;
clear;
clc;

%%
% get terminal constraints
loadConstraints = true;
saveConstraints = true;

if not(loadConstraints)
    [K,P,alpha] = getTerminalConstraintsVdP(saveConstraints);
else
    load('terminalConstraintsVdP.mat');
end



%%
% dynamics of the Van der Pol Oscillator
x1dot = @(x2) x2;
x2dot = @(x1, x2, u) 0.5*(1-x1^2)*x2 - x1 + u;


% initial conditions and initialization
t0 = 0;
tend = 20;
tstep = 0.01;
time = t0:tstep:tend;


% initial states
x1(1) = 6;
x2(1) = 8;

x1ol(1) = x1(1);
x2ol(1) = x2(1);

x1lin(1) = x1(1);
x2lin(1) = x2(1);


%% Run NMPC
% output to console
fprintf('---------------------------------------------------\n');
fprintf('-- NMPC Ietartion for the Van der Pol Oscillator --\n');
fprintf('---------------------------------------------------\n\n');

% closed loop dynamics simulation
for i = 1 : length(time)-1
    if [x1(i);x2(i)]'*P*[x1(i);x2(i)] > alpha
        [u, flags(i)] = ocp_van_der_pol(x1(i), x2(i), P, alpha);
        if isnan(u)
            if i == 1
                % catch if initial problem is infeasible
                controlInput(i) = 0;
                % controlInput(i) = -K*[x1(i);x2(i)];
            else
                % Handling of infeasible solutions
                controlInput(i) = controlInput(i-1);
                % controlInput(i) = -controlInput(i-1);
                % controlInput(i) = 0;
                % controlInput(i) = -K*[x1(i);x2(i)];
            end
        else
            % normal MPC iteration
            controlInput(i) = u;
        end
    else
        controlInput(i) = -K*[x1(i);x2(i)];
        flags(i) = 0;
    end

    % Euler-Cauchy Method for explicit solution of the IVP
    x2(i+1) = x2(i) + tstep*x2dot(x1(i), x2(i), controlInput(i));
    x1(i+1) = x1(i) + tstep*x1dot(x2(i));

    x2ol(i+1) = x2ol(i) + tstep*x2dot(x1ol(i), x2ol(i), 0);
    x1ol(i+1) = x1ol(i) + tstep*x1dot(x2ol(i));

    controlInputLinear(i) = -K*[x1lin(i);x2lin(i)];
    x2lin(i+1) = x2lin(i) + tstep*x2dot(x1lin(i), x2lin(i), controlInputLinear(i));
    x1lin(i+1) = x1lin(i) + tstep*x1dot(x2lin(i));

    if mod(i/(length(time)-1)*100, 10) == 0
        done = i/(length(time)-1)*100;
        disp([num2str(done),'% done']);
    end
end


%% Postprocessing

x = mpvar('x', [2 1]);
V = x'*P*x;

% Define the grid of points where you want to plot the vectors
[x11, x22] = meshgrid(-10:0.5:10, -10:0.5:10);

% Evaluate the vector field at each point
v = zeros(size(x11));
w = zeros(size(x22));
w_closed = zeros(size(x22));
q = zeros(size(x11));

for i = 1:numel(x11)
    v(i) = x1dot(x22(i));
    w(i) = x2dot(x11(i), x22(i), q(i));
    w_closed(i) = x2dot(x11(i), x22(i), -K*[x11(i); x22(i)]);
end

% Normalize the vectors for better visualization
magnitude = sqrt(v.^2 + w.^2);
% magnitude_closed = sqrt(v.^2 + w_closed.^2);
magnitude_closed = 1;

v_normalized = v ./ magnitude;
w_normalized = w ./ magnitude;
w_closed_normalized = w_closed ./ magnitude_closed;

figure;
hold on;
title('State Trajectory');
plot(x1, x2);
% plot(x1ol,x2ol);
plot(x1lin,x2lin);
[C,h] = pcontour(V, alpha, 10*[-1 1 -1 1]);
% quiver(x11, x22, v_normalized, w_normalized, 'AutoScale', 'on', 'AutoScaleFactor', 0.5);
quiver(x11, x22, v_normalized, w_closed_normalized, 'AutoScale', 'on', 'AutoScaleFactor', 2);
xlabel('x1');
ylabel('x2');
grid on;
legend('NMPC trajectory', 'terminal region', 'vector field');
axis equal;
hold off;


figure;
tiledlayout(3, 1);
nexttile;
hold on;
title('x1')
plot(time, x1);
plot(time, x1lin);
xlabel('Time');
ylabel('x1');
grid on;
hold off;
    nexttile;
    hold on;
    title('x2')
    plot(time, x2);
    plot(time, x2lin);
    xlabel('Time');
    ylabel('x2');
    grid on;
    hold off;
        nexttile;
        hold on;
        title('u')
        plot(time(1:end-1), controlInput);
        plot(time(1:end-1), controlInputLinear);
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











