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
muvdp = 0.75;

x1dot = @(x2) x2;
x2dot = @(x1, x2, u) muvdp*(1-x1^2)*x2 - x1 + u;


% initial conditions and initialization
t0 = 0;
tend = 40;
tstep = 0.01;
time = t0:tstep:tend;


% initial states
x1(1) = 3;
x2(1) = 4;

x1ol(1) = x1(1);
x2ol(1) = x2(1);

x1lin(1) = x1(1);
x2lin(1) = x2(1);

enableNoise = false;


%% Run NMPC
% output to console
fprintf('---------------------------------------------------\n');
fprintf('-- NMPC Ietartion for the Van der Pol Oscillator --\n');
fprintf('---------------------------------------------------\n\n');

% closed loop dynamics simulation
for i = 1 : length(time)-1
    if [x1(i);x2(i)]'*P*[x1(i);x2(i)] > alpha
        [u, flags(i), uopt{i}, xopt{i}, timeopt{i}] = ocp_van_der_pol(x1(i), x2(i), P, alpha);
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

    x1Noise = 0;
    x2Noise = 0;

    % create noise
    if enableNoise
        x1Noise = 0.01*randn(1,1);
        x2Noise = 0.01*randn(1,1);
    end

    % Euler-Cauchy Method for explicit solution of the IVP
    x2(i+1) = x2(i) + tstep*x2dot(x1(i), x2(i), controlInput(i)) + x1Noise;
    x1(i+1) = x1(i) + tstep*x1dot(x2(i)) + x2Noise;

    x2ol(i+1) = x2ol(i) + tstep*x2dot(x1ol(i), x2ol(i), 0);
    x1ol(i+1) = x1ol(i) + tstep*x1dot(x2ol(i));

    controlInputLinear(i) = -K*[x1lin(i);x2lin(i)];
    x2lin(i+1) = x2lin(i) + tstep*x2dot(x1lin(i), x2lin(i), controlInputLinear(i)) + x1Noise;
    x1lin(i+1) = x1lin(i) + tstep*x1dot(x2lin(i)) + x2Noise;

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
plot(x1lin,x2lin);
plot(x1ol,x2ol)
[C,h] = pcontour(V, alpha, 10*[-1 1 -1 1]);
% quiver(x11, x22, v_normalized, w_normalized, 'AutoScale', 'on', 'AutoScaleFactor', 0.5);
quiver(x11, x22, v_normalized, w_closed_normalized, 'AutoScale', 'on', 'AutoScaleFactor', 2);
xlabel('x1');
ylabel('x2');
grid on;
legend('NMPC trajectory', 'linear controller', 'open loop', 'terminal region', 'vector field');
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

%% NMPC Trajectories over time

% figure;
% hold on;
% grid on;
% skip = 10;
% for i = 1 : skip : length(xopt)
%     colorValue = (i - 1) / (length(xopt) - 1);
%     % Use the color value to create a color between blue and red
%     color = [colorValue, 0, 1-colorValue];  % RGB values for blue to red transition
%     uplot = uopt{i};
%     plot(timeopt{i}(1:end-1)+i*tstep,uplot,'Color',color);
% end
% xlabel('time [s]','Interpreter','latex');
% ylabel('input $u$','Interpreter','latex');
% hold off;

% Set the desired figure width and height
figureWidth = 1200;  % in pixels
figureHeight = 400;  % in pixels

% Create a figure for x1 = r and x2 = r_dot
figure('Position', [100, 100, figureWidth, figureHeight]);
tiledlayout(1,2);

nexttile;
hold on;
grid on;
plot(x1,x2, 'LineWidth',1.);
plot(x1ol,x2ol,'LineWidth',1.);
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');
legend('controlled','uncontrolled','Interpreter','latex','Location','northwest');
hold off;


nexttile;
hold on;
grid on;
skip = 10;
for i = 1 : skip : length(xopt)
    colorValue = (i - 1) / (length(xopt) - 1);
    % Use the color value to create a color between blue and red
    color = [colorValue, 0, 1-colorValue];  % RGB values for blue to red transition
    uplot = uopt{i};
    plot(timeopt{i}(1:end-1)+i*tstep,uplot(1,:),'Color',color);
end
xlabel('time [s]','Interpreter','latex');
ylabel('$u$','Interpreter','latex');
hold off;









