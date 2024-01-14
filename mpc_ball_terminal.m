close all;
clear;
clc;
tic

% dynamics of the Van der Pol Oscillator
global m k g
m = 10;
k = 6;
g = 9.81;

x1dot = @(x2) x2;
x2dot = @(x1, x3, x4) (m*x1*x4^2 - m*g*sin(x3))/(k);
x3dot = @(x4) x4;
x4dot = @(u) u;

% initial conditions and initialization
t0 = 0;
tend = 10;
tstep = 0.01;
time = t0:tstep:tend;


% initial states
x1(1) = -8;
x2(1) = 15;
x3(1) = deg2rad(-30);
x4(1) = deg2rad(0);

% Feedback gain
K = [-1.41421356237310	-2.09568031754641	19.6070713670136	6.27009910081389];
% level
alpha = 3.96;
% matrix P
P = [156.137272667873	165.979604294874	-1029.57531099601	-135.990290148617
165.979604294874	200.491036358359	-1372.76821971225	-197.100041534732
-1029.57531099601	-1372.76821971225	10831.4364857583	1827.92041664604
-135.990290148617	-197.100041534732	1827.92041664604	365.781641903947];


% output to console
fprintf('--------------------------------------\n');
fprintf('-- NMPC Iteration for the Paper --\n');
fprintf('------------------------------------\n\n');

controlInput = zeros(length(time)-1,1);
flags = 0;
mpc_counter = 0;

% closed loop dynamics simulation
for i = 1 : length(time)-1
    if [x1(i);x2(i);x3(i);x4(i)]'*P*[x1(i);x2(i);x3(i);x4(i)] > alpha
        disp('>')
        mpc_counter = mpc_counter + 1;
        [u, flags(i)] = Ball_OCP_terminal(x1(i), x2(i),x3(i),x4(i),P,alpha);
        if isnan(u)
            if i == 1
                % catch if initial problem is infeasible
                %controlInput(i) = -(K(1)*x1(i) + K(2)*x2(i));
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
        controlInput(i) = -K*[x1(i);x2(i);x3(i);x4(i)];
    end

    % Euler-Cauchy Method for explicit solution of the IVP
    x4(i+1) = x4(i) + tstep*x4dot(controlInput(i));
    x3(i+1) = x3(i) + tstep*x3dot(x4(i));
    x2(i+1) = x2(i) + tstep*x2dot(x1(i), x3(i), x4(i));
    x1(i+1) = x1(i) + tstep*x1dot(x2(i));

    % Runge-Kutta update (RK4)
    % k1_x2 = tstep*x2dot(x1(i), x2(i), controlInput(i));
    % k1_x1 = tstep*x1dot(x1(i), x2(i), controlInput(i));
    % 
    % k2_x2 = tstep*x2dot(x1(i) + 0.5*k1_x1, x2(i) + 0.5*k1_x2, controlInput(i));
    % k2_x1 = tstep*x1dot(x1(i) + 0.5*k1_x1, x2(i) + 0.5*k1_x2, controlInput(i));
    % 
    % k3_x2 = tstep*x2dot(x1(i) + 0.5*k2_x1, x2(i) + 0.5*k2_x2, controlInput(i));
    % k3_x1 = tstep*x1dot(x1(i) + 0.5*k2_x1, x2(i) + 0.5*k2_x2, controlInput(i));
    % 
    % k4_x2 = tstep*x2dot(x1(i) + k3_x1, x2(i) + k3_x2, controlInput(i));
    % k4_x1 = tstep*x1dot(x1(i) + k3_x1, x2(i) + k3_x2, controlInput(i));
    % 
    % x2(i+1) = x2(i) + (1/6)*(k1_x2 + 2*k2_x2 + 2*k3_x2 + k4_x2);
    % x1(i+1) = x1(i) + (1/6)*(k1_x1 + 2*k2_x1 + 2*k3_x1 + k4_x1);

    if abs(x1(i+1)) > 10 || abs(x2(i+1)) > 20 || abs(x3(i+1)) > pi/2 || abs(x4(i+1)) > 10
        disp('States out of bounds')
        time = time(1:i+1);
        controlInput = controlInput(1:i);
        break
    end
    if mod(i/(length(time)-1)*100, 10) == 0
        done = i/(length(time)-1)*100;
        disp([num2str(done),'% done']);
        disp(time(i))
    end
end


%% Postprocessing

figure(1);
hold on;
title('State Trajectory');
plot(x1, x2);
xlabel('r');
ylabel('rdot');
grid on;

% Define the quadratic form expression x'Px
% Generate a grid of points
[x1plot, x2plot] = meshgrid(linspace(-2, 2, 100), linspace(-2, 2, 100));

% Evaluate the quadratic form for each point
quad_form = x1plot .* (P(1,1) * x1plot + P(1,2) * x2plot) + x2plot .* (P(2,1) * x1plot + P(2,2) * x2plot);

% Create a contour plot
contour(x1plot, x2plot, quad_form, [alpha, alpha], 'LineWidth', 2);

hold off;

figure(2);
subplot(2, 2, 1);
hold on;
title('x1')
plot(time, x1);
xlabel('Time');
ylabel('x1');
grid on;

subplot(2,2,2)
title('x2')
plot(time, x2);
xlabel('Time');
ylabel('x2');
grid on;

subplot(2,2,3)
hold on;
title('x3')
plot(time, x3);
xlabel('Time');
ylabel('x3');
grid on;

subplot(2,2,4)
hold on;
title('x4')
plot(time, x4);
xlabel('Time');
ylabel('x4');
grid on;

figure(3);
hold on;
title('OCP Flags');
plot(flags);
ylabel('Flag');
xlabel('Iteration');
grid on;
hold off;

figure(4)
hold on;
title('u')
plot(time(1:end-1), controlInput);
xlabel('Time');
ylabel('u');
grid on;

% %% plot vector field
% 
% x1dot = @(x1, x2) x2 - K*[x1;x2]*(mu + (1-mu)*x1);
% x2dot = @(x1, x2) x1 - K*[x1;x2]*(mu - 4*(1 - mu)*x2);
% 
% % Define the grid of points where you want to plot the vectors
% [x11, x22] = meshgrid(-2:0.1:2, -2:0.1:2);
% 
% % Evaluate the vector field at each point
% u = zeros(size(x11));
% v = zeros(size(x22));
% 
% for i = 1:numel(x11)
%     u(i) = x1dot(x11(i), x22(i));
%     v(i) = x2dot(x11(i), x22(i));
% end
% 
% % Normalize the vectors for better visualization
% magnitude = sqrt(u.^2 + v.^2);
% u_normalized = u ./ magnitude;
% v_normalized = v ./ magnitude;
% 
% % Plot the vector field using quiver
% hold on
% figure(1);
% hold on
% quiver(x11, x22, u, v, 'AutoScale', 'on', 'AutoScaleFactor', 0.9);
% 
% % Set plot properties
% xlabel('x1');
% ylabel('x2');
% grid on;
% 
toc

%%
beamLength = 20;
ballRadius = 0.5;

figure(5)
axis equal

% Plot the beam
beamX = [-beamLength/2*cos(x3(1)), 0.5*beamLength * cos(x3(1))];
beamY = [-beamLength/2*sin(x3(1)), beamLength/2 * sin(x3(1))];
beamLine = plot(beamX, beamY, 'LineWidth', 2, 'Color', 'b');

% Plot the ball
ballX = beamLength * cos(x3(1)) + x1(1);
ballY = beamLength * sin(x3(1));
ball = viscircles([ballX, ballY], ballRadius, 'EdgeColor', 'r');

% Update plot in a loop
for i = 2:mpc_counter+1/tstep
    % Update beam position
    beamX = [-beamLength/2*cos(x3(i)), 0.5*beamLength * cos(x3(i))];
    beamY = [-beamLength/2*sin(x3(i)), beamLength/2 * sin(x3(i))];
    set(beamLine, 'XData', beamX, 'YData', beamY);

    % Delete the previous ball
    delete(ball);

    % Plot the updated ball
    ballX = cos(x3(i)) * x1(i);
    ballY = x1(i) * sin(x3(i)) + ballRadius + 0.1;
    ball = viscircles([ballX, ballY], ballRadius, 'EdgeColor', 'r');
    
    % Pause for animation
    axis equal
    axis([-12, 12, -10, 10]); % Set the axis limits
    pause(0.005);
    drawnow;
end

hold off;
