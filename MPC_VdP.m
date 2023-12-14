close all;
clear;
clc;

% dynamics of the Van der Pol Oscillator
x1dot = @(x2) x2;
x2dot = @(x1, x2, u) 0.5*(1-x1^2)*x2 - x1 + u;


% initial conditions and initialization
t0 = 0;
tend = 30;
tstep = 0.1;
time = t0:tstep:tend;


% initial states
x1(1) = 0.1;
x2(1) = 0.1;


% test robustness agaionst random noise
enableNoise = false;


% output to console
fprintf('---------------------------------------------------\n');
fprintf('-- NMPC Ietartion for the Van der Pol Oscillator --\n');
fprintf('---------------------------------------------------\n\n');


% closed loop dynamics simulation
for i = 1 : length(time)-1
    [u, flags(i)] = VdP_OCP(x1(i), x2(i));
    if isnan(u)
        if i == 1
            % catch if initial problem is infeasible
            controlInput(i) = 0;
        else
            % Handling of infeasible solutions
            % controlInput(i) = controlInput(i-1);
            % controlInput(i) = -controlInput(i-1);
            controlInput(i) = 0;
        end
    else
        % normal MPC iteration
        controlInput(i) = u;
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











