function [simres] = mpc_ball_beam(simcase)


%%
% get terminal constraints
loadConstraints = true;
saveConstraints = true;

if not(loadConstraints)
    [K,P,alpha] = getTerminalConstraintsBallBeam(saveConstraints);
else
    load('terminalConstraintsBallBeam.mat');
end



%%
% dynamics of the Van der Pol Oscillator
global m g k
m = 10;
g = 9.81;
k = 6;

x1dot = @(x2) x2;
x2dot = @(x1, x3, x4) k^-1 * (m*x1*x4^2 - m*g*sin(x3));
x3dot = @(x4) x4;
x4dot = @(u) u;


% initial conditions and initialization
t0 = 0;
tend = 10;
tstep = 0.01;
time = t0:tstep:tend;


% initial states
x1(1) = -9;
x2(1) = 18;
x3(1) = deg2rad(-30);
x4(1) = deg2rad(0);



%% Run NMPC
% output to console
fprintf('------------------------------------------------\n');
fprintf('-- NMPC Ietartion for the Ball & Beam Example --\n');
fprintf('------------------------------------------------\n\n');

mpc_counter = 0;
linearControlActive = [];

% closed loop dynamics simulation
for i = 1 : length(time)-1
    switch simcase
        case 'QIH'
            threshold = alpha;
        case 'TR'
            threshold = alpha;
        case 'TC'
            threshold = 0;
        case 'NTC'
            threshold = 0;
    end

    if [x1(i);x2(i);x3(i);x4(i)]'*P*[x1(i);x2(i);x3(i);x4(i)] > threshold
            mpc_counter = mpc_counter + 1;
            [u, flags(i)] = ocp_ball_beam(simcase, x1(i), x2(i), x3(i), x4(i), P, alpha);
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
            linearControlActive(i) = 0;
        else
            controlInput(i) = -K*[x1(i);x2(i);x3(i);x4(i)];
            flags(i) = 0;
            linearControlActive(i) = 1;
    end


    % Euler-Cauchy Method for explicit solution of the IVP
    x1(i+1) = x1(i) + tstep*x1dot(x2(i));
    x2(i+1) = x2(i) + tstep*x2dot(x1(i), x3(i), x4(i));
    x3(i+1) = x3(i) + tstep*x3dot(x4(i));
    x4(i+1) = x4(i) + tstep*x4dot(controlInput(i));

    if abs(x1(i+1)) > 10 || abs(x2(i+1)) > 20 || abs(x3(i+1)) > pi/2 || abs(x4(i+1)) > 10
        disp('States out of bounds')
        time = time(1:i+1);
        controlInput = controlInput(1:i);
        break
    end

    if mod(i/(length(time)-1)*100, 5) == 0
        done = i/(length(time)-1)*100;
        disp([num2str(done),'% done']);
    end

end

%% Postprocessing
simres.x1 = x1;
simres.x2 = x2;
simres.x3 = x3;
simres.x4 = x4;
simres.u = controlInput;
simres.flags = flags;
simres.time = time;
simres.linearControlActive = linearControlActive;


%% Animation

% beamLength = 20;
% ballRadius = 0.5;
% 
% figure(5)
% hold on;
% grid on;
% title('Balancing Ball on Beam with overkill NMPC');
% motor = scatter(0,0,300,'black','filled');
% axis equal
% 
% % Plot the beam
% beamX = [-beamLength/2*cos(x3(1)), 0.5*beamLength * cos(x3(1))];
% beamY = [-beamLength/2*sin(x3(1)), beamLength/2 * sin(x3(1))];
% beamLine = plot(beamX, beamY, 'LineWidth', 2, 'Color', 'b');
% 
% % Plot the ball
% ballX = beamLength * cos(x3(1)) + x1(1);
% ballY = beamLength * sin(x3(1));
% ball = viscircles([ballX, ballY], ballRadius, 'EdgeColor', 'r');
% legend('motor','beam');
% % Update plot in a loop
% for i = 2:mpc_counter+1/tstep
%     % Update beam position
%     beamX = [-beamLength/2*cos(x3(i)), 0.5*beamLength * cos(x3(i))];
%     beamY = [-beamLength/2*sin(x3(i)), beamLength/2 * sin(x3(i))];
%     set(beamLine, 'XData', beamX, 'YData', beamY);
% 
%     % Delete the previous ball
%     delete(ball);
% 
%     % Plot the updated ball
%     ballX = cos(x3(i)) * x1(i);
%     ballY = x1(i) * sin(x3(i)) + ballRadius + 0.1;
%     ball = viscircles([ballX, ballY], ballRadius, 'EdgeColor', 'r');
% 
%     % Pause for animation
%     axis equal
%     axis([-12, 12, -10, 10]); % Set the axis limits
%     pause(0.005);
%     drawnow;
% end
% hold off;

end