close all;
clear;
clc;

%% select simulation case

% simcase = 'QIH';    % Quasi-Infinite Horizon
% simcase = 'TR';     % Just Terminal Region
% simcase = 'TC';     % classic Terminal Constraint
% simcase = 'NTC';    % no Terminal Condition

simcases = {'QIH', 'TR', 'TC', 'NTC'};
simres = cell(1,4);

for i = 1 : 4
    simcase = simcases{i};
    simres{i} = mpc_ball_beam(simcase);
end


%% pots-processsing

% Set the desired figure width and height
figureWidth = 1200;  % in pixels
figureHeight = 400;  % in pixels

% linestyles
ls = {'-','--',':','-.'};

% Create a figure for x1 = r and x2 = r_dot
figure('Position', [100, 100, figureWidth, figureHeight]);
tiledlayout(1,2);

    nexttile;
    hold on;
    grid on;
    for i = 1 : length(simres)
        plot(simres{i}.time, simres{i}.x1,ls{i},'Color','Black','LineWidth',1.125);
        active = find(simres{i}.linearControlActive);
        if ~isempty(active)
            firstActive = active(1);
            scatter(simres{i}.time(firstActive), simres{i}.x1(firstActive),30,'black','filled',HandleVisibility='off');
        end
    end
    ylabel('$x_1 = r$','Interpreter','latex');
    xlabel('time [s]','Interpreter','latex');
    legend('QIH', 'TR', 'TC', 'NTC');
    hold off;
    
    nexttile;
    hold on;
    grid on;
    for i = 1 : length(simres)
        plot(simres{i}.time, simres{i}.x2,ls{i},'Color','Black','LineWidth',1.125);
        active = find(simres{i}.linearControlActive);
        if ~isempty(active)
            firstActive = active(1);
            scatter(simres{i}.time(firstActive), simres{i}.x2(firstActive),30,'black','filled',HandleVisibility='off');
        end
    end
    ylabel('$x_2 = \dot{r}$','Interpreter','latex');
    xlabel('time [s]','Interpreter','latex');
    legend('QIH', 'TR', 'TC', 'NTC');
    hold off;


%%
% Create a figure for x3 = alpha and x4 = alpha_dot
figure('Position', [100, 100, figureWidth, figureHeight]);
tiledlayout(1,2);

    nexttile;
    hold on;
    grid on;
    for i = 1 : length(simres)
        plot(simres{i}.time, simres{i}.x3,ls{i},'Color','Black','LineWidth',1.125);
        active = find(simres{i}.linearControlActive);
        if ~isempty(active)
            firstActive = active(1);
            scatter(simres{i}.time(firstActive), simres{i}.x3(firstActive),30,'black','filled',HandleVisibility='off');
        end
    end
    ylabel('$x_3 = \alpha$','Interpreter','latex');
    xlabel('time [s]','Interpreter','latex');
    legend('QIH', 'TR', 'TC', 'NTC');
    hold off;
    
    nexttile;
    hold on;
    grid on;
    for i = 1 : length(simres)
        plot(simres{i}.time, simres{i}.x4,ls{i},'Color','Black','LineWidth',1.125);
        active = find(simres{i}.linearControlActive);
        if ~isempty(active)
            firstActive = active(1);
            scatter(simres{i}.time(firstActive), simres{i}.x4(firstActive),30,'black','filled',HandleVisibility='off');
        end
    end
    ylabel('$x_4 = \dot{\alpha}$','Interpreter','latex');
    xlabel('time [s]','Interpreter','latex');
    legend('QIH', 'TR', 'TC', 'NTC');
    hold off;

%%
% Create a figure for control input u
figure('Position', [100, 100, figureWidth/2, figureHeight]);
hold on;
grid on;
for i = 1 : length(simres)
    plot(simres{i}.time(1:end-1), simres{i}.u,ls{i},'Color','Black','LineWidth',1.125);
    active = find(simres{i}.linearControlActive);
    if ~isempty(active)
        firstActive = active(1);
        scatter(simres{i}.time(firstActive), simres{i}.u(firstActive),30,'black','filled',HandleVisibility='off');
    end
end
ylabel('$u = \ddot{\alpha}$','Interpreter','latex');
xlabel('time [s]','Interpreter','latex');
legend('QIH', 'TR', 'TC', 'NTC');
hold off;

%%
% % Create a figure for flags
% figure('Position', [100, 100, figureWidth/2, figureHeight]);
% hold on;
% grid on;
% for i = 1 : length(simres)
%     scatter(1:length(simres{i}.flags),simres{i}.flags);
% end
% ylabel('flag','Interpreter','latex');
% xlabel('index','Interpreter','latex');
% legend('QIH', 'TR', 'TC', 'NTC');
% hold off;


%%
figure('Position', [100, 100, figureWidth, figureHeight]);
tiledlayout(1,2)
nexttile;
hold on;
grid on;
for i = 1 : length(simres)
        plot(simres{i}.x1, simres{i}.x2,ls{i},'Color','Black','LineWidth',1.125);
        active = find(simres{i}.linearControlActive);
        if ~isempty(active)
            firstActive = active(1);
            scatter(simres{i}.x1(firstActive), simres{i}.x2(firstActive),30,'black','filled',HandleVisibility='off');
        end
end
ylabel('$x_1$','Interpreter','latex');
xlabel('$x_2$','Interpreter','latex');
legend('QIH', 'TR', 'TC', 'NTC');
hold off;

nexttile;
% figure('Position', [100, 100, figureWidth/2, figureHeight]);
hold on;
grid on;
for i = 1 : length(simres)
        plot(simres{i}.x3, simres{i}.x4,ls{i},'Color','Black','LineWidth',1.125);
        active = find(simres{i}.linearControlActive);
        if ~isempty(active)
            firstActive = active(1);
            scatter(simres{i}.x3(firstActive), simres{i}.x4(firstActive),30,'black','filled',HandleVisibility='off');
        end
end
ylabel('$x_3$','Interpreter','latex');
xlabel('$x_4$','Interpreter','latex');
legend('QIH', 'TR', 'TC', 'NTC');
hold off;

