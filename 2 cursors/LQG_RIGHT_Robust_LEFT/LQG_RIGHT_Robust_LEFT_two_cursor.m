clear all
close all
clc

addpath(genpath(pwd));

%-------------------------------------------------------------------------%
% based on the papers:
%
% 1) "A very fast time scale of human motor adaptation: Within movement 
% adjustments of internal representations during reaching", Crevecoeur
% 2020.
%
% 2) "Optimal Task-Dependent Changes of Bimanual Feedback Control and 
% Adaptation", Diedrichsen 2007.

%% %---PARAMETERS---% %% 

global m k tau delta theta alpha learning_rates coeffQ time stab nStep N ...
    right_perturbation left_perturbation numoftrials catch_trials;

m              = 2.5;  % [kg]
k              = 0.1;  % [Nsm^-1]
tau            = 0.1;  % [s]
delta          = 0.01; % [s]  
theta          = 15;   % [N/(m/s)] - coeff pert force (F = +-L*dy/dt)
alpha          = [1000 1000 20 20 0 0];% [PosX, PosY, VelX, VelY, Fx, Fy]
learning_rates = [.1 .1];% [right left]
coeffQ         = 1;      % increase or decrease Q matrix during trials
time           = 0.6;    % [s] - experiment 1 - reaching the target
stab           = 0.01;   % [s] - experiment 1 - stabilization 
nStep          = round((time+stab)/delta)-1;
N              = round(time/delta);


% Protocol parameters

right_perturbation = 'BASELINE';     % CCW, CW or BASELINE (no FFs)
left_perturbation  = 'CW';% CCW, CW or BASELINE (no FFs)
numoftrials        = 15;        % number of trials 
catch_trials       = 0;         % number of catch trials


%% %---SYSTEM CREATION---% %%

global xinit_LQG xfinal_LQG xinit_Robust xfinal_Robust

% LQG - STATE VECTOR REPRESENTATION:
% (1-6)  RIGHT HAND x, y, dx, dy, fx and fy
% (7-12) LEFT HAND  x, y, dx, dy, fx and fy
xinit_LQG  = [.06 0 0 0 0 0 -.06 0 0 0 0 0]';    % [right left]
xfinal_LQG = [.06 .15 0 0 0 0 -.06 .15 0 0 0 0]';% [right left] 

% Robust - STATE VECTOR REPRESENTATION:
% (1-6)  RIGHT HAND x, y, dx, dy, fx and fy
% (7-12) LEFT HAND  x, y, dx, dy, fx and fy
% (13-14) MID POINT x and y
xinit_Robust  = [.06 0 0 0 0 0 -.06 0 0 0 0 0]';    % [right left]
xfinal_Robust = [.06 .15 0 0 0 0 -.06 .15 0 0 0 0]';% [right left] 

%---System---%

simout_LQG    = system_creation_LQG();
simout_Robust = system_creation_Robust();

%% %---COST FUNCTION---% %%

costfunin_LQG    = cost_function_LQG(simout_LQG);
costfunin_Robust = cost_function_Robust(simout_Robust);

%% %---OPTIMAL FEEDBACK GAIN---% %%

costfunin_LQG.L = LQGsolver(simout_LQG.A,simout_LQG.B,costfunin_LQG.Q,costfunin_LQG.R);

%% %---SIMULATION---% %%

out = simulation_LQGandRobust(simout_LQG, simout_Robust,costfunin_LQG, costfunin_Robust);

%% %---GRAPHS---% %%

translate = repmat([.06 .15 0 0 0 0 -.06 .15 0 0 0 0]',1,N,numoftrials);
x         = out.x(:,1:N,:) + translate;
control   = out.control;
avControl = out.avControl;

for trial = 1:numoftrials
    
    figure(1)
    
    % Position
    subplot(2,2,[1,2])
    hold on;
    midx = 0.5*(x(1,1:N,trial)+x(7,1:N,trial));
    midy = 0.5*(x(2,1:N,trial)+x(8,1:N,trial));
    plot(x(1,1:N,trial), x(2,1:N,trial)); 
    plot(x(7,1:N,trial), x(8,1:N,trial));
    plot(midx, midy);
    plot(0,0,'ro','LineWidth',2);
    plot(0,.15,'ro','MarkerSize',10,'LineWidth',2);
    plot(0.06,0,'ro','LineWidth',2);
    plot(0.06,.15,'ro','MarkerSize',10,'LineWidth',2);
    plot(-0.06,0,'ro','LineWidth',2);
    plot(-0.06,.15,'ro','MarkerSize',10,'LineWidth',2);
    xlabel('x-coord [m]'); ylabel('y-coord [m]'); title(['LQG (R) and Robust (L) models - 2 cursor - trajectories'],'FontSize',14);
    axis([-(max(x(1,:,1)) + 0.04) (max(x(1,:,1)) + 0.04)  -0.01 0.16])

    % Control
    subplot(2,2,4)
    plot([.01:.01:(nStep)*.01],control(1,:,trial)), hold on;
    plot([.01:.01:(nStep)*.01],avControl(1,:),'k','Linewidth',2)
    xlabel('Time [s]'); ylabel('Control [Nm]'); title('Control Vector - Right','FontSize',14);
    %axis square

    subplot(2,2,3)
    plot([.01:.01:(nStep)*.01],control(3,:,trial)), hold on;
    plot([.01:.01:(nStep)*.01],avControl(3,:),'k','Linewidth',2)
    xlabel('Time [s]'); ylabel('Control [Nm]'); title('Control Vector - Left','FontSize',14);
    %axis square

    %input(' ');    
end

% Velocity profiles

figure(2)
subplot(131)
plot([.01:.01:(nStep)*.01], x(3,:,end), 'b');hold on;
plot([.01:.01:(nStep)*.01], x(9,:,end), 'r');hold off;
xlabel('Time [s]');
ylabel('X Velocity [m/s]');
legend('right', 'left');

subplot(132)
plot([.01:.01:(nStep)*.01], x(4,:,end), 'b');hold on;
plot([.01:.01:(nStep)*.01], x(10,:,end), 'r');hold off;
xlabel('Time [s]');
ylabel('Y Velocity [m/s]');
legend('right', 'left');

subplot(133)
plot([.01:.01:(nStep)*.01], sqrt(x(3,:,end).^2+x(4,:,end).^2), 'b');hold on;
plot([.01:.01:(nStep)*.01], sqrt(x(9,:,end).^2+x(10,:,end).^2), 'r');hold off;
xlabel('Time [s]');
ylabel('Velocity [m/s]');
legend('right', 'left');

% Force profiles

figure(3)
subplot(131)
plot([.01:.01:(nStep)*.01], x(5,:,end), 'b');hold on;
plot([.01:.01:(nStep)*.01], x(11,:,end), 'r');hold off;
xlabel('Time [s]');
ylabel('X Force [N]');
legend('right', 'left');

subplot(132)
plot([.01:.01:(nStep)*.01], x(6,:,end), 'b');hold on;
plot([.01:.01:(nStep)*.01], x(12,:,end), 'r');hold off;
xlabel('Time [s]');
ylabel('Y Force [N]');
legend('right', 'left');

subplot(133)
plot([.01:.01:(nStep)*.01], sqrt(x(5,:,end).^2+x(11,:,end).^2), 'b');hold on;
plot([.01:.01:(nStep)*.01], sqrt(x(6,:,end).^2+x(12,:,end).^2), 'r');hold off;
xlabel('Time [s]');
ylabel('Force [N]');
legend('right', 'left');

