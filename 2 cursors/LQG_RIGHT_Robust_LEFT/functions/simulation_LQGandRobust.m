function out = simulation_LQGandRobust(simin_LQG, simin_Robust,costfunin_LQG, costfunin_Robust)


%---PARAMETERS---%

global m delta theta xinit_LQG xfinal_LQG xinit_Robust xfinal_Robust...
    coeffQ nStep N right_perturbation left_perturbation numoftrials ...
    catch_trials learning_rates;

% Robust matrixes
A_Robust     = simin_Robust.A;
Ahat_Robust  = simin_Robust.A_hat;
B_Robust     = simin_Robust.B;
DA_Robust    = simin_Robust.DA;
H_Robust     = simin_Robust.H;
E_Robust     = simin_Robust.E;
D_Robust     = costfunin_Robust.D;
Q_Robust     = costfunin_Robust.Q;
M_Robust     = costfunin_Robust.M;
TM_Robust    = costfunin_Robust.TM;
gamma_Robust = costfunin_Robust.gamma;
Csdn_Robust  = costfunin_Robust.Csdn;

ns_Robust    = simin_Robust.ns;
nc_Robust    = simin_Robust.nc;
nStep_Robust = nStep;
N_Robust     = N;

% LQG matrixes

A_LQG     = simin_LQG.A;
Ahat_LQG  = simin_LQG.A_hat;
B_LQG     = simin_LQG.B;
Q_LQG     = costfunin_LQG.Q;
R_LQG     = costfunin_LQG.R;
L_LQG     = costfunin_LQG.L;

ns_LQG    = simin_LQG.ns;
nc_LQG    = simin_LQG.nc;
nf_LQG    = simin_LQG.nf;
nStep_LQG = nStep;
N_LQG     = N;

%--- SIMULATION ---%

switch right_perturbation
    case 'CCW'
        A_LQG(3,4) = -delta*(theta/m);
        A_Robust(3,4) = -delta*(theta/m);
        A_LQG(4,3) = delta*(theta/m);
        A_Robust(4,3) = delta*(theta/m);
    case 'CW'
        A_LQG(3,4) = delta*(theta/m);
        A_Robust(3,4) = delta*(theta/m);
        A_LQG(4,3) = -delta*(theta/m);
        A_Robust(4,3) = -delta*(theta/m);
    case 'BASELINE'
        A_LQG(3,4) = 0;
        A_Robust(3,4) = 0;
        A_LQG(4,3) = 0;
        A_Robust(4,3) = 0;
    case 'RANDOM'
        %
    otherwise
        error('The perturbation choice is incorrect !')
end

switch left_perturbation
    case 'CCW'
        A_LQG(9,10) = -delta*(theta/m);
        A_Robust(9,10) = -delta*(theta/m);
        A_LQG(10,9) = delta*(theta/m);
        A_Robust(10,9) = delta*(theta/m);
    case 'CW'
        A_LQG(9,10) = delta*(theta/m);
        A_Robust(9,10) = delta*(theta/m);
        A_LQG(10,9) = -delta*(theta/m);
        A_Robust(10,9) = -delta*(theta/m);
    case 'BASELINE'
        A_LQG(9,10) = 0;
        A_Robust(9,10) = 0;
        A_LQG(10,9) = 0;
        A_Robust(10,9) = 0;
    case 'RANDOM'
        %
    otherwise
        error('The perturbation choice is incorrect !')
end

% Initialization simulation

 % Robust
x = zeros(ns_Robust,nStep_Robust+1,numoftrials);
xhat = x;

 % LQG
z               = zeros(2*ns_LQG,nStep_LQG+1,numoftrials);% Initialize the state
z(1:ns_LQG,1,:) = repmat(xinit_LQG, 1, 1,numoftrials);% Initialize the estimated state
zhat        = z;                              % Initialize the state estiamte

control     = zeros(nc_LQG,nStep_LQG,numoftrials);    % Initialize control
avControl   = zeros(nc_LQG,nStep_LQG);                % Average Control variable


% Random indexes for catch trials
catch_trials_idx = [];

if catch_trials ~= 0
    while length(catch_trials_idx) ~= catch_trials
        random = randi(numoftrials, 1, 1); 
        catch_trials_idx = [catch_trials_idx random];
        catch_trials_idx = unique(catch_trials_idx);
    end
end

for p = 1:numoftrials

    % Robust
    x(:,1,p) = xinit_Robust - xfinal_Robust;
    xhat(:,1,p) = x(:,1,p);
    u_Robust = zeros(nStep_Robust-1,size(B_Robust,2)); % size(B,2) is the control dimension
    w = zeros(ns_Robust,1);
    Oxi = 0.001*B_Robust*B_Robust';
    Omega = eye(6)*Oxi(5,5);

    %Parameters for State Estimation
    Sigma = zeros(ns_Robust,ns_Robust,nStep_Robust);
    Sigma(:,:,1) = eye(ns_Robust)*10^-2;
    SigmaK = Sigma;

     % LQG
    z(ns_LQG+1:end,1,p)    = xfinal_LQG;
    zhat(ns_LQG+1:end,1,p) = xfinal_LQG;

    for i = 1:nStep_LQG-1
        
       %--- Catch trials ---%
       A_old_Robust = A_Robust;
        
        if (~isempty(catch_trials_idx)) & (i == catch_trials_idx)
            A_Robust(3,4)  = 0;
            A_Robust(9,10) = 0;
            A_Robust(4,3)  = 0;
            A_Robust(10,9) = 0;
        end
        
        A_old_LQG = A_LQG;
        
        if (~isempty(catch_trials_idx)) & (i == catch_trials_idx)
            A_LQG(3,4)  = 0;
            A_LQG(9,10) = 0;
            A_LQG(4,3)  = 0;
            A_LQG(10,9) = 0;
        end


        %--- Robust ---%
        
        sensoryNoise = mvnrnd(zeros(size(Omega,1),1),Omega)';
        sensoryNoise = [sensoryNoise; sensoryNoise];
        motorNoise_Robust = mvnrnd(zeros(size(Oxi,1),1),Oxi)';

        %MINMAX HINFTY CONTROL ------------------------------------------------
        %Riccati Equation for the State Estimator
        Sigma(:,:,i+1) = Ahat_Robust*(Sigma(:,:,i)^-1+H_Robust'*(E_Robust*E_Robust')^-1*H_Robust-gamma_Robust^-2*Q_Robust(:,:,i))^-1*Ahat_Robust'+D_Robust*D_Robust';

        %Feedback Eequation
        yx = H_Robust*x(:,i,p) + sensoryNoise;

        %Minmax Simulation with State Estimator
        u_Robust(i,:) = -B_Robust'*(M_Robust(:,:,i+1)^-1+B_Robust*B_Robust'-gamma_Robust^-2*(D_Robust*D_Robust'))^-1*Ahat_Robust*...   %Control
        (eye(ns_Robust)-gamma_Robust^-2*Sigma(:,:,i)*M_Robust(:,:,1))^-1*xhat(:,i,p);


        %Signal Dependent Noise - Robust Control
        sdn = 0;

        for isdn = 1:nc_Robust
            sdn = sdn + normrnd(0,1)*Csdn_Robust(:,:,isdn)*u_Robust(i,:)';
        end

        xhat(:,i+1,p) = Ahat_Robust*xhat(:,i,p) + B_Robust*u_Robust(i,:)'+...
        Ahat_Robust*(Sigma(:,:,i)^-1+H_Robust'*(E_Robust*E_Robust')^-1*H_Robust-gamma_Robust^-2*Q_Robust(:,:,i))^-1*(gamma_Robust^-2*Q_Robust(:,:,i)*xhat(:,i,p)+H_Robust'*(E_Robust*E_Robust')^-1*(yx-H_Robust*xhat(:,i,p)));

        % Minmax Simulation
        DA_Robust = A_Robust - Ahat_Robust;
        wx = DA_Robust*x(:,i,p); % Non zero if there is a model error. 
        x(:,i+1,p) = Ahat_Robust*x(:,i,p) + B_Robust*u_Robust(i,:)'+D_Robust*wx + motorNoise_Robust + sdn;

        %--- LQG ---%

        motorNoise_LQG   = mvnrnd(zeros(2*ns_LQG,1),(B_LQG*B_LQG'))'; % motor noise

        % Computation control vector and update optimal feedback gains
        u_LQG = -L_LQG(:,:,1)*z(:,i,p); % control variable
        Q_LQG = coeffQ*Q_LQG;
        L_LQG = LQGsolver(Ahat_LQG,B_LQG,Q_LQG(:,:,i+1:end),R_LQG(:,:,i+1:nStep_LQG));
        % Computation next state and next estimated state 
        zhat(:,i+1,p) = Ahat_LQG*z(:,i,p) + B_LQG*u_LQG;       % State Estimate
        z(:,i+1,p) = A_LQG*z(:,i,p) + B_LQG*u_LQG + motorNoise_LQG; % dynamics

        % UPDATE STATE VECTOR %
 
        % Right = LQG model z
        % Left  = Robust model x
        
        z(7:12,i+1,p)    = x(7:12,i+1,p)+ [-.06 .15 0 0 0 0]';
        zhat(7:12,i+1,p) = xhat(7:12,i+1,p) + [-.06 .15 0 0 0 0]';
        x(1:6,i+1,p)     = z(1:6,i+1,p) - [.06 .15 0 0 0 0]';
        xhat(1:6,i+1,p)  = zhat(1:6,i+1,p)- [.06 .15 0 0 0 0]';
        
        eps = z(1:6,i+1,p)- zhat(1:6,i+1,p);
        
        theta_up_R = Ahat_LQG(3,4);
        dzhat_dL = zeros(1,6);
        dzhat_dL(1,3) = zhat(4,i+1,p);
        theta_up_R = theta_up_R + learning_rates(1)*dzhat_dL*eps;
        Ahat_LQG(3,4) = theta_up_R;
        Ahat_Robust(3,4) = theta_up_R;
        
        theta_up_R = Ahat_LQG(4,3);
        dzhat_dL = zeros(1,6);
        dzhat_dL(1,4) = zhat(3,i+1,p);
        theta_up_R = theta_up_R + learning_rates(1)*dzhat_dL*eps;
        Ahat_LQG(4,3) = theta_up_R;
        Ahat_Robust(4,3) = theta_up_R;
        
        eps = x(7:12,i+1,p)- xhat(7:12,i+1,p);
        
        theta_up_L = Ahat_Robust(9,10);
        dzhat_dL = zeros(1,6);
        dzhat_dL(1,3) = xhat(10,i+1,p);
        theta_up_L  = theta_up_L + learning_rates(2)*dzhat_dL*eps;
        Ahat_Robust(9,10) = theta_up_L;
        Ahat_LQG(9,10) = theta_up_L;
        
        theta_up_L = Ahat_Robust(10,9);
        dzhat_dL = zeros(1,6);
        dzhat_dL(1,4) = xhat(9,i+1,p);
        theta_up_L  = theta_up_L + learning_rates(2)*dzhat_dL*eps;
        Ahat_Robust(10,9) = theta_up_L;
        Ahat_LQG(10,9) = theta_up_L;

        u = [u_LQG(1) u_LQG(2) u_Robust(i,3) u_Robust(i,4)]';
        control(:,i,p) = u;
        
        % Catch trials
        A_Robust = A_old_Robust;
        A_LQG    = A_old_LQG;
    end
    avControl = avControl + control(:,:,p)/numoftrials;
end

out.x            = x(:,1:N_Robust,:);
out.z            = z(:,1:N_LQG,:);
out.control      = control;
out.avControl    = avControl;

end

