function costfunout = cost_function_Robust(simin)

%--- Parameters ---%
global time nStep N delta

ns      = simin.ns;
nc      = simin.nc;
A       = simin.A;
A_hat   = simin.A_hat;
B       = simin.B;
H       = simin.H;
E       = simin.E;

%--- Cost function ---%
Q = zeros(size(A,1),size(A,2),nStep);
M = Q;
TM = Q;
Id = eye(ns);

runningalpha = zeros(ns,nStep); 
for i = 1:nStep
    fact = min(1,(i*delta/time))^6;
    temp = [fact*10^6 fact*10^6 fact*10^4 fact*10^4 1 1];
    runningalpha(:,i) = [temp temp]';
end

%Filling in the cost matrices
for j = 1:nStep
    for i = 1:ns
        Q(:,:,j) = Q(:,:,j) + runningalpha(i,j)*Id(:,i)*Id(:,i)';       
    end
end

%Filling in the cost matrices
for j = 1:nStep
    for i = 1:ns
        
        Q(:,:,j) = Q(:,:,j) + runningalpha(i,j)*Id(:,i)*Id(:,i)';
        
    end
end

%Signal Dependent Noise
nc = size(B,2);
Csdn = zeros(size(B,1),nc,nc);
for i = 1:nc
    Csdn(:,i,i) = .1*B(:,i);    
end

M = Q;
TM = Q;
D = eye(ns);

% Implementing the backwards recursions

M(:,:,end) = Q(:,:,end);
L = zeros(size(B,2),size(A,1),nStep-1);  % Optimal Minimax Gains
Lambda = zeros(size(A,1),size(A,2),nStep-1);

% Optimization of gamma
gamma = 50000;
minlambda = zeros(nStep-1,1);
gammaK = 0.5;
reduceStep = 1;
positive = false;
relGamma = 1;

while (relGamma > .001 || ~positive)

    for k = nStep-1:-1:1

        % Minimax Feedback Control
        TM(:,:,k) = gamma^2*eye(size(A))-D'*M(:,:,k+1)*D;
        minlambda(k) = min(eig(TM(:,:,k)));

        Lambda(:,:,k) = eye(size(A_hat))+(B*B'-gamma^-2*(D*D'))*M(:,:,k+1);
        M(:,:,k) = Q(:,:,k)+A_hat'*(M(:,:,k+1)^-1+B*B'-gamma^-2*D*D')^-1*A_hat;
        L(:,:,k) = B'*M(:,:,k+1)*Lambda(:,:,k)^-1*A_hat;

    end

    oldGamma = gamma;

    if min(real(minlambda)) >= 0

        gamma = (1-gammaK)*gamma;
        relGamma = (oldGamma-gamma)/oldGamma;
        positive = true;

    elseif min(real(minlambda)) < 0

        gamma = (1-gammaK)^-1*gamma;
        reduceStep = reduceStep + 0.5;
        relGamma = -(oldGamma-gamma)/oldGamma;
        gammaK = gammaK^reduceStep;
        positive = false;

    end
end

costfunout.Q     = Q;
costfunout.M     = M;
costfunout.TM    = TM;
costfunout.L     = L;
costfunout.D     = D;
costfunout.gamma = gamma;
costfunout.Csdn  = Csdn;

end

