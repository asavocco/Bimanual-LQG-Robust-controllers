function [simout] = system_creation_Robust()

%---Parameters---%
global m k tau delta theta

%---System---%
A = [0 0 1 0 0 0; 0 0 0 1 0 0;...
	 0 0 -k/m 0 1/m 0;...
	 0 0 0 -k/m 0 1/m; 0 0 0 0 -1/tau 0;...
	 0 0 0 0 0 -1/tau];
A = blkdiag(A,A);
A = [A zeros(12,2); zeros(2,14)];
A(13,:) = [0 0 0.5 0 0 0 0 0 0.5 0 0 0 0 0];
A(14,:) = [0 0 0 0.5  0 0 0 0 0 0.5 0 0 0 0];
B = zeros(6,2);
B(5,1) = 1/tau;
B(6,2) = 1/tau;
B = blkdiag(B,B);
B = [B; zeros(2,4)];

ns = size(A,1);
nc = size(B,2);

A_hat = A;
DA = (A-A_hat)*delta; % Used when there is a model error
A = eye(size(A))+delta*A;
A_hat = eye(size(A_hat))+delta*A_hat;
B = delta*B;

% Observability Matrix
H = eye(size(A,1));
E = eye(ns,1)';          %See Basar and Tamer, pp. 171

simout.ns     = ns;
simout.nc     = nc;
simout.A      = A;
simout.A_hat  = A_hat;
simout.B      = B;
simout.DA     = DA;
simout.H      = H;
simout.E      = E;

end

