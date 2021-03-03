function [simout] = system_creation_LQG()

%---Parameters---%
global m k tau delta xfinal_LQG

% State vector representation
% 1-6:  RIGHT HAND x, y, dx, dy, fx and fy
% 7-12: LEFT HAND x, y, dx, dy, fx and fy
% Initial and final state vector for both hand (right, then left)
%---System---%
A = [0 0 1 0 0 0; 0 0 0 1 0 0;...
	 0 0 -k/m 0 1/m 0;...
	 0 0 0 -k/m 0 1/m; 0 0 0 0 -1/tau 0;...
	 0 0 0 0 0 -1/tau];
A = blkdiag(A,A);
	 
B = zeros(6,2);
B(5,1) = 1/tau;
B(6,2) = 1/tau;
B = blkdiag(B,B);

% DT representation
ns = size(A,1);
nc = size(B,2);
nf = size(xfinal_LQG,1);

Ad = eye(ns) + delta * A; % Adiscrete = 1 + delta_t*Acontinuous 
Bd = delta * B;           % Bdiscrete = delta_t*Bcontinuous 

% Expend matrixes for the final target state vector
A = [Ad,zeros(ns,nf);zeros(nf,ns),eye(nf)];
B = [Bd;zeros(nf,nc)];
A_hat = A;

simout.ns    = ns;
simout.nc    = nc;
simout.nf    = nf;
simout.A     = A;
simout.A_hat = A_hat;
simout.B     = B;

end

