function costfunout = cost_function_LQG(simin)

%---Parameters---%

global nStep N alpha

ns     = simin.ns;
nc     = simin.nc;
nf     = simin.nf;

%---Cost-function---%

Q = zeros(2*ns,2*ns,nStep+1);
R = repmat(10^-5*eye(nc), 1, 1, nStep);

% Fill in the cost of the last target

% Position
p = zeros(2*ns,2*ns);
p(:,1) = [0.5 zeros(1,5) 0.5 zeros(1,5) -0.5 zeros(1,5) -0.5 zeros(1,5)]';
p(:,2) = [0 0.5 zeros(1,4) 0 0.5 zeros(1,4) 0 -0.5 zeros(1,4) 0 -0.5 zeros(1,4)]';
p(:,[7, 13, 19]) = p(:,1).*[ones(2*ns,1) -ones(2*ns,2)];
p(:,[8, 14, 20]) = p(:,2).*[ones(2*ns,1) -ones(2*ns,2)];

% Velocity
vx   = [0;0;1;0;0;0];
temp = vx*vx';
Vx   = blkdiag(temp,temp,zeros(ns,ns));
vy   = [0;0;0;1;0;0];
temp = vy*vy';
Vy   = blkdiag(temp,temp,zeros(ns,ns));

% Sum        
C = p + Vx + Vy;

% Cost
beta  = [alpha(1) alpha(2) 0 0 0 0];
alpha = [alpha alpha beta beta];

for j = N+1:nStep+1
    for i = 1:2*ns
        temp = zeros(2*ns,2*ns);
        temp(:,i) = C(:,i);
        Q(:,:,j) = Q(:,:,end) + alpha(i)*temp;
    end
end

for t = 1:N
    Q(:,:,t) = (t/N)^3*Q(:,:,end);
end


costfunout.Q = Q;
costfunout.R = R;

end

