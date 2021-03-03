function costfunout = cost_function_LQG(simin)

%---Parameters---%

global nStep N alpha

ns     = simin.ns;
nc     = simin.nc;
nf     = simin.nf;

%---Cost-function---%
Q = zeros(2*ns,2*ns,nStep+1);
R = repmat(10^-5*eye(nc), 1, 1, nStep);% Cost parameter for control action

In = eye(ns);
alpha = [alpha alpha];

for j = N+1:nStep+1
    for i = 1:ns
        ei = [In(:,i);-In(:,i)];
        Q(:,:,j) = Q(:,:,end) + alpha(i)*(ei*ei');
    end
end

for t = 1:N
    Q(:,:,t) = (t/N)^25*Q(:,:,end);
end

costfunout.Q = Q;
costfunout.R = R;

end

