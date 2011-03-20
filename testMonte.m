%%%
%matlabpool local 4 % only need to run once per session, parallel setup 
% remains initialized after you clear the workspace
% Rather unfortunately, this will actually throw an error if matlabpool is 
% already initialized.
N = 4^7; 
% mean and standard deviation of "noisy" parameter
sigma = 1;
mu = 2; % mean of 2, std. deviation of 1

out = zeros(N,1);

b = 5; % some invariant parameter

parfor i = 1:N
   r = randn(1)*sigma+mu;
   % It would be hard to think of a physical system that this would
   % describe.
   out(i) = r*abs(r)-r*b;

end

figure
hist(out,100)
mean(out)
std(out)

wander = mu*abs(mu) - mu*b; % would be the same as mean(out) if sigma = 0;
disp(wander)

%%%