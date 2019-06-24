% 2019-05-22 LLW script file to generate transition matrix for A


t_start = 0.0
dt      = 0.01
t_stop  = 10.0

% construct 144 element initial condition - concatenate
% 12 unit vectors, each all zeros except 1 in the ith element
% for computing transition matrix of bias + east estimator A(t) matrix
x0 = my_stack(eye(12));

% solve ODE
tic
[t x] = ode45('sys',t_start:dt:t_stop, x0);
toc
