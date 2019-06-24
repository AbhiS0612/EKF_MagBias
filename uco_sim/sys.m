function xdot = sys(t,x,samp,s0) % also needs current state

% 2019-05-22 LLW ode45 function for \dot{x}(t) = A(t) * x(t), runs 12
%                      copies at the same time
% for computing transition matrix of bias + east estimator A(t) matrix

  %s = [0.6946 0.1560 0.2994 0 0 0 1.0000 0 0 1.0000 0 1.0000]';
  
  A_mat = A(t,samp,s0); %needs current state

  xdot = kron(eye(12), A_mat) * x;

%  likely faster, but less transparent approach to compute same result
%  for i = 0:11
%    xdot(1+(i*12):12+(i*12),1) = A_mat * x(1+(i*12):12+(i*12),1);
%   end
  

    
