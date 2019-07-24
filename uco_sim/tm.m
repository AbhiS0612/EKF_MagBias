function PHI=tm(samp,t_start,t_stop,s0)

% 2019-05-22 LLW for computing transition matrix of bias + east
% estimator A(t) matrix
  
% handle special case 
if (t_start == t_stop)
   PHI = eye(6);
  %PHI = eye(12);
  %PHI = eye(15);
  return;
end

%default to integration timestep 
if( exist('dt') ~= 1)
  dt = 0.1;
end

% handle short intervals
if((t_stop-t_start) < dt)
  dt = (t_stop-t_start)/10;
end 

% construct 144 element initial condition - concatenate
% 12 unit vectors, each all zeros except 1 in the ith element
 x0 = my_stack(eye(6));
%x0 = my_stack(eye(12));  
%x0 = my_stack(eye(15));

% solve ODE
[~, x] = ode45(@(t,x) sys(t,x,samp,s0),t_start:dt:t_stop, x0);

PHI = my_unstack(x(end,:)',6);
%PHI = my_unstack(x(end,:)',12);
%PHI = my_unstack(x(end,:)',15);
%PHI = eye(15);