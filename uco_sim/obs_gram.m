function mat = obs_gram(samp,t0,t1,dt,s0)
% 2019-05-22 LLW for computing observability grammian  of bias + east
% estimator

  %default to integration timestep 
  if( exist('dt') ~= 1)
    dt = 0.1;
  end
  
  dt_s = 0.01; %time step for ODE to propagate the state.
  
  mat = zeros(12,12);
  s = s0;
 
  for tau = t0+(0.5*dt):dt:t1
     %get the current state
    [~, s_vec] = ode45(@(t,s_vec) sysS(t,s_vec,samp),t0:dt_s:tau,s);
    s = s_vec(end,:)';
    
     %numerical integration 
    mat = mat + (obs_gram_integration_argument(samp,t0,tau,s) * dt);
    t0 = tau;
  end
end

function arg = obs_gram_integration_argument(samp,t0,t1,s0)

  % construct C - constant, not time varying
  C = [eye(3) zeros(3,9)];  

  % compute transition matrix
  Phi = tm(samp,t0,t1,s0); %compute transition matrix (a function of state)
  arg = Phi' * C' * C * Phi;

end

