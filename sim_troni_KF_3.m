%% OUTDATED: Use exp_troni_KF_4 instead
 % Simulation of Giancarlo's algorithm for estimatating sensor bias without
 % knowledge of the vehicles true attitude (i.e using angular rates). This
 % code only works when there is no angular velocity bias.
 % Method: Extended Kalman Filter
 % 06-13-19:  CREATED by Abhi Shah
 % 07-01-19:  It seems that T can be observed to within a scale factor. One
 % of the diagonal elements must be known, the other 5 can be estimated as
 % relative to it.
 % 07-03-19: Bu term is always zero because u(t) is always zero, because Tp
 % is always in the kernel of D[T]
 
 %% Notes
  % System sometimes diverges when Q matrix is chosen larger than R
  % One of the diagonal elements must be known, since T can only be 
  % known upto a scale factor, i.e it will scale itself based on the known
  % element of T
  % The u(t) = f(x) - A(x).x = 0, at all times in this case
  % Still not entirely sure why the ODE version oscillates so much.
  
%% Add paths
addpath '/home/abhis/Matlab/DSCL/dscl_matlab-master';
addpath '/home/abhis/Matlab/DSCL/andrew-matlab'; 
addpath(genpath_nosvn_nogit_nohg('/home/abhis/Matlab/DSCL/dscl_matlab-master'));
addpath(genpath_nosvn_nogit_nohg('/home/abhis/Matlab/DSCL/andrew-matlab'));
%addpath(genpath_nosvn_nogit_nohg('/home/abhis/Matlab/DSCL/uco_sim'));
%ms = read_microstrain('/home/abhis/Matlab/DSCL/log/2019_06_12_18_26.MST');

%% Generate Simulated Data
lat = 39.33; % degrees
hz = 20;  % frequency of data generation  
rate_d = hz; % rate of discretization 
t_end = 300; % seconds
ts = 1/rate_d; %discretization time interval
w_max = [0.5,0.3,0.5]; %set max ang. vel. for excitation

%set biases 
bias.ang = [0;0;0];
bias.acc = [0;0;0];
bias.mag = [.1; 0.05; -.1];
% must be symmetric, also: assume that T is scaled such that T(1,1) = +1
% bias.T=[ 1.0    0.1   -0.2
%          0.1    1.1    0.05
%         -0.2    0.05    0.9];
%     
bias.T=[ 1.0    0.1  -0.15
         0.1    1.1   0.05
        -0.15   0.05   0.8];

ms = gen_samples(lat, hz, t_end, bias, w_max);

%% Read data
%ms = read_microstrain('/home/abhis/Matlab/DSCL/log/sim1.MST');
ms.t = ms.t-ms.t(1);
w = ms.ang';
z = ms.mag';
time = ms.t';
%phins = read_phins_imbat('/home/abhis/Matlab/DSCL/log/sim1.INS');
 
 %% System Set up
 runTime = t_end; % time in seconds
 %x = zeros(12, runTime*rate_d); %uncomment this if program runs slow
 x(:,1) = [z(:,1); 0;0;0; 1;0;0;1;0;1];  %initialize state 
 Ald = A_lin_dis(1,ms,x(:,1),ts); %linearized discrete time A
 C = [eye(3) zeros(3,9)];
 w_pos = zeros(size(w));
 
 % Kalman Filter Setup
 % Note that values for Q and R are dependent on the sampling rate.
 
% These values seem to work pretty well @ 20Hz, when using linearized prediction, choose slightly larger Q. 
% These values are for the discretized prediction model.
%  q1 = [1 1 1]*10^(-10);
%  q2 = [1 1 1]*2*10^(-11);
%  q3 = [0 1 1 1 1 1]*10^(-10);
%  Q = diag([q1 q2 q3]);
%  
%  stdev = [1 1 1]*2*10^(-4);
%  R = stdev.^2;
 
 q1 = [1 1 1]*10^(-10);
 q2 = [1 1 1]*2*10^(-11);
 q3 = [0 1 1 1 1 1]*10^(-10); % q3(1) must be set to 0
 Q = diag([q1 q2 q3]);
 
 stdev = [1 1 1]*2*10^(-4);
 R = stdev.^2; 

 err_mag = stdev;
 err_magb = [1 1 1]*10^(-3);
 err_T = [0 1 1 1 1 1]*10^(-3);
 err = [err_mag err_magb err_T]; %expected initial error - should the first 3 be the same as R?
 Sig = err'*err;
 
 %% Kalman Filter implematation
 tic
 
 for i = 2:(runTime*rate_d)  
   % prediction of next state: x(:,i) = Ald*x(:,i-1)+ Bld*u(:,t-1);
   % u turns out to be zero at all time (notes for details)
   % the system can therefore be simplified to
    x(:,i) = Ald*x(:,i-1);
    Sig = Ald*Sig*Ald' + Q;
    
% find state using ode solver: slower, and oscillates a lot, WHY??
%      tspan = [time(i-1) time(i)];
%      [t, x_vec] = ode45(@(t,x_vec) sysS(t,x_vec,ms), tspan, x(:,i-1));
%      x(:,i) = x_vec(end,:)';
%      
     
    %kalman gain calculation 
     K = (Sig*C')/(C*Sig*C' + R);
     
    % update equations 
     x(:,i) = x(:,i) + K*(z(:,i) - C*x(:,i));
     Sig = (eye(12) - K*C)*Sig;
    
     %recompute A matrix
     Ald = A_lin_dis(i,ms,x(:,i),ts); %linearized discrete time A
     
     % numericl integration to find angular position
     w_pos(:,i) = w_pos(:,i-1) + w(:,i-1)*ts;
  
     if ~mod(i,100)
        str = sprintf('Filter status: %i percent',floor((i*100)/(runTime*rate_d)));
        disp(str);
     end
 end
 toc
 
 
 % Plotting
%  figure(1)
% % subplot(2,1,1)
%  plot(time(1:length(x)), x(1,:), time(1:length(x)), x(2,:), time(1:length(x)), x(3,:));
%  legend ('x','y','z');
%  title 'State estimate';
%  grid minor;
%  subplot(2,1,2)
%  plot(time(1:length(x)),z(:,1:length(x)));
%  legend('x','y','z');
%  xlabel 'time (s)';
%  title 'Mag Measurements';
%  grid minor;
%  figure 
%  plot(ms.t,(w_pos).*180/pi);
%  title 'Orientation of the vehicle';
%  legend ('x','y','z');
%  ylabel 'angular displacement (degrees)'; 
%  xlabel 'time(s)';
%  grid minor;
% 


% figure(2)
%  plot(time(1:length(x)), x(4,:), time(1:length(x)), x(5,:), time(1:length(x)), x(6,:));
%  title 'Hard-Iron Effects';
%  ylabel 'magnetometer bias terms';
%  legend ('x','y','z');
%  xlabel 'time (s)';
%  title 'Biases';
%  grid minor;
  
% figure(3)
%  plot(time(1:length(x)), x(7:12,:));
%  title 'Soft-Iron Effects'
%  ylabel 'elememts of symmetric T matrix'
%  legend ('a','b','c','d','e','f');
% % ylim([0.9 1.2]);
%  xlabel 'time (s)';
%  grid minor;
 
%  
% %  %% Display Results
%  disp 'Computed mean biases:';
%  b_mean = mean(x(4:6,length(x)-100:length(x))')'
%  disp 'True bias:';
%  bias.mag
%  disp 'Computed mean T matrix:';
%  Ts = mean(x(7:12,floor(length(x)-100):length(x))');
%  T_mean = [Ts(1:3)' [Ts(2) Ts(4) Ts(5)]' [Ts(3) Ts(5) Ts(6)]']
%  disp 'True T matrix:';
%  bias.T
 
% %% function to evaluate f_x 
%  function f = f_x(t,samp, x)
%   we = samp.ang(t,:)';
%   xm = x(1:3,1);
%   b = x(4:6,1);
%   Ts = x(7:12,1)';
%   T = [Ts(1:3)' [Ts(2) Ts(4) Ts(5)]' [Ts(3) Ts(5) Ts(6)]'];
% 
%   f=  [-(T*J(we)/T)*(xm-b);
%               zeros(9,1)];
%  end
%  
%  % function to compute the discrete time B matrix
%  function B = B_dis(A, ts)
%    Ald = expm(A*ts);
%    arg_sum = zeros(12,12);
%    step = ts/10;
%    for i = 0:step:ts
%        arg_sum = arg_sum + expm(-A*i)*step;
%    end
%    
%    B = Ald*arg_sum;
%  end
 