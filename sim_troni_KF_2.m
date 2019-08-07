%% OUTDATED: Use exp_troni_KF_4 instead
 % Simulation of Giancarlo's algorithm for estimatating sensor bias without
 % knowledge of the vehicles true attitude (i.e using angular rates).  
  
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
t_end = 1000; % seconds
ts = 1/rate_d; %discretization time interval
w_max = [0.5,0.3,.5]; %set max ang. vel. for excitation

%set biases 
bias.ang = [3;2;1]*10^(-3);  %values chosen based on experiment @ 20Hz
bias.acc = [0;0;0];
bias.mag = [.1; 0.05; -.1];
% must be symmetric, also: assume that T is scaled such that T(1,1) = +1
Ts = [1,0.1,-0.15,1.1,0.05,0.81];
bias.T= [Ts(1:3)' [Ts(2) Ts(4) Ts(5)]' [Ts(3) Ts(5) Ts(6)]']; 

ms = gen_samples(lat, hz, t_end, bias, w_max);

% Read data
%ms = read_microstrain('/home/abhis/Matlab/DSCL/log/sim1.MST');
ms.t = ms.t-ms.t(1);
w = ms.ang(:,1:2)';
z = ms.mag(:,1:2)';
time = ms.t';
%phins = read_phins_imbat('/home/abhis/Matlab/DSCL/log/sim1.INS');
 
 %% System Set up
 runTime = t_end; % time in seconds
 %x = zeros(15, runTime*rate_d); %uncomment this if program runs slow
 x(:,1) = [z(:,1); 0;0;0; 1;0;0;1;0;1; 0;0;0];  %initialize state 
 Al = A_lin_15(1,ms,x(:,1)); %linearized discrete time A
 Ald = expm(Al*ts);
 Bld = B_dis(Al, ts); %linear discrete time B
 u(:,1) = f_x(1,ms,x(:,1)) - Al*x(:,1);
 C = [eye(3) zeros(3,12)];
 w_pos = zeros(size(w));
 
 % Kalman Filter Setup
 % Note that values for Q and R are dependent on the sampling rate.
 
% These values seem to work pretty well @ 20Hz, when using linearized prediction, choose slightly larger Q. 
% These values are for the discretized prediction model
% Decreasig Q results in les noisy final value, but takes longer to arrive
% at that value.
% Values that work pretyy well for process noise.
%  q1 = [1 1 1]*10^(-10);
%  q2 = [1 1 1]*2*10^(-10);
%  q3 = [0 1 1 1 1 1]*10^(-10); % q3(1) must be set to 0
%  q4 = [1 1 1]*10^(-12);
%  Q = diag([q1 q2 q3 q4]);
 
 q1 = [1 1 1]*10^(-10);
 q2 = [1 1 1]*2*10^(-10);
 q3 = [0 1 1 1 1 1]*10^(-10); % q3(1) must be set to 0
 q4 = [1 1 1]*5*10^(-12);
 Q = diag([q1 q2 q3 q4]);
 
 stdev = [1 1 1]*2*10^(-4);
 R = stdev.^2; 
 err_mag = stdev;
 err_magb = [1 1 1]*10^(-2);
 err_T = [0 1 1 1 1 1]*10^(-2);
 err_angb = [1 1 1]*10^(-3);
 err = [err_mag err_magb err_T err_angb]; %expected initial error - should the first 3 be the same as R?
 Sig = err'*err;
 
 %% Kalman Filter implematation
 tic
 
 for i = 2:(runTime*rate_d)  
   % prediction of next state: x(:,i) = Ald*x(:,i-1)+ Bld*u(:,t-1);
   % u turns out to be zero at all time (notes for details)
   % the system can therefore be simplified to
    x(:,i) = Ald*x(:,i-1) + Bld*u(:,i-1);  % could it be because the B term is missing again?
    Sig = Ald*Sig*Ald' + Q;
    
% find state using ode solver: slower, and oscillates a lot, WHY??
%      tspan = [time(i-1) time(i)];
%      [t, x_vec] = ode45(@(t,x_vec) sysS_15(t,x_vec,ms), tspan, x(:,i-1));
%      x(:,i) = x_vec(end,:)';  
     
    %kalman gain calculation 
     K = (Sig*C')/(C*Sig*C' + R);
     
    % update equations 
     x(:,i) = x(:,i) + K*(z(:,i) - C*x(:,i));
     Sig = (eye(15) - K*C)*Sig;
    
     %recompute A and B matrices
     Al  = A_lin_15(i,ms,x(:,i)); %linearized A 
     Ald = expm(Al*ts);
     Bld = B_dis(Al, ts); %linear discrete time B
     u(:,i) = f_x(i,ms,x(:,i)) - Al*x(:,i);
     
     % numericl integration to find angular position
     w_pos(:,i) = w_pos(:,i-1) + w(:,i-1)*ts;
  
     if ~mod(i,100)
        str = sprintf('Filter status: %i percent',floor((i*100)/(runTime*rate_d)));
        disp(str);
     end
 end
 toc
 
 
% Plotting
% Convert state in terms of error
x(4:6,:) = x(4:6,:) - ones(size(x(4:6,:))).*bias.mag;
x(7:12,:) = x(7:12,:) - ones(size(x(7:12,:))).*Ts';
x(13:15,:) = x(13:15,:) - ones(size(x(13:15,:))).*bias.ang;
%  figure(1)
%  hold on;
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

 figure 
 plot(ms.t,(w_pos).*180/pi);
 title 'Orientation of the vehicle';
 legend ('x','y','z');
 ylabel 'angular displacement (degrees)'; 
 xlabel 'time(s)';
 grid minor;


% figure
% plot(time(1:length(x)), x(4,:), time(1:length(x)), x(5,:), time(1:length(x)), x(6,:));
% title 'Error: Hard Iron Bias';
% ylabel 'magnetometer bias terms';
% legend ('x','y','z');
% xlabel 'time (s)';
% ylim([-0.1 0.1]);
% %  xlim([450 500]);
% title 'Biases';
% grid minor;
%   
% figure
% plot(time(1:length(x)), x(7:12,:));
% title 'Error: Soft-Iron Bias'
% ylabel 'elememts of symmetric T matrix'
% legend ('a','b','c','d','e','f');
% ylim([-0.1 0.1]);
% %xlim([450 500]);
% xlabel 'time (s)';
% grid minor;
 
figure
plot(time(1:length(x)), x(13,:), time(1:length(x)), x(14,:), time(1:length(x)), x(15,:));
title 'Error: Angular Velocity Bias';
legend ('x','y','z');
ylim([-0.01 0.01]*10^(-1));
title 'Ang Bias';
grid minor;
 
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

%% function to evaluate f_x 
 function f = f_x(t,samp,x)
  we = samp.ang(t,:)';
  wb = x(13:15,1);
  xm = x(1:3,1);
  b = x(4:6,1);
  Ts = x(7:12,1)';
  T = [Ts(1:3)' [Ts(2) Ts(4) Ts(5)]' [Ts(3) Ts(5) Ts(6)]'];

  f=  [(T*J(T\(xm-b)))*(we-wb);
              zeros(12,1)];
 end
 
 % function to compute the discrete time B matrix
 function B = B_dis(A,ts)
   Ald = expm(A*ts);
   arg_sum = zeros(15,15);
   step = ts/50; %going below 100 does add much value
   for i = 0:step:ts
       arg_sum = arg_sum + expm(-A*i)*step;
   end
   
   B = Ald*arg_sum;
 end
 