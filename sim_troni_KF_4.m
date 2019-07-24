%% An Extention of sim_troni_3
 %This script expands the state space to include the angular velocity bias 
 % term (i.e 15 states) 
 %
 % Method: Extended Kalman Filter
 % 06-27-19:  CREATED by Abhi Shah
 
 %% Notes
  % Soft and hard iron effects for magnetometer.
  % should correct error covaiance for angular measurements.
  % 
  
%% Add paths
addpath '/home/abhis/Matlab/DSCL/dscl_matlab-master';
addpath '/home/abhis/Matlab/DSCL/andrew-matlab'; 
addpath(genpath_nosvn_nogit_nohg('/home/abhis/Matlab/DSCL/dscl_matlab-master'));
addpath(genpath_nosvn_nogit_nohg('/home/abhis/Matlab/DSCL/andrew-matlab'));
addpath(genpath_nosvn_nogit_nohg('/home/abhis/Matlab/DSCL/uco_sim'));
%ms = read_microstrain('/home/abhis/Matlab/DSCL/log/2019_06_12_18_26.MST');

%% Generate Simulated Data
lat = 39.33; % degrees
hz = 50;  % frequency of data generation  
rate_d = hz; % rate of discretization 
t_end = 500; % seconds
ts = 1/rate_d; %discretization time interval
w_max = [0.5,0.3,0.5]; %set max ang. vel. for excitation

%set biases 
bias.ang = [0.01;0.01;0.01];
bias.acc = [0;0;0];
bias.mag = [.1; 0.05; -.1];
% must be symmetric
T_bias=[ 1.1    0.0    0.0
         0.0    1.0    0.0
         0.0    0.0    1.0];

ms = gen_samples(lat, hz, t_end, bias, T_bias, w_max);

%% Read data
%ms = read_microstrain('/home/abhis/Matlab/DSCL/log/sim1.MST');
ms.t = ms.t-ms.t(1);
w = ms.ang';
z = ms.mag';
time = ms.t';
%phins = read_phins_imbat('/home/abhis/Matlab/DSCL/log/sim1.INS');
 
 %% System Set up
 runTime = t_end; % time in seconds
 x = zeros(15, runTime*rate_d); %uncomment this if program runs slow
 x(:,1) = [z(:,1); 0;0;0; 1;0;0;1;0;1; 0;0;0];  %initialize state 
 Al = A_lin_dis_15(1,ms,x(:,1),ts); %linearized discrete time A
 C = [  eye(3)      zeros(3,12)];
 w_pos = zeros(size(w));
 % Kalman Filter Setup
 % Note that values for Q and R are dependent on the sampling rate.
 
% These values seem to work pretty well @ 20Hz 
%  q1 = [1 1 1]*10^(-9);
%  q2 = [1 1 1]*10^(-9);
%  q3 = [1 1 1 1 1 1]*10^(-10);
%  Q = diag([q1 q2 q3]);
%  
%  stdev = [1 1 1]*3*10^(-4);
%  R = stdev.^2;
 
 q1 = [1 1 1]*10^(-9);
 q2 = [1 1 1]*10^(-9);
 q3 = [1 1 1 1 1 1]*10^(-10);
 q4 = [1 1 1]*10^(-9);
 Q = diag([q1 q2 q3 q4]);
 
 mag_std = [1 1 1]*2*10^(-4);
 R = diag(mag_std.^2);

 errT = ones(1,9)*0.1;
 err = [0.01 0.01 0.01 0.1 0.1 0.1 errT]; %expected initial error.
 Sig = err'*err;
 
 %% Kalman Filter implematation
 tic
 for t = 2:(runTime*rate_d)  
     x(:,t) = Al*x(:,t-1);
     
%     % find state using ode solver: slower, but slightly more accurate
%       [~, x_vec] = ode45(@(tx,x_vec) sysS(tx,x_vec,ms),((t-1):.2:t)/rate_d,x(:,t-1));
%       x(:,t) = x_vec(end,:)';
     
     Sig = Al*Sig*Al' + Q;
     K = (Sig*C')/(C*Sig*C' + R);
     
     x(:,t) = x(:,t) + K*(z(:,t) - C*x(:,t));
     Sig = (eye(15) - K*C)*Sig;
    
     %recompute A matrix
     Al = A_lin_dis_15(t,ms,x(:,t),ts);
     
     % numericl integration to find angular position
     w_pos(:,t) = w_pos(:,t-1) + w(:,t-1)*ts;
  
     if ~mod(t,100)
        str = sprintf('Filter status: %i percent',floor((t*100)/(runTime*rate_d)));
        disp(str);
     end
 end
 toc
 
 
 %% Plotting
%  figure
%  subplot(2,1,1)
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

 figure
 plot(time(1:length(x)), x(4,:), time(1:length(x)), x(5,:), time(1:length(x)), x(6,:));
 title 'Hard-Iron Effects';
 ylabel 'magnetometer bias terms';
 ylim([-0.2 0.2]);
 legend ('x','y','z');
 xlabel 'time (s)';
 title 'Biases';
 grid minor;
  
 figure
 plot(time(1:length(x)), x(7:12,:));
 title 'Soft-Iron Effects'
 ylabel 'elememts of symmetric T matrix'
 legend ('a','b','c','d','e','f');
 xlabel 'time (s)';
 grid minor;
 
 figure
 plot(time(1:length(x)), x(13,:), time(1:length(x)), x(14,:), time(1:length(x)), x(15,:));
 title 'Angular Velocity Bias';
 ylabel 'bias terms';
 ylim([-0.2 0.2]);
 legend ('x','y','z');
 xlabel 'time (s)';
 title 'Biases';
 grid minor;
 
 
%  %% Display Results
%  disp 'Computed mean biases:';
%  b_mean = mean(x(4:6,length(x)-100:length(x))')'
%  disp 'True bias:';
%  bias.mag
%  disp 'Computed mean T matrix:';
%  Ts = mean(x(7:12,floor(length(x)-100):length(x))');
%  T_mean = [Ts(1:3)' [Ts(2) Ts(4) Ts(5)]' [Ts(3) Ts(5) Ts(6)]']
%  disp 'True T matrix:';
%  T_bias
 
 