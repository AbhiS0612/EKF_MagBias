%% Simulation of Giancarlo's algorithm for estimatating sensor bias without
 % knowledge of the vehicles true attitude (i.e using angular rates).
 % Method: Kalman Filter
 % 06-13-19:  CREATED by Abhi Shah
 
 %% Notes
  % Soft and hard iron effects for magnetometer.
  % should correct error covaiance for angular measurements.
  % System diverges if T is is slightly different from the identity
  
%% Add paths
addpath '/home/abhis/Matlab/DSCL/dscl_matlab-master';
addpath '/home/abhis/Matlab/DSCL/andrew-matlab'; 
addpath(genpath_nosvn_nogit_nohg('/home/abhis/Matlab/DSCL/dscl_matlab-master'));
addpath(genpath_nosvn_nogit_nohg('/home/abhis/Matlab/DSCL/andrew-matlab'));
addpath(genpath_nosvn_nogit_nohg('/home/abhis/Matlab/DSCL/uco_sim'));
%  ms = read_microstrain('/home/abhis/Matlab/DSCL/log/2019_06_12_18_26.MST');

%% Generate Simulated Data
lat = 39.33; % degrees
hz = 50;  % frequency of data generation  
rate_d = hz; % rate of discretization 
t_end = 100; % seconds
ts = 1/rate_d; %discretization time interval

%set biases 
bias.ang = [0;0;0];
bias.acc = [0;0;0];
bias.mag = [.1; 0.05; -.1];
%T_bias=[.95,0,0; 0,1.1,0;0,0,1.05];
T_bias=[.9,0,0; 0,1.1,0;0,0,1];
%T_bias=[1,0,0; 0,1,0; 0,0,1];

ms = gen_samples(lat, hz, t_end, bias, T_bias);

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
 Al = A_lin_dis(1,ms,x(:,1),ts); %linearized discrete time A
 C = [eye(3) zeros(3,9)];
 
 % Kalman Filter Setup
 q1 = [1 1 1]*0.0001;
 q2 = [1 1 1]*0.00001;
 q3 = [1 1 1 1 1 1]*0.00001;
 
 Q = diag([q1 q2 q3]);
 R = eye(3)*0.001;
 
 errT = ones(1,6)*0.01;
 err = [0.01 0.01 0.01 0.1 0.1 0.1 errT]; %expected initial error.
 Sig = err'*err;
 
 %% Kalman Filter implematation
 s(1) = norm(Sig);
 tic
 for t = 2:(runTime*rate_d)  
     %x(:,t) = Al*x(:,t-1);
     
     % find state using ode solver: faster way to do this?
     [~, x_vec] = ode45(@(tx,x_vec) sysS(tx,x_vec,ms),((t-1):.2:t)/rate_d,x(:,t-1));
     x(:,t) = x_vec(end,:)';
     
     Sig = Al*Sig*Al' + Q;
     K = (Sig*C')/(C*Sig*C' + R);
     
     x(:,t) = x(:,t) + K*(z(:,t) - C*x(:,t));
     Sig = (eye(12) - K*C)*Sig;
     s(t) = norm(Sig);
    
     %recompute A matrix
     Al = A_lin_dis(t,ms,x(:,t),ts);
  
     if ~mod(t,100)
        str = sprintf('Filter status: %i percent',floor((t*100)/(runTime*rate_d)));
        disp(str);
     end
 end
 toc
 
 %% Plotting
 figure
 subplot(2,1,1)
 plot(time(1:length(x)), x(1,:), time(1:length(x)), x(2,:), time(1:length(x)), x(3,:));
 legend ('x','y','z');
 title 'State estimate';
 grid minor;
 subplot(2,1,2)
 plot(time(1:length(x)),z(:,1:length(x)));
 legend('x','y','z');
 xlabel 'time (s)';
 title 'Mag Measurements';
 grid minor;
 
 figure
 plot(time(1:length(x)), x(4,:), time(1:length(x)), x(5,:), time(1:length(x)), x(6,:));
 legend ('x','y','z');
 title 'Biases';
 hold on;
  
 figure
 plot(time(1:length(x)), x(7:12,:));
 legend ('a','b','c','d','e','f');
 
 %% Display Results
 disp 'Computed mean biases:';
 b_mean = mean(x(4:6,length(x)-1000:length(x))')'
 disp 'True bias:';
 bias.mag
 disp 'Computed mean T matrix:';
 Ts = mean(x(7:12,floor(.5*length(x)):length(x))');
 T_mean = [Ts(1:3)' [Ts(2) Ts(4) Ts(5)]' [Ts(3) Ts(5) Ts(6)]']
 disp 'True T matrix:';
 T_bias
 
 