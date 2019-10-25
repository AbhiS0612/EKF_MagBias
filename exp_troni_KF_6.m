%% Experimental Evalution of sim_troni_KF_4
%  08-23-19:  CREATED by Abhi Shah
%% Notes
 % Trying to make a 17 states version that exclude T(1,1) from the state
 % entirely.

%% Add paths 
addpath '/home/abhis/Matlab/DSCL/dscl_matlab-master';
addpath '/home/abhis/Matlab/DSCL/andrew-matlab'; 
addpath(genpath_nosvn_nogit_nohg('/home/abhis/Matlab/DSCL/dscl_matlab-master'));
addpath(genpath_nosvn_nogit_nohg('/home/abhis/Matlab/DSCL/andrew-matlab'));
% log_file_path = '/log/microstrain/2019_08_09_18-19.MST';
% log_file_path = '/log/kvh/2019_08_02_16_27.KVH';
% log_file_path_p = '/log/phins/2019_08_09_18-19.INS';

%% Generate Simulated Data (for comparison / testing)
clear all;
lat = 39.33; % degrees
hz = 20.1;  % frequency of data generation (be sure to choose appropriate noise specs)  
t_end = 2000; % seconds
w_max = [0,0.2,.2]; %for 17 states to work do you need motion about all 3 axes??

% set biases 
bias.ang = [0;2;1]*10^(-3); 
bias.acc = [0;0;0];
bias.mag = [0; 0.2; -0.2];
% must be symmetric, also: assume that T is scaled such that T(1,1) = +1
Ts = [1,0.01,0.03,0.95,0.01,1.4];

%bias.T= [Ts(1:3)' [Ts(2) Ts(4) Ts(5)]' [Ts(3) Ts(5) Ts(6)]']; 
bias.T = diag([1 1.2 0.6]);
%bias.T = eye(3);
ms = gen_samples(lat, hz, t_end, bias, w_max);

%% Read data
% ms = read_microstrain(log_file_path);
% phins = read_phins_imbat(log_file_path_p);
% 
% % Post Processing (optional)
%  ms1 = read_kvh(log_file_path);
%  ds_factor = 50; % downSample factor
%  ms1 = downSample(ms1, ds_factor);

% %  trunate data range (for large data samples)
% sp = 5;  %start percentage
% ep = 80;  %end percentage
% n = length(ms1.t);
% ms.mag = ms1.mag(floor(sp*n/100):floor(ep*n/100),:);
% ms.ang = ms1.ang(floor(sp*n/100):floor(ep*n/100),:);
% ms.t = ms1.t(floor(sp*n/100):floor(ep*n/100),:);

start_time = ms.t(1);  %record unix start time
ms.t = ms.t-ms.t(1);  %start time set to 0
w = ms.ang';  % ang vel. measurement
z = ms.mag';  % mag measurement
time = ms.t'; % time measurement

%% System and Filter Set up
 filter_rate = 10;   % (Hz)  
 ts = 1/filter_rate; % (s)
 filter_time = 0:ts:floor(time(end))-ts;  
 filter_points = filter_rate * floor(time(end)); %number of data records
 ft = filter_time(1);
 x(:,1) = [z(:,1); 0;0;0; 0;0;1;0;1; 0;0;0; w(:,1)];  %initialize state 
 Al = A_lin_17_exp(ft,ms,x(:,1)); %linearized A
 Ald = expm(Al*ts);  %discretized Al
 Bld = B_dis(Al, ts); %linear discrete time B
 u(:,1) = f_x(ft,ms,x(:,1)) - Al*x(:,1);  %psedo-control input
 C = [eye(3) zeros(3,14)
      zeros(3,14) eye(3)];
 
% KVH @20hz
 q1 = [1 1 1]*10^(-9);
 q2 = [1 1 1]*10^(-10);
 q3 = [1 1 1 1 1]*10^(-10);
 q4 = [1 1 1]*10^(-13);
 q5 = [1 1 1]*10^(-6);
 Q = diag([q1 q2 q3 q4 q5]);
 
% Noise Covariance Matrix
 stdev_mag = [1.7 1.7 3.0]*10^(-4);  %MST @20hz
 stdev_ang = [1 1 1]*3*10^(-5);  %MST @20hz
 R = diag([stdev_mag.^2 stdev_ang.^2]);
 
% Initial Error Covariance Matrix 
 err_mag = stdev_mag;
 err_magb = [1 1 1]*10^(-4);
 err_T = [1 1 1 1 1]*10^(-4);
 err_angb = [1 1 1]*10^(-4);
 err_ang = stdev_ang;
 err = [err_mag err_magb err_T err_angb err_ang];
 Sig = err'*err;
 
 %starting with zero error covariance
 %Sig = zeros(17);
 s(:,1) = diag(Sig);
%%  Kalman Filter implematation
 tic
 for i = 2:filter_points
    %prediction of next state: x(:,i) = Ald*x(:,i-1)+ Bld*u(:,i-1);
    x(:,i) = Ald*x(:,i-1) + Bld*u(:,i-1);
    Sig = Ald*Sig*Ald' + Q;
    
%   find state using ode solver: way slower, and oscillates a lot, WHY??
%      tspan = [time(i-1) time(i)];
%      [t, x_vec] = ode45(@(t,x_vec) sysS_15(t,x_vec,ms), tspan, x(:,i-1));
%      x(:,i) = x_vec(end,:)';  
     
    %kalman gain calculation 
     K = (Sig*C')/(C*Sig*C' + R);
     
    %update equations 
     %get interpolated measurement at current time and mag measurement
     ft = filter_time(i);
     z_m = interp1(time, z', ft)';
     w_m = interp1(time, w', ft)';
     zw_m = [z_m; w_m];
     x(:,i) = x(:,i) + K*(zw_m - C*x(:,i));
     Sig = (eye(17) - K*C)*Sig;
     s(:,i) = diag(Sig);
     %recompute A, B and pseudo-control input
     Al  = A_lin_17_exp(ft,ms,x(:,i)); %linearized A 
     Ald = expm(Al*ts); %discretized Al
     Bld = B_dis(Al, ts); %linear discrete time B
     u(:,i) = f_x(ft,ms,x(:,i)) - Al*x(:,i);
  
     %Display status of the filter 
     if ~mod(i,100)
        str = sprintf('Filter status: %i percent',floor((i*100)/(filter_points)));
        disp(str);
     end
 end
 toc

%% Convert state in terms of error
% x(4:6,:) = x(4:6,:) - ones(size(x(4:6,:))).*bias.mag;
% x(7:11,:) = x(7:11,:) - ones(size(x(7:11,:))).*Ts';
% x(12:14,:) = x(12:14,:) - ones(size(x(12:14,:))).*bias.ang;

%% Plot Estimated Biases
figure
plot(filter_time, x(1:3,:), ms.t, ms.mag');
title ('mag state v/s measurement');
ylabel 'gauss';
legend ('x','y','z');
%ylim([-0.1 0.3]);
%xlim([0 2500]);
grid minor;

figure;
plot(filter_time, s');
title 'Diagonal Elements of Sigma';
xlabel 'time(s)';
legend;
grid minor;

figure
subplot(3,1,1)
plot(filter_time, x(4:6,:));
title 'm_b';
ylabel 'gauss';
legend ('x','y','z');
%ylim([-0.1 0.3]);
%xlim([0 2500]);
grid minor;

subplot(3,1,2)
plot(filter_time, x(7:11,:));
title 't_p'
ylabel '';
legend ('b','c','d','e','f');
%ylim([-0.2 1.6]);
%xlim([0 2500]);
grid minor;
 
subplot(3,1,3)
plot(filter_time, (x(12:14,:))*(180/pi));
title 'w_b';
legend ('x','y','z');
ylabel 'deg/s';
xlabel 'time (s)';
%ylim([-0.5 0.5]);
%xlim([0 2500]);
grid minor;

%% Print Results
%  disp 'Computed biases:';
%  b = mean(x(4:6,length(x)-100:length(x))')
%  disp 'Computed T matrix:';
%  Ts = mean(x(7:12,floor(length(x)-100):length(x))');
%  T = [Ts(1:3)' [Ts(2) Ts(4) Ts(5)]' [Ts(3) Ts(5) Ts(6)]']
%  wm = mean(x(13:15,length(x)-200:length(x))')
 
%% Plot Attitude Error
%  mst.mag = (T\(ms.mag' - b'));
%  mst.acc = ms.acc';
%  mst.t = ms.t + start_time;
%  att = calc_heading(mst.t, mst, phins);
%  disp 'RMS Error: ';
%  rms(att.att_error-mean(att.att_error,2),2)'

%% Helper functions (could put these in seperate files)
%function to evatuate f(x) at the current state
 function f = f_x(t,samp,x)
  %wm = interp1(samp.t,samp.ang,t)';
  wb = x(12:14,1);
  xm = x(1:3,1);
  b = x(4:6,1);
  Ts = [1 x(7:11,1)'];
  wm = x(15:17,1);
  T = [Ts(1:3)' [Ts(2) Ts(4) Ts(5)]' [Ts(3) Ts(5) Ts(6)]'];

  f=  [(T*J(T\(xm-b)))*(wm-wb);
              zeros(14,1)];
 end
 
 %function to compute the discrete time B matrix
 function B = B_dis(A,ts)
   Ald = expm(A*ts);
   arg_sum = zeros(17);
   step = ts/50; %going below 100 does add any significant value
   for i = 0:step:ts
       arg_sum = arg_sum + expm(-A*i)*step;
   end
   B = Ald*arg_sum;
 end
 