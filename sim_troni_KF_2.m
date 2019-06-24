%% Simulation of Giancarlo's algorithm for estimatating sensor bias without
 % knowledge of the vehicles true attitude (i.e using angular rates).
 % Method: Kalman Filter (hard iron and soft-iron effects
 % 06-13-19:  CREATED by Abhi Shah
 
 %% Notes
  % Attempt at both hard iron and soft iron effects for magnetometer.
  
  
%% Add paths
addpath '/home/abhis/Matlab/DSCL/dscl_matlab-master';
addpath '/home/abhis/Matlab/DSCL/andrew-matlab'; 
addpath(genpath_nosvn_nogit_nohg('/home/abhis/Matlab/DSCL/dscl_matlab-master'));
addpath(genpath_nosvn_nogit_nohg('/home/abhis/Matlab/DSCL/andrew-matlab'));
%  ms = read_microstrain('/home/abhis/Matlab/DSCL/log/2019_06_12_18_26.MST');

%% Generate Simulated Data
lat = 39.33; % degrees
hz = 20;  % frequency of data generation  
rate_d = hz; % rate of discretization 

t_end = 300; % seconds
bias.ang = [0;0;0];
bias.acc = [0;0;0];
bias.mag = [.1; 0.05; -.1];

%ms = gen_samples(lat, hz, t_end, bias);
%ms = read_microstrain('/home/abhis/Matlab/DSCL/log/sim1.MST');
w = ms.ang';
z = ms.mag;
time = ms.t;
%phins = read_phins_imbat('/home/abhis/Matlab/DSCL/log/sim1.INS');

% %% For using experimental data
%  rate_orig = 1000; %Hz
%  rate_d = 10; % Hz 
%  rate_req = rate_d;
%  [time, w, z'] = downSample(ms.t, ms.ang, ms.mag, rate_orig, rate_req);
 
 %% System Set up
 T = eye(3);
 A = [-T*J(w(1,:))/T   T*J(w(1,:))/T   zeros(3,6);
        zeros(9,12)];
 B = zeros(12,1);
 C = [eye(3) zeros(3,9)];
 D = zeros(3,1);
 sys = ss(A,B,C,D);
 ts = 1/rate_d;
 sysd = c2d(sys, ts);
 A = sysd.A;
 
 % Kalman Filter Setup
 q1 = [1 1 1]*0.01;
 q2 = [1 1 1]*0.01;
 q3 = [1 1 1 1 1 1]*0.00;
 Q = diag([q1 q2 q3]);
 R = eye(3)*0.001; % how to choose this, and also the Q
 errT = ones(1,6)*0.00;
 err = [0.001 0.001 0.001 0.3 0.3 0.3 errT]; %expected initial error.
 Sig = err'*err;
 
 %% Kalman Filter implematation
 runTime = t_end; % time in seconds
 x = zeros(12, runTime*rate_d);
 x(:,1) = [z(:,1); zeros(3,1); [1 0 0 1 0 1]'];  %initialize state 
 s(1) = norm(Sig);
 for t = 2:(runTime*rate_d)
    
     x(:,t) = A*x(:,t-1);
     Sig = A*Sig*A' + Q;
     K = (Sig*C')/(C*Sig*C' + R);
     
     x(:,t) = x(:,t) + K*(z(:,t) - C*x(:,t));
     Sig = (eye(12) - K*C)*Sig;
    % s(t) = norm(Sig);
     
    % Recompute T 
     xt = x(7:12,t)';
     T = [xt(1:3)' [xt(2) xt(4) xt(5)]' [xt(3) xt(5) xt(6)]']; 
     
     %recompute A matrix
    A = [-J(w(1,:))   J(w(1,:))   zeros(3,6);
        zeros(9,12)];
    sys = ss(A,B,C,D);
    sysd = c2d(sys, ts);
    A = sysd.A;
     
     if ~mod(t,100)
        str = sprintf('Filter status: %i percent',floor((t*100)/(runTime*rate_d)));
        disp(str);
     end
 end
 
 %% Plotting
%  figure(1)
%  subplot(2,1,1)
%  plot(time(1:length(x)), x(1,:), time(1:length(x)), x(2,:), time(1:length(x)), x(3,:));
%  legend ('x','y','z');
%  grid minor;
%  subplot(2,1,2)
%  plot(time,z);
%  legend('x','y','z');
%  title 'Mag stuff';
%  grid minor;
 
 figure(2)
 plot(time(1:length(x)), x(4,:), time(1:length(x)), x(5,:), time(1:length(x)), x(6,:));
 legend ('x','y','z');
 title 'Biases';
 
 figure(3)
 plot(time(1:length(x)), x(7:12,:));
 legend ('a','b','c','d','e','f');
 
 %figure(4)
 %plot(time(1:length(s)), s);
 %% Hat operation
 function h = J(w)
 h = [0  -w(3)  w(2)
      w(3)  0  -w(1)
      -w(2) w(1)  0];
 end
 
 %% Down sample readings; quick bad method
 function [t w z] = downSample(t0, w0, z0, r0, r1)
 step = floor(r0/r1);
 j = 1;
 sum_t = 0;
 sum_w = 0;
 sum_z = 0;
 
 for i = 1:length(t0)
        sum_t = sum_t + t0(i,:);
        sum_w = sum_w + w0(i,:);
        sum_z = sum_z + z0(i,:);
        
        if(mod(i,step) == 0)
          t(j,:) = sum_t/step;
          w(j,:) = sum_w/step;
          z(j,:) = sum_z/step;
          j = j+1;
          sum_t = 0;
          sum_w = 0;
          sum_z = 0;
        end
 end
 end
 