%% Simulation of Giancarlo's algorithm for estimatating sensor bias without
 % knowledge of the vehicles true attitude (i.e using angular rates).
 % Method: Kalman Filter
 % 06-13-19:  CREATED by Abhi Shah
 
 %% Notes
  % Only does hard iron effects for magnetometer.
  % should correct error covaiance for angular measurements and try this
  % this again.
  
  
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
ts = 1/rate_d; 

t_end = 300; % seconds
bias.ang = [0;0;0];
bias.acc = [0;0;0];
bias.mag = [.1; 0.05; -.1];

ms = gen_samples(lat, hz, t_end, bias);
%ms = read_microstrain('/home/abhis/Matlab/DSCL/log/sim1.MST');
w = ms.ang';
z = ms.mag';
time = ms.t';
%phins = read_phins_imbat('/home/abhis/Matlab/DSCL/log/sim1.INS');


% %% For using experimental data
%  rate_orig = 1000; %Hz
%  rate_d = 10; % Hz 
%  rate_req = rate_d;
%  [time, w, z'] = downSample(ms.t, ms.ang, ms.mag, rate_orig, rate_req);
 
 %% System Set up
 Ac = [  -J(w(:,1))   J(w(:,1))
          zeros(3)    zeros(3)];
 A = expm(Ac*ts);
 C = [eye(3) zeros(3)];
 % Kalman Filter Setup
 Q = eye(6)*0.01;
 R = eye(3)*0.001; % how to choose this, and also the Q
 err = [0.001 0.001 0.001 0.3 0.3 0.3]; %expected initial error.
 Sig = err'*err;
 
 %% Kalman Filter implematation
 runTime = t_end; % time in seconds
 x = zeros(6, runTime*rate_d);
 x(:,1) = [z(:,1); 0;0;0];  %initialize state 
 s(1) = norm(Sig);
 for t = 2:(runTime*rate_d)
    
     x(:,t) = A*x(:,t-1);
     Sig = A*Sig*A' + Q;
     K = (Sig*C')/(C*Sig*C' + R);
     
     x(:,t) = x(:,t) + K*(z(:,t) - C*x(:,t));
     Sig = (eye(6) - K*C)*Sig;
     %s(t) = norm(Sig);
     
     %recompute A matrix
     Ac = [-J(w(:,t))    J(w(:,t))
           zeros(3)    zeros(3)];
     A = expm(Ac*ts);
     
     if ~mod(t,100)
        str = sprintf('Filter status: %i percent',floor((t*100)/(runTime*rate_d)));
        disp(str);
     end
 end
 
 %% Plotting
%  figure
%  subplot(2,1,1)
%  plot(time(1:length(x)), x(1,:), time(1:length(x)), x(2,:), time(1:length(x)), x(3,:));
%  legend ('x','y','z');
%  grid minor;
%  subplot(2,1,2)
%  plot(time,z);
%  legend('x','y','z');
%  title 'Mag stuff';
%  grid minor;
 
 figure
 plot(time(1:length(x)), x(4,:), time(1:length(x)), x(5,:), time(1:length(x)), x(6,:));
 legend ('x','y','z');
 title 'Biases';
 hold on;
 
 %figure(3)
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
 