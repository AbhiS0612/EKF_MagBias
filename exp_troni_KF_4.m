%% Experimental Evalution of sim_troni_KF_4
 % Minor changes have been made from sim_troni_KF_4. 
 % 07-24-19:  CREATED by Abhi Shah

 %% Notes
  % The MST data doesnt come in at the sepcified rate - need to account for
  % that.
  % This code has been written specifically for the MST
  % 
%% Add paths
addpath '/home/abhis/Matlab/DSCL/dscl_matlab-master';
addpath '/home/abhis/Matlab/DSCL/andrew-matlab'; 
addpath(genpath_nosvn_nogit_nohg('/home/abhis/Matlab/DSCL/dscl_matlab-master'));
addpath(genpath_nosvn_nogit_nohg('/home/abhis/Matlab/DSCL/andrew-matlab'));
log_file_path = '/log/microstrain/2019_07_23_16_01.MST';

%% Read data
ms = read_microstrain(log_file_path);
%One might wish to downsample the data. Do that here, using downSample.m
ds_factor = 5; % downSample factor
ms = downSample(ms, ds_factor);
ms.t = ms.t-ms.t(1);
w = ms.ang';
z = ms.mag';
time = ms.t';
 
 %% System Set up
 data_recs = length(time); %number of data records
 runTime = time(end); %total runtime in seconds
 ts = data_recs/runTime;  %frequency of data recs (Hz) 
 %x = zeros(15, data_recs); %uncomment this if program runs slow
 x(:,1) = [z(:,1); 0;0;0; 1;0;0;1;0;1; 0;0;0];  %initialize state 
 Al = A_lin_15(1,ms,x(:,1)); %linearized discrete time A
 Ald = expm(Al*ts);
 Bld = B_dis(Al, ts); %linear discrete time B
 u(:,1) = f_x(1,ms,x(:,1)) - Al*x(:,1);
 C = [eye(3) zeros(3,12)];
 w_pos = zeros(size(w));
 
 % Kalman Filter Setup
 q1 = [1 1 1]*10^(-8);
 q2 = [1 1 1]*2*10^(-13);
 q3 = [0 1 1 1 1 1]*10^(-12); % q3(1) must be set to 0
 q4 = [1 1 1]*2*10^(-12);
 Q = diag([q1 q2 q3 q4]);
 
 stdev = [1 1 1]*2*10^(-4);
 R = stdev.^2; 
 err_mag = stdev;
 err_magb = [1 1 1]*10^(-1);
 err_T = [0 1 1 1 1 1]*10^(-1);
 err_angb = [1 1 1]*10^(-2);
 err = [err_mag err_magb err_T err_angb]; %expected initial error - should the first 3 be the same as R?
 Sig = err'*err;
 %% Kalman Filter implematation
 tic
 for i = 2:(data_recs)/2 
   % prediction of next state: x(:,i) = Ald*x(:,i-1)+ Bld*u(:,t-1);
   % this step can also be done using ode45. 
    x(:,i) = Ald*x(:,i-1) + Bld*u(:,i-1);
    Sig = Ald*Sig*Ald' + Q;
   
    %kalman gain calculation 
     disp (i);
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
        str = sprintf('Filter status: %i percent',floor((i*100)/(data_recs)));
        disp(str);
     end
 end
 toc
 
% %%  figure(1)
% %  hold on;
% % % subplot(2,1,1)
% %  plot(time(1:length(x)), x(1,:), time(1:length(x)), x(2,:), time(1:length(x)), x(3,:));
% %  legend ('x','y','z');
% %  title 'State estimate';
% %  grid minor;
% %  subplot(2,1,2)
% %  plot(time(1:length(x)),z(:,1:length(x)));
% %  legend('x','y','z');
% %  xlabel 'time (s)';
% %  title 'Mag Measurements';
% 
% figure 
% plot(ms.t,(w_pos).*180/pi);
% title 'Orientation of the vehicle';
% legend ('x','y','z');
% ylabel 'angular displacement (degrees)'; 
% xlabel 'time(s)';
% grid minor;
% 
% 
% figure
% plot(time(1:length(x)), x(4,:), time(1:length(x)), x(5,:), time(1:length(x)), x(6,:));
% title 'Error: Hard Iron Bias';
% ylabel 'bias terms';
% legend ('x','y','z');
% xlabel 'time (s)';
% %ylim([-0.1 0.1]);
% %xlim([400 500]);
% grid minor;
%   
% figure
% plot(time(1:length(x)), x(7:12,:));
% title 'Error: Soft-Iron Bias'
% ylabel 'elememts of symmetric T matrix'
% legend ('a','b','c','d','e','f');
% ylim([-0.1 0.1]);
% %xlim([400 500]);
% xlabel 'time (s)';
% grid minor;
%  
% figure
% plot(time(1:length(x)), x(13,:), time(1:length(x)), x(14,:), time(1:length(x)), x(15,:));
% title 'Error: Angular Velocity Bias';
% legend ('x','y','z');
% ylim([-0.01 0.01]);%*10^(-1));
% %xlim([400 500]);
% ylabel 'bias terms';
% xlabel 'time (s)';
% grid minor;
%  
% %  
% % %  %% Display Results
% %  disp 'Computed mean biases:';
% %  b_mean = mean(x(4:6,length(x)-100:length(x))')'
% %  disp 'True bias:';
% %  bias.mag
% %  disp 'Computed mean T matrix:';
% %  Ts = mean(x(7:12,floor(length(x)-100):length(x))');
% %  T_mean = [Ts(1:3)' [Ts(2) Ts(4) Ts(5)]' [Ts(3) Ts(5) Ts(6)]']
% %  disp 'True T matrix:';
% %  bias.T



%% Helper functions
% function to compute 
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
 