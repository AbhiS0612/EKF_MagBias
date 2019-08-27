%% Script to generate all the EKF Figures


%% Self Validation - Data Set B
load 'Workspaces/2019_08_27_B_MST.mat';
B.mb = b;
B.T = T;
B.wb = w_bias;
B.err = ans;

figure
subplot(3,1,1)
plot(filter_time, x(4:6,:));
title ('\bf{Magnetometer hard-iron bias}','FontSize', 10);
ylabel('$m_b$ (G)', 'FontSize', 11 );
legend ('x','y','z');
ylim([-0.1 0.3]);
xlim([0 3000]);
grid minor;

subplot(3,1,2)
plot(filter_time, x(7:12,:));
title ('\bf{Magnetometer soft-iron bias}','FontSize', 10);
ylabel ('$t_p$' ,'FontSize', 11);
legend ('a','b','c','d','e','f');
ylim([-0.2 1.7]);
xlim([0 3000]);
grid minor;
 
subplot(3,1,3)
plot(filter_time, (x(13:15,:))*(180/pi));
title ('\bf{Angular velocity bias}','FontSize', 10);
legend ('x','y','z');
ylabel ('$w_b$ ($^\circ$/s)','FontSize', 11);
xlabel 'time (s)';
ylim([-.8 .8]);
xlim([0 3000]);
grid minor;

% Plot Attitude Error
 mst.mag = (T\(ms.mag' - b'));
 mst.acc = ms.acc';
 mst.t = ms.t + start_time;
 att = calc_heading(mst.t, mst, phins);
 disp 'RMS Error: ';
 rms(att.att_error-mean(att.att_error,2),2)'


%% Self Validation - Data Set C
load 'Workspaces/2019_08_27_C_MST.mat';
Cc.mb = b;
Cc.T = T;
Cc.wb = w_bias;
Cc.err = ans;

figure
subplot(3,1,1)
plot(filter_time, x(4:6,:));
title ('\bf{Magnetometer hard-iron bias}','FontSize', 10);
ylabel('$m_b$ (G)', 'FontSize', 11 );
legend ('x','y','z');
ylim([-0.1 0.3]);
xlim([0 3000]);
grid minor;

subplot(3,1,2)
plot(filter_time, x(7:12,:));
title ('\bf{Magnetometer soft-iron bias}','FontSize', 10);
ylabel ('$t_p$' ,'FontSize', 11);
legend ('a','b','c','d','e','f');
ylim([-0.2 1.7]);
xlim([0 3000]);
grid minor;
 
subplot(3,1,3)
plot(filter_time, (x(13:15,:))*(180/pi));
title ('\bf{Angular velocity bias}','FontSize', 10);
legend ('x','y','z');
ylabel ('$w_b$ ($^\circ$/s)','FontSize', 11);
xlabel 'time (s)';
ylim([-.8 .8]);
xlim([0 3000]);
grid minor;

% Plot Attitude Error
 mst.mag = (T\(ms.mag' - b'));
 mst.acc = ms.acc';
 mst.t = ms.t + start_time;
 att = calc_heading(mst.t, mst, phins);
 disp 'RMS Error: ';
 rms(att.att_error-mean(att.att_error,2),2)'
 
 %% Cross Validation - Data Set C, biases from B
 mst.mag = (B.T\(ms.mag' - B.mb'));
 mst.acc = ms.acc';
 mst.t = ms.t + start_time;
 att = calc_heading(mst.t, mst, phins);
 disp 'RMS Error: ';
 rms(att.att_error-mean(att.att_error,2),2)'
 
 %% Cross Validation - Data Set B, biases from C
load 'Workspaces/2019_08_27_B_MST.mat';
 mst.mag = (Cc.T\(ms.mag' - Cc.mb'));
 mst.acc = ms.acc';
 mst.t = ms.t + start_time;
 att = calc_heading(mst.t, mst, phins);
 disp 'RMS Error: ';
 rms(att.att_error-mean(att.att_error,2),2)'
 
 %% Comparison Graph using data Set C
  %comment out attitude graph from calc_heading, add a figure number and add hold on;
 clear all;
 
 %Raw
 load 'Workspaces/2019_08_27_C_MST.mat';
 mst.mag = ms.mag';
 mst.acc = ms.acc';
 mst.t = ms.t + start_time;
 att = calc_heading(mst.t, mst, phins);
 disp 'RMS Error: ';
 rms(att.att_error-mean(att.att_error,2),2)'
 
 %Only Hard Iron
 mst.mag = ms.mag' - b';
 att = calc_heading(mst.t, mst, phins);
 disp 'RMS Error: ';
 rms(att.att_error-mean(att.att_error,2),2)'
 
 %Hard and Soft Iron
 mst.mag = (T\(ms.mag' - b'));
 att = calc_heading(mst.t, mst, phins);
 disp 'RMS Error: ';
 rms(att.att_error-mean(att.att_error,2),2)'
 
 legend ('Raw', 'Hard-iron', 'Hard and soft-iron');
