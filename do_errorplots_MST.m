%% Microstrain Error Testing
addpath '/home/abhis/Matlab/DSCL/dscl_matlab-master';
addpath(genpath_nosvn_nogit_nohg('/home/abhis/Matlab/DSCL/dscl_matlab-master'));
ms = read_microstrain('/log/microstrain/2019_06_25_17_26.MST');

%% Plot Histograms 
bins = 100;

% plot_title = 'MST Accelerometer Data from 2019\_06\_25\_14\_09.MST, sampled @ 100Hz';
% histplots(ms.acc, bins, plot_title);

plot_title = 'MST Angular Velocity Data from 2019\_06\_25\_14\_09.MST, sampled @ 100Hz';
histplots(ms.ang, bins, plot_title);

plot_title = 'MST Magnetometer Data from 2019\_06\_25\_14\_09.MST, sampled @ 100Hz';
histplots(ms.mag, bins, plot_title);

% %% Plot Raw Data
% ms.t = ms.t - ms.t(1);
% 
% figure
% subplot(311)
% title 'MST Data';
% plot(ms.t, ms.acc)
% ylabel 'Acceleration'
% legend ('a_x','a_y','a_z');
% 
% subplot(312)
% plot(ms.t, ms.ang)
% ylabel 'Angular Velocity';
% legend ('w_x', 'w_y','w_z');
% 
% subplot(313)
% plot(ms.t, ms.mag);
% ylabel 'Magnetometer'
% legend ('m_x', 'm_y','m_z');
% xlabel 'time (s)';
