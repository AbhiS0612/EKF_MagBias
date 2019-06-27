%% Microstrain Error Testing
addpath '/home/abhis/Matlab/DSCL/dscl_matlab-master';
addpath(genpath_nosvn_nogit_nohg('/home/abhis/Matlab/DSCL/dscl_matlab-master'));
log_file_path = '/log/microstrain/2019_06_27_15_06.MST';
ms = read_microstrain(log_file_path);
%ms = read_microstrain('/home/abhis/Matlab/DSCL/log/sim1.MST');

%% Plot Histograms 
i = 7;
bins = 30;
log_file = '2019\_06\_27\_15\_06.MST';
secs = 104.915; % duration of trial
Hz = 1000.255; % run section above to print precise Hz

plot_title = sprintf('MST Accelerometer Data from %s sampled at %.3fHz (%.3fsecs)', log_file, Hz, secs);
spec_acc{i} = histplots(ms.acc, bins, plot_title);
print('-fillpage','hist_acc_1000','-dpdf')

plot_title = sprintf('MST Angular Velocity Data from %s sampled at %.3fHz (%.3fsecs)', log_file, Hz, secs);
spec_ang{i} = histplots(ms.ang, bins, plot_title);
print('-fillpage','hist_ang_1000','-dpdf')

plot_title = sprintf('MST Magnetometer Data from %s sampled at %.3fHz (%.3fsecs)', log_file, Hz, secs);
spec_mag{i} = histplots(ms.mag, bins, plot_title);
print('-fillpage','hist_mag_1000','-dpdf')

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
