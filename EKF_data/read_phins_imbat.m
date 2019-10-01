function phins = read_phins_imbat(filename)
% MODIFICATION HISTORY
% 2018-08-05 LLW Created and written
% 2018-08-27 LLW Revised convert accels from g to m/s^2 

fprintf(1,'Reading Phins IMBAT ASCII data from file %s \n', filename);

[   phins.string,...
    phins.date(:,1),...
    phins.date(:,2),...
    phins.date(:,3),...
    phins.time(:,1),...
    phins.time(:,2),...
    phins.time(:,3),...
    phins.rov_time,...
    phins.ros_time,...
    phins.phins_data_timestamp,...
    phins.ang(:,1),...  % deg/s - will be converted to rad/s
    phins.ang(:,2),...  % deg/s - will be converted to rad/s
    phins.ang(:,3),...  % deg/s - will be converted to rad/s
    phins.acc(:,1),...  % g  - will be converted to m/s^2
    phins.acc(:,2),...  % g  - will be converted to m/s^2 
    phins.acc(:,3),...  % g  - will be converted to m/s^2
    phins.att(:,1),...  % roll  in deg
    phins.att(:,2),...  % pitch in deg
    phins.att(:,3),...  % yaw   in deg
    phins.heave,...     % m
    phins.status] = textread(filename,'%3s %f/%f/%f %f:%f:%f %f %f %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f');

    phins.unix_time = ymdhms_to_sec(phins.date(:,1), phins.date(:,2), phins.date(:,3),...
				 phins.time(:,1), phins.time(:,2), phins.time(:,3));

    phins.t = phins.rov_time;

    % convert phins angular velocities from deg/sec to rad/sec
    phins.ang = phins.ang *(pi/180.0);

    % convert phins acc from g to m/s^2
    phins.acc = phins.acc * 9.81;

    dt = phins.t(end) - phins.t(1);
    n  = length(phins.t);
    hz = n/dt;
    fprintf(1,'Phins IMBAT ASCII file %s has %f sec of data, %d data records, %f Hz\n', filename, dt, n, hz);    

return    

