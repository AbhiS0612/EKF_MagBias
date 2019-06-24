function ms = read_microstrain(filename)
% MODIFICATION HISTORY
% 2018-08-05 LLW created and written
% 2018-08-27 LLW comments

fprintf(1,'Reading Microstrain IMU 3DM_GX5 ASCII data from from file %s \n', filename);

[   ms.string,...
    ms.date(:,1),...
    ms.date(:,2),...
    ms.date(:,3),...
    ms.time(:,1),...
    ms.time(:,2),...
    ms.time(:,3),...
    ms.rov_time,...
    ms.ros_time,...
    ms.ms_header_timestamp,...
    ms.ang(:,1),...  % rad/s
    ms.ang(:,2),...  % rad/s
    ms.ang(:,3),...  % rad/s
    ms.acc(:,1),...  % m/s^2
    ms.acc(:,2),...  % m/s^2
    ms.acc(:,3),...  % m/s^2
    ms.mag(:,1),...  % Gauss
    ms.mag(:,2),...  % Gauss
    ms.mag(:,3),...  % Gauss
    ms.pressure ] = textread(filename,'%3s %f/%f/%f %f:%f:%f %f %f %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f');

    % compute unix_time from ymdhms
    ms.unix_time = ymdhms_to_sec(ms.date(:,1), ms.date(:,2), ms.date(:,3),...
				 ms.time(:,1), ms.time(:,2), ms.time(:,3));

    % use rov_time as the best time available
    ms.t = ms.rov_time;

    % newsy
    dt = ms.t(end) - ms.t(1);
    n  = length(ms.t);
    hz = n/dt;
    fprintf(1,'Microstrain (MST) file %s has %f sec of data, %d data records, %f Hz\n', filename, dt, n, hz);

return;
