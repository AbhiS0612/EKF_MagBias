function kvh = read_kvh(filename)
% MODIFICATION HISTORY
% 2018-08-04 LLW created and written
% 2018-08-27 LLW document units


fprintf(1,'Reading KVH 1775 ASCII data from from file %s \n', filename);

[   kvh.string,...
    kvh.date(:,1),...
    kvh.date(:,2),...
    kvh.date(:,3),...
    kvh.time(:,1),...
    kvh.time(:,2),...
    kvh.time(:,3),...
    kvh.rov_time,...
    kvh.ros_time,...
    kvh.ang(:,1),...  % rad/s
    kvh.ang(:,2),...  % rad/s
    kvh.ang(:,3),...  % rad/s
    kvh.acc(:,1),...  % g - will be converted to m/s^2
    kvh.acc(:,2),...  % g - will be converted to m/s^2
    kvh.acc(:,3),...  % g - will be converted to m/s^2
    kvh.mag(:,1),...  % Gauss
    kvh.mag(:,2),...  % Gauss
    kvh.mag(:,3),...  % Gauss
    kvh.temp,...      % C
    kvh.seq_num,...
    kvh.timestamp,...
    kvh.comp_timestamp,...
    kvh.status(:,1),...
    kvh.status(:,2),...
    kvh.status(:,3),...
    kvh.status(:,4),...
    kvh.status(:,5),...
    kvh.status(:,6) ] = textread(filename,'%3s %f/%f/%f %f:%f:%f %f %f %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%f,%f,%d,%d,%d,%d,%d,%d');

    kvh.unix_time = ymdhms_to_sec(kvh.date(:,1), kvh.date(:,2), kvh.date(:,3),...
    kvh.time(:,1), kvh.time(:,2), kvh.time(:,3));

    % scale acc to m/s^2
    kvh.acc = kvh.acc * 9.81;

    kvh.t = kvh.rov_time;

    dt = kvh.t(end) - kvh.t(1);
    n  = length(kvh.t);
    hz = n/dt;
    fprintf(1,'KVH ASCII file %s has %f sec of data, %d data records, %f Hz\n', filename, dt, n, hz);    
    
