%% Down sample data from microstrain or kvh by averaging. 
% Since the downsampled rate must be a diviser of the original sampling rate,
% this function simply requires the downsampling factor, rather than the
% required rate.
% the function returns a struct, with average time, acc, ang, mag values
% function struct1 = downSample(struct0, req_factor)
%  Input
%   struct0:    struct of microstrain or kvh data
%   req_factor: required averaging factor (must be a positive integer)
%  Output
%   struct1:    struct with fields t, ang, acc, mag containing the
%   downsampled data
%
% e.g. if the data in struct0 is at 99.5hz, and req_factor is 10, stuct1
% will contain data at 9.95hz

 function st1 = downSample(st0, rf)
 k = 1;
  for i = 1:rf:(length(st0.t)-rf)
    st1.t(k,:) = mean(st0.t(i:(i+rf-1)));
    st1.acc(k,:) = mean(st0.acc(i:(i+rf-1),:));
    st1.ang(k,:) = mean(st0.ang(i:(i+rf-1),:));
    st1.mag(k,:) = mean(st0.mag(i:(i+rf-1),:));
    k = k+1;
  end
  
  % Print info
    n  = (length(st0.t) - mod(length(st0.t),rf));
    dt = st0.t(n) - st0.t(1);
    n = n/rf;
    hz = n/dt;
    fprintf(1,' Down sampled data has %f sec of data, %d data records, %f Hz\n', dt, n, hz);
 end