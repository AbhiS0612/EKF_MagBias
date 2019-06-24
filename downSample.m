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