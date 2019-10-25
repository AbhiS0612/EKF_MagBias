function newx = resample2(t,x,newt,interp_method)
%
% function newx = resample(t,x,newt,interp_method)
% 
% given a vector pair t and x containing data points x(i) at times t(i)
% and vector newt containing a new set of sample times newt(i), returns
% a resampled vector newx of interpolatated points newx(i) at newt(i)
% 
% if newt range exceeds t range, pads with x(i) value closest in time
% 
% t         time sequence of original data points
% x         original data points           
% newt      time sequence for new interpolated points
%
% interp_method = 'linear'  works as above
% interp_method = 'nearest' works by matching closest time
%
%
% June  1998  Louis Whitcomb  Created and written
% Oct 2 1998  Louis Whitcomb  Added nearest time match method

% default to linear interpolation
if ~exist('interp_method')
  interp_method = 'linear';
end

% msgtext = ['RESAMPLE.M: interpolating method is --> ' interp_method ' <--\n'];
fprintf(1, ['RESAMPLE.M: interpolating method is --> ' interp_method ' <--\n']);

% find time range of original data
rti   = t(1);
rtf   = t(length(t));

% find time range of new time sequence
dti   = newt(1);
dtf   = newt(length(newt));

if( (rti > dti ) | (rtf < dtf))
   fprintf(1,'RESAMPLE.M: original  data from t=%f to t=%f \n',rti,rtf);
   fprintf(1,'RESAMPLE.M: resampled data from t=%f to t=%f \n',dti,dtf);
   fprintf(1,'RESAMPLE.M: Padding missing X data points with adjacent values.\n');
   if( (abs(rti - dti) < 1) & (abs(rtf - dtf) < 1))
       fprintf(1,'RESAMPLE.M: Missing less than 1 second of X data, so no problem.\n');    
   else
       fprintf(1,'RESAMPLE.M: MISSING MORE THAN 1 SECOND OF X data!!!! PROBLEM!!!!.\n');    
   end
   if( rti > dti ) 
     t   = [ dti;     t  ];
     x   = [ x(1,:); x];
   end
   if( rtf < dtf ) 
     t   = [ t  ; dtf];
     x   = [ x; x(length(x),:)];
   end
end

if 1 == strcmp(interp_method,'linear') 
  newx = interp1(t,x,newt,'linear');
elseif 1 == strcmp(interp_method,'nearest')
   dt = [0; diff(t);];
%   t  = t + (1e-10 * ((dt == 0) .* (1:length(t))'));
   newx = interp1(t,x,newt,'nearest');
   % newx = nearest(t,x,newt);
else
  fprintf(1,'RASAMPLE: ERROR: Bad interpolation method.\n');
  return;
end

return;  



% ----------------------------------------------------------------------
function newx = nearest(t,x,newt)
%
% Oct  4 1998  Louis Whitcomb  Created and written

rti   = t(1);
rtf   = t(length(t));

% find time range of new time sequence
dti   = newt(1);
dtf   = newt(length(newt));

bestj     = 1;
bestjabs  = abs(t(1) - newt(1));

newx = zeros(length(newt),length(x(1,:)));

for i=1:length(newt)
   for j = bestj:length(t)
      if abs(newt(i) - t(j)) < bestjabs
        bestj    = j;
        bestjabs = abs(newt(i) - t(j));
      end   
   end
   newx(i,:) = x(bestj); 
end

return

  
  



