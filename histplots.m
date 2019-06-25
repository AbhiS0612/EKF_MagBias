%% Histplot function
 % Description: The function plots and fits histograms, for visualizing
 % data
 % Inputs: nx3 vector of data, 
 %         no. of bins (optional), 
 %         chart title (optional).
 % Output: 3x1 histogram plots, also displays mean and std of each
 %         distribution.
 % Author: Abhi Shah (06/15/19)
 
function histplots(x, bins, str)

% defualts to 100 bins
if( exist('bins','var') ~= 1)
  bins = 100;
end

% default title
if( exist('str','var') ~= 1)
  str = 'Histogram Plots';
end

figure
subplot(311)
histfit(x(:,1), bins);
title (str);
[v,h] = hist(x(:,1),bins);
yp = 0.8*max(v);
xp = h(floor(0.8*bins));
str = sprintf('\\mu: %f',mean(x(:,1)));
text(xp,yp, str);
yp = 0.6*max(v);
str = sprintf('\\sigma: %f',std(x(:,1)));
text(xp,yp, str);

subplot(312)
histfit(x(:,2), bins);
[v,h] = hist(x(:,2),bins);
yp = 0.8*max(v);
xp = h(floor(0.8*bins));
str = sprintf('\\mu: %f',mean(x(:,2)));
text(xp,yp, str);
yp = 0.6*max(v);
str = sprintf('\\sigma: %f',std(x(:,2)));
text(xp,yp, str);

subplot(313)
histfit(x(:,3), bins);
[v,h] = hist(x(:,3),bins);
yp = 0.8*max(v);
xp = h(floor(0.8*bins));
str = sprintf('\\mu: %f',mean(x(:,3)));
text(xp,yp, str);
yp = 0.6*max(v);
str = sprintf('\\sigma: %f',std(x(:,3)));
text(xp,yp, str);
xlabel 'Data Values';
end