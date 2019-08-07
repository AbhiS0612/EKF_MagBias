%% Histplot function
 % Description: The function plots and fits histograms, for visualizing
 % data
 % Inputs: nx3 vector of data, 
 %         no. of bins (optional), 
 %         chart title (optional).
 % Output: Returns the mean and std of each n-vector.
 %         Displays 3x1 histogram plots withtheir mean and std
 % Author: Abhi Shah (06/15/19)
 
function specs = histplots(x, bins, str)

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
histogram(x(:,1), bins);
title (str);
[v,h] = hist(x(:,1),bins);
yp = 0.8*max(v);
xp = h(floor(0.8*bins));
str = sprintf('\\mu: %f',mean(x(:,1)));
text(xp,yp, str);
yp = 0.6*max(v);
str = sprintf('\\sigma: %f',std(x(:,1)));
text(xp,yp, str);
ylabel 'x';

subplot(312)
histogram(x(:,2), bins);
[v,h] = hist(x(:,2),bins);
yp = 0.8*max(v);
xp = h(floor(0.8*bins));
str = sprintf('\\mu: %f',mean(x(:,2)));
text(xp,yp, str);
yp = 0.6*max(v);
str = sprintf('\\sigma: %f',std(x(:,2)));
text(xp,yp, str);
ylabel 'y';

subplot(313)
histogram(x(:,3), bins);
[v,h] = hist(x(:,3),bins);
yp = 0.8*max(v);
xp = h(floor(0.8*bins));
str = sprintf('\\mu: %f',mean(x(:,3)));
text(xp,yp, str);
yp = 0.6*max(v);
str = sprintf('\\sigma: %f',std(x(:,3)));
text(xp,yp, str);
ylabel 'z';
xlabel 'Data Values';

specs.me = mean(x);
specs.st = std(x);
end