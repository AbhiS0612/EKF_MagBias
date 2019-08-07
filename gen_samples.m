%% function samp = gen_samples(lat,hz,t_end,bias,w_max)
 % Function to generate simulated data for MST, KVH, PHINS
 % Author: Andrew Spielvogel
 % 
 % 06/25/19: Modified by Abhi Shah, updated std devs. for MST.
 %           Added the ability control the max angular velocity.
 
function samp = gen_samples(lat,hz,t_end,bias,w_max)
dt = 1/hz;
t= 0:dt:(t_end);
num = length(t);

lat = lat*pi/180;

Rni{num} = eye(3);
Rni{1} = eye(3);

ang = zeros(3,num);
acc = zeros(3,num);
samp.E = zeros(3,num);
T = bias.T;

%w_sig = 6.32 * 10^(-3)*pi/180;  % measured 1775, units are rad/sec
%a_sig = 0.0037;            % measured 1775, units are g, not m/s^2
a_sig = 8 * 10^(-4);  % MST @ 20Hz
w_sig = 3 * 10^(-4);  % MST @ 20Hz
%m_sig = 2 * 10^(-4);  % MST @ 20Hz

%most recent m_sig - seems to have different cov for different axes
m_sig = [1.7;1.7;3] * 10^(-4);  %housing MS @20hz, with filter
%m_sig = 2 * 10^(-3);  % random testing

% generate a_n
r = 6371*1000; %earth radius in m
Ren = [-sin(lat),0,-cos(lat);0,1,0;cos(lat),0,-sin(lat)];
a_e = [cos(lat);0;sin(lat)] - (15.04*pi/180/3600)^2*cos(lat)*[r;0;0]/9.81;
a_n = Ren'*a_e;

m_n = [0.205796;-0.040654;0.468785];
%m_n = [1;0;0];

w_E_e = [0;0;1]*15.04*pi/180/3600; % multiply by 0 to set no earth rate
w_E_n = Ren'*w_E_e*0;

%fileID = fopen('/home/spiels/log/sim/kvh/sim5.KVH','w');
fileMST = fopen('/home/abhis/Matlab/DSCL/log/sim1.MST','w');
%filePhins = fopen('/home/abhis/Matlab/DSCL/log/sim1.INS','w');

%T=[.95,0,0; 0,1.1,0;0,0,1.05]; %ijrr, soft iron T

for i=1:num
    % get w_v, Rsn
    w_veh = get_w(t(i),w_max);

    if i~=num
        Rni{i + 1} = Rni{i}*expm(J(w_veh)*dt);
    end

    ang(:,i) = w_veh + Rni{i}'*w_E_n + normrnd(bias.ang, w_sig,[3,1]);
    acc(:,i) = Rni{i}'*a_n + bias.acc + normrnd(0,a_sig,[3,1]);%+J(w_veh)*[0.1;0;0];
    samp.acc_nv(:,i) = acc(:,i);
    samp.E(:,i) = Rni{i}'*[0;1;0]*norm(J(w_E_n)*a_n);
    samp.D(:,i) = Rni{i}'*a_n;
    samp.att(:,i) = rot2rph(Rni{i});
    samp.w_E(:,i) = Rni{i}'*w_E_n;
    samp.w_E_n(:,i) = Rni{i}'*[w_E_n(1);0;0];
    %T=eye(3);
    samp.mag(:,i) = T*Rni{i}'*m_n + normrnd(0,m_sig,[3,1]) + bias.mag;
    samp.mag_true(:,i) = Rni{i}'*m_n;
    %samp.mag(:,i) = Rni{i}'*m_n + m_sig*randn(3,1) + bias.mag;
    
    % print progress
    if ~mod(t(i),30)
        str = sprintf('Made %i:%i0 of data at %i hz',floor(t(i)/60),mod(t(i),60)/10,hz);
        disp(str);
    end
    R = Rni{i};
    rpy_phins = samp.att(:,i)*180/pi;%rot2rph(R);
    %fprintf(fileID,'IMU_RAW, %.40f,%.40f,%.40f, %.35f,%.35f,%.35f,%f,%f,%f, 0, 0, %.30f,0,1,1,1,1,1,1, %f,%f,%f,%f,%f,%f,%f,%f,%f,%.40f,%.40f,%.40f \n',ang(1,i),ang(2,i),ang(3,i),acc(1,i),acc(2,i),acc(3,i),samp.mag(1,i),samp.mag(2,i),samp.mag(3,i),t(i),R(1,1),R(1,2),R(1,3),R(2,1),R(2,2),R(2,3),R(3,1),R(3,2),R(3,3),rpy_phins(1),rpy_phins(2),rpy_phins(3));
    %fprintf(fileID,'KVH 0/0/0 0:0:0 %f %f %.40f,%.40f,%.40f, %.35f,%.35f,%.35f,%f,%f,%f,0,%d,%f,%f,1,1,1,1,1,1 \n',t(i),t(i),ang(1,i),ang(2,i),ang(3,i),acc(1,i),acc(2,i),acc(3,i),samp.mag(1,i),samp.mag(2,i),samp.mag(3,i),mod(i,128),t(i),t(i));
   % fprintf(filePhins,'INS 0/0/0 0:0:0 %f %f %f,0,0,0,0,0,0,%f,%f,%f,0,0 \n',t(i),t(i),t(i),rpy_phins(1),rpy_phins(2),rpy_phins(3));
    fprintf(fileMST,'MST 0/0/0 0:0:0 %f %f %f,%f,%f,%f,%f,%f,%f,%f,%f,%f,1 \n',t(i),t(i),t(i),ang(1,i),ang(2,i),ang(3,i),acc(1,i),acc(2,i),acc(3,i),samp.mag(1,i),samp.mag(2,i),samp.mag(3,i));

end

samp.t = t;
samp.acc = acc';
samp.ang = ang';
samp.mag = samp.mag';
samp.Rni = Rni;
samp.bias = bias;


function w = get_w(t, w_max)
%w = zeros(3,1);    %no movement, used to check noise
%w = [0;0;cos(t/11)*2];   %not PE, can't converge on z

% sampled that aim to replicate data from 'log/microstrain/2019_08_02_16_27.MST'
w = [cos(2*t)*w_max(1); cos(t) * w_max(2); -cos(t/20) * w_max(3)];
end

end