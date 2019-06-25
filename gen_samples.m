function samp = gen_samples(lat,hz,t_end,bias,T)

dt = 1/hz;
t= 0:dt:t_end;
num = size(t,2);

lat = lat*pi/180;

Rni{num} = eye(3);
Rni{1} = eye(3);

ang = zeros(3,num);
acc = zeros(3,num);
samp.E = zeros(3,num);


%w_sig = 6.32 * 10^(-3)*pi/180;  % measured 1775, units are rad/sec
a_sig = 0.0037;            % measured 1775, units are g, not m/s^2

w_sig = 5 * 10^(-4); % roughly what was measured from experiments with MST, sketchy run though
m_sig = 0.001;       % measured from experiments, sketchy run though

% generate a_n
r = 6371*1000;
Ren = [-sin(lat),0,-cos(lat);0,1,0;cos(lat),0,-sin(lat)];
a_e = [cos(lat);0;sin(lat)] - (15.04*pi/180/3600)^2*cos(lat)*[r;0;0]/9.81;
a_n = Ren'*a_e;

m_n = [0.205796;-0.040654;0.468785];
%m_n = [1;0;0];

w_E_e = [0;0;1]*15.04*pi/180/3600;
w_E_n = Ren'*w_E_e;

%fileID = fopen('/home/spiels/log/sim/kvh/sim5.KVH','w');
fileMST = fopen('/home/abhis/Matlab/DSCL/log/sim1.MST','w');
%filePhins = fopen('/home/abhis/Matlab/DSCL/log/sim1.INS','w');

%T=[.95,0,0; 0,1.1,0;0,0,1.05]; %ijrr, soft iron T

%T=eye(3);

for i=1:num

    % get w_v, Rsn
    w_veh = get_w(t(i));

    if i~=num
        Rni{i + 1} = Rni{i}*expm(skew(w_veh)*dt);
    end

    ang(:,i) = w_veh + Rni{i}'*w_E_n + bias.ang + normrnd(0,w_sig,[3,1]);
    acc(:,i) = Rni{i}'*a_n + bias.acc + normrnd(0,a_sig,[3,1]);%+skew(w_veh)*[0.1;0;0];
    samp.acc_nv(:,i) = acc(:,i);
    samp.E(:,i) = Rni{i}'*[0;1;0]*norm(skew(w_E_n)*a_n);
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


function w = get_w(t)

%%%% IROS2018

% w = [cos(t/2)/20;sin(t/5)/15;-cos(t/30)/10]*0;
% 
% w = [cos(t/50)/25;-sin(t/9)/7*0;cos(t/10)/6]; %exp1
% w = [cos(t/50)/25*0;-sin(t/9)/7*0;cos(t/10)/6]; %exp2
% w = [cos(t/50)/100;-sin(t/9)/7*0;cos(t/10)/6]; %exp3
% w = [cos(t/50)/25;-sin(t/25)/10;cos(t/10)/5]; %exp4 sim2
% w = [cos(t/20)/55*0;-sin(t/5)/30*0;cos(t/30)/5]; %exp5 10hz optimization

%w = [cos(t/20)/20;sin(t/50)/32;cos(t/60)/10];

%w = [0;0;cos(t/50)/15]; % 6(10hz) sim0
%w = [0;0;0]; % (1000hz) sim1
%w = [0;0;cos(t/60)/10]; % (1000hz) sim2
%w = [cos(t/5)/50;sin(t/10)/40;cos(t/60)/10]; % (1000hz) sim3
%w = [sin(t/50)/8+cos(t/10)/5;sin(t/10)/40+cos(t/30)/20;cos(t/60)/10]; % (1000hz) sim3

%w = [sin(t/20)/20+cos(t/10)/15;sin(t/10)/40+cos(t/30)/25;cos(t/60)/10]; % (1000hz) sim4 (sim3 in IJRR paper)

%w = [(cos(t/7.2)*2 +sin(t/3));  (-cos(t/5.1));  (cos(t/11))*2];
w = [0;0;cos(t/11)*2];


if t>2500
    %w=0*[0;0;sin(t/20)/10];
end