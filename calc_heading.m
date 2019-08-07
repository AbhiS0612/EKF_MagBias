function out = calc_heading(t,samp,phins)


    %samp.mag = rph2R([0;0;pi])*samp.mag; %kvh
    %R = rph2rot([-90.00;0.13;89.98]*pi/180); %kvh
    R = rph2rot([0.48;0.07;-179.86]*pi/180); %mst


    samp.acc = R*samp.acc;
    samp.mag = R*samp.mag;


    
    %samp.acc = my_lowpass(samp.acc',50,1,1)';
    %samp.mag = my_lowpass(samp.mag',50,1,1)';

    roll  = atan2(-samp.acc(2,:),-samp.acc(3,:));
    pitch = atan2(samp.acc(1,:),sqrt(samp.acc(2,:).*samp.acc(2,:) + samp.acc(3,:).*samp.acc(3,:)));
    
    mag = zeros(size(samp.mag));
    
    %phins_att = resample2(phins.t,phins.att,t)';
    %roll = phins_att(1,:)*pi/180;
    %pitch = phins_att(2,:)*pi/180;
    
    for i=1:size(samp.t,1)
       
        mag(:,i) =  rph2R([roll(i),pitch(i),0])*samp.mag(:,i);

    end

        
    heading = atan2(-mag(2,:),mag(1,:));

    att = [roll',pitch',heading'];
    
    att = unwrap(att);

    att = resample2(t,att,phins.t)';
    phins.att = unwrap(phins.att*pi/180)*180/pi;
    
    out.att = att;
    out.att_error = wrapTo180(att'*180/pi-phins.att)';
    out.mag_ned = mag;
    %figure;plot(out.mag_ned(1,:),out.mag_ned(2,:),'.');grid on;axis equal;xlabel('x (gaus)');ylabel('y (gaus)');title('Cor Mag - Local Level Plane');
    pt = taxis(phins.t);
    figure;plot(taxis(phins.t),out.att_error,'.');grid on;xlabel(tlabel(phins.t));ylabel('degrees');title('Att Error');xlim([pt(1),pt(end)]);ylim([-50,50]);
    figure;plot(taxis(phins.t),phins.att,'.',taxis(phins.t),out.att*180/pi,'.');grid on;legend('px','py','pz','x','y','z');xlabel(tlabel(phins.t));title('Attitude');xlim([pt(1),pt(end)]);ylim([-200,200]);
    %figure;plot(taxis(samp.t),samp.b);grid on;xlabel(tlabel(samp.t));legend('x','y','z');title('bias');ylabel('gaus');xlim([pt(1),pt(end)]);
    %figure;plot(taxis(samp.t),samp.theta);grid on;legend('a','b','c','d','e','f','a1','a2','a3','btb');xlabel(tlabel(samp.t));title('theta');xlim([pt(1),pt(end)]);

        
        