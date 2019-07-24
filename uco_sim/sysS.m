function xdot = sysS(t,x,samp)
  we = interp1(samp.t,samp.ang,t)';
 % wb = x(13:15,1);
  xm = x(1:3,1);
%   b = x(4:6,1);
%   Ts = x(7:12,1)';
%   T = [Ts(1:3)' [Ts(2) Ts(4) Ts(5)]' [Ts(3) Ts(5) Ts(6)]'];

    b = [0.1, 0.05, -0.1]';
    Td = x(4:6,1)';
    T = [Td(1)   0.1   -0.2;
         0.1    Td(2)  0.05;
        -0.2    0.05  Td(3)];

  xdot = [-(T*J(we)/T)*(xm-b);
              zeros(3,1)];
end