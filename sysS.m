function xdot = sysS(t,x,samp)
  we = interp1(samp.t,samp.ang,t)';
  xm = x(1:3,1);
  b = x(4:6,1);
  Ts = x(7:12,1)';
  T = [Ts(1:3)' [Ts(2) Ts(4) Ts(5)]' [Ts(3) Ts(5) Ts(6)]'];

  xdot = [-(T*J(we)/T)*(xm-b);
              zeros(9,1)];
end