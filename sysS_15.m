%% function sysS_15(t,x,samp): ode for 15 state model
function xdot = sysS_15(t,x,samp)
  we = interp1(samp.t,samp.ang,t)';
  xm = x(1:3,1);
  b = x(4:6,1);
  Ts = x(7:12,1)';
  T = [Ts(1:3)' [Ts(2) Ts(4) Ts(5)]' [Ts(3) Ts(5) Ts(6)]'];
  wb = x(13:15,1);

  xdot = [-(T*J(we-wb)/T)*(xm-b);
              zeros(12,1)];
end