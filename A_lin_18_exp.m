 %% Function to return linearized A matrix in continuous time, using interpolation
  % This function is to be used for with exp_troni_KF_5 (18 states)
  % Modified from (2019-05-22 LLW) 
  
 function A_mat=A_lin_18_exp(t,samp,s0) % needs current state
  %ang_m = interp1(samp.t,samp.ang,t)';
  mag_m = s0(1:3,1);
  mag_b = s0(4:6,1);
  Ts = s0(7:12,1)';
  T = [Ts(1:3)' [Ts(2) Ts(4) Ts(5)]' [Ts(3) Ts(5) Ts(6)]']; 
  ang_b = s0(13:15,1);
  ang_m = s0(16:18,1);
  
  Z = zeros(3,3);
 % T = eye(3);
  
  %jacobian of T-stack:  dim(9,6)
  D_Ts = [1 0 0 0 0 0;
          0 1 0 0 0 0;
          0 0 1 0 0 0;
          0 1 0 0 0 0;
          0 0 0 1 0 0;
          0 0 0 0 1 0;
          0 0 1 0 0 0;
          0 0 0 0 1 0;
          0 0 0 0 0 1];
  %jacobian of inv(T)-stack:  dim(9,6)
  D_Tsi = -kron(inv(T)', inv(T))*D_Ts;
        
  A_mat = [-T*J(ang_m-ang_b)/T   T*J(ang_m-ang_b)/T   -(kron(((J(ang_m-ang_b)/T)*(mag_m-mag_b))', eye(3))*D_Ts +(kron((mag_m-mag_b)',(T*J(ang_m-ang_b)))*D_Tsi))  -T*J( T\(mag_m-mag_b))   T*J( T\(mag_m-mag_b));          
              zeros(15,18)];
          
 end
 