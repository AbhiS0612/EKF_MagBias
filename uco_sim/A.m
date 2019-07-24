function A_mat=A(t,samp,s0) % needs current state
% 2019-05-22 LLW create time varying A(t) 
   
  ang_m = interp1(samp.t,samp.ang,t)';
  mag_m = s0(1:3,1);
%   mag_b = s0(4:6,1);
%   Ts = s0(7:12,1)';
%   T = [Ts(1:3)' [Ts(2) Ts(4) Ts(5)]' [Ts(3) Ts(5) Ts(6)]']; 
  %ang_mag_b = s0(13:15,1);
  
   mag_b = [0.1, 0.05, -0.1]';
   Td = s0(4:6,1)';
    T = [Td(1)   0.1   -0.2;
         0.1    Td(2)  0.05;
        -0.2    0.05  Td(3)];

  
 % Z = zeros(3,3);
 % T = eye(3);
  
  %jacomag_bian of T-stack:  dim(9,6)
%   D_Ts = [1 0 0 0 0 0;
%           0 1 0 0 0 0;
%           0 0 1 0 0 0;
%           0 1 0 0 0 0;
%           0 0 0 1 0 0;
%           0 0 0 0 1 0;
%           0 0 1 0 0 0;
%           0 0 0 0 1 0;
%           0 0 0 0 0 1];
%   %jacomag_bian of inv(T)-stack:  dim(9,6)
%   D_Tsi = -kron(inv(T)', inv(T))*D_Ts;
%   
%   % jacomag_bian of J(w_mag_b)-stack: dim(9,6)
%   D_Js = [0  0  0;
%           0  0  1;
%           0 -1  0;
%           0  0 -1;
%           0  0  0;
%           1  0  0;
%           0  1  0;
%          -1  0  0;
%           0  0  0]; 
      
  D_Td = [1 0 0
          0 0 0
          0 0 0
          0 0 0
          0 1 0
          0 0 0
          0 0 0
          0 0 0
          0 0 1];
      
  D_Tdi = -kron(inv(T)', inv(T))*D_Td;     
          
  
    A_mat = [-T*J(ang_m)/T, -(kron((J(ang_m)*(T\(mag_m-mag_b)))', eye(3))*D_Td + (kron((mag_m-mag_b)',(T*J(ang_m)))*D_Tdi));
              zeros(3,6)];
      
%   A_mat = [-T*J(ang_m)/T,   T*J(ang_m)/T,   -(kron((J(ang_m)*(T\(mag_m-mag_b)))', eye(3))*D_Ts + (kron((mag_m-mag_b)',(T*J(ang_m)))*D_Tsi));
%               zeros(9,12)];    
%       
%   A_mat = [-T*J(ang_m-ang_mag_b)/T,   T*J(ang_m-ang_mag_b)/T,   -(kron(((J(ang_m-ang_mag_b)/T)*(mag_m-mag_mag_b))',eye(3))*D_Ts + (kron((mag_m-mag_mag_b)',(T*J(ang_m-ang_mag_b)))*D_Tsi)),  -T*J(T\(mag_m-mag_mag_b));          
%               zeros(12,15)];
                      
 end