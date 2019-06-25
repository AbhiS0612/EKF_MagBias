 %% Function to return linearized A matrix in discrete time
 function A_mat=A_lin_dis(t,samp,s0,ts) % needs current state
% 2019-05-22 LLW create time varying A(t) 
    
  we = samp.ang(t,:)';
  xm = s0(1:3,1);
  b = s0(4:6,1);
  Ts = s0(7:12,1)';
  T = [Ts(1:3)' [Ts(2) Ts(4) Ts(5)]' [Ts(3) Ts(5) Ts(6)]']; 
  Z = zeros(3,3);
 % T = eye(3);
  
  %jacobina of T-stack:  dim(9,6)
  D_Ts = [1 0 0 0 0 0;
          0 1 0 0 0 0;
          0 0 1 0 0 0;
          0 1 0 0 0 0;
          0 0 0 1 0 0;
          0 0 0 0 1 0;
          0 0 1 0 0 0;
          0 0 0 0 1 0;
          0 0 0 0 0 1];
  
  D_Tsi = -kron(inv(T)', inv(T))*D_Ts;

%   A_mat = [-J(we),  J(we),  eye(3),  eye(3);...
%             Z,       Z,       Z,      Z;...
%             Z,       Z,       Z,      Z;...
%             Z,       Z,       Z,      Z];
        
  A_mat_c = [-T*J(we)/T   T*J(we)/T   -(kron(((J(we)/T)*(xm-b))', eye(3))*D_Ts +(kron((xm-b)',(T*J(we)))*D_Tsi));
              zeros(9,12)];
          
  A_mat = expm(A_mat_c * ts); 
 end
 