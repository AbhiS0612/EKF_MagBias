function w_e_vec = w_e(t)
% 2019-05-23 LLW create time varying w_e(t), called by A(t)   

w_e_vec = 1.0 * [ 0.0 * sin(1.2*t); 0.2 * cos(2.3*t); 0.5 * sin(1.5 * t)];  
