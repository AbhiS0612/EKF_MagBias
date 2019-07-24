% 2019-05-23 LLW script file to make obs grammian plots

%% Generate Simulated Data
samp = read_microstrain('/home/abhis/Matlab/DSCL/log/sim1.MST');
samp.t = samp.t-samp.t(1);
t0 = 0;
t1 = 60;
dt = 1;

t0 = samp.t(1);
t1 = t0+t1;

% row vector
t = t0:dt:t1;

%initial condition, column vector
 s0 = [samp.mag(1,:) [1 1 1]]';
%s0 = [samp.mag(1,:) [0 0 0] [1 0 0 1 0 1]]';

% function is vectorized
w_e_mat = interp1(samp.t,samp.ang,t);

% function is vectorized
a_e_mat = interp1(samp.t,samp.acc,t);

figure;
subplot(211)
plot(samp.t,samp.ang);
grid on;
xlabel('Time (s)');
ylabel('\omega_m (rad)');
legend('\omega_m_1', '\omega_m_2', '\omega_m_3');
title('\omega_m vs time');

subplot(212)
plot(samp.t, samp.acc);
grid on;
xlabel('Time (s)');
ylabel('a_m (m/s^2)');
legend('a_m_1', 'a_m_2', 'a_m_3'); 
title('a_m vs time');


% plot w_e and a_e vs time
figure

subplot(211)
plot(t, w_e_mat);
grid on;
xlabel('Time (s)');
ylabel('\omega_e (rad)');
legend('\omega_e_1', '\omega_e_2', '\omega_e_3');
title('\omega_e vs time');

subplot(212)
plot(t, a_e_mat);
grid on;
xlabel('Time (s)');
ylabel('a_e (m/s^2)');
legend('a_e_1', 'a_e_2', 'a_e_3'); 
title('a_e vs time');

% compute observability grammian
fprintf(1,'Computing Observability Grammian from t0=%.1f to t1=%.1f with dt=%.02f:\n', t0, t1, dt);
OG = obs_gram(samp,t0, t1, dt, s0);

fprintf(1,'Rank of Observability Grammian=%d\n', rank(OG));

fprintf(1,'Eigenvalues of Observability Grammian:\n');
[vec, val] = eig(OG);
diag(val)

fprintf(1,'Determinant of Observability Grammian=%e\n', det(OG)');