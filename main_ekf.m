rng default
R_E = s3_constants('R_EARTH');
mu = s3_constants('GM_EARTH');
r = R_E + 500;
sat_data = [r,0,0,100;
    0,r,0,150;
    0,0,r,300;
    -r,0,0,50;
    0,-r,0,175;
    0,0,-r,250];
num_gps = size(sat_data,1);

eci_c = prop_data.cart_chief;
eci_d = prop_data.cart_dep;
t = prop_data.time_vec;
n_steps = length(t);

% true initial state
b_c_true = 110;
b_d_true = 178;
x_init = [eci_c(:,1); eci_d(:,1); b_c_true; b_d_true];
num_states = length(x_init);

x = x_init + mvnrnd(zeros(1,num_states), diag([50 50 50 10 10 10 50 50 50 10 10 10 5 5]))'; % add some noise
Sigma = 10 * eye(num_states,num_states);
R = diag([5 5 5 100 100 100 5 5 5 100 100 100 5 5]);
Q = 5 * eye(2*num_gps,2*num_gps);

% store x and sigma
x_hist = zeros(num_states,n_steps);
Sigma_hist = zeros(num_states,num_states,n_steps); 
x_hist(:,1) = x;
Sigma_hist(:,:,1) = Sigma;

for i = 2:n_steps
    x_true = [eci_c(:,i);eci_d(:,i); b_c_true; b_d_true];
    meas_true = meas_pred(x_true,sat_data);
    meas = meas_true + mvnrnd(zeros(1,2*num_gps), diag(3*ones(1,2*num_gps)))'; % add some noise
    dt = t(i) - t(i-1);
   
    [x, Sigma] = ekf(x, Sigma, meas, dt, mu, sat_data, Q, R);
    x_hist(:,i) = x;
    Sigma_hist(:,:,i) = Sigma;
end

rel_eci = x_hist(7:12,:) - x_hist(1:6,:);
rtn_est = zeros(6,n_steps);
rel_cov_rtn = zeros(6,6,n_steps);
for i = 1:n_steps
    x_c = x_hist(1:6,i);
    x_d = x_hist(7:12,i);
    rel_eci_i = rel_eci(:,i);
    Rot =  R_RTN_to_ECI(x_c)';
    R_full = [Rot zeros(3,3); zeros(3,3) Rot];
    rel_rtn = R_full * rel_eci_i;
    rel_cov_rtn(:,:,i) = R_full * Sigma_hist(1:6,1:6,i) * R_full' + R_full * Sigma_hist(7:12,7:12,i) * R_full';
    rtn_est(:,i) = rel_rtn;
end

chi_vec = prop_data.chi_vec;
rr_true = chi_vec(1,:);
rt_true = chi_vec(2,:);
rn_true = chi_vec(3,:);

rr_est = rtn_est(1,:);
rt_est = rtn_est(2,:);
rn_est = rtn_est(3,:);
figure
% subplot(1,2,1)
hold on
plot(rt_true,rr_true)
plot(rt_est,rr_est,'.')
theta = linspace(0,2*pi);
for i = 1:n_steps
    rel_pos_cov_rtn_i = rel_cov_rtn(1:3,1:3,i);
    sig_rr = rel_pos_cov_rtn_i(1,1);
    sig_rt = rel_pos_cov_rtn_i(1,2);
    sig_tt = rel_pos_cov_rtn_i(2,2);
    cov = [sig_tt, sig_rt; sig_rt, sig_rr];
    xy_ellipse = cov_ellipse([rt_est(i); rr_est(i)],cov,5.991,theta);
    plot(xy_ellipse(1,:), xy_ellipse(2,:),'g')
end
axis equal
grid on
xlabel('T Rel Position (m)')
ylabel('R Rel Position (m)')
legend('True Position', 'Mean Position Estimate', 'Estimate 95% Confidence Interval')
% subplot(1,2,2)
figure
hold on
plot(rn_true,rr_true)
plot(rn_est,rr_est,'.')
for i = 1:n_steps
    rel_pos_cov_rtn_i = rel_cov_rtn(1:3,1:3,i);
    sig_rr = rel_pos_cov_rtn_i(1,1);
    sig_rn = rel_pos_cov_rtn_i(1,3);
    sig_nn = rel_pos_cov_rtn_i(3,3);
    cov = [sig_nn, sig_rn; sig_rn, sig_rr];
    xy_ellipse = cov_ellipse([rn_est(i); rr_est(i)],cov,5.991,theta);
    plot(xy_ellipse(1,:), xy_ellipse(2,:),'g')
end
axis equal
grid on
xlabel('N Rel Position (m)')
ylabel('R Rel Position (m)')
%% for each state estimate, get the error between real and estimated absolute states and relative state
r_eci_c_error = zeros(1,n_steps);
r_eci_d_error = zeros(1,n_steps);
v_eci_c_error = zeros(1,n_steps);
v_eci_d_error = zeros(1,n_steps);
r_rtn_rel_error = zeros(1,n_steps);
v_rtn_rel_error = zeros(1,n_steps);
for i = 1:n_steps
    r_eci_c_true = eci_c(1:3,i);
    r_eci_d_true = eci_d(1:3,i);
    v_eci_c_true = eci_c(4:6,i);
    v_eci_d_true = eci_d(4:6,i);
    
    r_eci_c_est = x_hist(1:3,i);
    r_eci_d_est = x_hist(7:9,i);
    v_eci_c_est = x_hist(4:6,i);
    v_eci_d_est = x_hist(10:12,i);
    
    r_rtn_rel_true = chi_vec(1:3,i);
    v_rtn_rel_true = chi_vec(4:6,i);
    r_rtn_rel_est = rtn_est(1:3,i);
    v_rtn_rel_est = rtn_est(4:6,i);
    
    r_eci_c_error(i) = norm(r_eci_c_true-r_eci_c_est);
    r_eci_d_error(i) = norm(r_eci_d_true-r_eci_d_est);
    v_eci_c_error(i) = norm(v_eci_c_true-v_eci_c_est);
    v_eci_d_error(i) = norm(v_eci_d_true-v_eci_d_est);
    r_rtn_rel_error(i) = norm(r_rtn_rel_true - r_rtn_rel_est);
    v_rtn_rel_error(i) = norm(v_rtn_rel_true - v_rtn_rel_est);
end
figure
hold on
plot(t, r_eci_c_error)
plot(t, r_eci_d_error)
plot(t, r_rtn_rel_error)
grid on
title('Position Mean Estimate Error')
legend('Chief Error','Deputy Error','Relative Error')
xlabel('Time (s)')
ylabel('Norm of Position Error from Truth (m)')

figure
hold on
plot(t, v_eci_c_error)
plot(t, v_eci_d_error)
plot(t, v_rtn_rel_error)
grid on
title('Velocity Mean Estimate Error')
legend('Chief Error','Deputy Error','Relative Error')
xlabel('Time (s)')
ylabel('Norm of Velocity Error from Truth (m/s)')
%% at each estimate, propagate the orbit and check minimum separation
dt = 20;
t_max = 100;
n_samples = floor(t_max / dt);
confidence = 7.815; % 95%
min_seps = zeros(1,n_samples);
figure
for i = 1:n_steps
   x = x_hist(:,i);
   Sigma = Sigma_hist(:,:,1);
   
   x_prop = x;
   min_d = 1000; % high number to start
   rtn_prop_hist = zeros(6,n_samples);
   Sigma_rtn_hist = zeros(6,6,n_samples);
   x_prop_hist = zeros(14,n_samples);
   for j = 1:n_samples
       F = dynamics_jac(x_prop,mu,dt);
       x_prop = F*x_prop;
       Sigma = F * Sigma * F';
       
       Rot =  R_RTN_to_ECI(x_prop(1:6))';
       R_full = [Rot zeros(3,3); zeros(3,3) Rot];
       rel_rtn = R_full * (-x_prop(1:6) + x_prop(7:12));
       rel_cov_rtn = R_full * Sigma(1:6,1:6) * R_full' + R_full * Sigma(7:12,7:12) * R_full';
       d = dist2cov(rel_rtn(1:3),rel_cov_rtn(1:3,1:3),confidence);
       min_d = min([min_d,d]);
       
       rtn_prop_hist(:,j) = rel_rtn;
       Sigma_rtn_hist(:,:,j) = rel_cov_rtn;
       x_prop_hist(:,j) = x_prop;
   end
   
   
   w = sqrt(mu/(r^3));
   STM_rtn = [4-3*cos(w*dt), 0, 0, 1/w*sin(w*dt), 2/w*(1-cos(w*dt)), 0;
    6*sin(w*dt)-6*w*dt, 1, 0, 2/w*(cos(w*dt)-1), 4/w*sin(w*dt)-3*dt, 0;
    0, 0, cos(w*dt), 0, 0, 1/w*sin(w*dt);
    3*w*sin(w*dt), 0, 0, cos(w*dt), 2*sin(w*dt), 0;
    6*w*cos(w*dt)-6*w, 0, 0, -2*sin(w*dt), 4*cos(w*dt)-3, 0;
    0, 0, -w*sin(w*dt), 0, 0, cos(w*dt)];
    Sigma = Sigma_hist(:,:,1);
    Rot = R_RTN_to_ECI(x(1:6))';
    R_full = [Rot zeros(3,3); zeros(3,3) Rot];
    rel_rtn = R_full * (-x(1:6) + x(7:12));
    hcw_rtn = zeros(6,n_samples);
    Sigma_rtn = R_full * Sigma(1:6,1:6) * R_full' + R_full * Sigma(7:12,7:12) * R_full';
    Sigma_hcw_hist = zeros(6,6,n_samples);
   for j = 1:n_samples
       
       rel_rtn = STM_rtn * rel_rtn;
       Sigma_rtn = STM_rtn * Sigma_rtn * STM_rtn'; 
%        F = dynamics_jac(x_prop,mu,dt);
%        x_prop = F*x_prop;
%        Sigma = F * Sigma * F';
%        
%        Rot =  R_RTN_to_ECI(x_prop(1:6))';
%        R_full = [Rot zeros(3,3); zeros(3,3) Rot];
%        rel_rtn = R_full * (x_prop(1:6) - x_prop(7:12));
%        rel_cov_rtn = R_full * Sigma(1:6,1:6) * R_full' + R_full * Sigma(7:12,7:12) * R_full';
%        d = dist2cov(rel_rtn(1:3),rel_cov_rtn(1:3,1:3),confidence);
%        min_d = min([min_d,d]);
       
       hcw_rtn(:,j) = rel_rtn;
       Sigma_hcw_hist(:,:,j) = Sigma_rtn;
%        rtn_prop_hist(:,j) = rel_rtn;
%        Sigma_rtn_hist(:,:,j) = rel_cov_rtn;
%        x_prop_hist(:,j) = x_prop;
   end
   
   
   
   
   min_seps(i) = min_d
   i
   
   subplot(1,2,1)
   plot(x_prop_hist(1,:),x_prop_hist(2,:))
   subplot(1,2,2)
   plot(x_prop_hist(1,:),x_prop_hist(3,:))
   
   clf
   subplot(1,2,1)
   plot(hcw_rtn(2,:),hcw_rtn(1,:))
   axis equal
   subplot(1,2,2)
   plot(hcw_rtn(3,:),hcw_rtn(1,:))
   axis equal
   
   clf
   subplot(1,2,1)
   hold on
   plot(rtn_prop_hist(2,:),rtn_prop_hist(1,:))
   theta = linspace(0,2*pi);
   for j = 1:n_samples
%        rel_pos_cov_rtn_i = Sigma_rtn_hist(1:3,1:3,j);
       rel_pos_cov_rtn_i = Sigma_hcw_hist(1:3,1:3,j); %careful, mixing up methods here

       sig_rr = rel_pos_cov_rtn_i(1,1);
       sig_rt = rel_pos_cov_rtn_i(1,2);
       sig_tt = rel_pos_cov_rtn_i(2,2);
       cov = [sig_tt, sig_rt; sig_rt, sig_rr];
       xy_ellipse = cov_ellipse([rtn_prop_hist(2,j);rtn_prop_hist(1,j)],cov,5.991,theta);
       plot(xy_ellipse(1,:), xy_ellipse(2,:),'g')
   end
   axis equal
   grid on
   xlabel('T Rel Position (m)')
   ylabel('R Rel Position (m)')
   
   subplot(1,2,2)
   hold on
   plot(rtn_prop_hist(3,:),rtn_prop_hist(1,:))
   for j = 1:n_samples
%        rel_pos_cov_rtn_i = Sigma_rtn_hist(1:3,1:3,j);
       rel_pos_cov_rtn_i = Sigma_hcw_hist(1:3,1:3,j);

       sig_rr = rel_pos_cov_rtn_i(1,1);
       sig_rn = rel_pos_cov_rtn_i(1,3);
       sig_nn = rel_pos_cov_rtn_i(3,3);
       cov = [sig_nn, sig_rn; sig_rn, sig_rr];
       xy_ellipse = cov_ellipse([rtn_prop_hist(3,j);rtn_prop_hist(1,j)],cov,5.991,theta);
       plot(xy_ellipse(1,:), xy_ellipse(2,:),'g')
   end
   axis equal
   grid on
   xlabel('N Rel Position (m)')
   ylabel('R Rel Position (m)')
   
   
   clf
   
   
   
end

% %% at each estimate, propagate the orbit and check minimum separation
% dt = 30;
% t_max = 100;
% n_samples = floor(t_max / dt);
% confidence = 7.815; % 95%
% min_seps = zeros(1,n_samples);
% figure
% for i = 1:n_steps
%    x = x_hist(:,i);
%    Sigma = Sigma_hist(:,:,1);
%    Rot = R_RTN_to_ECI(x(1:6))';
%    R_full = [Rot zeros(3,3); zeros(3,3) Rot];
%    rel_rtn = R_full * (-x(1:6) + x(7:12));
%    A = rtn_to_roe_matrix(x(1:6), mu);
%    roe = A * rel_rtn;
%    Sigma_roe = A * Sigma * A';
%    
%    
%    
%    
%    
% end
