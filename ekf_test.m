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

b_c_true = 110;
b_d_true = 178;
x_init = [eci_c(:,1); eci_d(:,1); b_c_true; b_d_true];
num_states = length(x_init);

x = x_init + mvnrnd(zeros(1,num_states), diag(50*ones(1,num_states)))'; % add some noise
Sigma = 10 * eye(num_states,num_states);
R = 10 * eye(num_states,num_states);
Q = 10 * eye(2*num_gps,2*num_gps);

x_hist = zeros(num_states,n_steps);
Sigma_hist = zeros(num_states,num_states,n_steps); 

x_hist(:,1) = x;
Sigma_hist(:,:,1) = Sigma;

for i = 2:n_steps
    x_true = [eci_c(:,i);eci_d(:,i); b_c_true; b_d_true];
    meas_true = meas_pred(x_true,sat_data);
    
    meas = meas_true + mvnrnd(zeros(1,2*num_gps), diag(10*ones(1,2*num_gps)))'; % add some noise
    dt = t(i) - t(i-1);
    
%     x_bar = dynamics(x,mu,dt);
%     F = dynamics_jac(x,mu,dt);
%     H = meas_jac(x,sat_data);
%     Sigma_bar = F * Sigma * F' + R;
% %     K = Sigma_bar * H' * inv(H * Sigma_bar * H' + Q);
%     K = Sigma_bar * H' / (H * Sigma_bar * H' + Q);
%     x = x_bar + K * (meas - meas_pred(x_bar,sat_data));
%     Sigma = (eye(num_states) - K*H) * Sigma_bar;

    [x, Sigma] = ekf(x, Sigma, meas, dt, mu, sat_data, Q, R);
    x_hist(:,i) = x;
    Sigma_hist(:,:,i) = Sigma;
end

% figure
% plot(eci_c(1,:),eci_c(2,:));
% hold on
% plot(x_hist(1,:), x_hist(2,:),'.')
% legend('True','Estimated')
% grid on
% axis equal
% title('X-Y Chief')
% 
% figure
% plot(eci_c(1,:),eci_c(3,:));
% hold on
% plot(x_hist(1,:), x_hist(3,:),'.')
% legend('True','Estimated')
% grid on
% axis equal
% title('X-Z Chief')
% 
% figure
% plot(eci_d(1,:),eci_d(2,:));
% hold on
% plot(x_hist(7,:), x_hist(8,:),'.')
% legend('True','Estimated')
% grid on
% axis equal
% title('X-Y Dep')
% 
% figure
% plot(eci_d(1,:),eci_d(3,:));
% hold on
% plot(x_hist(7,:), x_hist(9,:),'.')
% legend('True','Estimated')
% grid on
% axis equal
% title('X-Z Dep')

rel_eci = x_hist(7:12,:) - x_hist(1:6,:);
rr_est = zeros(1,n_steps);
rt_est = zeros(1,n_steps);
rn_est = zeros(1,n_steps);
rel_pos_cov_rtn = zeros(3,3,n_steps);
for i = 1:n_steps
    x_c = x_hist(1:6,i);
    x_d = x_hist(7:12,i);
    rel_eci_i = rel_eci(:,i);
    R =  R_RTN_to_ECI(x_c)';
    rel_rtn = R * rel_eci(1:3,i);
    rr_est(i) = rel_rtn(1);
    rt_est(i) = rel_rtn(2);
    rn_est(i) = rel_rtn(3);
    rel_pos_cov_rtn(:,:,i) = R * (Sigma(1:3,1:3) + Sigma(7:9,7:9)) * R';
end

chi_vec = prop_data.chi_vec;
rr_true = chi_vec(1,:);
rt_true = chi_vec(2,:);
rn_true = chi_vec(3,:);

figure
subplot(1,2,1)
hold on
plot(rt_true,rr_true)
plot(rt_est,rr_est,'.')
theta = linspace(0,2*pi);
for i = 1:n_steps
    rel_pos_cov_rtn_i = rel_pos_cov_rtn(:,:,i);
    sig_rr = rel_pos_cov_rtn_i(1,1);
    sig_rt = rel_pos_cov_rtn_i(1,2);
    sig_tt = rel_pos_cov_rtn_i(2,2);
    cov = [sig_tt, sig_rt; sig_rt, sig_rr];
    xy_ellipse = cov_ellipse([rt_est(i); rr_est(i)],cov,5.991,theta);
    plot(xy_ellipse(1,:), xy_ellipse(2,:),'g')
end
axis equal
grid on

subplot(1,2,2)
hold on
plot(rn_true,rr_true)
plot(rn_est,rr_est,'.')
for i = 1:n_steps
    rel_pos_cov_rtn_i = rel_pos_cov_rtn(:,:,i);
    sig_rr = rel_pos_cov_rtn_i(1,1);
    sig_rn = rel_pos_cov_rtn_i(1,3);
    sig_nn = rel_pos_cov_rtn_i(3,3);
    cov = [sig_nn, sig_rn; sig_rn, sig_rr];
    xy_ellipse = cov_ellipse([rn_est(i); rr_est(i)],cov,5.991,theta);
    plot(xy_ellipse(1,:), xy_ellipse(2,:),'g')
end
axis equal
grid on