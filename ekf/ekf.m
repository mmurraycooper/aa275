function [x, Sigma] = ekf(x, Sigma, meas, dt, mu, sat_data, Q, R)
x_bar = dynamics(x,mu,dt);
F = dynamics_jac(x,mu,dt);
H = meas_jac(x,sat_data);
Sigma_bar = F * Sigma * F' + R;
K = Sigma_bar * H' / (H * Sigma_bar * H' + Q);
x = x_bar + K * (meas - meas_pred(x_bar,sat_data));
Sigma = (eye(length(x)) - K*H) * Sigma_bar;
end