function H = meas_jac(x,sat_data)
% c = 299792458; # just used time*c as the quantity estimated
x_c = x(1:3);
x_d = x(7:9);
num_gps = size(sat_data,1);

diff_xc = x_c' - sat_data(:,1:3);
diff_xd = x_d' - sat_data(:,1:3);

norm_xc = vecnorm(diff_xc,2,2);
norm_xd = vecnorm(diff_xd,2,2);

H = [ diff_xc ./ norm_xc, zeros(num_gps,3), zeros(num_gps,6), ones(num_gps,1), zeros(num_gps,1);
    zeros(num_gps,6), diff_xd ./ norm_xd, zeros(num_gps,3), zeros(num_gps,1), ones(num_gps,1)];
end