function meas = meas_pred(x,sat_data)
% c = 299792458; # just used time*c as the quantity estimated
x_c = x(1:3);
x_d = x(7:9);
b_c = x(13);
b_d = x(14);
 
meas = [vecnorm(sat_data(:,1:3) - x_c',2,2) + (b_c - sat_data(:,end));
    vecnorm(sat_data(:,1:3) - x_d',2,2) + (b_d - sat_data(:,end))];
end