function R_rot = R_RTN_to_ECI(x_c)

r = x_c(1:3);
v = x_c(4:6);

R = r / norm(r);
N = cross(r,v) / norm(cross(r,v));
T = cross(N,R);

R_rot = [R T N];
end
