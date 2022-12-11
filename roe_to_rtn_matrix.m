function A = roe_to_rtn_matrix(x_c_eci,mu)
r = norm(x_c_eci(1:3));
v = norm(x_c_eci(4:6));
h = cross(r, v);
hh = h / norm(h);
u = atan2(r(3), -r(1) * hh(2) + r(2) * hh(1));
a =  1 / (2 / norm(r) - norm(v) * norm(v) / mu);

rr = [1., 0., -cos(u), -sin(u), 0., 0.];
rt = [0., 1., 2. * sin(u), -2. * cos(u), 0., 0.];
rn = [0., 0., 0., 0., sin(u), -cos(u)];
vr = [0., 0., v / a * sin(u), -v / a * cos(u), 0., 0.];
vt = [-3. / 2. * v / a,    0., 2. * v / a * cos(u),2. * v / a * sin(u), 0., 0.];
vn = [0., 0., 0., 0., v / a * cos(u), v / a * sin(u)];
A = [rr; rt; rn; vr; vt; vn];
end
