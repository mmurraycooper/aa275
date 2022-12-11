% // Linearized mapping matrix from RTN to quasi-nonsingular ROE
function A = rtn_to_roe_matrix(x_c_eci, mu)
r = x_c_eci(1:3);
v = x_c_eci(4:6);
h = cross(r, v);
hh = h / norm(h);
u = atan2(r(3), -r(1) * hh(2) + r(2) * hh(1));
a =  1 / (2 / norm(r) - norm(v) * norm(v) / mu);
ada = [4., 0., 0., 0., 2., 0.];
adlam = [0., 1., 0., -2., 0., 0.];
adex = [3. * cos(u), 0., 0., sin(u), 2. * cos(u), 0.];
adey = [3. * sin(u), 0., 0., -cos(u), 2. * sin(u), 0.];
adix = [0., 0., sin(u), 0., 0., cos(u)];
adiy = [0., 0., -cos(u), 0., 0., sin(u)];
A = [ada; adlam; adex; adey; adix; adiy] * diag([1., 1., 1., a / v, a / v, a / v]);
end