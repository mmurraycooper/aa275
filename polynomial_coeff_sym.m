syms x1_sq x2_sq x3_sq sig_11 sig_22 sig_33 h
d1 = (sig_11+h)^2;
d2 = (sig_22+h)^2;
d3 = (sig_33+h)^2;
eq = d1*d2*d3 - (x1_sq * d2 * d3 + x2_sq * d1 * d3 + x3_sq * d1 * d2);
c = coeffs(eq, h);
c = c([end:-1:1])