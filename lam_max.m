function lam = lam_max(x1_sq, x2_sq, x3_sq, sig_11, sig_22, sig_33)
coeffs  = [1, 2*sig_11 + 2*sig_22 + 2*sig_33, 4*sig_11*sig_22 - x2_sq - x3_sq - x1_sq + sig_11^2 + sig_22^2 + sig_33^2 + 2*sig_33*(2*sig_11 + 2*sig_22), 2*sig_33*(sig_11^2 + 4*sig_11*sig_22 + sig_22^2) + sig_33^2*(2*sig_11 + 2*sig_22) - 2*sig_11*x2_sq - 2*sig_11*x3_sq - 2*sig_22*x1_sq - 2*sig_22*x3_sq - 2*sig_33*x1_sq - 2*sig_33*x2_sq + 2*sig_11*sig_22^2 + 2*sig_11^2*sig_22, sig_11^2*sig_22^2 + 2*sig_33*(2*sig_11^2*sig_22 + 2*sig_11*sig_22^2) + sig_33^2*(sig_11^2 + 4*sig_11*sig_22 + sig_22^2) - sig_11^2*x2_sq - sig_11^2*x3_sq - sig_22^2*x1_sq - sig_22^2*x3_sq - sig_33^2*x1_sq - sig_33^2*x2_sq - 4*sig_11*sig_22*x3_sq - 4*sig_11*sig_33*x2_sq - 4*sig_22*sig_33*x1_sq, sig_33^2*(2*sig_11^2*sig_22 + 2*sig_11*sig_22^2) - 2*sig_11*sig_22^2*x3_sq - 2*sig_11^2*sig_22*x3_sq - 2*sig_11*sig_33^2*x2_sq - 2*sig_11^2*sig_33*x2_sq - 2*sig_22*sig_33^2*x1_sq - 2*sig_22^2*sig_33*x1_sq + 2*sig_11^2*sig_22^2*sig_33, sig_11^2*sig_22^2*sig_33^2 - x3_sq*sig_11^2*sig_22^2 - x2_sq*sig_11^2*sig_33^2 - x1_sq*sig_22^2*sig_33^2];
lams = roots(coeffs);
lams = lams(imag(lams)==0);
lam = max(lams);
end
