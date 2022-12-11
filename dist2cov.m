function d = dist2cov(mu,Sigma,confidence)
x_p = [0;0;0];
x_hat = mu;
[U,S,V] = svd(Sigma*(confidence));
S_sqrt = sqrtm(S);
x_tilde = S_sqrt * V' * (x_p - x_hat);
lam = lam_max(x_tilde(1)^2, x_tilde(2)^2, x_tilde(3)^2, S(1,1), S(2,2), S(3,3));
x = x_hat + U * S_sqrt / (S + lam*eye(3,3)) * S_sqrt * V' * (x_p - x_hat);
d = norm(x);
end