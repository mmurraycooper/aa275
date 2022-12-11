function A = STM(x,mu)
scale = -mu / (norm(x(1:3))^3);
v_block = diag(scale * ones(1,3));
A = [zeros(3,3), eye(3,3); v_block, zeros(3,3)]; 
end