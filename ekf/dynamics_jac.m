function F = dynamics_jac(x,mu,dt)
x_c = x(1:6); 
x_d = x(7:12);
F = [eye(6,6) + dt * STM(x_c,mu), zeros(6,6), zeros(6,1), zeros(6,1);
    zeros(6,6), eye(6,6) + dt * STM(x_d,mu), zeros(6,1), zeros(6,1);
    zeros(1,6), zeros(1,6), 1, 0;
     zeros(1,6), zeros(1,6), 0, 1;];
end