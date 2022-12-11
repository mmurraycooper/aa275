function x_next = dynamics(x,mu,dt)
x_next = dynamics_jac(x,mu,dt) * x;
end
