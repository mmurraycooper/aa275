function rtn = HCW_dynamics(rtn_init,w,t)
STM = [4-3*cos(w*t), 0, 0, 1/w*sin(w*t), 2/w*(1-cos(w*t)), 0;
    6*sin(w*t)-6*w*t, 1, 0, 2/w*(cos(w*t)-1), 4/w*sin(w*t)-3*t, 0;
    0, 0, cos(w*t), 0, 0, 1/w*sin(w*t);
    3*w*sin(w*t), 0, 0, cos(w*t), 2*sin(w*t), 0;
    6*w*cos(w*t)-6*w, 0, 0, -2*sin(w*t), 4*cos(w*t)-3, 0;
    0, 0, -w*sin(w*t), 0, 0, cos(w*t)];
rtn = (STM*rtn_init')';
end
