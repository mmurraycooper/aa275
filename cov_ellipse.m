function  xy_ellipse = cov_ellipse(mu,Sigma,confidence,theta)
%confidence comes from 2dof chi squared 
[v,e] = eig(Sigma);
[~,max_ind] = max(max(abs(e)));
max_vec = v(:,max_ind);
alpha = atan2(max_vec(2),max_vec(1));
major_r = sqrt(confidence*max(max(e)));
minor_r = sqrt(confidence*min(max(e)));
xy_ellipse = [major_r * cos(theta); minor_r * sin(theta)];
xy_ellipse = ([cos(alpha) -sin(alpha); sin(alpha) cos(alpha)] * xy_ellipse)  + mu;
end