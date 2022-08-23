function [rc_touch,theta_touch] = Position_block_wall_(rc_center,n_vec,r_circle,n_wall,rc_wall)
% calcualat the block position of sphere rolling on cirlce by plane at rc_wall facing n_wall
% by xiachj 20191229
%% INPUT
% rn_center: center of circle
% n_vec: normal vector of the circle
% r_circle:radius of the circle
% n_wall: normal vector of wall facing inward
% rc_all: one point on the wall
%% OUTPUT
% rc_touch: position where the sphere touch rc_n
% theta_touch: theta where the sphere touch rc_n

%% name of variables:
% w: wall
% c: circle center
% p: projection of center to plane
% r: ring


%%
g_vec=[0 0 -1]';

%%
x_vec=n_vec;% set x axis to be the normal vector of circle
y_vec=cross(x_vec,g_vec);% set y in horizontal direction
y_vec=y_vec/norm(y_vec);
z_vec=cross(x_vec,y_vec);

theta_p=atan(dot(y_vec,n_wall)/dot(z_vec,n_wall));% equation obtained based on vector geometry
if isnan(theta_p)
    rc_touch=[];
    theta_touch=[];
    return
end
if theta_p>=0% modify atan range value
    theta_p=[theta_p theta_p+pi];% second solution of atan
else
    theta_p=[theta_p+pi theta_p+2*pi];
end

% dist_r=dot(y_vec*sin(theta_p)+z_vec*cos(theta_p),repmat(n_wall,1,2));
dist_r=dot(y_vec*sin(theta_p)+z_vec*cos(theta_p),[n_wall n_wall]);% faster
theta_p=theta_p(dist_r==min(dist_r));

dist_pc=abs(dot(rc_center-rc_wall,n_wall));
theta_cw=acos(abs(dot(n_vec,n_wall)));
diff_theta=acos((dist_pc-0.5)/sin(theta_cw)/r_circle);

theta_touch=[theta_p-diff_theta theta_p+diff_theta];
% theta_touch(theta_touch<0)=theta_touch(theta_touch<0)+2*pi;% limit theta value between [0 2*pi]
% theta_touch(theta_touch>=2*pi)=theta_touch(theta_touch>=2*pi)-2*pi;
theta_touch=mod(real(theta_touch),2*pi);% faster; use real for error in rare cases

rc_touch=y_vec*r_circle*sin(theta_touch)+z_vec*r_circle*cos(theta_touch)+rc_center;% upper and lower touch psotion

%% modify rc_touch to eliminiate finite accuracy problem
rc_touch(n_wall~=0,:)=rc_wall(n_wall~=0)+n_wall(n_wall~=0)*0.5;

end

