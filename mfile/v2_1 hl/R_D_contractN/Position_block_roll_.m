function [rc_touch,theta_touch] = Position_block_roll_(n_vec,r_circle,rc_n)
% calcualat the block position of sphere rolling on cirlce by one sphere at rc_n
% by xiachj 20191229
%% INPUT
% n_vec: normal vector of the circle
% r_circle:radius of the circle
% rc_n: RELATIVE position from center of the circle to a close sphere
%% OUTPUT
% rc_touch: position where the sphere touch rc_n
% theta_touch: theta where the sphere touch rc_n

%% name of variables:
% c: center of circle
% t: touch position
% n: rc_n
% p: projection of n to circle plane


%%
g_vec=[0 0 -1]';

%%
x_vec=n_vec;% set x axis to be the normal vector of circle
y_vec=cross(x_vec,g_vec);% set y in horizontal direction
y_vec=y_vec/norm(y_vec);
z_vec=cross(x_vec,y_vec);

x_coor=dot(x_vec,rc_n);% relative coordinates
y_coor=dot(y_vec,rc_n);
z_coor=dot(z_vec,rc_n);
[~,theta_p,dist_pc]=cart2sph(0,y_coor,z_coor);% projection theta and distance
theta_p=pi/2-theta_p;
if y_coor<0
    theta_p=2*pi-theta_p;
end

dist_pt=sqrt(1-x_coor^2);% distance projection to touch position
cos_delta_theta=(r_circle^2+dist_pc^2-dist_pt^2)/2/r_circle/dist_pc;% geometry
delta_theta=acos(cos_delta_theta);
theta_touch=[theta_p+delta_theta theta_p-delta_theta];% upper and lower touch theta

% theta_touch(theta_touch<0)=theta_touch(theta_touch<0)+2*pi;% limit theta value between [0 2*pi]
% theta_touch(theta_touch>2*pi)=theta_touch(theta_touch>2*pi)-2*pi;
theta_touch=mod(real(theta_touch),2*pi);% faster; use real to avoid error

rc_touch=y_vec*r_circle*sin(theta_touch)+z_vec*r_circle*cos(theta_touch);% upper and lower touch psotion

end

