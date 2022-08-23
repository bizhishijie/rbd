function dist_2circle = Distance_plane2circle_(rc_center,n_vec,r_circle,n_wall,rc_wall)
%% distance between plane at rc_wall facing n_wall and circle at rc_center facing n_vec
% rc_center: position of the circle center
% n_vec: normal vector of the circle
% r_circle: radius of the circle
% n_wall: normal vector of the wall
% rc_wall: one position on the wall

%% variable name:
% w: wall
% c: circle center
% p: projection of center to plane

%%
dist_pc=abs(dot(rc_center-rc_wall,n_wall));
dist_r=sqrt(1-dot(n_vec,n_wall)^2)*r_circle;% sin(theta) = sin(acos(cos(nw*nc)))
dist_2circle=dist_pc-dist_r;

end

