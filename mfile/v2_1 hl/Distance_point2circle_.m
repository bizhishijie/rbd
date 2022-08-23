function dist_2circle = Distance_point2circle_(rc_center,n_vec,r_circle,rn_all)
% rc_center: position of the circle center
% n_vec: normal vector of the circle
% rn_all: positions of other particles

%%
%dist_nc=pdist2(rc_center',rn_all');% slower
dist_nc=sqrt((rn_all(1,:)-rc_center(1)).^2+(rn_all(2,:)-rc_center(2)).^2+(rn_all(3,:)-rc_center(3)).^2);% faster
% dist_np=dot(rn_all-repmat(rc_center,1,size(rn_all,2)),repmat(n_vec,1,size(rn_all,2)));
dist_np=(rn_all(1,:)-rc_center(1))*n_vec(1)+(rn_all(2,:)-rc_center(2))*n_vec(2)+(rn_all(3,:)-rc_center(3))*n_vec(3);% faster
dist_2circle=sqrt((sqrt(dist_nc.^2-dist_np.^2)-r_circle).^2+dist_np.^2);

end

