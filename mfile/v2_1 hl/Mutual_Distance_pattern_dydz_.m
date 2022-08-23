function delta = Mutual_Distance_pattern_dydz_ (d,Rc1,Rc2)

Rc1(3,:)=Rc1(3,:)+d;
dz_tmp=repmat(Rc1(3,:)',1,size(Rc2,2))-repmat(Rc2(3,:),size(Rc1,2),1);
dy_tmp=repmat(Rc1(2,:)',1,size(Rc2,2))-repmat(Rc2(2,:),size(Rc1,2),1);
dist_yz_tmp=sqrt(dz_tmp.^2+dy_tmp.^2);
% dz_tmp=abs(dz_tmp);
delta=sqrt(mean(min(dist_yz_tmp,[],2).^2));
end