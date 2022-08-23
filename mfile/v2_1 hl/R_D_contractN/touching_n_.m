function [p] = touching_n_(Rc_sp,Rc_posi,tole_degree)
pp=zeros(1,6);
idx_nowall=abs(Rc_sp(2,:))<((max(Rc_sp(2,:))-min(Rc_sp(2,:)))/2)-2*eps;
Rc_sp=Rc_sp(:,idx_nowall);
for ii=1:size(Rc_sp,2)
    r=Rc_sp(:,ii);
    Rc_tmp=Rc_posi(:,Rc_posi(3,:)<r(3)+3&Rc_posi(3,:)>r(3)-3);
%     Rc_tmp(:,ii)=[];
    d=pdist2(r',Rc_tmp');
    pt=length(find(d<=1+tole_degree*eps))-1;%remove itself
%     pt=length(find(d<=1+tole_degree*eps));
    pp(pt)=pp(pt)+1;
end
% p=pp./sum(pp);
p=pp;
end