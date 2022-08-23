function [Rc_final,idx_touch_new,LOGIC_stable] = RD_mono_3d_wall_ (Rc_initial,idx_touch_old,Rc_old,wall_range,tole_degree)
% drop or roll one sphere from initial position to final position along the existing packing
% RD: Random Deposition, and also, Roll and Drop
% by xiachj 20191229

%% input:
% Rc_initial: current position of the dropping sphere
% Rc_old: positions of existing spheres
% idx_touch_old: index of Rc_old that in contact with Rc_initial
% tole_degree: degree of tolerance of finite accuracy, in unit of eps
%% output:
% Rc_final: new position of the dropping sphere
% idx_touch_new: index of Rc_old that in contact with Rc_final
% LOGIC_stable: to indicate current position is stable or not

%%
g_vec=[0 0 -1]';% gravity vector
LOGIC_stable=false;% defaulf value, not stable

wall_x=wall_range(1);
wall_y=wall_range(2);
LOGIC_wall_x=false;
LOGIC_wall_y=false;
if abs(abs(Rc_initial(1))-(wall_x-0.5))/eps<tole_degree
    LOGIC_wall_x=true;
end
if abs(abs(Rc_initial(2))-(wall_y-0.5))/eps<tole_degree
    LOGIC_wall_y=true;
end
idx_wall_xy=1*LOGIC_wall_x+2*LOGIC_wall_y;% 1 for x 2 for y; 0 for no touching wall

%% drop to bottom
if isempty(Rc_old)% if spheres exist
    Rc_final=[Rc_initial([1 2]);0.5];% drop to the bottom and stop
    LOGIC_stable=true;
    idx_touch_new=[];
    return
end

%% find first touch below
if isempty(idx_touch_old)% if sphere touches no one, it drop freely
    dist_xy=sqrt(sum((Rc_old([1 2],:)-repmat(Rc_initial([1 2]),1,size(Rc_old,2))).^2,1));% relative horizontal distance between all sphere and dropping sphere
    dist_z=Rc_old(3,:)-Rc_initial(3);% relative horizontal distance between all sphere and dropping sphere
    idx_below=find(dist_xy<1&dist_z<0);% possible spheres that may touch the dropping sphere from below
    if isempty(idx_below)% if no spheres exist below
        Rc_final=[Rc_initial([1 2]);0.5];% drop to the bottom and stop
        LOGIC_stable=true;
        idx_touch_new=[];
        return
    else% if one or more sphere may touch the dropping sphere from below
        z_touch_candidate=Rc_old(3,idx_below)+sqrt(1-dist_xy(idx_below).^2);% height of touching point of dropping sphere with spheres below
        idx_touch_old=idx_below(z_touch_candidate==max(z_touch_candidate));% find the highest
        %         idx_touch_old=idx_below(z_touch_possible>=max(z_touch_possible));
        Rc_initial=[Rc_initial([1 2]);max(z_touch_candidate)];% drop the sphere to the highest touch point; use max() because may have more than one value
    end
end

if LOGIC_wall_x&&LOGIC_wall_y% if at the cornor of box; and hold by a below sphere
    Rc_final=Rc_initial;
    idx_touch_new=idx_touch_old;
    LOGIC_stable=true;
    return
end

%% calculate roll direction
v_vec1_all=zeros(3,length(idx_touch_old));% direction rolling down on ONE particle

for ii=1:length(idx_touch_old)% roll direction
    v_vec1_tmp=Direction_roll_wall_ (Rc_initial,Rc_old(:,idx_touch_old(ii)),idx_wall_xy);% rolling direction
    v_vec1_all(:,ii)=v_vec1_tmp;
end

if length(idx_touch_old)<2% if only one base
    idx_touch_pair=[];
    v_vec2_all=[];
else
    idx_touch_pair=nchoosek(idx_touch_old,2);% each pair of base sphere
    v_vec2_all=zeros(3,size(idx_touch_pair,1));
    for ii=1:size(idx_touch_pair,1)
        v_vec2_tmp=Direction_roll_wall_ (Rc_initial,Rc_old(:,idx_touch_pair(ii,:)),0);% rolling direction
        v_vec2_all(:,ii)=v_vec2_tmp;
    end
end


%% determine feasible roll direction
r_vec_zi=repmat(Rc_initial,1,length(idx_touch_old))-Rc_old(:,idx_touch_old);% vector from all base to initial position

LOGIC_overlap1_all=true(1,size(v_vec1_all,2));% one means rolling down will cause immediate overlap
for ii=1:size(v_vec1_all,2)
    diff_dist_tmp=dot(r_vec_zi,repmat(v_vec1_all(:,ii),1,length(idx_touch_old)));% differential distance to all base if roll towards this direction
    diff_dist=diff_dist_tmp(1:size(v_vec1_all,2)~=ii);% remove the current base
    if sum(diff_dist<0)==0% if no overlap; diff_dist<0 means will cause overlap
        LOGIC_overlap1_all(ii)=false;
    end
end

LOGIC_overlap2_all=true(1,size(v_vec2_all,2));
for ii=1:size(idx_touch_pair,1)
    if LOGIC_wall_x||LOGIC_wall_y
        if v_vec2_all(idx_wall_xy,ii)*Rc_initial(idx_wall_xy)>=0% if will overlap with wall; >= is important, for wall_x=0.5 special case
            LOGIC_overlap2_all(ii)=true;% set this direction overlap true (unnecessary: already true)
            continue% NECESSARY: don't judge other spheres, since already touch wall
        end
    end

    diff_dist_tmp=dot(r_vec_zi,repmat(v_vec2_all(:,ii),1,length(idx_touch_old)));% differential distance to all base if roll towards this direction
    [~,idx_tmp_overlap]=ismember(idx_touch_pair(ii,:),idx_touch_old);% remove the current base
    diff_dist=diff_dist_tmp(~ismember(1:length(idx_touch_old),idx_tmp_overlap));
    if isempty(diff_dist)% no third base (unnecessary)
        LOGIC_overlap2_all(ii)=false;% no overlap
        continue
    end
    if sum(diff_dist<0)==0% if no overlap; diff_dist<0 means will cause overlap
        LOGIC_overlap2_all(ii)=false;
    end
end

%% determine steepest roll direction
v_vec_all=[v_vec1_all v_vec2_all];
LOGIC_overlap_all=[LOGIC_overlap1_all LOGIC_overlap2_all];
idx_touch_all=[repmat(idx_touch_old,2,1) idx_touch_pair'];
v_vec_all(3,LOGIC_overlap_all==true)=inf;
if sum(~isinf(v_vec_all(3,:)))==0% no feasible roll direction
    Rc_final=Rc_initial;
    idx_touch_new=idx_touch_old;
    LOGIC_stable=true;
    return
end
idx_tmp_roll=find(v_vec_all(3,:)==min(v_vec_all(3,:)));
idx_touch=unique(idx_touch_all(:,idx_tmp_roll(1))');% (1) means if roll from one and two base are same, use one base


%% define the rolling circle, and find all possible sphere to block the rolling
theta_max=pi/2;
if length(idx_touch)==1% roll from one sphere
    if LOGIC_wall_x||LOGIC_wall_y
        n_vec=zeros(3,1);
        n_vec(idx_wall_xy)=1;
        rc_center=Rc_old(:,idx_touch);
        rc_center(idx_wall_xy)=Rc_initial(idx_wall_xy);
        r_circle=norm(rc_center-Rc_initial);
        theta_initial=acos((Rc_initial(3)-Rc_old(3,idx_touch))/r_circle);

        y_vec=cross(n_vec,g_vec);
        y_coor=dot(y_vec,Rc_initial-rc_center);
        if y_coor<0
            n_vec=-n_vec;
        end
    else
        theta_initial=acos(Rc_initial(3)-Rc_old(3,idx_touch));% inital theta before rolling
        rc_center=Rc_old(:,idx_touch);%center of the circle over which the sphere roll
        n_vec=cross(g_vec,Rc_initial-rc_center);%normal vector of the circle over which the sphere roll
        n_vec=n_vec/norm(n_vec);
        r_circle=1;%radius of the circle over which the sphere roll
    end

    dist_2circle=Distance_point2circle_(rc_center,n_vec,r_circle,Rc_old);% distance from all point to circle
    idx_block_all=find(dist_2circle<1);
    idx_block_all=idx_block_all(idx_block_all~=idx_touch);% remove the base (unnecessary ?)
    if ~(LOGIC_wall_x||LOGIC_wall_y)
        idx_block_all=idx_block_all(~ismember(idx_block_all,idx_touch_old));% if one leaves a top base and roll down on the lower base, should not touch the top base again % (right ? yes for no wall, no for wall)
    end

elseif length(idx_touch)==2
    %%
    rc_center=mean(Rc_old(:,idx_touch),2);%center of the circle over which the sphere roll
    n_vec=Rc_old(:,idx_touch(1))-Rc_old(:,idx_touch(2));%normal vector of the circle over which the sphere roll
    n_vec=n_vec/norm(n_vec);
    dist_base12=norm(Rc_old(:,idx_touch(1))-Rc_old(:,idx_touch(2)));
    r_circle=sqrt(1-(dist_base12/2)^2);%radius of the circle over which the sphere roll

    x_vec=n_vec;
    y_vec=cross(x_vec,g_vec);
    y_vec=y_vec/norm(y_vec);
    y_coor=dot(y_vec,Rc_initial-rc_center);
    if y_coor<0% set initial position at positive y direction
        n_vec=-n_vec;x_vec=-x_vec;y_vec=-y_vec;
    elseif y_coor==0&&(LOGIC_wall_x||LOGIC_wall_y)
        if dot(y_vec,Rc_initial)>0% if defined to roll outward
            n_vec=-n_vec;x_vec=-x_vec;y_vec=-y_vec;
        end
    end
    z_vec=cross(x_vec,y_vec);
    [~,theta_initial,~]=cart2sph(dot(x_vec,Rc_initial-rc_center),dot(y_vec,Rc_initial-rc_center),dot(z_vec,Rc_initial-rc_center));
    theta_initial=pi/2-theta_initial;% inital theta before rolling

    dist_2circle=Distance_point2circle_(mean(Rc_old(:,idx_touch),2),n_vec,r_circle,Rc_old);
    idx_block_all=find(dist_2circle<1);
    idx_block_all=idx_block_all(~ismember(idx_block_all,idx_touch));% remove base (unnecessary)

end

%% block position by wall
wall_vec_loop=[1 0 0;-1 0 0;0 1 0;0 -1 0]';
wall_pos_loop=[-wall_x 0 0;wall_x 0 0 ;0 -wall_y 0 ;0 wall_y 0]';
dist_2plane=zeros(1,4);
theta_touch_plane=inf(2,4);
for ii=1:4
    n_wall=wall_vec_loop(:,ii);
    rc_wall=wall_pos_loop(:,ii);
    dist_2plane_tmp=Distance_plane2circle_(rc_center,n_vec,r_circle,n_wall,rc_wall);
    dist_2plane(:,ii)=dist_2plane_tmp;
    if dist_2plane_tmp<0.5
        [~,theta_touch_tmp]=Position_block_wall_(rc_center,n_vec,r_circle,n_wall,rc_wall);
        theta_touch_plane(:,ii)=theta_touch_tmp';
    end
end
theta_touch_plane(theta_touch_plane==0)=inf;%% zero degree is unstable; if v_vec allow the sphere to roll
theta_touch_plane(theta_touch_plane<=theta_initial-eps*tole_degree|theta_touch_plane>theta_max)=inf;% remove unphysical theta % <= very important; % tolerence: some sphere may actually block the rolling but not included due to finite accuracy
LOGIC_block_wall=true;
if sum(~isinf(theta_touch_plane(:)))==0% if no wall can block movement on the ring
    LOGIC_block_wall=false;
else
    [idx_top_bottom,idx_tmp_block]=find(theta_touch_plane==min(theta_touch_plane(:)));% find the top touch position
    [rc_block_tmp,theta_touch_tmp]=Position_block_wall_(rc_center,n_vec,r_circle,wall_vec_loop(:,idx_tmp_block(1)),wall_pos_loop(:,idx_tmp_block(1)));% idx_block(1) as different idx_block() correspond to same block position
    Rc_final_wall=rc_block_tmp(:,idx_top_bottom(1));% block position . (1) is for rare case that two equivalent contact formes
    theta_final_wall=theta_touch_tmp(idx_top_bottom(1));
end


%% roll to final position
if isempty(idx_block_all)% if no sphere can block movement on the ring
    if LOGIC_block_wall==false
        Position_separate_base_% move the sphere to the spearation position with the base
    else
        Rc_final=Rc_final_wall;
        idx_touch_new=idx_touch;
    end
    return
end

theta_touch=zeros(2,length(idx_block_all));% all theta touching
for ii=1:length(idx_block_all)
    [~,theta_touch_tmp]=Position_block_roll_(n_vec,r_circle,Rc_old(:,idx_block_all(ii))-rc_center);
    theta_touch(:,ii)=theta_touch_tmp';
end
theta_touch(theta_touch<=theta_initial-eps*tole_degree|theta_touch>theta_max)=inf;% remove unphysical theta % <= very important; % tolerence: some sphere may actually block the rolling but not included due to finite accuracy

if sum(ismember(idx_block_all,idx_touch_old))~=0% if original touch sphere may block the rolling; % only useful for rolling from two bases
    idx_tmp_block=ismember(idx_block_all,idx_touch_old)==1;% get the idex
    theta_touch_tmp=theta_touch(:,idx_tmp_block);% get the theta
    theta_touch_tmp(abs((theta_touch_tmp-theta_initial))/eps<tole_degree)=inf;% remove the same initial point (upper one); let it roll to lower block position
    theta_touch(:,idx_tmp_block)=theta_touch_tmp;
end

if sum(~isinf(theta_touch(:)))==0% if no sphere can block movement on the ring
    if LOGIC_block_wall
        Rc_final=Rc_final_wall;
        idx_touch_new=idx_touch;
    else
        Position_separate_base_% move the sphere to the spearation position with the base
    end
else
    %     Position_block_top_% idx_touch_new defined here
    [idx_top_bottom,idx_tmp_block]=find(theta_touch==min(theta_touch(:)));% find the top touch position
    idx_block=idx_block_all(idx_tmp_block);
    [rc_block_tmp,theta_touch_tmp]=Position_block_roll_(n_vec,r_circle,Rc_old(:,idx_block(1))-rc_center);% idx_block(1) as different idx_block() correspond to same block position
    Rc_final=rc_block_tmp(:,idx_top_bottom(1))+rc_center;% block position . (1) is for rare case that two equivalent contact formes
    theta_final=theta_touch_tmp(idx_top_bottom);
    %     idx_touch_new=unique([idx_touch idx_block]);% touch index; original touch ones plus the new block one; unique: (uncessary)
    idx_touch_new=[idx_touch idx_block];% faster;

    if LOGIC_block_wall
        if theta_final_wall<=theta_final+eps*tole_degree% if theta of rolling down a sphere is very close to theta rolling along a surface; set the value to surface value which is more accurate
            Rc_final=Rc_final_wall;
            idx_touch_new=idx_touch;
        end
    end

    if abs(theta_final-theta_initial)/eps<tole_degree% if theta change very little, means it does not truely left the oiginal contact spheres
        idx_touch_new=unique([idx_touch_new idx_touch_old]);
    end
end


end

