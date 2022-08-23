function v_vec = Direction_roll_wall_ (Rc_initial,Rc_base,idx_wall_xy)
% calculate the rolling direction of Rc_initial from Rc_base
%% INPUT
% Rc_initial: initial position
% Rc_base:  positions of base
% idx_wall_xy: 1 for touch x wall; 2 for y; (x wall means wall at +/- wall_x)
%% OUTPUT
% v_vec: vector of rolling direction from Rc_initial

g_vec=[0 0 -1]';

if min(Rc_base(3,:))>=Rc_initial(3)% SPECIAL case: if the base is higher than initial position
    v_vec=g_vec;
    return
end

%%
if idx_wall_xy==0% not influenced by wall
    base_num=size(Rc_base,2);
   
    switch base_num
        case 1 % roll from one base
            %%
            r_vec_1i=Rc_initial-Rc_base;% vector from base 1 to initial position
            if sum(r_vec_1i([1 2])~=0)==0% SPECIAL case: if Rc_initial is excatly on top of base
                th_roll_tmp=rand(1)*2*pi;% roll towards a random direction
                v_vec=[cos(th_roll_tmp) sin(th_roll_tmp) 0];
                return
            end
            %         if v_vec(3)==0
            %             th_roll_tmp=rand(1)*2*pi;
            %             v_vec=[cos(th_roll_tmp) sin(th_roll_tmp) 0];
            %         end
            v_vec=cross(r_vec_1i,cross(g_vec,r_vec_1i));
            v_vec=v_vec/norm(v_vec);% in common case, should already be right here
            if v_vec(3)>0% set to roll downwards (should be uncessary)
                v_vec=-v_vec;
            end
            
        case 2 % roll from two bases
            %%
            r_vec_1i=Rc_initial-Rc_base(:,1);% vector from base 1 to initial position
            r_vec_2i=Rc_initial-Rc_base(:,2);% vector from base 1 to initial position
            v_vec=cross(r_vec_1i,r_vec_2i);
            v_vec=v_vec/norm(v_vec);% normalization
            if v_vec(3)==0% % SPECIAL case: if balanced on the base; with wall, it usually happens
                v_vec=v_vec*(2*(rand(1)>0.5)-1);% set a random direction
                %%%%%%%%%%%%%%%%%%%%%%
            end
            if v_vec(3)>0% set to roll downwards (NECESSARY ! )
                v_vec=-v_vec;
            end
    end
    
else %influenced by wall
    %%
    idx_dir_xy=3-idx_wall_xy;% 1->2 ; 2->1
    
    r_vec_1i=Rc_initial-Rc_base;% vector from base 1 to initial position
    if r_vec_1i(idx_dir_xy)==0%% SPECIAL case: if leaning on the wall in the center direction; undetimined which direction to roll
        v_vec=zeros(3,1);
        v_vec(idx_dir_xy)=2*(rand(1)>0.5)-1;% set a random direction
        %%%%%%%%%%%%%%%%%%%%%%%
        return
    end
    if sum(r_vec_1i([1 2])~=0)==0%% SPECIAL case: if exactly on top of base
        LOGIC_away_wall=false;
        while LOGIC_away_wall==false
            th_roll_tmp=rand(1)*2*pi;% roll towards a random direction
            v_vec=[cos(th_roll_tmp) sin(th_roll_tmp) 0];
            if v_vec(idx_wall_xy)*Rc_initial(idx_wall_xy)<0% random direction move the sphere away from the wall
                LOGIC_away_wall=true;
            end
        end
        return
    end
    
    r_vec_1i(idx_wall_xy)=0;% set the direction parrallel to wall
    wall_vec=zeros(3,1);
    wall_vec(idx_wall_xy)=1;% direction of the wall
    v_vec=cross(wall_vec,r_vec_1i);% roll along the walll direction
    
    v_vec=v_vec/norm(v_vec);
    if v_vec(3)>0% set to roll downwards (should be uncessary)
        v_vec=-v_vec;
    end
end

end

