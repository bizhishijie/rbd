%% calculate the final position of rolling sphere separate with the base(s)
% by xiachj 20191229

%%
base_num=length(idx_touch);
switch base_num
    case 1% roll from one base
        v_vec_tmp=Rc_initial-rc_center;
        v_vec_tmp(3)=0;% horizontal direction
        Rc_final=rc_center+v_vec_tmp/norm(v_vec_tmp)*r_circle;% new position
        idx_touch_new=[];% separate from the base
    case 2% roll from two base
        Rc_final=rc_center+y_vec*r_circle;% new position
        idx_touch_new=idx_touch(Rc_old(3,idx_touch)==min(Rc_old(3,idx_touch)));% keep touch with lower base
        if length(idx_touch_new)==2% if the two base are horizontal
            idx_touch_new=[];% separate from both bases
        end
end
