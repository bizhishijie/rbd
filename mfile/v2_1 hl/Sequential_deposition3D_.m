function [Rc_old,tole_degree_loop]=Sequential_deposition3D_(wx,wy,N_total)
%to create a single N number deposition

wall_x=wx/2;%0.515
wall_y=wy/2;% 3
Rc_initial_all=[(rand(1,N_total)-0.5)*(2*wall_x-1);(rand(1,N_total)-0.5)*(2*wall_y-1);inf(1,N_total)];

N_top_est=ceil((wall_y*2-1)+eps)*ceil((wall_x*2-1)+eps)*4;
N_top_est=N_top_est*2;
% N_top_est=N_total;

%%

% clc
Rc_old=[Rc_initial_all([1 2],1);0.5];
tole_degree_loop=zeros(1,size(Rc_initial_all,2));
for nn=2:N_total% loop, to drop spheres one after another
    Rc_initial=Rc_initial_all(:,nn);% initial position to drop
    idx_touch_old=[];% index of touching sphere
    tole_degree_tmp=[];% to record tolerance
    LOGIC_stable=false;% to indicate whether sphere nn is stable or not
    %     Rc_trace=Rc_initial;
    Rc_old_tmp=Rc_old(:,max(1,size(Rc_old,2)-N_top_est):end);
    
    while LOGIC_stable~=true% loop to find a stable point for sphere nn
        %tole_degree=1;% degree of tolerance due to finite accuracy of computer
        tole_degree=ceil(nn/10)*2;% use larger tole_degree to accelerate computing
        LOGIC_step_end=false;% to indicate whether a current step end correctly
        while LOGIC_step_end==false% loop until end correctly, if not, relax tolerance
            [Rc_final,idx_touch_new,LOGIC_stable] = RD_mono_3d_wall_ ...
                (Rc_initial,idx_touch_old,Rc_old_tmp,[wall_x wall_y],tole_degree);
            %             dist_nn=pdist2(Rc_final',Rc_old_tmp');% distance between sphere in and others
            dist_fn=sqrt((Rc_old_tmp(1,:)-Rc_final(1)).^2+(Rc_old_tmp(2,:)-Rc_final(2)).^2+(Rc_old_tmp(3,:)-Rc_final(3)).^2);% faster
            dist_wn=max(max(abs(Rc_old(1,:)))-(wall_x-0.5),max(abs(Rc_old(2,:)))-(wall_y-0.5));

            if (1-min(dist_fn))/eps>tole_degree||dist_wn/eps>tole_degree% if this step cause overlap; finite accuracy may cause overlap
                % tole_degree=tole_degree+1;% relax tolerance by one eps
                tole_degree=tole_degree*2;% relax tolerance by one eps
            else% if no overlap
                LOGIC_step_end=true;% for most time, the current loop end here
            end
            
            if mean(abs(Rc_initial-Rc_final)/eps)<=tole_degree&&...
                    isequal(sort(idx_touch_new),sort(idx_touch_old))&&...
                    LOGIC_stable~=1% if get stuck in a dead loop (output = input) due to finite accuracy,
                LOGIC_step_end=false;% continue the loop
                % tole_degree=tole_degree+1;% relax tolarance by one eps
                tole_degree=tole_degree*2;% relax tolerance by one eps
            end

        end
        
        tole_degree_tmp=[tole_degree_tmp tole_degree];% when this step end, record the final tolerance
        Rc_initial=Rc_final;% set nex initial position to be the final position in this step
        idx_touch_old=idx_touch_new;% set index of new touch spheres
    end
    
    tole_degree_loop(nn)=max(tole_degree_tmp);% when stable position found, record the maximum tolerance
    Rc_old=[Rc_old Rc_final];% add sphere nn to the packing
end

end
