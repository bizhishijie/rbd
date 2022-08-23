%%
%find the transient length
tl30t=[];tl42t=[];tl56t=[];tn30t=[];tn42t=[];tn56t=[];
[~,idx_s]=sort(Rc_old(3,:));% NECESSARY to sort the particle for calculate period
Rc_tmp=Rc_old(:,idx_s);
p_posible=[30 42 56];
delta_thres=0.05;
p_minlenthre=5;%the min length
height_period_est=(2*wall_y-1)*sqrt(3)*1.5;% estimated period length
idx_pattern_old=find(Rc_tmp(3,:)>=0.5&Rc_tmp(3,:)<0.5+height_period_est);% get a structure to compare with others
height_range=[0.5 0.5+3*height_period_est];
pstart=0;
while pstart==0&&max(idx_pattern_old)<1300
    Rc_pattern=Rc_tmp(:,idx_pattern_old);% get a structure to compare with others
    idx_block=find(Rc_tmp(3,:)>=height_range(1)&Rc_tmp(3,:)<height_range(2));% get a block of structure to be compared
    Rc_block=Rc_tmp(:,idx_block);
    length_pattern=max(Rc_pattern(3,:))-min(Rc_pattern(3,:));
    length_block=max(Rc_block(3,:))-min(Rc_block(3,:));
    
    dist_loop=p_minlenthre:0.05:length_block-length_pattern;% move the pattern from bottom to top of the block
    delta_loop=zeros(size(dist_loop));% average distance difference between spheres
    delta_max_loop=zeros(size(dist_loop));% max distance difference between spheres
    Rc_pattern_tmp=Rc_pattern;
    for ii=1:length(dist_loop)% loop for different moving distances
        Rc_pattern_tmp(3,:)=Rc_pattern(3,:)+dist_loop(ii);% move the pattern in z direction
        dz_tmp=repmat(Rc_pattern_tmp(3,:)',1,size(Rc_block,2))-repmat(Rc_block(3,:),size(Rc_pattern,2),1);
        dy_tmp=repmat(Rc_pattern_tmp(2,:)',1,size(Rc_block,2))-repmat(Rc_block(2,:),size(Rc_pattern,2),1);
        dist_yz_tmp=sqrt(dz_tmp.^2+dy_tmp.^2);
        delta_loop(ii)=sqrt(mean(min(dist_yz_tmp,[],2).^2));% get the shortest distance, and their average
        delta_max_loop(ii)=max(min(dist_yz_tmp,[],2));% % get the shortest distance, and their max
    end
    delta_loop(delta_max_loop>delta_thres)=1;% remove similar but mis-matched peaks
    %         delta_loop(dist_loop<p_minlenthre)=1;
    
    [height_peak,idx_peak]=findpeaks(-delta_loop);% find peak
    idx_peak=idx_peak(height_peak>-delta_thres);% find peak higher than threshold
    height_peak=height_peak(idx_peak~=0);% remove first (unnecessary)
    idx_peak=idx_peak(idx_peak~=0);
    
    if ~isempty(idx_peak)% if find peak, means repeated structures may found
        idx_peak=min(idx_peak);
        dist_coarse=dist_loop(idx_peak);% get the moving distance of matching
        
        [dist_fine,delta_min]=fminsearch(@(dz)Mutual_Distance_pattern_dydz_(dz,Rc_pattern,Rc_block),dist_coarse);% find the closest match
        Rc_pattern_tmp(3,:)=Rc_pattern(3,:)+dist_fine;
        dz_tmp=repmat(Rc_pattern_tmp(3,:)',1,size(Rc_block,2))-repmat(Rc_block(3,:),size(Rc_pattern,2),1);
        dy_tmp=repmat(Rc_pattern_tmp(2,:)',1,size(Rc_block,2))-repmat(Rc_block(2,:),size(Rc_pattern,2),1);
        dist_yz_tmp=sqrt(dz_tmp.^2+dy_tmp.^2);
        [~,idx_min]=min(dist_yz_tmp,[],2);
        gap_tmp=mode(idx_block(idx_min)-idx_pattern_old);% calculate different between sphere index
        judge_tmp=p_posible-gap_tmp;
        idx_correct=abs(judge_tmp)<=2;
        if sum(idx_correct)==0
            gap_tmp=gap_tmp/2;%double
            judge_tmp=p_posible-gap_tmp;
            idx_correct=abs(judge_tmp)<=2;
        end
        if sum(idx_correct)~=0
            pstart=1;
            p_tmp=p_posible(idx_correct);
            transient_ball_ordinal=min(idx_pattern_old);
            transient_length=Rc_tmp(3,transient_ball_ordinal);
            if p_tmp==30
                tl30t=[tl30t transient_length];;tn30t=[tn30t transient_ball_ordinal];
            elseif p_tmp==42
                tl42t=[tl42t transient_length];;tn42t=[tn42t transient_ball_ordinal];
            elseif p_tmp==56
                tl56t=[tl56t transient_length];;tn56t=[tn56t transient_ball_ordinal];
            end
        end
    end
    idx_pattern_old=idx_pattern_old+1;
end
