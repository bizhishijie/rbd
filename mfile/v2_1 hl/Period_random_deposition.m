%the last version by haol 20200213
%is gonna be replaced by a function
%just for the wy==5

p_posible=[30 42 56];
p_minlenthre=5;%the min length
[~,idx_s]=sort(Rc_old(3,:));% NECESSARY to sort the particle for calculate period
Rc=Rc_old(:,idx_s);
% Rc(3,:)=max(Rc(3,:))-Rc(3,:);% turn the packing up-side-down to make algorithm easy to understand
clear Rc_old

delta_thres=0.05;% threshold to determine period structure, in unit of sphere diameter
multiple=3;

%%
x_save=[];dx_save=[];
de42t=[];dz42t=[];dy42t=[];de30t=[];dz30t=[];dy30t=[];de56t=[];dz56t=[];dy56t=[];
pll30t=[];pll42t=[];pll56t=[];pln30t=[];pln42t=[];pln56t=[];
pl30mt=[];pl30at=[];pl42mt=[];pl42at=[];pl56mt=[];pl56at=[];

y_range=-wall_y:1:wall_y;
y_int=floor(Rc(2,:));
y_loop=min(y_int):max(y_int);
height_y=zeros(size(y_loop));
for iy=1:length(y_loop)
    height_y(iy)=max(Rc(3,y_int==y_loop(iy)));
end
height_findend=min(height_y);% start from this height; above which the packing has not fully been formed

height_start=0.5;
height_period_est=(2*wall_y-1)*sqrt(3)*1.5;% estimated period length
idx_pattern_old=find(Rc(3,:)>=height_start&Rc(3,:)<height_start+height_period_est);% get a structure to compare with others
height_range=[height_start height_start+3*height_period_est];

%%
LOGIC_period=false;
tl30t=[];tl42t=[];tl56t=[];tn30t=[];tn42t=[];tn56t=[];
period_num=0;
idx_period=zeros(2,0);
delta_all=[];dy_all=[];dz_all=[];dx_all=[];x_all=[];pl_all=[];

while height_range(2)<height_findend-12% loop for all the packing (a -12 for avoiding unexpected break)
    if isempty(delta_all)%no period before
        Rc_pattern=Rc(:,idx_pattern_old);% get a structure to compare with others
        idx_block=find(Rc(3,:)>=height_range(1)&Rc(3,:)<height_range(2));% get a block of structure to be compared
        Rc_block=Rc(:,idx_block);
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
            if sum(idx_correct)==0%%wrongly found
                idx_pattern_old=idx_pattern_old+1;
            else
                period_num=period_num+1;% a new period strucutre is found
                idx_period=[idx_period [min(idx_pattern_old) 0]'];% record the starting index
                idx_start=min(idx_block(idx_min));
                p_tmp=p_posible(idx_correct);
                idx_pattern_old=idx_start:(idx_start+p_tmp-1);
                height_old=dist_fine;
                delta_all=[delta_all delta_min];%core too stand for a structure is founded
                if period_num==1%the first period
                    transient_ball_ordinal=idx_period(1,1);
                    transient_length=Rc(3,transient_ball_ordinal);
                    if p_tmp==30
                        tl30t=[tl30t transient_length];;tn30t=[tn30t transient_ball_ordinal];
                    elseif p_tmp==42
                        tl42t=[tl42t transient_length];;tn42t=[tn42t transient_ball_ordinal];
                    elseif p_tmp==56
                        tl56t=[tl56t transient_length];;tn56t=[tn56t transient_ball_ordinal];
                    end
                end
            end
        else% if no repeated structures found
            idx_pattern_old=idx_pattern_old+1;% define new pattern downwards
        end
        height_start=min(Rc(3,idx_pattern_old));% define new start height of block
        height_range=[height_start height_start+3*height_period_est];
        
    else%do has period before
        Rc_pattern=Rc(:,idx_pattern_old);% get a structure to compare with others
        length_pattern=max(Rc_pattern(3,:))-min(Rc_pattern(3,:));
        idx_block=find(Rc(3,:)>=max(Rc_pattern(3,:))-1&Rc(3,:)<max(Rc_pattern(3,:))+height_old+1);% get a block of structure to be compared
        Rc_block=Rc(:,idx_block);
        length_block=max(Rc_block(3,:))-min(Rc_block(3,:));
        
        delz=1/2*(length_block-length_pattern);
        dist_loop=[0:0.05:length_block-length_pattern]+height_old-delz;% move the pattern from bottom to top of the block
        delta_loop=zeros(size(dist_loop));% average distance difference between spheres
        Rc_pattern_tmp=Rc_pattern;
        for ii=1:length(dist_loop)% loop for different moving distances
            Rc_pattern_tmp(3,:)=Rc_pattern(3,:)+dist_loop(ii);% move the pattern in z direction
            dz_tmp=repmat(Rc_pattern_tmp(3,:)',1,size(Rc_block,2))-repmat(Rc_block(3,:),size(Rc_pattern,2),1);
            dy_tmp=repmat(Rc_pattern_tmp(2,:)',1,size(Rc_block,2))-repmat(Rc_block(2,:),size(Rc_pattern,2),1);
            dist_yz_tmp=sqrt(dz_tmp.^2+dy_tmp.^2);
            delta_loop(ii)=sqrt(mean(min(dist_yz_tmp,[],2).^2));% get the shortest distance, and their average
        end
        
        [~,idx_peak]=min(delta_loop);% find min
        idx_peak=min(idx_peak);
        dist_coarse=dist_loop(idx_peak);% get the moving distance of matching
        %         gap_tmp=mode(idx_block(idx_min)-idx_pattern_old);% calculate different between sphere index
        [dist_fine,delta_min]=fminsearch(@(dz)Mutual_Distance_pattern_dydz_(dz,Rc_pattern,Rc_block),dist_coarse);% find the closest match
        Rc_pattern_tmp(3,:)=Rc_pattern(3,:)+dist_fine;
        dz_tmp=repmat(Rc_pattern_tmp(3,:)',1,size(Rc_block,2))-repmat(Rc_block(3,:),size(Rc_pattern,2),1);
        dy_tmp=repmat(Rc_pattern_tmp(2,:)',1,size(Rc_block,2))-repmat(Rc_block(2,:),size(Rc_pattern,2),1);
        dist_yz_tmp=sqrt(dz_tmp.^2+dy_tmp.^2);
        [~,idx_min]=min(dist_yz_tmp,[],2);
        gap_tmp=mode(idx_block(idx_min)-idx_pattern_old);% calculate different between sphere index
        
        if delta_min<=mean(delta_all)*multiple&&abs(gap_tmp-p_tmp)<=2%core! period will continue
            period_end=false;
            %             pend_thre=mean(delta_all)*multiple;
        else %period will end
            period_end=true;
        end
        if period_end
            %find the exact endding point
            idx_max_tmp=min(idx_pattern_old);
            pend=0;idx_pattern_old=idx_pattern_old-p_tmp;%close to the top
            while pend==0
                idx_pattern_old=idx_pattern_old+1;
                Rc_pattern=Rc(:,idx_pattern_old);% get a structure to compare with others
                length_pattern=max(Rc_pattern(3,:))-min(Rc_pattern(3,:));
                idx_block=find(Rc(3,:)>=max(Rc_pattern(3,:))-1&Rc(3,:)<max(Rc_pattern(3,:))+height_old+1);% get a block of structure to be compared
                Rc_block=Rc(:,idx_block);
                length_block=max(Rc_block(3,:))-min(Rc_block(3,:));
                
                delz=1/2*(length_block-length_pattern);
                dist_loop=[0:0.05:length_block-length_pattern]+height_old-delz;% move the pattern from bottom to top of the block
                delta_loop=zeros(size(dist_loop));% average distance difference between spheres
                Rc_pattern_tmp=Rc_pattern;
                for ii=1:length(dist_loop)% loop for different moving distances
                    Rc_pattern_tmp(3,:)=Rc_pattern(3,:)+dist_loop(ii);% move the pattern in z direction
                    dz_tmp=repmat(Rc_pattern_tmp(3,:)',1,size(Rc_block,2))-repmat(Rc_block(3,:),size(Rc_pattern,2),1);
                    dy_tmp=repmat(Rc_pattern_tmp(2,:)',1,size(Rc_block,2))-repmat(Rc_block(2,:),size(Rc_pattern,2),1);
                    dist_yz_tmp=sqrt(dz_tmp.^2+dy_tmp.^2);
                    delta_loop(ii)=sqrt(mean(min(dist_yz_tmp,[],2).^2));% get the shortest distance, and their average
                end
                [~,idx_peak]=min(delta_loop);% find min
                idx_peak=min(idx_peak);
                dist_coarse=dist_loop(idx_peak);% get the moving distance of matching
                [dist_fine,delta_min]=fminsearch(@(dz)Mutual_Distance_pattern_dydz_(dz,Rc_pattern,Rc_block),dist_coarse);% find the closest match
                Rc_pattern_tmp(3,:)=Rc_pattern(3,:)+dist_fine;
                dz_tmp=repmat(Rc_pattern_tmp(3,:)',1,size(Rc_block,2))-repmat(Rc_block(3,:),size(Rc_pattern,2),1);
                dy_tmp=repmat(Rc_pattern_tmp(2,:)',1,size(Rc_block,2))-repmat(Rc_block(2,:),size(Rc_pattern,2),1);
                dist_yz_tmp=sqrt(dz_tmp.^2+dy_tmp.^2);
                [~,idx_min]=min(dist_yz_tmp,[],2);
                gap_tmp=mode(idx_block(idx_min)-idx_pattern_old);% calculate different between sphere index
                if delta_min<=mean(delta_all)*multiple&&abs(gap_tmp-p_tmp)<=2% find min
                    pend=0;
                else
                    pend=1;
                end
                if min(idx_pattern_old)==idx_max_tmp
                    pend=1;
                end
            end
            
            idx_period(2,end)=max(idx_pattern_old)+p_tmp-1;
            if p_tmp==30
                de30t=[de30t delta_all];
                x_save=[x_save x_all];
                dx_save=[dx_save dx_all];
                dy30t=[dy30t dy_all];
                dz30t=[dz30t dz_all];
                pl30mt=[pl30mt mean(pl_all)];
                pl30at=[pl30at pl_all];
                pll30t=[pll30t Rc(3,idx_period(2,end))-Rc(3,idx_period(1,end))];
                pln30t=[pln30t idx_period(2,end)-idx_period(1,end)];
            elseif p_tmp==42
                de42t=[de42t delta_all];
                x_save=[x_save x_all];
                dx_save=[dx_save dx_all];
                dy42t=[dy42t dy_all];
                dz42t=[dz42t dz_all];
                pl42mt=[pl42mt mean(pl_all)];
                pl42at=[pl42at pl_all];
                pll42t=[pll42t Rc(3,idx_period(2,end))-Rc(3,idx_period(1,end))];
                pln42t=[pln42t idx_period(2,end)-idx_period(1,end)];
            elseif p_tmp==56
                de56t=[de56t delta_all];
                x_save=[x_save x_all];
                dx_save=[dx_save dx_all];
                dy56t=[dy56t dy_all];
                dz56t=[dz56t dz_all];
                pl56mt=[pl56mt mean(pl_all)];
                pl56at=[pl56at pl_all];
                pll56t=[pll56t Rc(3,idx_period(2,end))-Rc(3,idx_period(1,end))];
                pln56t=[pln56t idx_period(2,end)-idx_period(1,end)];
            end
            
            %             delta_save=[delta_save delta_all];
            delta_all=[];dy_all=[];dz_all=[];dx_all=[];x_all=[];pl_all=[];
            height_pend=Rc(3,max(idx_pattern_old)+p_tmp);
            idx_pattern_old=find(Rc(3,:)>=height_pend&Rc(3,:)<height_pend+height_period_est);
            
        else
            dz_tmp=repmat(Rc_pattern_tmp(3,:)',1,size(Rc_block,2))-repmat(Rc_block(3,:),size(Rc_pattern,2),1);
            dy_tmp=repmat(Rc_pattern_tmp(2,:)',1,size(Rc_block,2))-repmat(Rc_block(2,:),size(Rc_pattern,2),1);
            dist_yz_tmp=sqrt(dz_tmp.^2+dy_tmp.^2);
            [~,idx_min]=min(dist_yz_tmp,[],2);
            idx_pattern_new=idx_block(idx_min);% use the matched structure as next pattern
            dz_all=[dz_all Rc_pattern_tmp(3,:)-Rc_block(3,idx_min)];
            dy_all=[dy_all Rc_pattern_tmp(2,:)-Rc_block(2,idx_min)];
            dx_all=[dx_all Rc_pattern_tmp(1,:)-Rc_block(1,idx_min)];
            x_all=[x_all Rc_block(1,idx_min)];pl_all=[pl_all dist_fine];
            
            delta_all=[delta_all delta_min];
%             gap_tmp=mode(idx_pattern_new-idx_pattern_old);% calculate different between sphere index
            idx_pattern_old=idx_pattern_new;
            height_old=dist_fine;
        end
        
        height_start=min(Rc(3,idx_pattern_old));% define new start height of block
        height_range=[height_start height_start+3*height_period_est];
    end
end

%%
% p_num=1:max(period_loop)-1;% index of different period structure; remove the last structure, calculation not be finished
% P=zeros(1,length(p_num));% period
% for ii=1:length(p_num)
%     P(ii)=mode(gap_loop(period_loop==ii));
% end
% N=idx_period(2,1:end-1)-idx_period(1,1:end-1);% number of particles in one period
% L=Rc(3,idx_period(2,1:end-1))-Rc(3,idx_period(1,1:end-1));% length of one period

% if idx_period(2,end)==0
%     idx_period(:,end)=[];
% end


% if idx_period(2,end)==0
%     idx_period(2,end)=size(Rc,2);
% end
% % plot different periods
% figure(2);clf;hold on
% for ii=1:size(idx_period,2)
%     idx_tmp=idx_period(:,ii);
%     plot(Rc(2,idx_tmp(1):idx_tmp(2)),Rc(3,idx_tmp(1):idx_tmp(2)),'.','Color','b')
%     if ii>1
%         plot(Rc(2,idx_tmpl:idx_tmp(1)),Rc(3,idx_tmpl:idx_tmp(1)),'.','Color','r')
%     end
%     idx_tmpl=idx_period(2,ii);
% end
% 
% Rc(3,idx_period).*10