clear
load('1.1_5_100000_4.mat')

[~,idx_s]=sort(Rc(3,:));
Rc=Rc(:,idx_s);

%%
Nmax=2*ceil(wy);
C_IDX=cell(1,size(Rc,2));
for ii=1:size(Rc,2)
    idx_range=[max(1,ii-Nmax):1:min(size(Rc,2),ii+Nmax)];
    dis_tmp=sqrt(sum((Rc(:,idx_range)-Rc(:,ii)).^2,1));
    C_IDX{ii}=idx_range(dis_tmp<=1+tole_degree_loop(ii)*eps&idx_range~=ii);
end

idx_L=find(Rc(2,:)==-((wy-1)/2));
idx_R=find(Rc(2,:)==((wy-1)/2));
zigzag_IDX=cell(1,length(idx_L));
zigzag_loop=zeros(1,length(idx_L));
for ii=1:length(idx_L)
    idx_start=idx_L(ii);
    idx_zigzag=idx_start;
    up_or_down=1;
    zigzag_fail=0;
    while 1==1
        c_idx=C_IDX{idx_start};
        delta_y=Rc(2,c_idx)-Rc(2,idx_start);
        delta_z=Rc(3,c_idx)-Rc(3,idx_start);
        c_idx=c_idx(delta_y>0);
        delta_z=delta_z(delta_y>0);
        delta_z=delta_z*up_or_down;
        idx_connect=c_idx(delta_z==max(delta_z)&delta_z>0);
        if isempty(idx_connect)
            zigzag_fail=1;
            break
        end
        up_or_down=-1*up_or_down;
        idx_zigzag=[idx_zigzag idx_connect];
        idx_start=idx_connect;
        if ismember(idx_connect,idx_R)
            break
        end
    end
    %     clc
    %     zigzag_fail
    %     plot3(Rc(1,idx_zigzag),Rc(2,idx_zigzag),Rc(3,idx_zigzag),'o-')
    %     ii=ii+1
    
    if zigzag_fail==0
        zigzag_IDX{ii}=idx_zigzag;
        zigzag_loop(ii)=1;
    end
end

%%
N_zigzag=cellfun('length',zigzag_IDX);
idx_successive=N_zigzag(2:end-1)==N_zigzag(1:end-2)&...
    N_zigzag(2:end-1)==N_zigzag(3:end)&...
    zigzag_loop(2:end-1)~=0;
idx_successive=[0 idx_successive 0];
L_idx_successive=bwlabel(idx_successive);

IDX_period=cell(1,max(L_idx_successive));
for ii=1:max(L_idx_successive)
    idx=find(L_idx_successive==ii);
    idx=[idx(1)-1 idx idx(end)+1];
    n=N_zigzag(idx(1));
    
    IDX_period{ii}=cell2mat(zigzag_IDX(idx)');
end


%%

dy2_all=zeros(0,0);
idx_add=1;

for ii=1:length(IDX_period)
    idx_period=IDX_period{ii};
    n=size(idx_period,2);
    
    t_loop_tmp=1:floor(size(idx_period,1)/(n-1))-1;
    
    dy2_ALL=cell(n-1,n-2);
    for jj=2:n-1
        for kk=1:6
            idx_tmp=idx_period(kk:6:end,jj);
            y_tmp=Rc(2,idx_tmp);
            d2=(y_tmp-y_tmp').^2;
            d2=d2(:);
            t=(1:length(y_tmp))-(1:length(y_tmp))';
            t=t(:);
            d2=d2(t>0);
            t=t(t>0);
            
            dy2_loop_tmp=zeros(1,length(t_loop_tmp));
            for tt=1:length(t_loop_tmp)
                dy2_loop_tmp(tt)=mean(d2(t==tt));
            end
            dy2_ALL{kk,jj-1}=dy2_loop_tmp;
        end
    end
    
    dy2_mean=mean(cell2mat(dy2_ALL(:)),1);
    %     loglog(1:size(dy2_all,2),mean(dy2_all,1),'o-')
    
    if n==7
        dy2_all(idx_add,1:length(dy2_mean))=dy2_mean;
        idx_add=idx_add+1;
    end
end

dy2_all(dy2_all==0)=nan;

figure(1);hold on
loglog(1:size(dy2_all,2),nanmean(dy2_all,1),'o-')

