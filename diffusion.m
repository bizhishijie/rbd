clear
wx=1.1;wy=5;N=100000;
n=6;%一个block有42个颗粒
block=n*(n+1);
DELTA=zeros(100,floor(N/block)-1);
tw=100;%起始位置
for ii =1:100
    clear Rc tole_degree_loop
    fileName=sprintf('disappear/%.1f_%d_%d_%d.mat',wx,wy,N,ii);
    disp(fileName);
    load(fileName);

    delta=zeros(1,floor(N/block)-1);
    Rc=sortRc(Rc);%对Rc 按照zig zag编号
    PQ=Rc(:,(tw-1)*block+1:tw*block);%最靠下侧的
    for jj=tw+1:floor(N/block)-1
        P=Rc(:,jj*block+1:(jj+1)*block);%靠上的任意一个
        z0=(sum(P(3,:))-sum(PQ(3,:)))/block;
        P(3,:)=P(3,:)-z0;
        [k,d] = dsearchn(P',PQ'); %k是数据点的索引，d是最近点的距离
        if any(diff(diff(sort(k))))
            break
        end
        %         d=sum(P-PQ);
        delta(jj)=1/block*sum(d.^2);
    end
    DELTA(ii,:)=delta;
end
