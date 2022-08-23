clear
wx=1.1;wy=5;N=100000;
n=6;%一个block有42个颗粒
block=n*(n+1);
DELTA=zeros(100);
tw=1;%

Rc_cell=cell(1,100);
delta=zeros(1,floor(N/block)-1);
for ii=1:100
    clear Rc tole_degree_loop
    fileName=sprintf('disappear/%.1f_%d_%d_%d.mat',wx,wy,N,ii);
    disp(fileName);
    load(fileName);
    Rc=sortRc(Rc);
    Rc_cell{ii}=Rc;
end
%之后需要拆分block
for ii=1:100
    for jj=ii:100
        delta=mean(sum(sum((Rc_cell{ii}-Rc_cell{jj}).^2)))/block;
        DELTA(ii,jj)=delta;
        DELTA(jj,ii)=delta;
    end
end
