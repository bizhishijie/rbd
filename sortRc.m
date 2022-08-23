function Rc=sortRc(Rc)
Rc=sortrows(Rc',3)';
RcD=diff(Rc(3,:));%区分不同层
[~,locs] =findpeaks(RcD,'MinPeakDistance',2);
locsD=diff(locs);

for ii=length(locsD):1
    if locsD(ii)==6%ii+1 - ii
        locs=[locs(1:ii) locs(ii)+3 locs(ii+1:end)];
    elseif locsD(ii)==7 && locsD(ii-1)==3 %ii+1 - ii
        locs=[locs(1:ii) locs(ii)+4 locs(ii+1:end)];
    elseif locsD(ii)==7 && locsD(ii-1)==4 %ii+1 - ii
        locs=[locs(1:ii) locs(ii)+3 locs(ii+1:end)];
    elseif locsD(ii)==8
        locs=[locs(1:ii) locs(ii)+4 locs(ii+1:end)];
    end
end

locs=[0,locs]+1;
for  ii=1:length(locs)-1
    Rc(:,locs(ii):locs(ii+1)-1)=sortrows(Rc(:,locs(ii):locs(ii+1)-1)',2)';%按照某一方向排序
end
end