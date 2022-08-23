function sortByZigzag(Rc,threWall,threNear)
threWall=0.1;
threNear=8;
sizeRc=size(Rc,2);

Rc=sortrows(Rc',3)';%按照z坐标预排序
yMin=min(Rc(2,:));%所有的小球最靠近左侧墙壁的
yMax=max(Rc(2,:));%右侧的

idMin=find(Rc(2,:)<yMin+threWall);%靠近左侧墙壁的小球的坐标
idMax=find(Rc(2,:)>yMax-threWall);%靠近右侧的

for ii=1:length(idMin)
    RcTemp=Rc(max(idMin(ii)-threNear,1)...
        ,min(idMin+threNear,sizeRc));%减少筛选范围
    [k,d] = dsearchn(RcTemp',Rc(:,idMin(ii))');
end
end