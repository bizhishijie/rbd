clear
wx = 1.1; wy = 5; N = 100000; ii = 2;
fileName = sprintf('%.1f_%d_%d_%d.mat', wx, wy, N, ii);
load(fileName, 'Rc');

threNear = 10; %找前后多少个球
threSide = 0.01; %;临近小球的阈值
threWall = 0.00; %多靠近墙认为是靠近墙

Rc = sortrows(Rc', 3)'; %按照z坐标预排序
RcNew = zeros(size(Rc));
yMin = min(Rc(2, :)); %所有的小球最靠近左侧墙壁的
yMax = max(Rc(2, :)); %右侧的

idMin = sort(find(Rc(2, :) <= yMin + threWall)); %靠近左侧墙壁的小球的坐标
idMax = sort(find(Rc(2, :) >= yMax - threWall)); %靠近右侧的
cnt = 1;

for ii = 1:length(idMin) %对zigzag遍历
    sgn = 1; %默认第一个向上找小球，1代表向上，-1代表向下，初始化
    idRcTemp = idMin(ii); %使用左侧的小球初始化
    RcNew(:, cnt) = Rc(:, idRcTemp);
    cnt = cnt + 1;

    while ~ismember(idRcTemp, idMax) %只要没找到另外一端
        idNearList = RcNear(Rc, idRcTemp, threNear, threSide); %找临近的小球的编号
        RcNearList = Rc(:, idNearList); %找临近的小球的坐标
        delta = RcNearList(2:3, :) - repmat(Rc(2:3, idRcTemp), 1, size(RcNearList, 2)); %计算临近的小球与当前小球的坐标差值
        theta = atan(delta(2, :) ./ delta(1, :) .* (delta(1, :) > 0)); %利用坐标差值计算得到角度
        theta(theta == 0) = nan; %只要方向不对（向右找）就干掉
        idRcTemp = idNearList(dsearchn(theta', sgn * pi / 2)); %最接近某个角度的
        RcNew(:, cnt) = Rc(:, idRcTemp); %将找到的小球坐标加到RcNew里
        sgn = -sgn; %转换向上找还是向下找
        cnt = cnt + 1;
        fprintf('%d        %d\n', cnt, idRcTemp)
    end

end

figure(1)
showRc(Rc, 1, 40)
figure(2)
showRc(RcNew, 1, 20)
