function idList = RcNear(Rc, id0, threNear, threSide)
    Rc0 = Rc(:, id0); %待寻找的小球的位置
    sizeRc = size(Rc, 2); %一共有的小球数目
    RcTemp = Rc(:, max(id0 - threNear, 1):min(id0 + threNear, sizeRc)); %构建待寻找的小球附近的列表
    Rc1 = repmat(Rc0, 1, size(RcTemp, 2)); %格式化数据
    id = find(abs(sum((RcTemp - Rc1).^2) - 1) < threSide); %寻找距离为1的小球
    idList = -find(sum(RcTemp == Rc(:, id0)) == 3) + id0 + id; %还原
end
