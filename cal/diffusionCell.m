clear
wx = 1.1; wy = 5; N = 100000;
n = 6; %一个block有42个颗粒
block = n * (n + 1);
DELTA = zeros(100, floor(N / block) - 1);
blockStart = 10; %起始位置
threNear = 10; %找前后多少个球
tw = 3;
threSide = 0.01; %;临近小球的阈值
threWall = 0.00; %多靠近墙认为是靠近墙

for ii = 1:100
    clear Rc tole_degree_loop
    fileName = sprintf('../disappear/%.1f_%d_%d_%d.mat', wx, wy, N, ii);
    disp(fileName);
    load(fileName);
    delta = zeros(1, floor(N / block) - 1);
    Rc = sortByZigzagCell(Rc, threNear, threSide, threWall);

    for jj = blockStart:floor(N / block) - tw - 1
        P = Rc{jj}; %当前的block
        Q = Rc{jj+1}; %靠上的block
        z0 = (sum(P(3, :)) - sum(Q(3, :))) / block;
        P(3, :) = P(3, :) - z0;
        d = P - Q;
        delta(jj) = 1 / block * sum(sum(d.^2));
    end

    DELTA(ii, :) = delta;
end
