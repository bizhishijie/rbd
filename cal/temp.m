clear
wx = 1.1; wy = 5; N = 100000;
n = 6; %一个block有42个颗粒
block = n * (n + 1);
DELTA = zeros(100, floor(N / block) - 1);
blockStart = 10; %起始位置
threNear = 10; %找前后多少个球
threSide = 0.01; %;临近小球的阈值
threWall = 0.00; %多靠近墙认为是靠近墙
ii = 1;
fileName = sprintf('../disappear/%.1f_%d_%d_%d.mat', wx, wy, N, ii);
disp(fileName);
load(fileName, "Rc");
P = Rc(:, blockStart * block + 4:(blockStart + 1) * block + 3);
Rc = sortByZigzag(Rc, threNear, threSide, threWall);

for jj = blockStart:1000
    Q = Rc(:, (jj + blockStart) * block + 4:(jj + blockStart + 1) * block + 3); %靠上的block
    z0 = (sum(Q(3, :)) - sum(P(3, :))) / block;
    Q(3, :) = Q(3, :) - z0;
    showRc(Q, true)
    title(jj)
    drawnow
    saveas(gcf, sprintf('./fig/%d.jpg', jj))
    clf
end
