for i = 10:100
    str = strcat(sprintf('./fig/%d.jpg', i));
    disp(str)
    A = imread(str);
    [I, map] = rgb2ind(A, 256);

    if (i == 10)
        imwrite(I, map, 'movefig.gif', 'DelayTime', 0.1, 'LoopCount', Inf)
    else
        imwrite(I, map, 'movefig.gif', 'WriteMode', 'append', 'DelayTime', 0.1)
    end

end
