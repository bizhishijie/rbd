function showRc(Rc, num1, num2, isBall)
    r = 0.5;

    if nargin == 1
        plot3(Rc(1, :), Rc(2, :), Rc(3, :), 'r.');

        for ii = 1:length(Rc)
            text(Rc(1, ii), Rc(2, ii), Rc(3, ii), sprintf('%d', ii));
        end

    elseif nargin == 2 && num1
        plot3(Rc(1, :), Rc(2, :), Rc(3, :), 'r.');

        for ii = 1:length(Rc)
            hold on
            [x, y, z] = sphere(30);
            X = x * r + Rc(1, ii);
            Y = y * r + Rc(2, ii);
            Z = z * r + Rc(3, ii);
            mesh(X, Y, Z);
        end

    elseif nargin == 3
        plot3(Rc(1, num1:num2), Rc(2, num1:num2), Rc(3, num1:num2), 'r.');

        for ii = num1:num2
            text(Rc(1, ii), Rc(2, ii), Rc(3, ii), sprintf('%d', ii));
        end

    elseif nargin == 4 && isBall
        plot3(Rc(1, num1:num2), Rc(2, num1:num2), Rc(3, num1:num2), 'r.');

        for ii = num1:num2
            hold on
            [x, y, z] = sphere(30);
            X = x * r + Rc(1, ii);
            Y = y * r + Rc(2, ii);
            Z = z * r + Rc(3, ii);
            mesh(X, Y, Z);
        end

    end

    axis equal
    grid on
end
