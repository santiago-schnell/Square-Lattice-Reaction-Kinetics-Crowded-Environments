function [x,y] = boundarycheck(x,y,N)

if x > N, x = 1;
else if x < 1, x = N;
    end
end

if y > N, y = 1;
else if y < 1, y = N;
    end
end
