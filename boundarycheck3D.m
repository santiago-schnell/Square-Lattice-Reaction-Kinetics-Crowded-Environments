function [x,y,z] = boundarycheck3D(x,y,z,N)

if x > N, x = 1;
else if x < 1, x = N;
    end
end

if y > N, y = 1;
else if y < 1, y = N;
    end
end

if z > N, z = 1;
else if z < 1, z = N;
    end
end
