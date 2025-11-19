function [lon, lat] = v_projection(x, y, z)
    R = 6378.137;
    lat = asind(z / R);
    lon = atan2d(x, y);  %
end
