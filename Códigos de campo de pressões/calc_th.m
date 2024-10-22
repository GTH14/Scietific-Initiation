function [theta] = calc_th(x0, y0, z0, x, y, z, normal, r)
    theta = acos((normal(1).*(x-x0) + normal(2).*(y-y0) + normal(3).*(z-z0))./r);
end