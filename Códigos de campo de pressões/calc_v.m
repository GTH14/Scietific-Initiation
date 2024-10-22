function [v] = calc_v(k, a, theta)
    v = k*a*sin(theta)+1e-12;
end