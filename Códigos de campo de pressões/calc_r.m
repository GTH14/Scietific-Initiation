function [r] = calc_r(x0,y0,z0,x,y,z)
    r=sqrt((x-x0).^2 + (y-y0).^2 + (z-z0).^2); 
end