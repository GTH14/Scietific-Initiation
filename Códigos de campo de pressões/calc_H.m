function [H] = calc_H(v)
    H = 2*besselj(1,v)./v;
end