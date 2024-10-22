function [p] = calc_p(Pax, H, omega, k, t, r, phase)
    p = Pax.*H.*exp(j*(omega*t + k*r - phase));
end