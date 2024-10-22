function [Fd] = calc_Drag(v, bolinha, transdutor)
nu = transdutor.nu; rho = transdutor.rho0; R = bolinha.R;
mod_v = sqrt(v'*v);
Re = rho*mod_v*2*R/nu+eps;
Cd = 24/Re*sqrt(1+3/16*Re) ;
Fd = -1/2*Cd*(pi/4)*(2*R)^2*rho*mod_v*v*0.85;