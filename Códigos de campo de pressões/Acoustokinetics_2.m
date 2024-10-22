function [F] = Acoustokinetics_2(x, y, z, transdutor)
    disc = transdutor.disc;
    k = transdutor.k;
    c = transdutor.c;
    R = transdutor.R;
    rho = transdutor.rho0;
    rho_e = 30;
    c_e = 900;
    
    f0 = 1 - rho*c^2/(rho_e*c_e^2);
    f1 = 2*(rho_e - rho)/(2*rho_e + rho);
    alpha_a = 4*pi*R^3/(3*rho*c^2)*f0*(-1 + 1i/3*(f0 + f1)*(k*R)^3);
    beta_a = 2*pi*R^3/(rho*c^2)*f1*(1 + 1i/6*f1*(k*R)^3);
    
    p = Calc_pressao_mex(x, y, z, transdutor);
    
    p_x_mais = Calc_pressao_mex(x+disc/2, y, z, transdutor);
    p_x_menos = Calc_pressao_mex(x-disc/2, y, z, transdutor);
    p_y_mais = Calc_pressao_mex(x, y+disc/2, z, transdutor);
    p_y_menos = Calc_pressao_mex(x, y-disc/2, z, transdutor);
    p_z_mais = Calc_pressao_mex(x, y, z+disc/2, transdutor);
    p_z_menos = Calc_pressao_mex(x, y, z-disc/2, transdutor);
    
    p_x_mais_y_mais = Calc_pressao_mex(x+disc/2, y+disc/2, z, transdutor);
    p_x_mais_y_menos = Calc_pressao_mex(x+disc/2, y-disc/2, z, transdutor);
    p_x_menos_y_mais = Calc_pressao_mex(x-disc/2, y+disc/2, z, transdutor);
    p_x_menos_y_menos = Calc_pressao_mex(x-disc/2, y-disc/2, z, transdutor);
    
    p_x_mais_z_mais = Calc_pressao_mex(x+disc/2, y, z+disc/2, transdutor);
    p_x_mais_z_menos = Calc_pressao_mex(x+disc/2, y, z-disc/2, transdutor);
    p_x_menos_z_mais = Calc_pressao_mex(x-disc/2, y, z+disc/2, transdutor);
    p_x_menos_z_menos = Calc_pressao_mex(x-disc/2, y, z-disc/2, transdutor);
    
    p_y_mais_z_mais = Calc_pressao_mex(x, y+disc/2, z+disc/2, transdutor);
    p_y_mais_z_menos = Calc_pressao_mex(x, y+disc/2, z-disc/2, transdutor);
    p_y_menos_z_mais = Calc_pressao_mex(x, y-disc/2, z+disc/2, transdutor);
    p_y_menos_z_menos = Calc_pressao_mex(x, y-disc/2, z-disc/2, transdutor);
    
    grad_p_x = (p_x_mais - p_x_menos)/disc;
    grad_p_y = (p_y_mais - p_y_menos)/disc;
    grad_p_z = (p_z_mais - p_z_menos)/disc;
    
    grad_pconj_x = conj(grad_p_x);
    grad_pconj_y = conj(grad_p_y);
    grad_pconj_z = conj(grad_p_z);
    
    grad_pconj_x2 = conj((p_x_mais - 2*p + p_x_menos)/(disc/2)^2);
    grad_pconj_y2 = conj((p_y_mais - 2*p + p_y_menos)/(disc/2)^2);
    grad_pconj_z2 = conj((p_z_mais - 2*p + p_z_menos)/(disc/2)^2);
    grad_pconj_xy = conj((p_x_mais_y_mais - p_x_mais_y_menos - p_x_menos_y_mais...
        + p_x_menos_y_menos)/(disc)^2);
    
    grad_pconj_xz = conj((p_x_mais_z_mais - p_x_mais_z_menos - p_x_menos_z_mais...
        + p_x_menos_z_menos)/(disc)^2);
    grad_pconj_yz = conj((p_y_mais_z_mais - p_y_mais_z_menos - p_y_menos_z_mais...
        + p_y_menos_z_menos)/(disc)^2);
  
    Fx1 = alpha_a*p*grad_pconj_x;
    Fy1 = alpha_a*p*grad_pconj_y;
    Fz1 = alpha_a*p*grad_pconj_z;

    Fx2 = beta_a*k^(-2)*(grad_p_x*grad_pconj_x2 + grad_p_y*grad_pconj_xy + grad_p_z*grad_pconj_xz);
    Fy2 = beta_a*k^(-2)*(grad_p_x*grad_pconj_xy + grad_p_y*grad_pconj_y2 + grad_p_z*grad_pconj_yz);
    Fz2 = beta_a*k^(-2)*(grad_p_x*grad_pconj_xz + grad_p_y*grad_pconj_yz + grad_p_z*grad_pconj_z2);

    Fx = 1/2*real(Fx1 + Fx2);
    Fy = 1/2*real(Fy1 + Fy2);
    Fz = 1/2*real(Fz1 + Fz2);

    F = [Fx; Fy; Fz];

end