function [U] = Calc_Gorkov_Simulacao_Dinamica(x, y, z, transdutor)
    R = transdutor.R;
    rho0 = transdutor.rho0;
    c = transdutor.c;
    p_total = 0+0*j;
    
    p_total = Calc_pressao_mex(x, y, z, transdutor);
    [vx, vy, vz] = Calc_vel1(x, y, z, transdutor);
    m_p = 1/2*abs(p_total).^2;
    m_v = 1/2*(vx.*conj(vx) + vy.*conj(vy) + vz.*conj(vz));
    U = 2*pi*R^3*(m_p/(3*rho0*c^2) - rho0/2*(m_v));
end   