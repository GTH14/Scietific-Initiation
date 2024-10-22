function [vx, vy, vz] = Calc_vel2(x, y, z, transdutor)
    disc = transdutor.disc;
    rho0 = transdutor.rho0;
    omega = transdutor.omega;
    vx = -(Calc_pressao_2_mex(x+disc/2, y, z, transdutor)-Calc_pressao_2_mex(x-disc/2, y, z, transdutor))/(j*omega*rho0*disc);
    vy = -(Calc_pressao_2_mex(x, y+disc/2, z, transdutor)-Calc_pressao_2_mex(x, y-disc/2, z, transdutor))/(j*omega*rho0*disc);
    vz = -(Calc_pressao_2_mex(x, y, z+disc/2, transdutor)-Calc_pressao_2_mex(x, y, z-disc/2, transdutor))/(j*omega*rho0*disc);
end