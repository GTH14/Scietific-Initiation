function [x,y,z,vx,vy,vz] = f_Simulacao_dinamica_v3(transdutor,bolinha)
    n = length(bolinha.t_v);
    h = bolinha.h;
    x = zeros(1,n)+zeros(1,n)+j*zeros(1,n); y = zeros(1,n)+j*zeros(1,n); z = zeros(1,n)+j*zeros(1,n);

    vx = zeros(1,n)+j*zeros(1,n); vy = zeros(1,n)+j*zeros(1,n); vz = zeros(1,n)+j*zeros(1,n);

    i = 1;

    x(i) = bolinha.pos_i(1); y(i) = bolinha.pos_i(2); z(i) = bolinha.pos_i(3);

    vx(i) = bolinha.v_i(1); vy(i) = bolinha.v_i(2); vz(i) = bolinha.v_i(3);

    disc = transdutor.disc;
    m = bolinha.rho*4/3*pi*bolinha.R^3;
    % Teste inicial com integração de Euler Explicito
    while i <= n-2
        x = bolinha.x(i); y = bolinha.y(i); z = bolinha.z(i);
        vx = bolinha.vx(i); vy = bolinha.vy(i); vz = bolinha.vz(i);
        Fd = calc_Drag(v, bolinha, transdutor);
        Fx = -(Calc_Gorkov(x+disc/2, y, z, transdutor) - Calc_Gorkov(x-disc/2, y, z, transdutor))/disc;
        Fy = -(Calc_Gorkov(x, y+disc/2, z, transdutor) - Calc_Gorkov(x, y-disc/2, z, transdutor))/disc;
        Fz = -(Calc_Gorkov(x, y, z+disc/2, transdutor) - Calc_Gorkov(x, y, z-disc/2, transdutor))/disc;
        
        vxi = vx + Fx/m*h; vyi = vy + Fy/m*h; vzi = vz + Fz/m*h;

        xi = x + vxi*h; yi = y + vyi*h; zi = z + vzi*h;

        x(i+1) = xi; y(i+1) = yi; z(i+1) = zi;

        vx(i+1) = vxi; vy(i+1) = vyi; vz(i+1) = vzi;
        i = i+1;
    end

end