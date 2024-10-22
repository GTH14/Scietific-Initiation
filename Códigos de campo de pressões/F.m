function [pdot] = F(t,p)
    pdot = [0;0;0;0;0;0];
    x = p(1); y = p(2); z = p(3);
    v = p(4:6);
    Fd = calc_Drag(v, bolinha, transdutor);
    Fgx = -(Calc_Gorkov(x+disc/2, y, z, transdutor) - Calc_Gorkov(x-disc/2, y, z, transdutor))/disc;
    Fgy = -(Calc_Gorkov(x, y+disc/2, z, transdutor) - Calc_Gorkov(x, y-disc/2, z, transdutor))/disc;
    Fgz = -(Calc_Gorkov(x, y, z+disc/2, transdutor) - Calc_Gorkov(x, y, z-disc/2, transdutor))/disc;
    pdot(1) = p(4); pdot(2) = p(5); pdot(3) = p(6);
    pdot(4) = (Fgx+Fd(1))/m; pdot(5) = (Fgy+Fd(2))/m; pdot(6) = (Fgz+Fd(3))/m;
end