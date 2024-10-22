for theta = 0:0.1:6*pi
    figure(4)
    surf(x*1e3,z*1e3,real(P1*exp(j*theta)));
    colormap(jet);
    shading interp;
    view(2);
    a = colorbar;
    caxis([-1300 1300])
    ylabel(a,'Pressão acústica [Pa]');
    xlabel('Eixo x [mm]');
    xlim([-50e-3,50e-3]*1e3);
    ylim([-71e-3, 71e-3]*1e3);
    ylabel('Eixo z [mm]');
    pause(.1);
end