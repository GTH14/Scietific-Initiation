freq=1e3; %frequencia [Hz]
c = 340; % velocidade do som [m/s]
A=1; %Amplitude [Pa]
rho0 = 1.2922;
c = 340;
omega=2*pi*freq;
k=omega/c;
[x,y]=meshgrid(-1:.05:1,-1:.05:1);

r=sqrt(x.^2 + y.^2);
rx = x./r;
ry = y./r;

for t=0:.1e-4:5e-3

    p=(A./r).*exp(1j*(omega*t-k*r));
    ux = rx.*(1- 1j./(k*r)).*p./(rho0*c);
    uy = ry.*(1- 1j./(k*r)).*p./(rho0*c);
    figure(1)
%     surf(x,y,real(p))
%     colormap(jet)
%     shading interp
%     view(2)
%     caxis([-10 10])
%     axis equal
%     colorbar
%     hold on
    q = quiver(x,y,real(ux),real(uy));
    
    pause(.05)
end