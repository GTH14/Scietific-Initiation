
freq=1e3; %frequencia [Hz]
c=340; % velocidade do som [m/s]
A=1; %Amplitude [Pa]
R = 0.1;
rho0 =1.2;
omega=2*pi*freq;
k=omega/c;
[x,y]=meshgrid(-1:.005:1,-1:.005:1);

r=sqrt(x.^2 + y.^2);

% for t=0:.5e-4:5e-3
% t = 0
    p=(A./r).*exp(j*(omega*t-k*r));
    media_p = 1/2*abs(p).^2;
    [dpx, dpy] = gradient(p, 0.005);
    vx = dpx/(j*omega*rho0);
    vy = dpy/(j*omega*rho0);
    media_v = (vx.*conj(vx)+vy.*conj(vy))/2;
    U_Gorkov = 2*pi*R^3*(media_p/(3*rho0*c^2) - rho0*media_v/2);
    i = find(sqrt(x.^2+y.^2)<0.05);
    U_Gorkov(i) = NaN;
    surf(x,y,U_Gorkov)
    
%     figure(1)
%     surf(x,y,real(p))
%     title(num2str(t))
    colormap(jet)
    shading interp
    view(2)
%     caxis([-10 10])
    axis equal
    colorbar
%     figure(2)
%     vx(i) = NaN;
%     vy(i) = NaN;
%     q = quiver(x,y,real(vx), real(vy));
%     q.AutoScaleFactor = 2;
%     title(num2str(t))
%     axis equal
%     pause(.05)
% end