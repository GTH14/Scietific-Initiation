freq=2e3; %frequencia [Hz]
c=340; % velocidade do som [m/s]
A=1; %Amplitude [Pa]
rho0= 1.2922;
c = 340;

%emissores
x1=.25;
y1=0;
x2=-.25;
y2=0;
x3 = .15;
y3 = 0;
x4 = .05;
y4 = 0;
x5 = -.05;
y5 = 0;
x6 = -.15;
y6 = 0;
omega=2*pi*freq;
k=omega/c;
[x,y]=meshgrid(-2:.05:2,-2:.05:2);


r1=sqrt((x-x1).^2 + (y-y1).^2);
r2=sqrt((x-x2).^2 + (y-y2).^2);
r3=sqrt((x-x3).^2 + (y-y3).^2);
r4=sqrt((x-x4).^2 + (y-y4).^2);
r5=sqrt((x-x5).^2 + (y-y5).^2);
r6=sqrt((x-x6).^2 + (y-y6).^2);
rx1 = (x-x1)./r1; ry1 = (y-y1)./r1;
rx2 = (x-x2)./r2; ry2 = (y-y2)./r2;
rx3 = (x-x3)./r3; ry3 = (y-y3)./r3;
rx4 = (x-x4)./r4; ry4 = (y-y4)./r4;
rx5 = (x-x5)./r5; ry5 = (y-y5)./r5;
rx6 = (x-x6)./r6; ry6 = (y-y6)./r6;
for t=0:.05e-4:1e-3
% t = .2e-4;
    p1=(A./r1).*exp(j*(omega*t-k*r1));
    p2=(A./r2).*exp(j*(omega*t-k*r2));
    p3=(A./r3).*exp(j*(omega*t-k*r3));
    p4=(A./r4).*exp(j*(omega*t-k*r4));
    p5=(A./r5).*exp(j*(omega*t-k*r5));
    p6=(A./r6).*exp(j*(omega*t-k*r6));
    [ux1,uy1] = calc_u(rx1, ry1, k, r1, p1, rho0, c);
    [ux2,uy2] = calc_u(rx2, ry2, k, r2, p2, rho0, c);
    [ux3,uy3] = calc_u(rx3, ry3, k, r3, p3, rho0, c);
    [ux4,uy4] = calc_u(rx4, ry4, k, r4, p4, rho0, c);
    [ux5,uy5] = calc_u(rx5, ry5, k, r5, p5, rho0, c);
    [ux6,uy6] = calc_u(rx6, ry6, k, r6, p6, rho0, c);
    ux = ux1+ux2+ux3+ux4+ux5+ux6;
    uy = uy1+uy2+uy3+uy4+uy5+uy6;
    p= p1 + p2 + p3 + p4 + p5 + p6;
    
%     figure(1)
%     surf(x,y,real(p))
%     colormap(jet)
%     shading interp
%     view(2)
%     caxis([-10 10])
%     axis equal
%     colorbar
    figure(2)
    u = sqrt(real(ux).^2+real(uy).^2);
    i = find(u>0.5E-1);
    ux(i) = NaN;
    uy(i) = NaN;
    q = quiver(x, y, real(ux), real(uy));
    pause(.05)
end
function [ux,uy] = calc_u(rx, ry, k, r, p, rho0, c)
    ux = rx.*(1- 1j./(k*r)).*p./(rho0*c);
    uy = ry.*(1- 1j./(k*r)).*p./(rho0*c);
end