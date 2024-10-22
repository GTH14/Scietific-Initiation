freq=.5e3; %frequencia [Hz]
c=340; % velocidade do som [m/s]
A=100; %Amplitude [Pa]
omega=2*pi*freq;
k = omega/c;
a = 10/k;
[x,y]=meshgrid(-20:.1:20,-20:.1:20);
r=sqrt(x.^2 + y.^2);

theta = atan(y./x)+eps;
Pax = (A./r);
v = k*a*sin(theta);
H = abs(2*besselj(1,v)./v);
% R = 5;

for t=0:.2e-4:5e-3
    pressao = Pax.*H.*exp(exp(j*(omega*t-k*r)));
%     i = find(r<R);
%     pressao(i) = NaN;
    surf(x,y,real(pressao))
    colormap(jet)
    shading interp
    view(2)
    caxis([-10 10])
%     axis equal
    colorbar
    pause(.05)
end
