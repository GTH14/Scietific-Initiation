freq=2e3; %frequencia [Hz]
c=340; % velocidade do som [m/s]
A=1; %Amplitude [Pa]
% rho0 = 	1.2922;
% U0 = 1;
% a = 0.01;
% L = 0.5; %Comprimento do emissor
omega=2*pi*freq;
k=omega/c;
L = 24/k;
[x,y]=meshgrid(-2:.005:2,-2:.005:2);
y0 = 0;
x0 = 0;
dx = 0.01;
x0_v =  (-L/2+x0):dx:(L/2+x0);
s = size(x);
for t=0:.2e-4:5e-3
% t = 0;
pressao = zeros(s);
    for x0_i = x0_v
            r=sqrt((x-x0_i).^2 + (y-y0).^2);

        p_i= (A./r).*exp(j*(omega*t-k*r));
        pressao = p_i+pressao;

    end

    surf(x,y,real(pressao))
    colormap(jet)
    shading interp
    view(2)
    caxis([-10 10])
%     axis equal
    colorbar
    pause(.01)
end