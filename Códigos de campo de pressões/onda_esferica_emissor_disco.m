freq=2e3; %frequencia [Hz]
c=340; % velocidade do som [m/s]
A=1; %Amplitude [Pa]
% rho0 = 1.2922;
% U0 = 1;
% a = 0.01;
% L = 0.5; %Comprimento do emissor
omega=2*pi*freq;
k=omega/c;
L = 16*pi/k;
[x,y]=meshgrid(-2:.01:2,-2:.01:2);
% R = 200;
% theta = linspace(0,2*pi,1000);
% x = R*sin(theta);
% y = R*cos(theta);
y0 = 0;
x0 = 0;
z0 = 0;
dx = 0.05;
x0_v =  (-L/2+x0):dx:(L/2+x0);
s = size(x);
for t=0:.5e-4:5e-3
% t = 0;
    pressao = zeros(s);
    for x0_i = x0_v
        Lz = sqrt((L/2)^2-(x0_i-x0)^2);
        z0_v = (-Lz+z0):dx:(Lz+z0);
        for z0_i = z0_v
            r=sqrt((x-x0_i).^2 + (y-y0).^2 + (z0_i-z0)^2);

            p_i= (A./r).*exp(j*(omega*t-k*r));
            pressao = p_i+pressao;
        end
    end
    surf(x,y,abs(pressao))
    colormap(jet)
    shading interp
    view(2)
%     caxis([-0 100])
% %     axis equal
%     colorbar
%     pause(.1)
end