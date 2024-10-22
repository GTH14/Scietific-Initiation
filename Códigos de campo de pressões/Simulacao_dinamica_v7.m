%% Parâmetros de entrada
clear
clc
close
% Do transdutor
transdutor.freq = 4e4; %frequencia [Hz]
transdutor.c = 340; % velocidade do som [m/s]
transdutor.A = 7*.462; %Amplitude [Pa]
transdutor.rho0 = 1.18; %Massa específica do ar [kg/m^3]
transdutor.lambda = transdutor.c/transdutor.freq; %Comprimento de onda [m]
transdutor.omega = 2*pi*transdutor.freq; %frequência angular [rad/s]
transdutor.k = transdutor.omega/transdutor.c; %número de onda [1/m]
transdutor.a = 4.9e-3; %raio do pistão circular [m]
transdutor.disc = transdutor.lambda/10000*2; %Discretização [m]

% Posição e direção dos transdutores
[x0_v, y0_v] = meshgrid(-10.16*3.5:10.16:10.16*3.5,-10.16*3.5:10.16:10.16*3.5);
x0_v = 1e-3*reshape(x0_v,[64,1]);
y0_v = 1e-3*reshape(y0_v,[64,1]);
z0_v = zeros(64,1);
transdutor.posx = x0_v;
transdutor.posy = y0_v;
transdutor.posz = z0_v;

transdutor.normal = [0 0 1];
transdutor.foco = [0e-3; 0; 32e-3]; %posição focal [m]
% Da esfera
bolinha.rho = 30; %Densidade da bolinha [kg/m^3]
bolinha.R = 1.35e-3; %Raio da esfera [m]
pos_i = [0; 0; 30.5785e-3]; %Posição inicial da bolinha [m]
v_i = [0; 0; 0]; %Velocidade inicial da bolinha [m/s]
% Da simulação dinâmica
h = 0.1e-2; %Discretização do tempo (passo de integração)
ti = 0; 
tf = 1; %Tempo final
t_v = ti:h:tf;

%% Variáveis calculadas
dist = sqrt((x0_v-transdutor.foco(1)).^2 + (y0_v-transdutor.foco(2)).^2 + (z0_v-transdutor.foco(3)).^2);
phase_v = 2*pi*dist/transdutor.lambda;
f = find(x0_v>0);
phase_v(f) = phase_v(f)+pi;
transdutor.phase = phase_v;

clear dist phase_v f x0_v y0_v z0_v

%% Simulação dinâmica
n = length(t_v);
bolinha.x = zeros(1,n); bolinha.y = zeros(1,n); bolinha.z = zeros(1,n);

bolinha.vx = zeros(1,n); bolinha.vy = zeros(1,n); bolinha.vz = zeros(1,n);

i = 1;

bolinha.x(i) = pos_i(1); bolinha.y(i) = pos_i(2); bolinha.z(i) = pos_i(3);

bolinha.vx(i) = v_i(1); bolinha.vy(i) = v_i(2); bolinha.vz(i) = v_i(3);

disc = transdutor.disc;
m = bolinha.rho*4/3*pi*bolinha.R^3;
transdutor.R = bolinha.R;
transdutor.nu = 17.2e-6;
t1 =0;
g = 9.81;
Peso = -m*g;
b1 = 3e-6;
b2 = 3e-6;
% Teste inicial com integração de Euler Explicito
while i <= n-1
%     if i==1000
%         transdutor.foco = [0e-3; 0; 30e-3]; 
%         x0_v = transdutor.posx;
%         y0_v = transdutor.posy;
%         z0_v = transdutor.posz;
%         dist = sqrt((x0_v-transdutor.foco(1)).^2 + (y0_v-transdutor.foco(2)).^2 + (z0_v-transdutor.foco(3)).^2);
%         phase_v = 2*pi*dist/transdutor.lambda;
%         f = find(x0_v>0);
%         phase_v(f) = phase_v(f)+pi;
%         transdutor.phase = phase_v;
%         clear dist phase_v f x0_v y0_v z0_v
%     end
    x = bolinha.x(i); y = bolinha.y(i); z = bolinha.z(i);
    vx = bolinha.vx(i); vy = bolinha.vy(i); vz = bolinha.vz(i);
    v = [vx; vy; vz];
%     Fd = [0;0;0];
%     Fgx = -(Calc_Gorkov_Simulacao_Dinamica(x+disc/2, y, z, transdutor)...
%         - Calc_Gorkov_Simulacao_Dinamica(x-disc/2, y, z, transdutor))/disc;
%     Fgy = -(Calc_Gorkov_Simulacao_Dinamica(x, y+disc/2, z, transdutor)...
%         - Calc_Gorkov_Simulacao_Dinamica(x, y-disc/2, z, transdutor))/disc;
%     Fgz = -(Calc_Gorkov_Simulacao_Dinamica(x, y, z+disc/2, transdutor)...
%         - Calc_Gorkov_Simulacao_Dinamica(x, y, z-disc/2, transdutor))/disc;
    [F] = Acoustokinetics_2(x, y, z, transdutor);
%     Fd = calc_Drag(v, bolinha, transdutor);
    Fd = -b2*v*sqrt(v'*v) - b1*v;
    Fgx = F(1); Fgy = F(2); Fgz = F(3);
    k1x = h*(Fgx+Fd(1))/m; k1y = h*(Fgy+Fd(2))/m; k1z = h*(Fgz+Fd(3)+Peso)/m;
    l1x = h*vx; l1y = h*vy; l1z = h*vz;
    
    v = v + 1/2*[k1x; k1y; k1z];
    Fd = calc_Drag(v, bolinha, transdutor);
    F = Acoustokinetics_2(x+1/2*l1x, y+1/2*l1y, z+1/2*l1z, transdutor);
    Fgx = F(1); Fgy = F(2); Fgz = F(3);
    k2x = h*(Fgx+Fd(1))/m; k2y = h*(Fgy+Fd(2))/m; k2z = h*(Fgz+Fd(3)+Peso)/m;
    l2x = h*(vx+1/2*k1x); l2y = h*(vy+1/2*k1y); l2z = h*(vz+1/2*k1z);
    
    v = v + 1/2*[k2x; k2y; k2z];
    Fd = calc_Drag(v, bolinha, transdutor);
    [F] = Acoustokinetics_2(x+1/2*l2x, y+1/2*l2y, z+1/2*l2z, transdutor);
    Fgx = F(1); Fgy = F(2); Fgz = F(3);
    k3x = h*(Fgx+Fd(1))/m; k3y = h*(Fgy+Fd(2))/m; k3z = h*(Fgz+Fd(3)+Peso)/m;
    l3x = h*(vx+1/2*k2x); l3y = h*(vy+1/2*k2y); l3z = h*(vz+1/2*k2z);
    
    v = v + [k3x; k3y; k3z];
    Fd = calc_Drag(v, bolinha, transdutor);
    [F] = Acoustokinetics_2(x+l3x, y+l3y, z+l3z, transdutor);
    Fgx = F(1); Fgy = F(2); Fgz = F(3);
    k4x = h*(Fgx+Fd(1))/m; k4y = h*(Fgy+Fd(2))/m; k4z = h*(Fgz+Fd(3)+Peso)/m;
    l4x = h*(vx+k3x); l4y = h*(vy+k3y); l4z = h*(vz+k3z);
    
%     vxi = vx + (Fgx+Fd(1))/m*h; vyi = vy + (Fgy+Fd(2))/m*h; vzi = vz + (Fgz+Fd(3)+Peso)/m*h;
    vxi = vx+1/6*(k1x+2*k2x+2*k3x+k4x);vyi = vy+1/6*(k1y+2*k2y+2*k3y+k4y);
    vzi = vz+1/6*(k1z+2*k2z+2*k3z+k4z);
%     bolinha.Fgz(i) = Fgz;
%     xi = x + vxi*h; yi = y + vyi*h; zi = z + vzi*h;
    xi = x+1/6*(l1x+2*l2x+2*l3x+l4x); yi = y+1/6*(l1y+2*l2y+2*l3y+l4y);
    zi = z+1/6*(l1z+2*l2z+2*l3z+l4z);
    
    bolinha.x(i+1) = xi; bolinha.y(i+1) = yi; bolinha.z(i+1) = zi;
    
    bolinha.vx(i+1) = vxi; bolinha.vy(i+1) = vyi; bolinha.vz(i+1) = vzi;
    i = i+1;
end
% [xs,zs] = meshgrid(-10.16*3.5:10.16:10.16*3.5,0:0.1e-3:40e-3);
% U = Calc_pressao(xs,zs,transdutor);
% figure(1)
% surf(xs,zs, U);
% colormap(jet)
% shading interp
% view(2)
% % caxis([-2000 2000])
% axis equal
% colorbar
figure(2)
plot(t_v,bolinha.z);
grid on;
figure(3)
plot(t_v,bolinha.y);
grid on;
figure(4)
plot(t_v,bolinha.x);
grid on;

bolinha.t_v = t_v;
bolinha.h = h;
% 
