%% Parâmetros de entrada
clear all
clc
load('resultados.mat')
% Do transdutor
transdutor.freq = 4e4; %frequencia [Hz]
transdutor.c = 340; % velocidade do som [m/s]
transdutor.A = 7*.462; %Amplitude [Pa], .462 valor experimental.
transdutor.rho0 = 1.1839; %Massa específica do ar [kg/m^3]
transdutor.lambda = transdutor.c/transdutor.freq; %Comprimento de onda [m]
transdutor.omega = 2*pi*transdutor.freq; %frequência angular [rad/s]
transdutor.k = transdutor.omega/transdutor.c; %número de onda [1/m]
transdutor.a = 4.9e-3; %raio do pistão circular [m]
transdutor.disc = 1e-5; %Discretização [m]
disc = transdutor.disc;
[x,z]= meshgrid(-50e-3:disc:50e-3,-71.5e-3:disc:71.5e-3);
% Posição e direção dos transdutores
[x0_v, y0_v] = meshgrid(-10.16*3.5:10.16:10.16*3.5,-10.16*3.5:10.16:10.16*3.5);
x0_v = 1e-3*reshape(x0_v,[64,1]);
x0_v = [x0_v; x0_v];
y0_v = 1e-3*reshape(y0_v,[64,1]);
y0_v = [y0_v; y0_v];
z_dist = 94.5e-3;
z0_v1 = z_dist/2*ones(64,1);
z0_v2 = -z_dist/2*ones(64,1);
z0_v = [z0_v1; z0_v2];
transdutor.posx = x0_v;
transdutor.posy = y0_v;
transdutor.posz = z0_v;

transdutor.normal = [0 0 1];
% transdutor.normal2 = [0 0 -1];
transdutor.foco = [0;0 ; 0]; %posição focal [m]
% transdutor.foco = [0; 0; 40e-3]; %posição focal [m]
% Da esfera
bolinha.rho = 1000; %Densidade da bolinha [kg/m^3]
bolinha.R =.43e-3; %Raio da esfera [m]
pos_i = [0e-3;0; 0.e-3 ]; %Posição inicial da bolinha [m]
v_i = [0; 0; 0]; %Velocidade inicial da bolinha [m/s]
% Da simulação dinâmica
h = 5e-4; %Discretização do tempo (passo de integração)
ti = 0;
tf = 12; %Tempo final
t_v = ti:h:tf;
% transdutor.R = bolinha.R;
%% Variáveis calculadas
dist = sqrt((x0_v-transdutor.foco(1)).^2 + (y0_v-transdutor.foco(2)).^2 + (z0_v-transdutor.foco(3)).^2);
phase_v = 2*pi*dist/transdutor.lambda;
f = find(z0_v>0);
phase_v(f) = phase_v(f)+pi;
% phase_v = zeros(128,1);
transdutor.phase = phase_v;

% clear dist phase_v f x0_v y0_v z0_v

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
b1 = 0.3e-6;
b2 = 0.3e-6;
% Teste inicial com integração de Euler Explicito
while i <= n-1
    if i==1850*2
        transdutor.foco = [0; 0; 0.5e-3]; 
        x0_v = transdutor.posx;
        y0_v = transdutor.posy;
        z0_v = transdutor.posz;
        dist = sqrt((x0_v-transdutor.foco(1)).^2 + (y0_v-transdutor.foco(2)).^2 + (z0_v-transdutor.foco(3)).^2);
        phase_v = 2*pi*dist/transdutor.lambda;
        f = find(z0_v>0);
        phase_v(f) = phase_v(f)+pi;
        transdutor.phase = phase_v;
        clear dist phase_v f x0_v y0_v z0_v
    end
    x = bolinha.x(i); y = bolinha.y(i); z = bolinha.z(i);
    vx = bolinha.vx(i); vy = bolinha.vy(i); vz = bolinha.vz(i);
    v = [vx; vy; vz];
%     Fd = calc_Drag(v, bolinha, transdutor);
    Fd = -b2*v*sqrt(v'*v) - b1*v;
%     Fd = [0;0;0];
    Fgx = -(Calc_Gorkov_estacionario(x+disc/2, y, z, transdutor) -...
        Calc_Gorkov_estacionario(x-disc/2, y, z, transdutor))/disc;
    Fgy = -(Calc_Gorkov_estacionario(x, y+disc/2, z, transdutor) -...
        Calc_Gorkov_estacionario(x, y-disc/2, z, transdutor))/disc;
    Fgz = -(Calc_Gorkov_estacionario(x, y, z+disc/2, transdutor) -...
        Calc_Gorkov_estacionario(x, y, z-disc/2, transdutor))/disc;
    
    vxi = vx + (Fgx+Fd(1))/m*h; vyi = vy + (Fgy+Fd(2))/m*h; vzi = vz + (Fgz+Fd(3))/m*h;
    
    xi = x + vxi*h; yi = y + vyi*h; zi = z + vzi*h;
    
    bolinha.x(i+1) = xi; bolinha.y(i+1) = yi; bolinha.z(i+1) = zi;
    
    bolinha.vx(i+1) = vxi; bolinha.vy(i+1) = vyi; bolinha.vz(i+1) = vzi;
    i = i+1;
end
% U = Calc_Gorkov(bolinha.x,bolinha.y,bolinha.z,transdutor);
% figure(1)
% surf(xs,zs, U);
% colormap(jet)
% shading interp
% view(2)
% caxis([-2000 2000])
% axis equal
% colorbar
% figure(2)
% plot(t_v,bolinha.z);
% figure(3)
% plot(t_v,bolinha.y);
% figure(4)
% plot(t_v,bolinha.x);
% bolinha.t_v = t_v;
% bolinha.h = h;
% % 
x_CM = x_CM/26.6-17.55;
y_CM = 13.2415-y_CM/26.6;
t = 0:4e-3:(length(x_CM)-1)*4e-3;
figure(2)
plot(t, y_CM);
hold on;
plot(t_v,bolinha.z*1e3);
xlim([0, 12])
grid on;
xlabel('tempo [s]');
ylabel('z [mm]');
