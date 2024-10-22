%% Parâmetros de entrada
clear all
close all
clc
% Do transdutor
transdutor.freq = 4e4; %frequencia [Hz]
transdutor.c = 340; % velocidade do som [m/s]
transdutor.A = 2*.462; %Amplitude [Pa], .462 valor experimental.
transdutor.rho0 = 1.225; %Massa específica do ar [kg/m^3]
transdutor.lambda = transdutor.c/transdutor.freq; %Comprimento de onda [m]
transdutor.omega = 2*pi*transdutor.freq; %frequência angular [rad/s]
transdutor.k = transdutor.omega/transdutor.c; %número de onda [1/m]
transdutor.a = 4.9e-3; %raio do pistão circular [m]
transdutor.disc = 5e-4; %Discretização [m]
disc = transdutor.disc;
[x,y]= meshgrid(-30e-3:disc:30e-3,-30e-3:disc:30.5e-3);
% Posição e direção dos transdutores
[x0_v, y0_v] = meshgrid(-10.16*3.5:10.16:10.16*3.5,-10.16*3.5:10.16:10.16*3.5);
x0_v = 1e-3*reshape(x0_v,[64,1]);
y0_v = 1e-3*reshape(y0_v,[64,1]);
z0_v = zeros(64,1);
transdutor.posx = x0_v;
transdutor.posy = y0_v;
transdutor.posz = z0_v;

transdutor.normal = [0 0 1];
transdutor.foco = [0; 0; 40e-3]; %posição focal [m]
% Da esfera
bolinha.rho = 30; %Densidade da bolinha [kg/m^3]
bolinha.R = 0.75e-3; %Raio da esfera [m]
pos_i = [0; 0; 0.040 ]; %Posição inicial da bolinha [m]
v_i = [0; 0; 0]; %Velocidade inicial da bolinha [m/s]
% Da simulação dinâmica
h = 1e-4; %Discretização do tempo (passo de integração)
ti = 0;
tf = 0; %Tempo final
t_v = ti:h:tf;
transdutor.R = bolinha.R;
%% Variáveis calculadas
dist = sqrt((x0_v-transdutor.foco(1)).^2 + (y0_v-transdutor.foco(2)).^2 + (z0_v-transdutor.foco(3)).^2);
phase_v = 2*pi*dist/transdutor.lambda;
f = find(x0_v>0);
phase_v(f) = phase_v(f)+pi;
transdutor.phase = phase_v;

clear dist phase_v f x0_v y0_v z0_v

%% Simulação dinâmica


P1 = Calc_pressao(x,y,40e-3,transdutor);
U = Calc_Gorkov(x,y,40e-3,transdutor);
figure(1)
surf(x*1e3,y*1e3,abs(P1));
colormap(jet);
shading interp;
view(2);
a = colorbar;
ylabel(a,'Pressão acústica [Pa]');
xlabel('Eixo x [mm]');
xlim([-30e-3,30e-3]*1e3);
ylim([-30e-3, 30e-3]*1e3);
ylabel('Eixo y [mm]');

figure(2)
surf(x*1e3,y*1e3,U);
colormap(jet);
shading interp;
view(2);
a = colorbar;
ylabel(a,'Potencial de Gor''kov [J]');
xlabel('Eixo x [mm]');
xlim([-30e-3,30e-3]*1e3);
ylim([-30e-3, 30e-3]*1e3);
ylabel('Eixo y [mm]');