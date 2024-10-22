%% Parâmetros de entrada
clear all
clc
% Do transdutor
transdutor.freq = 4e4; %frequencia [Hz]
transdutor.c = 340; % velocidade do som [m/s]
transdutor.A = 2; %Amplitude [Pa]
transdutor.rho0 = 1.2; %Massa específica do ar [kg/m^3]
transdutor.lambda = transdutor.c/transdutor.freq; %Comprimento de onda [m]
transdutor.omega = 2*pi*transdutor.freq; %frequência angular [rad/s]
transdutor.k = transdutor.omega/transdutor.c; %número de onda [1/m]
transdutor.a = 4.9e-3; %raio do pistão circular [m]
transdutor.disc = 1e-3; %Discretização [m]

% Posição e direção dos transdutores
[x0_v, y0_v] = meshgrid(-10.16*3.5:10.16:10.16*3.5,-10.16*3.5:10.16:10.16*3.5);
x0_v = 1e-3*reshape(x0_v,[64,1]);
y0_v = 1e-3*reshape(y0_v,[64,1]);
z0_v = zeros(64,1);
transdutor.posx = x0_v;
transdutor.posy = y0_v;
transdutor.posz = z0_v;

transdutor.normal = [0 0 1];
% Da esfera
bolinha.rho = 30; %Densidade da bolinha [kg/m^3]
bolinha.R = 0.75e-3; %Raio da esfera [m]
pos_i = [0; 0; 0.046 ]; %Posição inicial da bolinha [m]
v_i = [0; 0; 0]; %Velocidade inicial da bolinha [m/s]
% Da simulação dinâmica
h = 1e-3; %Discretização do tempo (passo de integração)
ti = 0;
tf = 0.5; %Tempo final
t_v = ti:h:tf;

%% Variáveis calculadas
% dist = sqrt((x0_v-transdutor.foco(1)).^2 + (y0_v-transdutor.foco(2)).^2 + (z0_v-transdutor.foco(3)).^2);
% phase_v = 2*pi*dist/transdutor.lambda;
% f = find(y0_v>0);
% phase_v(f) = phase_v(f)+pi;
% transdutor.phase = phase_v;
% 
% clear dist phase_v f x0_v y0_v z0_v

transdutor.R = bolinha.R;
transdutor.nu = 17.2e-6;
transdutor.foco = [0; 0; 0.05]; %posição focal [m]
dist = sqrt((0-transdutor.foco(1)).^2 + (0-transdutor.foco(2)).^2 + (0-transdutor.foco(3)).^2);
phase_v = 2*pi*dist/transdutor.lambda;
f = find(y0_v>0);
% phase_v(f) = phase_v(f)+pi;
transdutor.phase = phase_v;
[xs,zs] = meshgrid(-0.05:transdutor.disc:0.05,0:transdutor.disc:0.1);
pos_x1 = 0;
pos_y1 = 0;
pos_z1 = 0;
r = calc_r(pos_x1, pos_y1, pos_z1, xs, 0,zs);
theta = calc_th(pos_x1, pos_y1, pos_z1,  xs, 0,zs , transdutor.normal, r);
Pax = calc_Pax(transdutor.A, r);
v = calc_v(transdutor.k, transdutor.a, theta);
H = calc_H(v);
p = calc_p(Pax, H, transdutor.omega, transdutor.k, 0, r, transdutor.phase(1));
surf(xs,zs,abs(p))
shading interp
colorbar
view(2)