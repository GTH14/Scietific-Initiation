%% Parâmetros de entrada

% Do transdutor
transdutor.freq=2e4; %frequencia [Hz]
transdutor.c=340; % velocidade do som [m/s]
transdutor.A=10; %Amplitude [Pa]
transdutor.rho0 = 1.2; %Massa específica do ar [kg/m^3]
transdutor.lambda = transdutor.c/transdutor.freq; %Comprimento de onda [m]
transdutor.omega=2*pi*transdutor.freq; %frequência angular [rad/s]
transdutor.k = transdutor.omega/transdutor.c; %número de onda [1/m]
transdutor.a = 4.9e-3; %raio do pistão circular [m]

% Posição e direção dos transdutores
[x0_v, y0_v] = meshgrid(-10.16*3.5:10.16:10.16*3.5,-10.16*3.5:10.16:10.16*3.5);
x0_v = 1e-3*reshape(x0_v,[64,1]);
y0_v = 1e-3*reshape(y0_v,[64,1]);
z0_v = zeros(64,1);
transdutor.posx = x0_v;
transdutor.posy = y0_v;
transdutor.posz = z0_v;

transdutor.normal = [1 0 0];
transdutor.foco = [0; 0; 0.2]; %posição focal [m]
% Da esfera
bolinha.rho = 30; %Densidade da bolinha [kg/m^3]
bolinha.R = 0.25e-3; %Raio da esfera [m]
pos_i = [0; 0.0005; 0.2 ]; %Posição inicial da bolinha [m]
v_i = [0; 0; 0]; %Velocidade inicial da bolinha [m/s]
% Da simulação dinâmica


%% Variáveis calculadas
dist = sqrt((x0_v-transdutor.foco(1)).^2 + (y0_v-transdutor.foco(2)).^2 + (z0_v-transdutor.foco(3)).^2);
phase_v = 2*pi*dist/transdutor.lambda;
f = find(y0_v>0);
phase_v(f) = phase_v(f)+pi;
transdutor.phase = phase_v;
iter1 = 1;
max1 = 1e3;

clear dist phase_v f x0_v y0_v z0_v
%% Simulação dinâmica
err = 1;
while or(err > 0.001, iter1 == max1)
    
    transdutor.disc = 1e-5; %Discretização [m]
%     for iter2 = 1:100:1e3
    h = 1e-3/iter1; %Discretização do tempo (passo de integração)
    ti = 0;
    tf = .50; %Tempo final
    t_v = ti:h:tf;
    n = length(t_v);
    bolinha.x = zeros(1,n); bolinha.y = zeros(1,n); bolinha.z = zeros(1,n);

    bolinha.vx = zeros(1,n); bolinha.vy = zeros(1,n); bolinha.vz = zeros(1,n);

    i = 1;

    bolinha.x(i) = pos_i(1); bolinha.y(i) = pos_i(2); bolinha.z(i) = pos_i(3);

    bolinha.vx(i) = v_i(1); bolinha.vy(i) = v_i(2); bolinha.vz(i) = v_i(3);

    disc = transdutor.disc;
    m = bolinha.rho*4/3*pi*bolinha.R^3;
    transdutor.R = bolinha.R;
    t1 =0;

    % Teste inicial com integração de Euler Explicito
    while i<=n-1
        x = bolinha.x(i); y = bolinha.y(i); z = bolinha.z(i);
        vx = bolinha.vx(i); vy = bolinha.vy(i); vz = bolinha.vz(i);
        Fx = (Calc_Gorkov(x+disc/2, y, z, transdutor) - Calc_Gorkov(x-disc/2, y, z, transdutor))/disc;
        Fy = (Calc_Gorkov(x, y+disc/2, z, transdutor) - Calc_Gorkov(x, y-disc/2, z, transdutor))/disc;
        Fz = (Calc_Gorkov(x, y, z+disc/2, transdutor) - Calc_Gorkov(x, y, z-disc/2, transdutor))/disc;

        vxi = vx + Fx/m*h; vyi = vy + Fy/m*h; vzi = vz + Fz/m*h;

        xi = x + vxi*h; yi = y + vyi*h; zi = z + vzi*h;

        bolinha.x(i+1) = xi; bolinha.y(i+1) = yi; bolinha.z(i+1) = zi;

        bolinha.vx(i+1) = vxi; bolinha.vy(i+1) = vyi; bolinha.vz(i+1) = vzi;

        clear vxi vyi vzi xi yi zi
        i = i+1;
    end
    if iter1 ~= 1
        err = abs(bolinha.z(end)-oldz(end))/bolinha.z(end) ;
    end
    oldx = bolinha.x; oldy = bolinha.y; oldz = bolinha.z;
%     end
    iter1 = iter1+10;
    figure(1)
    plot(t_v, bolinha.z)
    hold on
end