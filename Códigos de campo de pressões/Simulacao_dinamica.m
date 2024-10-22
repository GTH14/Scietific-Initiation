%% Parâmetros de entrada

% Do transdutor
transdutor.freq = 4e4; %frequencia [Hz]
transdutor.c = 340; % velocidade do som [m/s]
transdutor.A = 4; %Amplitude [Pa]
transdutor.rho0 = 1.2; %Massa específica do ar [kg/m^3]
transdutor.lambda = transdutor.c/transdutor.freq; %Comprimento de onda [m]
transdutor.omega = 2*pi*transdutor.freq; %frequência angular [rad/s]
transdutor.k = transdutor.omega/transdutor.c; %número de onda [1/m]
transdutor.a = 4.9e-3; %raio do pistão circular [m]
transdutor.disc = 5e-4; %Discretização [m]
% Posição e direção dos transdutores
[x0_v, y0_v] = meshgrid(-10.16*3.5:10.16:10.16*3.5,-10.16*3.5:10.16:10.16*3.5);
x0_v = 1e-3*reshape(x0_v,[64,1]);
y0_v = 1e-3*reshape(y0_v,[64,1]);
z0_v = zeros(64,1);
transdutor.posx = x0_v;
transdutor.posy = y0_v;
transdutor.posz = z0_v;
clear x0_v y0_v z0_v

transdutor.normal_v = [1 0 0];
pos_f = [0; .0; 0.5]; %posição focal [m]
% Da esfera
rho_e = 30; transdutor.rho_e = rho_e; %Densidade da bolinha [kg/m^3]
R = 0.25e-3; transdutor.R = R; %Raio da esfera [m]
pos_i = [0; 0.0005; 0.4 ]; %Posição inicial da bolinha [m]
v_i = [0; 0; 0]; %Velocidade inicial da bolinha [m/s]

% Da simulação dinâmica
h = 0.0001; %Discretização do tempo (passo de integração)
ti = 0;
tf = .50; %Tempo final
t_v = ti:h:tf;
% phase_v = [0 0 0 0 0 0 0 0];

%% Variáveis calculadas
dist = sqrt((x0_v-pos_f(1)).^2 + (y0_v-pos_f(2)).^2 + (z0_v-pos_f(3)).^2);
phase_v = 2*pi*dist/lambda;
f = find(y0_v>0);
phase_v(f) = phase_v(f)+pi;
clear dist

%% Simulação dinâmica
n = length(t_v);
xb = zeros(1,n);
yb = zeros(1,n);
zb = zeros(1,n);

vxb = zeros(1,n);
vyb = zeros(1,n);
vzb = zeros(1,n);

index = 1;

xb(1,index) = pos_i(1);
yb(1,index) = pos_i(2);
zb(1,index) = pos_i(3);

vxb(1,index) = v_i(1);
vyb(1,index) = v_i(2);
vzb(1,index) = v_i(3);

m = rho_e*4/3*pi*R^3;
t1 =0;
while index<=n
    [x_m, y_m] = meshgrid((xb(index)-.1e-3):disc:(xb(index)+.1e-3), (yb(index)-.1e-3):disc:(yb(index)+.1e-3));
    s = size(x_m);
    n1 = length(x0_v);
    r = zeros(s(1), s(2),n1);
    theta = zeros(s(1), s(2),n1);
    Pax = zeros(s(1), s(2),n1);
    v = zeros(s(1), s(2),n1);
    H = zeros(s(1), s(2),n1);
    p = zeros(s(1), s(2),n1);
    for i = 1:n1

        r(:,:,i)=calc_r(x0_v(i), y0_v(i), z0_v(i), x_m, y_m,0);
        theta(:,:,i) = calc_th(x0_v(i), y0_v(i), z0_v(i), x_m, y_m,0, normal_v, r(:,:,i));
        Pax(:,:,i) = calc_Pax(A, r(:,:,i));
        v(:,:,i) = calc_v(k, a, theta(:,:,i));
        H(:,:,i) = calc_H(v(:,:,i));
    end
    clear theta v
    for i = 1:n1
        p(:,:,i) = calc_p(Pax(:,:,i), H(:,:,i), omega, k, t1, r(:,:,i), phase_v(i));
    end
    p_total = sum(p,3);
    potential_velocity = p_total./(j*omega*rho0);
    [vx, vy] = gradient(potential_velocity);
    m_p = 1/2*abs(p_total).^2;
    m_v = 1/2*(vx.*conj(vx)+vy.*conj(vy));
    U = 2*pi*R^3*(m_p/(3*rho0*c^2)-rho0/2*(m_v));
    [Fx, Fy] = gradient(U,disc);
    Fx = -Fx;
    Fy = -Fy;
    clear potential velocity vx vy m_p m_v
    Fq_x1 = interp2(x_m, y_m, Fx, xb(index), yb(index), 'cubic');
    Fq_y1 = interp2(x_m, y_m, Fy, xb(index), yb(index), 'cubic');
    a_x1 = Fq_x1/m;
    a_y1 = Fq_y1/m;
    
    vx1 = vxb(index) + h*a_x1;
    vy1 = vyb(index) + h*a_y1;
    xb1 = xb(index) + h*(vxb(1,i)+vx1)/2;
    yb1 = yb(index) + h*(vyb(1,i)+vy1)/2;

    
    Fq_x2 = interp2(x_m, y_m, Fx, xb1, yb1, 'cubic');
    Fq_y2 = interp2(x_m, y_m, Fy, xb1, yb1, 'cubic');
    
    a_x2 = Fq_x2/m;
    a_y2 = Fq_y2/m;
    vx2 = vxb(1,index) + h*(a_x2+a_x1)/2;
    vy2 = vyb(1,index) + h*(a_y2+a_y1)/2;
    
    vxb(1,index+1) = vx2;
    vyb(1,index+1) = vy2;
    xb(1,index+1) = xb(1,index) + h*(vxb(1,index)+vxb(1,index+1))/2;
    yb(1,index+1) = yb(1,index) + h*(vyb(1,index)+vyb(1,index+1))/2;
    index = index+1;
end
hold on
figure(3)
plot(xb(1,:), yb(1,:));
%% FUNÇÕES
function [r] = calc_r(x0,y0,z0,x,y,z)
    r=sqrt((x-x0).^2 + (y-y0).^2 + (z-z0).^2); 
end
% function [theta] = calc_th(x0, y0, z0, x, y)
%     theta = atan(sqrt((y-y0).^2 + z0^2)./(x-x0));
% end
function [theta] = calc_th(x0, y0, z0, x, y, z, normal, r)
    theta = acos((normal(1).*(x-x0) + normal(2).*(y-y0) + normal(3).*(z-z0))./r);
end
function [Pax] = calc_Pax(A, r)
    Pax = A./r;
end
function [v] = calc_v(k, a, theta)
    v = k*a*sin(theta)+eps;
end
function [H] = calc_H(v)
    H = 2*besselj(1,v)./v;
end
function [p] = calc_p(Pax, H, omega, k, t, r, phase)
    p = Pax.*H.*exp(j*(omega*t - k*r + phase));
end