%% Parâmetros de entrada

% Do transdutor
freq=2e4; %frequencia [Hz]
c=340; % velocidade do som [m/s]
A=10; %Amplitude [Pa]
rho0 = 1.2; %Massa específica do ar [kg/m^3]
lambda = c/freq; %Comprimento de onda [m]
omega=2*pi*freq; %frequência angular [rad/s]
k = omega/c; %número de onda [1/m]
a = 4.9e-3; %raio do pistão circular [m]


% Do espaço trabalhado
disc = 0.005; %Discretização [m]
limxy = .1; %Limite [m]
% limxz = 1;
[x1,y] = meshgrid(-0.45:disc:-0.30,-0.16:disc:0.16); %plano em que se calcula a pressão e o potencial de Gorkov
% [x2,z] = meshgrid(-limxz:disc:limxz,-limxz:disc:limxz);
pos_f = [-.2; .0; 0]; %posição focal [m]

% Posição e direção dos transdutores
x0_v = -0.5*ones(1,12);
y0_v = [-0.055 -.045 -.035 -.025 -.015 -.005 .005 .015 .025 .035 .045 .055];
z0_v = [0 0 0 0 0 0 0 0 0 0 0 0];
normal_v = [1 0 0];

% Da esfera
rho_e = 30; %Densidade da bolinha [kg/m^3]
R = 0.25e-3; %Raio da esfera [m]
pos_i = [-.3; 0]; %Posição inicial da bolinha [m]
v_i = [0; 0]; %Velocidade inicial da bolinha [m/s]

% Da simulação dinâmica
h = 0.001; %Discretização do tempo (passo de integração)
ti = 0;
tf = 30; %Tempo final
t_v = ti:h:tf;
% phase_v = [0 0 0 0 0 0 0 0];

%% Variáveis calculadas
dist_xy = sqrt((x0_v-pos_f(1)).^2 + (y0_v-pos_f(2)).^2 + (z0_v-pos_f(3)).^2);
phase_v = 2*pi*dist_xy/lambda;
clear dist_xy
% i = find(y0_v>0);
% phase_v(i) = phase_v(i)+pi;
s_xy = size(x1);
% s_xz = size(x2);
n = length(x0_v);
r_xy = zeros(s_xy(1), s_xy(2),n);
theta_xy = zeros(s_xy(1), s_xy(2),n);
Pax_xy = zeros(s_xy(1), s_xy(2),n);
v_xy = zeros(s_xy(1), s_xy(2),n);
H_xy = zeros(s_xy(1), s_xy(2),n);
p_xy = zeros(s_xy(1), s_xy(2),n);
for i = 1:n
    
    r_xy(:,:,i)=calc_r_xy(x0_v(i), y0_v(i), z0_v(i), x1, y,0);
%     theta_xy1(:,:,i) = calc_th(x0_v(i), y0_v(i), z0_v(i), x1, y);
    theta_xy(:,:,i) = calc_th1(x0_v(i), y0_v(i), z0_v(i), x1, y, 0 , normal_v, r_xy(:,:,i));
    Pax_xy(:,:,i) = calc_Pax(A, r_xy(:,:,i));
    v_xy(:,:,i) = calc_v(k, a, theta_xy(:,:,i));
    H_xy(:,:,i) = calc_H(v_xy(:,:,i));
   
end
clear theta v_xy
% for i = 1:n
%     
%     r_xz(:,:,i)= calc_r_xy(x0_v(i), y0_v(i), z0_v(i), x2, y_0, z);
%     theta_xz(:,:,i) = calc_th(x0_v(i), y0_v(i), z0_v(i), x2, y, z);
%     Pax_xz(:,:,i) = calc_Pax(A, r_xz(:,:,i));
%     v_xz(:,:,i) = calc_v(k, a, theta_xz(:,:,i));
%     H_xz(:,:,i) = calc_H(v_xz(:,:,i));
%    
% end

%% SIMULAÇÃO ESTÁTICA
% for t=0:.2e-4:5e-3
t = 0;
    for i = 1:n
        p_xy(:,:,i) = calc_p(Pax_xy(:,:,i), H_xy(:,:,i), omega, k, t, r_xy(:,:,i), phase_v(i));
    end
    p_total = sum(p_xy,3);
    potential_velocity = p_total./(j*omega*rho0);
    [vx, vy] = gradient(potential_velocity, 0.05);
    m_p = 1/2*abs(p_total).^2;
    m_v = 1/2*(vx.*conj(vx)+vy.*conj(vy));
    U = 2*pi*R^3*(m_p/(3*rho0*c^2)-rho0/2*(m_v));
    [Fx, Fy] = gradient(U); 
    Fx = -Fx;
    Fy = -Fy;
    figure(1)
    i = find(y == 0);
    plot(x1(33,:), abs(p_total(33,:)))
%     surf(x1,y,abs(p_total))
%     colormap(jet)
%     shading interp
%     view(2)
%     caxis([-2000 2000])
% %     axis equal
%     colorbar
%     
%     figure(2)
%     surf(x1, y, real(U))
%     colormap(jet)
%     shading interp
%     view(2)
%     caxis([-1e-10 1e-10])
% %     axis equal
%     colorbar
%     pause(.05)
% end
% clear U p_total potential_velocity m_p m_v p_xy x1 y
% %% SIMULAÇÃO DINÂMICA
% n = length(t_v);
% xb = zeros(1,n);
% yb = zeros(1,n);
% vxb = zeros(1,n);
% vyb = zeros(1,n);
% i = 1;
% xb(1,i) = pos_i(1);
% yb(1,i) = pos_i(2);
% vxb(1,i) = v_i(1);
% vyb(1,i) = v_i(2);
% 
% m = rho_e*4/3*pi*R^3;
% for i = 1:(n-1)
% %     k = find(x1>-(xb(1,i)+0.05) & x1<(xb(1,i)+0.05) & y>-(yb(1,i)+0.05) & y<(xb(1,i)+0.05));
%     
%     Fq_x1 = interp2(x1, y, Fx, xb(1,i), yb(1,i), 'cubic');
%     Fq_y1 = interp2(x1, y, Fy, xb(1,i), yb(1,i), 'cubic');
%     a_x1 = Fq_x1/m;
%     a_y1 = Fq_y1/m;
%     
%     vx1 = vxb(1,i) + h*a_x1;
%     vy1 = vyb(1,i) + h*a_y1;
%     xb1 = xb(1,i) + h*(vxb(1,i)+vx1)/2;
%     yb1 = yb(1,i) + h*(vyb(1,i)+vy1)/2;
% 
%     
%     Fq_x2 = interp2(x1, y, Fx, xb1, yb1, 'cubic');
%     Fq_y2 = interp2(x1, y, Fy, xb1, yb1, 'cubic');
%     
%     a_x2 = Fq_x2/m;
%     a_y2 = Fq_y2/m;
%     vx2 = vxb(1,i) + h*(a_x2+a_x1)/2;
%     vy2 = vyb(1,i) + h*(a_y2+a_y1)/2;
%     
%     vxb(1,i+1) = vx2;
%     vyb(1,i+1) = vy2;
%     xb(1,i+1) = xb(1,i) + h*(vxb(1,i)+vxb(1,i+1))/2;
%     yb(1,i+1) = yb(1,i) + h*(vyb(1,i)+vyb(1,i+1))/2;
%     
% end
% figure(3)
% plot(xb(1,:), yb(1,:));

%% FUNÇÕES
function [r] = calc_r_xy(x0,y0,z0,x,y,z)
    r=sqrt((x-x0).^2 + (y-y0).^2 + (z-z0).^2); 
end
function [theta] = calc_th(x0, y0, z0, x, y)
    theta = atan(sqrt((y-y0).^2 + z0^2)./(x-x0));
end
function [theta] = calc_th1(x0, y0, z0, x, y, z, normal, r)
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
