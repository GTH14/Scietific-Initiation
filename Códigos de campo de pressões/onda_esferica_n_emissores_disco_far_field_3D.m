%% Parâmetros de entrada
freq=2e4; %frequencia [Hz]
c=340; % velocidade do som [m/s]
A=10; %Amplitude [Pa]
rho0 = 1.2; %Massa específica do ar [kg/m^3]
lambda = c/freq; %Comprimento de onda [m]
omega=2*pi*freq; %frequência angular [rad/s]
k = omega/c; %número de onda [1/m]
a = 4.9e-3; %raio do pistão circular [m]
R = 0.25e-3; %Raio da esfera [m]
dx = 0.005; %Discretização [m]
limxy = .5; %Limite [m]
% limxz = 1;
[x,y,z] = meshgrid(-limxy:dx:limxy,-limxy:dx:limxy,-limxy:dx:limxy ); %plano em que se calcula a pressão e o potencial de Gorkov
pos_f = [-.4; .0; 0]; %posição focal
pos_inicial = [-.4; 0; 0];
x0_v = -limxy*ones(1,12);
y0_v = [-0.055 -.045 -.035 -.025 -.015 -.005 .005 .015 .025 .035 .045 .055];
z0_v = [0 0 0 0 0 0 0 0 0 0 0 0];
normal_v = [1 0 0];
% phase_v = [0 0 0 0 0 0 0 0];

%% Variáveis calculadas
dist_xy = sqrt((x0_v-pos_f(1)).^2 + (y0_v-pos_f(2)).^2 + (z0_v-pos_f(3)).^2);
phase_v = 2*pi*dist_xy/lambda;
clear dist_xy
i = find(y0_v>0);
phase_v(i) = phase_v(i)+pi;
s_xyz = size(x);
% s_xz = size(x2);
n = length(x0_v);
r = zeros(s_xyz(1), s_xyz(2),s_xyz(3),n);
theta = zeros(s_xyz(1), s_xyz(2),s_xyz(3),n);
Pax = zeros(s_xyz(1), s_xyz(2),s_xyz(3),n);
v_xy = zeros(s_xyz(1), s_xyz(2),s_xyz(3),n);
H_xy = zeros(s_xyz(1), s_xyz(2),s_xyz(3),n);
p = zeros(s_xyz(1), s_xyz(2),s_xyz(3),n);
for i = 1:n
    
    r(:,:,:,i) = calc_r_xy(x0_v(i), y0_v(i), z0_v(i), x, y, z);
    theta(:,:,:,i) = calc_th(x0_v(i), y0_v(i), z0_v(i), x, y, z, normal_v, r(:,:,:,i));
%     theta_alt(:,:,:,i) = calc_th1(x0_v(i), y0_v(i), z0_v(i), x, y, z);
    Pax(:,:,:,i) = calc_Pax(A, r(:,:,:,i));
    v_xy(:,:,:,i) = calc_v(k, a, theta(:,:,:,i));
    H_xy(:,:,:,i) = calc_H(v_xy(:,:,:,i));
   
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
        p(:,:,:,i) = calc_p(Pax(:,:,:,i), H_xy(:,:,:,i), omega, k, t, r(:,:,:,i), phase_v(i));
        
    end
    clear H_xy Pax r
    p_total = sum(p,4);
    potential_velocity = p_total./(j*omega*rho0);
    [vx, vy, vz] = gradient(potential_velocity, 0.05);
    m_p = 1/2*abs(p_total).^2;
    m_v = 1/2*(vx.*conj(vx)+vy.*conj(vy));
    U = 2*pi*R^3*(m_p/(3*rho0*c^2)-rho0/2*(m_v));
    clear m_v m_p
    [Fx, Fy] = gradient(U);
    Fx = -Fx;
    Fy = -Fy;
    figure(1)
    slice(x,y,z,abs(p_total),-0.4,0,0)
    colormap(jet)
    shading interp
%     view(2)
    caxis([-2000 2000])
% %     axis equal
    colorbar
%     
%     figure(2)
%     surf(x, y, real(U))
%     colormap(jet)
%     shading interp
%     view(2)
%     caxis([-1e-10 1e-10])
% %     axis equal
%     colorbar
%     pause(.05)
% end

%% SIMULAÇÃO DINÂMICA
% for t=0:dt:tf
%     
%     
% end


%% FUNÇÕES
function [r] = calc_r_xy(x0,y0,z0,x,y,z)
    r=sqrt((x-x0).^2 + (y-y0).^2 + (z-z0).^2); 
end
function [theta] = calc_th(x0, y0, z0, x, y, z, normal, r)
    theta = acos((normal(1).*(x-x0) + normal(2).*(y-y0) + normal(3).*(z-z0))./r);
end
function [theta] = calc_th1(x0, y0, z0, x, y, z)
    theta = atan(sqrt((y-y0).^2 + (z-z0).^2)./(x-x0));
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