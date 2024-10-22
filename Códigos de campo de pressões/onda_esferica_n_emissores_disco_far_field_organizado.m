freq=1e2; %frequencia [Hz]
c=340; % velocidade do som [m/s]
A=100; %Amplitude [Pa]
rho0 = 1.2;
omega=2*pi*freq;
k = omega/c;
a = 10/k;
R = 0.1;
disc = 0.05;
[x,y]= meshgrid(-20:disc:20,-20:disc:20);
% [xz,z] = meshgrid(-20:disc:20,-20:disc:20);
x0_v = [20 20 20 20 -20 -20 -20 -20];
y0_v = [-1 -.5 .5 1 -1 -.5 .5 1];
z0_v = [0 0 0 0 0 0 0 0];
phase_v = [0 0 0 0 0 0 0 0];

s = size(x);
n = length(x0_v);
r = zeros(s(1), s(2),n);
theta = zeros(s(1), s(2),n);
Pax = zeros(s(1), s(2),n);
v = zeros(s(1), s(2),n);
H = zeros(s(1), s(2),n);
p = zeros(s(1), s(2),n);
for i = 1:n
    
    r(:,:,i)=calc_r(x0_v(i), y0_v(i), z0_v(i), x, y);
    theta(:,:,i) = calc_th(x0_v(i), y0_v(i), z0_v(i), x, y);
    Pax(:,:,i) = calc_Pax(A, r(:,:,i));
    v(:,:,i) = calc_v(k, a, theta(:,:,i));
    H(:,:,i) = calc_H(v(:,:,i));
   
end


% for t=0:.2e-4:5e-3
t = 0;
    for i = 1:n
        p(:,:,i) = calc_p(Pax(:,:,i), H(:,:,i), omega, k, t, r(:,:,i), phase_v(i));
    end
    p_total = sum(p,3);
    potential_velocity = p_total./(j*omega*rho0);
    [vx, vy] = gradient(potential_velocity, 0.05);
    m_p = 1/2*abs(p_total).^2;
    m_v = 1/2*(vx.*conj(vx)+vy.*conj(vy));
    U = 2*pi*R^3*(m_p/(3*rho0*c^2)-rho0/2*(m_v));
    
    figure(1)
    surf(x,y,real(p_total))
    colormap(jet)
    shading interp
    view(2)
    caxis([-100 100])
%     axis equal
    colorbar
    
    figure(2)
    surf(x, y, real(U))
    colormap(jet)
    shading interp
    view(2)
    caxis([-1e-4 1e-4])
%     axis equal
    colorbar
%     pause(.05)
% end
%     m_p = 1/2*abs(pressao).^2;
%     m_v = 1/2*(vx.*conj(vx)+vy.*conj(vy));
%     U = 2*pi*R^3*(m_p/(3*rho0*c^2)-rho0/2*(m_v));
% end

function [r] = calc_r(x0,y0,z0,x,y)
    r=sqrt((x-x0).^2 + (y-y0).^2 + z0^2); 
end
function [theta] = calc_th(x0, y0, z0, x, y)
    theta = atan((y-y0)./sqrt((x-x0).^2 + z0^2));
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