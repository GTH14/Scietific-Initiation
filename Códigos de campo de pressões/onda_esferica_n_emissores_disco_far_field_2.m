freq=1e3; %frequencia [Hz]
c=340; % velocidade do som [m/s]
A=1000; %Amplitude [Pa]
rho0 = 1.2;
omega=2*pi*freq;
k = omega/c;
a = 10/k;
R = 0.1;
[x,y]=meshgrid(-20:0.05:20,-20:0.05:20);
x1 = 0; y1 = .5; z1 = 0;
x2 = 0; y2 = 1; z2 = 0;
x3 = 0; y3 = 1.5; z3 = 0;
x4 = 0; y4 = 2; z4 = 0;
x5 = 0; y5 = -.5; z5 = 0;
x6 = 0; y6 = -1; z6 = 0;
x7 = 0; y7 = -1.5; z7 = 0;
x8 = 0; y8 = -2; z8 = 0;
% x9 = 0; y9 = 2; z9 = 0;
phase1 = 0;
phase2 = 0;
phase3 = 0;
phase4 = 0;
phase5 = pi;
phase6 = pi;
phase7 = pi;
phase8 = pi;

r1=calc_r(x1,y1,z1,x,y);
r2=calc_r(x2,y2,z2,x,y);
r3=calc_r(x3,y3,z3,x,y);
r4=calc_r(x4,y4,z4,x,y);
r5=calc_r(x5,y5,z5,x,y);
r6=calc_r(x6,y6,z6,x,y);
r7=calc_r(x7,y7,z7,x,y);
r8=calc_r(x8,y8,z8,x,y);
% r9=calc_r(x9,y9,z9,x,y);
theta1 = calc_th(x1, y1, z1, x, y);
theta2 = calc_th(x2, y2, z2, x, y);
theta3 = calc_th(x3, y3, z3, x, y);
theta4 = calc_th(x4, y4, z4, x, y);
theta5 = calc_th(x5, y5, z5, x, y);
theta6 = calc_th(x6, y6, z6, x, y);
theta7 = calc_th(x7, y7, z7, x, y);
theta8 = calc_th(x8, y8, z8, x, y);
% theta9 = calc_th(x9, y9, z9, x, y);
Pax1 = calc_Pax(A, r1);
Pax2 = calc_Pax(A, r2);
Pax3 = calc_Pax(A, r3);
Pax4 = calc_Pax(A, r4);
Pax5 = calc_Pax(A, r5);
Pax6 = calc_Pax(A, r6);
Pax7 = calc_Pax(A, r7);
Pax8 = calc_Pax(A, r8);
% Pax9 = calc_Pax(A, r9);
v1 = calc_v(k, a, theta1);
v2 = calc_v(k, a, theta2);
v3 = calc_v(k, a, theta3);
v4 = calc_v(k, a, theta4);
v5 = calc_v(k, a, theta5);
v6 = calc_v(k, a, theta6);
v7 = calc_v(k, a, theta7);
v8 = calc_v(k, a, theta8);
% v9 = calc_v(k, a, theta9);
H1 = calc_H(v1);
H2 = calc_H(v2);
H3 = calc_H(v3);
H4 = calc_H(v4);
H5 = calc_H(v5);
H6 = calc_H(v6);
H7 = calc_H(v7);
H8 = calc_H(v8);
% H9 = calc_H(v9);
% R = 5;

% for t=0:.2e-4:5e-3
    t = 0;
    pressao1 = calc_p(Pax1, H1, omega, k, t, r1, phase1);
    pressao2 = calc_p(Pax2, H2, omega, k, t, r2, phase2);
    pressao3 = calc_p(Pax3, H3, omega, k, t, r3, phase3);
    pressao4 = calc_p(Pax4, H4, omega, k, t, r4, phase4);
    pressao5 = calc_p(Pax5, H5, omega, k, t, r5, phase5);
    pressao6 = calc_p(Pax6, H6, omega, k, t, r6, phase6);
    pressao7 = calc_p(Pax7, H7, omega, k, t, r7, phase7);
    pressao8 = calc_p(Pax8, H8, omega, k, t, r8, phase8);
%     pressao9 = calc_p(Pax9, H9, omega, k, t, r9);
    pressao = pressao1+pressao2+pressao3+pressao4+pressao5+pressao6+pressao7+pressao8;%+pressao9;
    potential_velocity = pressao./(j*omega*rho0);
    [vx, vy] = gradient(potential_velocity, 0.05);
    m_p = 1/2*abs(pressao).^2;
    m_v = 1/2*(vx.*conj(vx)+vy.*conj(vy));
    U = 2*pi*R^3*(m_p/(3*rho0*c^2)-rho0/2*(m_v));
%     i = find(r<R);
%     pressao(i) = NaN;
%     surf(x,y,abs(pressao))
%     colormap(jet)
%     shading interp
%     view(2)
%     caxis([-1e-6 1e-6])
% %     axis equal
%     colorbar
%     pause(.05)
    surf(x, y, real(U))
    colormap(jet)
    shading interp
    view(2)
    caxis([-1e-3 1e-3])
%     axis equal
    colorbar
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