q = 1; % Carga total Unidade [C]
e0 = 8.85E-12; % Unidade [F/m]
R = 0.5; % Raio do disco [m]
mi = q/(pi*R^2);
dx = 0.005; % Elemento infinitesimal de comprimento [m]
dq = mi*dx^2;
x0 = 0;
z0 = 0;
x0_v = (-R+x0):dx:(R+x0);
% s_x0 = length(x0_v);
% dq = q/s_x0;
% for xi = x0_v
%     Lz = 2*sqrt(R^2-(xi-x0)^2);
%     z0_i = (-Lz+z0):dx:
% end
y0 = 0;
[x, y] = meshgrid(-1:0.02:1, -1:0.02:1);
s = size(x);
E_x = zeros(s);
E_y = zeros(s);
for x0_1 = x0_v
    Lz = sqrt(R^2-(x0_1 - x0)^2);
    z0_v = (-Lz+z0):dx:(Lz+z0);
    for z0_1 = z0_v
        r = sqrt((x - x0_1).^2+(y - y0).^2 + (z0_1-z0)^2);
        r_v_x = (x - x0_1)./r;
        r_v_y = (y - y0)./r;
        E_x_n = dq./(4*pi*e0*r.^2).*r_v_x;
        E_y_n = dq./(4*pi*e0*r.^2).*r_v_y;
        E_x = E_x + E_x_n;
        E_y = E_y + E_y_n;
    end
    %q.AutoScaleFactor = 2;

end
E = sqrt(E_x.^2 + E_y.^2);
i = find(E > 1E12);
E_x(i) = NaN;
E_y(i) = NaN;
figure(1)
qui = quiver(x, y, E_x, E_y);
axis equal
i_y = find(y>0);
y_1 = y(i_y);
E_y_1 = E_y(i_y);
x_1 = x(i_y);
% i_x = find(x == 0);
% E_x(i_x) = NaN;
% E_y(i_x) = NaN;
% E_analitico = q*y(i_x)./(4*pi*e0*(y(i_x).^2 + R^2).^(3/2));
E_analitico = mi/(2*e0)*(1-(y_1(x_1==0)./sqrt(y_1(x_1==0).^2 + R^2)));
figure(2)
% qui_2 = quiver(x, y, E_x, E_y);
plot(y_1(x_1==0), E_y_1(x_1==0))
hold on
plot(y_1(x_1==0), E_analitico, 'r.')