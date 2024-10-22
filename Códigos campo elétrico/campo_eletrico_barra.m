q = 1; % Carga total Unidade [C]
e0 = 8.85E-12; % Unidade [F/m]
L = 0.5; % Raio do disco [m]
dx = 0.01; % Elemento infinitesimal de comprimento [m]
x0 = 0;
z0 = 0;
x0_v = (-L/2+x0):dx:(L/2+x0);
s_x0 = length(x0_v);
dq = q/s_x0;
y0 = 0;
[x, y] = meshgrid(-1:0.02:1, -1:0.02:1);
s = size(x);
E_x = zeros(s);
E_y = zeros(s);
for x0_1 = x0_v
    r = sqrt((x - x0_1).^2+(y - y0).^2);
    r_v_x = (x - x0_1)./r;
    r_v_y = (y - y0)./r;
    E_x_n = dq./(4*pi*e0*r.^2).*r_v_x;
    E_y_n = dq./(4*pi*e0*r.^2).*r_v_y;
    E_x = E_x + E_x_n;
    E_y = E_y + E_y_n;
end
    %q.AutoScaleFactor = 2;

E = sqrt(E_x.^2 + E_y.^2);
i = find(E > 1E13);
E_x(i) = NaN;
E_y(i) = NaN;
qui = quiver(x, y, E_x, E_y);
axis equal