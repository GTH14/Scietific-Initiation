q = 1; % Carga total Unidade [C]
e0 = 8.85E-12; % Unidade [F/m]
n_pontos = 1000;
dq = q/n_pontos;
R = 0.5; % Raio do anel [m]
y0 = 0;
anel_x = R*cos(linspace(0,2*pi-2*pi/n_pontos,n_pontos));
anel_z = R*sin(linspace(0,2*pi-2*pi/n_pontos,n_pontos));
[x, y] = meshgrid(-1:0.02:1, -1:0.02:1);
s = size(x);
E_x = zeros(s);
E_y = zeros(s);

for i = 1:n_pontos
    r = sqrt((x - anel_x(i)).^2+(y - y0).^2 + anel_z(i)^2);
    r_v_x = (x - anel_x(i))./r;
    r_v_y = (y - y0)./r;
    E_x_n = dq./(4*pi*e0*r.^2).*r_v_x;
    E_y_n = dq./(4*pi*e0*r.^2).*r_v_y;
    E_x = E_x + E_x_n;
    E_y = E_y + E_y_n;
end
E = sqrt(E_x.^2 + E_y.^2);
i = find(E > 1E12);
E_x(i) = NaN;
E_y(i) = NaN;
figure(1)
qui = quiver(x, y, E_x, E_y);
axis equal
i_x = find(x == 0);
% E_x(i_x) = NaN;
% E_y(i_x) = NaN;
E_analitico = q*y(i_x)./(4*pi*e0*(y(i_x).^2 + R^2).^(3/2));
figure(2)
qui_2 = quiver(x, y, E_x, E_y);
plot(y(i_x), E_y(i_x))
hold on
plot(y(i_x), E_analitico, 'r.')