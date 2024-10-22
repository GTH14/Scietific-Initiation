q = 1; % Carga total Unidade [C]
e0 = 8.85E-12; % Unidade [F/m]
[x, y] = meshgrid(-1:0.02:1, -1:0.02:1);
n_pontos_t = 10000;
dq = q/n_pontos_t;
R = 0.5; % Raio do disco [m]
y0 = 0;
n = 2000;
Raios = linspace(R/n,R,n);
cte = n_pontos_t/sum(1./Raios);
disco_x = [];
disco_z = [];
for i = 1:n
    Raio = Raios(i);
    n_pontos = round(cte/Raio);
    disco_x_i = Raio*cos(linspace(0,2*pi-2*pi/n_pontos,n_pontos));
    disco_x = [disco_x disco_x_i];
    disco_z_i = Raio*sin(linspace(0,2*pi-2*pi/n_pontos,n_pontos));
    disco_z = [disco_z disco_z_i];
end
% anel_x = R*cos(linspace(0,2*pi-2*pi/n_pontos,n_pontos));
% anel_z = R*sin(linspace(0,2*pi-2*pi/n_pontos,n_pontos));
s = size(x);
E_x = zeros(s);
E_y = zeros(s);

for i = 1:n_pontos_t
    r = sqrt((x - disco_x(i)).^2+(y - y0).^2 + disco_z(i)^2);
    r_v_x = (x - disco_x(i))./r;
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