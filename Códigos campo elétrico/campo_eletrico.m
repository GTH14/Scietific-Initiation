q = 1; % Unidade [C]
e0 = 8.85E-12; % Unidade [F/m]
x0 = 0;
y0 = 0;
[x, y] = meshgrid(-1:0.05:1, -1:0.05:1);
r = sqrt((x - x0).^2+(y - y0).^2);
r_v_x = (x - x0)./r;
r_v_y = (y - y0)./r;
E_x = q./(4*pi*e0*r.^2).*r_v_x;
E_y = q./(4*pi*e0*r.^2).*r_v_y;
i = find(r < 0.2);
E_x(i) = NaN;
E_y(i) = NaN;
q = quiver(x, y, E_x, E_y);
axis equal
% q.AutoScaleFactor = 2;