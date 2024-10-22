%% Parâmetros de entrada
function [x] = Plot_U(transdutor)
[xs,zs] = meshgrid(-0.05:transdutor.disc:0.05,0:transdutor.disc:0.1);
% zs = 0:1e-4:.3;
% n = length(zs);
s = size(xs);
U = zeros(s)+zeros(s)*i;
p = zeros(s)+zeros(s)*i;

% for n = 1:s(1)
%     for j = 1:s(2)
U = Calc_Gorkov(xs,0,zs,transdutor);
p = Calc_pressao(xs,0,zs,transdutor);
%     end
% end
% max(max(abs(p)))
figure(1)

surf(xs,zs, U);
colormap(jet)
shading interp
view(2)
caxis([-10e-11 0])
axis equal
colorbar

figure(2)
surf(xs,zs,abs(p))
colormap(jet)
shading interp
view(2)
% caxis([0 2000])
axis equal
colorbar
x = 1;
