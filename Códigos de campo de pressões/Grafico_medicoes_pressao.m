load('campo_pressao_plano_xy.mat');
close all
P = tensao_RMS*sqrt(2)/(1e-3*1.6); % Valor de conversão 1e-3 verificar!
x = -30:0.5:30;
y = -30.5:0.5:31;
figure(1);
% surf(gradient(P))
surf(y,x,P);
colormap(jet);
shading interp;
view(2);
a = colorbar;
ylabel(a,'Pressão acústica [Pa]');
xlabel('Eixo x [mm]');
xlim([-30,30]);
ylim([-30,30]);
ylabel('Eixo y [mm]');
% caxis([0,1100])