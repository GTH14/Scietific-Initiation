

%Add path
addpath(genpath(pwd));

%Create array
array = Class_Array

%Transducers - position%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[x,y]=meshgrid(-10.16*3.5:10.16:10.16*3.5,-10.16*3.5:10.16:10.16*3.5);
x=x*1e-3; %from [mm] to [m]
y=y*1e-3; %from [mm] to [m]
array.xt=[ reshape(x,64,1)];
array.yt=[reshape(y,64,1)];
array.zt=[zeros(64,1)];
%Surface normal vector
array.nx(64,1)=0;
array.ny(64,1)=0;
array.nz=ones(64,1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




c=array.c0;
freq=array.freq;
%focal point
xf=0;
yf=0;
zf=40e-3;
lambda=c/freq;
array.phases=2*pi*sqrt( (array.xt-xf).^2 + (array.yt-yf).^2 + (array.zt-zf).^2)/lambda;

%TWIN TRAP %%%%%%%%%%%%%%
I=find(array.xt>0);
array.phases(I)=array.phases(I)-pi;
%%%%%%%%%%%%%%%%%%%%%%%%%

% % VORTEX BEAM %%%%%%%%%%%
% Zcomplex= array.xt+i*array.yt;
% angulo=angle(Zcomplex);
% array.phases=array.phases-angulo;
% %%%%%%%%%%%%%%%%%%%%%%%%%


[xp,zp]=meshgrid(-70e-3:.5e-3:70e-3,10e-3:.5e-3:190e-3);
[M,N]=size(xp);
xp=reshape(xp,M*N,1);
zp=reshape(zp,M*N,1);
yp=0*ones(M*N,1);
pressure=array.Calculate_pressure(xp,yp,zp);
hold on
surf(1000*reshape(xp,M,N),1000*reshape(yp,M,N),1000*reshape(zp,M,N),abs(reshape(pressure,M,N)))
shading interp
colormap(jet)
xlabel('x [mm]','fontsize',14)
ylabel('y [mm]','fontsize',14)
zlabel('z [mm]','fontsize',14)
h=colorbar
ylabel(h,'Pressure amplitude [Pa]','fontsize',14)
set(gcf','position',[680   232   902   746])
box on

hold on

[xp,yp]=meshgrid(-50e-3:.5e-3:50e-3,-50e-3:.5e-3:50e-3);
[M,N]=size(xp);
xp=reshape(xp,M*N,1);
yp=reshape(yp,M*N,1);
zp=100e-3*ones(M*N,1);
pressure=array.Calculate_pressure(xp,yp,zp);
hold on
surf(1000*reshape(xp,M,N),1000*reshape(yp,M,N),1000*reshape(zp,M,N),abs(reshape(pressure,M,N)))
shading interp
colormap(jet)

%Plot array
array.Draw_transducers
axis equal
view(-150,30)

% % % 
% % % U=array.Calculate_Gorkov_potential(xp,yp,zp);
% % % figure
% % % surf(1000*reshape(xp,M,N),1000*reshape(yp,M,N),1000*reshape(zp,M,N),(reshape(U,M,N)))
% % % shading interp
% % % colormap(jet)
% % % array.Draw_transducers
% % % 
% % % 
% % % %axial force %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % zp=50e-3:.2e-3:150e-3;
% % % xp=zeros(1,length(zp));
% % % yp=zeros(1,length(zp));
% % % 
% % % [Fx_Gorkov,Fy_Gorkov,Fz_Gorkov]=array.Calculate_radiation_force_Gorkov(xp,yp,zp);
% % % [Fx_Acoustokinetics,Fy_Acoustokinetics,Fz_Acoustokinetics]=array.Calculate_radiation_force_Acoustokinetics(xp,yp,zp);
% % % figure(5)
% % % plot(1000*zp,Fz_Gorkov)
% % % hold on
% % % plot(1000*zp,Fz_Acoustokinetics,'r--')
% % % xlabel('z [mm]','fontsize',14)
% % % ylabel('Fz [N]','fontsize',14)
% % % legend('Gorkov','Acoustokinetics')
% % % title('Fz (x = 0 mm, y = 0 mm)')
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % % force Fx %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % xp=-50e-3:.02e-3:50e-3;
% % % yp=zeros(1,length(xp));
% % % zp=100e-3*ones(1,length(xp));
% % % 
% % % [Fx_Gorkov,Fy_Gorkov,Fz_Gorkov]=array.Calculate_radiation_force_Gorkov(xp,yp,zp);
% % % [Fx_Acoustokinetics,Fy_Acoustokinetics,Fz_Acoustokinetics]=array.Calculate_radiation_force_Acoustokinetics(xp,yp,zp);
% % % figure(6)
% % % plot(1000*xp,Fx_Gorkov)
% % % hold on
% % % plot(1000*xp,Fx_Acoustokinetics,'r--')
% % % xlabel('x [mm]','fontsize',14)
% % % ylabel('Fx [N]','fontsize',14)
% % % legend('Gorkov','Acoustokinetics')
% % % title('Fx (y = 0 mm, z = 100 mm)')
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % % force Fy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % yp=-50e-3:.02e-3:50e-3;
% % % xp=zeros(1,length(yp));
% % % zp=100e-3*ones(1,length(xp));
% % % 
% % % [Fx_Gorkov,Fy_Gorkov,Fz_Gorkov]=array.Calculate_radiation_force_Gorkov(xp,yp,zp);
% % % [Fx_Acoustokinetics,Fy_Acoustokinetics,Fz_Acoustokinetics]=array.Calculate_radiation_force_Acoustokinetics(xp,yp,zp);
% % % figure(7)
% % % plot(1000*yp,Fy_Gorkov)
% % % hold on
% % % plot(1000*yp,Fy_Acoustokinetics,'r--')
% % % xlabel('y [mm]','fontsize',14)
% % % ylabel('Fy [N]','fontsize',14)
% % % legend('Gorkov','Acoustokinetics')
% % % title('Fy (x = 0 mm, z = 100 mm)')
% % % %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % 
% % % 
% % % % [Fx,Fy,Fz]=array.Calculate_radiation_force_Acoustokinetics(xp,yp,zp);
% % % % 
% % % % figure
% % % % quiver3(1000*reshape(xp,M,N),1000*reshape(yp,M,N),1000*reshape(zp,M,N),reshape(Fx,M,N),reshape(Fy,M,N),reshape(Fz,M,N))
% % % % view(2)
% % % % xlim([-10 10])
% % % % ylim([-10 10])
% % % 



