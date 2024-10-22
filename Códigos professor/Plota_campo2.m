

load COOR_CURR.dat
x_wire=COOR_CURR(:,1);
y_wire=COOR_CURR(:,2);
z_wire=COOR_CURR(:,3);
lx_wire=COOR_CURR(:,4);
ly_wire=COOR_CURR(:,5);
lz_wire=COOR_CURR(:,6);

load BFIELD.dat
x_field=BFIELD(:,1);
y_field=BFIELD(:,2);
z_field=BFIELD(:,3);
Bx_field=BFIELD(:,4);
By_field=BFIELD(:,5);
Bz_field=BFIELD(:,6);

B_magnitude=sqrt(Bx_field.^2 + By_field.^2 + Bz_field.^2);

plot3([(x_wire-lx_wire) (x_wire+lx_wire)]',[(y_wire-ly_wire) (y_wire+ly_wire)]',[(z_wire-lz_wire) (z_wire+lz_wire)]','r','linewidth',3) 
axis equal




hold on
q=quiver3(x_field',y_field',z_field',Bx_field'./B_magnitude',By_field'./B_magnitude',Bz_field'./B_magnitude',.5)
q.Color=[0 0 0]

grid on
xlabel('x','fontsize',14)
ylabel('y','fontsize',14)
zlabel('z','fontsize',14)
set(gca,'fontsize',14)
set(gcf,'color','white')

Nx=100;
Ny=10;
Nz=400;
[xx,zz]=meshgrid(linspace(min(x_field),max(x_field),Nx),linspace(min(z_field),max(z_field),Nz));
B_plot = griddata(x_field,z_field,B_magnitude,xx,zz);


B_max=.2
yy=zeros(Nz,Nx);
surf(xx,yy,zz,B_plot);
shading interp
caxis([0 B_max])
colorbar
colormap(hot)


