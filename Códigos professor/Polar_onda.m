A = 1;
R = 200;
theta = linspace(0,2*pi,1000);
kL = 24;
Pax = A*kL/R;
v = 1/2*kL*sin(theta);
H = abs(sin(v)./v);
P = Pax*H;
figure(1)
polar(theta, log10(100000*P));
figure(2)
plot(180*theta/pi, log10(100000*P))