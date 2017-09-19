clear all

fin = 80;%how many points are used for the meshgrid

iso = 3e-2;

x = [-12 12];
y = [-12 12];
z = [-8 8];

[X,Y,Z] = meshgrid(linspace(x(1),x(2),fin),linspace(y(1),y(2),fin),linspace(z(1),z(2),fin));
a = 2
b = 2
c = 2 
orbital('Cu', 2, a, b, c, 'imag')