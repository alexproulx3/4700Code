nx = 50;
ny = nx*3/2;
x = linspace(0,50);
y = linspace(0,75);
[xx,yy] = meshgrid(x,y);
nSum = 100;

z = zeros(length(xx),length(yy));
for n = 1:2:nSum
    z = z + (4/pi)*(1/n)*(cosh(n*pi*xx/ny)/cosh(n*pi*nx/ny)).*sin(n*pi*yy/ny);
end
z = z+fliplr(z);

figure
surf(xx,yy,z)
view(225,30)
title('Numerical solution for Voltage')
xlabel('Y')
ylabel('X')
zlabel('Voltage (V)')