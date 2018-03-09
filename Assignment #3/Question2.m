nx = 50;
ny = nx*3/2;
% nx = 200;
% ny = 100;
sigma0 = 1;
sigma1 = 1e-2;
dim = 1e-6;
Box = [nx*0.4 nx*0.6 ny*0.4 ny*0.6];

G = sparse(nx+(ny-1)*ny);
B = zeros(1,nx+(ny-1)*ny);

%sigma plot
sigma = zeros(nx,ny);
for i = 1:nx
    for j = 1:ny
        if i>Box(1)&&i<Box(2)&&(j<Box(3)||j>Box(4))
            sigma(i,j) = sigma1;
        else
            sigma(i,j) = sigma0;
        end
    end
end

figure
surf(sigma)
view(90,90)
axis 'tight'

for i = 1:nx
    for j = 1:ny
        n = i + (j - 1) * ny;
        if i==1
            G(n,:) = 0;
            G(n,n) = 1;
            B(n) = 1;
        elseif i==nx
            G(n,:) = 0;
            G(n,n) = 1;
        elseif j==1
            G(n,n) = -sigma(i+1,j)-sigma(i-1,j)-sigma(i,j+1);
            G(n,n+1) = sigma(i+1,j);
            G(n,n-1) = sigma(i-1,j);
            G(n,n+ny) = sigma(i,j+1);
        elseif j==ny
            G(n,n) = -sigma(i+1,j)-sigma(i-1,j)-sigma(i,j-1);
            G(n,n+1) = sigma(i+1,j);
            G(n,n-1) = sigma(i-1,j);
            G(n,n-ny) = sigma(i,j-1);
        else
            G(n,n) = -sigma(i+1,j)-sigma(i-1,j)-sigma(i,j+1)-sigma(i,j-1);
            G(n,n+1) = sigma(i+1,j);
            G(n,n-1) = sigma(i-1,j);
            G(n,n+ny) = sigma(i,j+1);
            G(n,n-ny) = sigma(i,j-1);
        end
    end
end

figure
title('G matrix')
spy(G)

V = G\B';

Vmap = double(zeros(nx,ny));
for i = 1:nx
    for j = 1:ny
        n = i + (j - 1) * ny;
        Vmap(i, j) = V(n);
    end
end

%figure('units','normalized','outerposition',[0 0 1 1])
figure
surf(Vmap)
axis 'tight'
view(225,45)
title('Voltage in two contact region')
xlabel('X (um)')
ylabel('Y (um)')
zlabel('Voltage (V)')

%electric field
[Ex, Ey] = gradient(Vmap);
E = (Ex.^2+Ey.^2).^(0.5);

figure
quiver(-Ex,-Ey)
axis tight
title('Electric Field in a two contact region')
xlabel('X (um)')
ylabel('Y (um)')