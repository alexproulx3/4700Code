%function [] = Iter1()

% Initial variables
nx = 100;
ny = nx;
iter = 1000;

%V = zeros(nx,ny);
V = rand(nx);
%E = zeros(nx,ny);
q = 0;%1.602e-19;
k = 0;%1e-12;
dx = 0;%1e-3;
[x,y] = meshgrid(1:nx,1:ny);

%Plot settings
figure('units','normalized','outerposition',[0 0 1 1])
subplot(1,2,1)
title('V')
subplot(1,2,2)
title('E')
axis([0 nx 0 ny])

for n=1:iter
    for i = 1:nx
        for j = 1:ny
%             if (i==1 && j==1)||(i==1 && j==ny)||(i==nx && j==1)||(i==nx && j==ny)
%                 V(i,j) = 0;
            if i == 1
                V(i,j) = 0;
            elseif i == nx
                V(i,j) = 1;
            elseif j == 1
                V(i,j)=(V(i+1,j)+V(i-1,j)+V(i,j+1))/3;
            elseif j == ny
                V(i,j)=(V(i+1,j)+V(i-1,j)+V(i,j-1))/3;
            else
                V(i,j)=(V(i+1,j)+V(i-1,j)+V(i,j+1)+V(i,j-1))/4;
            end
        end
    end

    [Ex,Ey] = gradient(V);
    subplot(1,2,1)
    surf(V)
    subplot(1,2,2)
    quiver(x,y,-Ex,-Ey)
    pause(0.0001)
end


%end

