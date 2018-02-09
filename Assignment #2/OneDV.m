nx = 50;
ny = nx*3/2;
sigma = 1e-3;

G = sparse(nx*ny);
B = zeros(1,nx*ny);

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
            G(n,n) = -3;
            G(n,n+1) = 1;
            G(n,n-1) = 1;
            G(n,n+ny) = 1;
        elseif j==ny
            G(n,n) = -3;
            G(n,n+1) = 1;
            G(n,n-1) = 1;
            G(n,n-ny) = 1;
        else
            G(n,n) = -4;
            G(n,n+1) = 1;
            G(n,n-1) = 1;
            G(n,n+ny) = 1;
            G(n,n-ny) = 1;
        end 
    end
end

figure
spy(G)

E = G\B';
    
Emap = double(zeros(nx,ny));
for i = 1:nx
    for j = 1:ny
        n = i + (j - 1) * ny;
        Emap(i, j) = E(n);
    end
end

figure('units','normalized','outerposition',[0 0 1 1])
%figure
surf(Emap)
axis 'tight'