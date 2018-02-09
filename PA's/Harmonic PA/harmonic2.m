nx = 50;
ny = nx;
soln = 9;

G = sparse(nx*ny);

for i = 1:nx
    for j = 1:ny
        n = i + (j - 1) * ny;
        if i==1||i==nx||j==1||j==ny
            G(n,:) = 0;
            G(n,n) = 1;           
        else
            G(n,n) = -4;
            G(n,n+1) = 1;
            G(n,n-1) = 1;
            G(n,n+ny) = 1;
            G(n,n-ny) = 1;   
        end
%         if i>10||i<20||j<10||j>20
%             G(n,n) = -2;
%         end
    end
end

figure
spy(G)

[E,D] = eigs(G,soln,'SM');
[~,Pr] = sort(diag(D),'descend');
D = D(Pr,Pr);
E = E(:,Pr);
    
Emap = double(zeros(nx,ny,soln));
for z = 1:soln
    for i = 1:nx
        for j = 1:ny
            n = i + (j - 1) * ny;
            Emap(i, j, z) = E(n,z);
        end
    end
end

figure('units','normalized','outerposition',[0 0 1 1])
for p = 1:soln
    subplot(3,3,p)
    surf(Emap(:,:,p))
end