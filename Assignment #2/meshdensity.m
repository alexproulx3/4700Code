amax = 10;
curr = zeros(1,amax);
nxx = zeros(1,amax);

for a = 1:amax
    nx = 10*a;
    ny = nx*3/2;
    sigma0 = 1;
    sigma1 = 1e-2;
    dim = 1e-6;
    Box = [nx*0.4 nx*0.6 ny*0.4 ny*0.6];
    
    G = sparse(nx+(ny-1)*ny);
    B = zeros(1,nx+(ny-1)*ny);
    
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
                if i>Box(1)&&i<Box(2)
                    G(n,n) = -3;
                    G(n,n+1) = sigma1;
                    G(n,n-1) = sigma1;
                    G(n,n+ny) = sigma1;
                else
                    G(n,n) = -3;
                    G(n,n+1) = sigma0;
                    G(n,n-1) = sigma0;
                    G(n,n+ny) = sigma0;
                end
            elseif j==ny
                if i>Box(1)&&i<Box(2)
                    G(n,n) = -3;
                    G(n,n+1) = sigma1;
                    G(n,n-1) = sigma1;
                    G(n,n-ny) = sigma1;
                else
                    G(n,n) = -3;
                    G(n,n+1) = sigma0;
                    G(n,n-1) = sigma0;
                    G(n,n-ny) = sigma0;
                end
            else
                if i>Box(1)&&i<Box(2)&&(j<Box(3)||j>Box(4))
                    G(n,n) = -4;
                    G(n,n+1) = sigma1;
                    G(n,n-1) = sigma1;
                    G(n,n+ny) = sigma1;
                    G(n,n-ny) = sigma1;
                else
                    G(n,n) = -4;
                    G(n,n+1) = sigma0;
                    G(n,n-1) = sigma0;
                    G(n,n+ny) = sigma0;
                    G(n,n-ny) = sigma0;
                end
            end
        end
    end
    
    V = G\B';
    
    Vmap = double(zeros(nx,ny));
    for i = 1:nx
        for j = 1:ny
            n = i + (j - 1) * ny;
            Vmap(i, j) = V(n);
        end
    end
    
    %sigma
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
    
    %electric field
    [Ex, Ey] = gradient(Vmap);
    E = sqrt(Ex.^2+Ey.^2);
    
    %current density
    J = sigma.*E;
    
    I = 0;
    for p = Box(1):Box(2)
        for q = Box(3):Box(4)
            I = I + J(p,q);
        end
    end
    I=I*(Box(2)-Box(1))*(Box(4)-Box(3))*dim^2;
    curr(a) = I;
    nxx(a) = nx;
end

figure
plot(nxx,curr)
axis 'tight'
title('Current vs Mesh Size')
xlabel('X dimension size (um)')
ylabel('Current (A)')