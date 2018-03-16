close all;
clear;

global C
%Defining constants
C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665;                      % metres (32.1740 ft) per s^2
C.am = 1.66053892e-27;              % atomic mass

%Defining variables
Temp = 300; %K
nPart = 10; %Number of particles
dt = 1e-15; %nx/vth/100%
Iter = 1000;
tStop = Iter * dt;
tmn = 0.2e-12;
Pscat = 1-exp(-dt/tmn);
goPlot = 1;
nx = 100;
ny = nx*3/2;
sigma0 = 1;
sigma1 = 1e-2;
dim = 1e-9;
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

V = G\B';

Vmap = double(zeros(nx,ny));
for i = 1:nx
    for j = 1:ny
        n = i + (j - 1) * ny;
        Vmap(i, j) = V(n);
    end
end

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

%Redefine position
nx = 100*dim;
ny = nx*3/2;
Box = Box*dim;

%Assign position
fprintf('Assigning particle position \n')
x(1:nPart) = nx*rand(1,nPart);
y(1:nPart) = ny*rand(1,nPart);
xnot(1:nPart) = (x(1:nPart) > Box(1)) & (x(1:nPart) < Box(2));
ynot(1:nPart) = (y(1:nPart) < Box(3)) | (y(1:nPart) > Box(4));
flag = xnot & ynot;
while sum(flag) ~= 0
    x(flag) = nx*rand(sum(flag),1);
    y(flag) = ny*rand(sum(flag),1);
    
    xnot(1:nPart) = (x(1:nPart) > Box(1)) & (x(1:nPart) < Box(2));
    ynot(1:nPart) = (y(1:nPart) < Box(3)) | (y(1:nPart) > Box(4));

    flag = xnot & ynot;
    %plot(x,y,'*');
end
fprintf('Particle position Complete\n')

%Assign velocity
std0 = sqrt(2 * C.kb * Temp / C.m_0); %Thermal Velocity
Vx(1:nPart) = std0 * randn(1, nPart);
Vy(1:nPart) = std0 * randn(1, nPart);

%Electric potential acceleration
Vxv = (Vmap - (C.q_0/C.m_0)*Ey*dt);
Vyv = (Vmap - (C.q_0/C.m_0)*Ex*dt);
    
%Calculate mean free path
mfp = zeros(0,Iter);

%Calculating thermal velocity
vth = zeros(0,Iter);

%Calculating temperature
T = zeros(0,Iter);

%Calculating Current
It = zeros(0,Iter);

%Initialize loop variables
xp(1:nPart) = x - dt * Vx;
yp(1:nPart) = y - dt * Vy;
t = 0;
c = 1;
time(c) = 0;

%indexing
Ix = ones(1,nPart);
Iy = ones(1,nPart);

if goPlot
    %set up figure
    figure('units','normalized','outerposition',[0 0 1 1])
    hold on
    xlim([0 nx])
    ylim([0 ny])
    line([Box(1), Box(2), Box(2), Box(1), Box(1)], [ny, ny, Box(4), Box(4), ny]);
    line([Box(1), Box(2), Box(2), Box(1), Box(1)], [0, 0, Box(3), Box(3), 0]);
    
    col = hsv(nPart);
end

%Main loop
while t < tStop
    %Changing particle positions
    for q = 1:nPart
        xq = ceil(x(q)/dim);
        yq = ceil(y(q)/dim);
        if xq <= 0
            xq = 1;
        end
        if xq >= nx/dim+1
            xq = nx/dim;
        end
        if yq <= 0
            yq = 1;
        end
        if yq >= ny/dim+1
            yq = ny/dim;
        end
        Vx = Vx + Vxv(xq,yq);
        Vy = Vy + Vyv(xq,yq);
    end
    
    x = xp + dt * Ix .* Vx; 
    y = yp + dt * Iy .* Vy; 
    
    % x-y boundary conditions
    i = find(x>nx);
    j = find(x<0);
    k = find(y>ny);
    l = find(y<0);
    if ~isempty(i)
        xp(1,i) = x(1,i) - nx;
        x(1,i) = 0;
    end
    if ~isempty(j)
        xp(1,j) = x(1,j) + nx;
        x(1,j) = nx;
    end
    if ~isempty(k)
        Iy(1,k) = -1*Iy(1,k);
        y(1,k) = yp(1,k) + 2 * dt * Iy(1,k) .* Vy(1,k);
    end
    if ~isempty(l)
        Iy(1,l) = -1*Iy(1,l);
        y(1,l) = yp(1,l) + 2 * dt * Iy(1,l) .* Vy(1,l);
    end
    
    % Box boundaries
    n = find(x>Box(1) & x<Box(2) & y>0 & y<Box(3));
    m = find(x>Box(1) & x<Box(2) & y>Box(4) & y<ny);
    if ~isempty(n)
        if sum((xp<Box(1)&x>Box(1))|(xp>Box(2)&x<Box(2)))
            Ix(1,n) = -1*Ix(1,n);
            x(1,n) = xp(1,n) + 2 * dt * Ix(1,n) .* Vx(1,n);
        end
        if sum((yp<Box(4)&y>Box(4))|(yp>Box(3)&y<Box(3)))
            Iy(1,n) = -1*Iy(1,n);
            y(1,n) = yp(1,n) + 2 * dt * Iy(1,n) .* Vy(1,n);
        end
    end
    if ~isempty(m)
        if sum((xp<Box(1)&x>Box(1))|(xp>Box(2)&x<Box(2)))
            Ix(1,m) = -1*Ix(1,m);
            x(1,m) = xp(1,m) + 2 * dt * Ix(1,m) .* Vx(1,m);
        end
        if sum((yp<Box(4)&y>Box(4))|(yp>Box(3)&y<Box(3)))
            Iy(1,m) = -1*Iy(1,m);
            y(1,m) = yp(1,m) + 2 * dt * Iy(1,m) .* Vy(1,m);
        end
    end
    
    %scattering
    rando = Pscat > rand(1,nPart);
    Vx(rando) = std0 * randn(sum(rando),1);
    Vy(rando) = std0 * randn(sum(rando),1);
    
    %Plot
    if goPlot
        for p = 1:nPart
            plot([x(p) xp(p)],[y(p) yp(p)],'color',col(p,:))
        end
    end
    
    %Calculating Temperature
    vth(1,c) = sum(Vx.^2+Vy.^2)/nPart;
    T(1,c) = (vth(1,c)*C.m_0)/(2*C.kb);
    
    %Calculate mean free path
    mfp(1,c) = tmn*vth(1,c);
    
    %Calculate current
    It(1,c) = vth(1,c)*nx*ny*C.q_0*10e15/1e-4;
    
    %Iterating loop parameters
    xp = x;
    yp = y;
    c = c + 1;
    t  = t + dt;
    time(c) = t;
    
    %Heartbeat
    fprintf('time: %g (%5.3g %%)\n', t, t / tStop * 100);
    pause(0.0001)
end
hold off

%Thermal velocity histogram
figure
hist(vth)
hold on
title('Thermal Velocity over simluation (1000 iterations)')
xlabel('Thermal Velocity (m/s)')
ylabel('Counts')
hold off

%Temperature plot
figure
plot(dt:dt:tStop,T)
hold on
title('System Temperature over simluation (1000 iterations)')
xlabel('Time (s)')
ylabel('Temperature (ºK)')
hold off

%Current plot
figure
plot(dt:dt:tStop,It)
hold on
title('System Current over simluation (1000 iterations)')
xlabel('Time (s)')
ylabel('Current (A)')
hold off

%electron density map
map = zeros(nx/dim,ny/dim);
X = ceil(nx/dim*x/nx);
Y = ceil(ny/dim*y/ny);

for p = 1:nPart
    if X(p) <= 0
        X(p) = 1;
    end
    if X(p) >= nx/dim+1
        X(p) = nx/dim;
    end
    if Y(p) <= 0
        Y(p) = 1;
    end
    if Y(p) >= ny/dim+1
        Y(p) = ny/dim;
    end
    map(X(p),Y(p)) = map(X(p),Y(p))+1;
end

figure
h1 = mesh(map);
%h1.EdgeColor = 'none';
hold on
axis tight
title('Electron Density Map')
xlabel('X (nm)')
ylabel('Y (nm)')
zlabel('Number of electrons')
hold off

%temperature map
vth = zeros(1,nPart);
T = vth;
for p = 1:nPart
    vth(p) = Vx(p).^2+Vy(p).^2;
    T(p) = (vth(p)*C.m_0)/(2*C.kb);
    map(Y(p),X(p)) = T(p);
end

figure
h2 = mesh(map);
axis tight
%h2.EdgeColor = 'none';
hold on
title('Temperature Map')
xlabel('X (nm)')
ylabel('Y (nm)')
zlabel('Temperature (ºK)')
hold off

%Output minor calculations
format long
%fprintf('Electric field strength is %f N/C\n ',V/nx)
%fprintf('Electron Force is %f N\n',C.q_0*V/nx)
%fprintf('Electron Acceleration is %f m/s^2\n',C.q_0*V/(C.m_0*nx))
fprintf('Mean free path was calculated to be %f m\n',sum(mfp)/c)
fprintf('Mean collision time is %e s\n',nx/(sum(vth)/c))