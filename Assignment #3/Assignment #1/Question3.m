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
dt = 5e-15; %nx/vth/100%
Iter = 1000;
tStop = Iter * dt;
tmn = 0.2e-12;
Pscat = 1-exp(-dt/tmn);
goPlot = 0;
nx = 200e-9;
ny = 100e-9;

%Define Box
Box = [nx*0.4 nx*0.6 ny*0.4 ny*0.6];

%Assign position
while true
    x(1:nPart) = nx*rand(1,nPart);
    y(1:nPart) = ny*rand(1,nPart);
    xnot = (x(:) > Box(1)) & (x(:) < Box(2));
    ynot = (y(:) > Box(3)) & (y(:) < Box(4));
    if xnot == 0
        if ynot == 0
            break;
        end
    end
end

%Assign velocity
std0 = sqrt(2 * C.kb * Temp / C.m_0); %Thermal Velocity
Vx(1:nPart) = std0 * randn(1, nPart);
Vy(1:nPart) = std0 * randn(1, nPart);

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
    %axis ([0 nx 0 ny])
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

%electron density map
map = zeros(10,20);
X = ceil(20*x/nx);
Y = ceil(10*y/ny);

for p = 1:nPart
    map(Y(p),X(p)) = map(Y(p),X(p))+1;
end

figure
surf(map)
hold on
title('Electron Density Map')
xlabel('X')
ylabel('Y')
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
surf(map)
hold on
title('Temperature Map')
xlabel('X')
ylabel('Y')
zlabel('Temperature (ºK)')
hold off