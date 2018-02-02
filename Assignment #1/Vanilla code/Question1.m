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
C.g = 9.80665;                      % metres (32.1740 ft) per s²
C.am = 1.66053892e-27;              % atomic mass

%Defining variables
Temp = 300; %K
nPart = 10; %Number of particles
dt = 1e-5;  %nx/vth/100%
tStop = 10000 * dt;
tmn = 0.2e-12;
Pscat = 0.005;
goPlot = 1;
nx = 200;
ny = 100;

%Assign position
x(1:nPart) = nx*rand(1,nPart);
y(1:nPart) = ny*rand(1,nPart);

%Assign velocity
std0 = sqrt(2 * C.kb * Temp / C.m_0); %Thermal Velocity
Vx(1:nPart) = std0 * randn(1, nPart);
Vy(1:nPart) = std0 * randn(1, nPart);

%Calculate mean free path
mfp = tmn .* sqrt(Vx.^2+Vy.^2);

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