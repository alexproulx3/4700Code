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
nPart = 30000; %Number of particles
dt = 5e-15; %nx/vth/100%
Iter = 1000;
tStop = Iter * dt;
tmn = 0.2e-12;
Pscat = 1-exp(-dt/tmn);
goPlot = 0;
nx = 200e-9;
ny = 100e-9;
V = 1;

%Assign position
x(1:nPart) = nx*rand(1,nPart);
y(1:nPart) = ny*rand(1,nPart);

%Assign velocity
std0 = sqrt(2 * C.kb * Temp / C.m_0); %Thermal Velocity
Vx(1:nPart) = std0 * randn(1, nPart);
Vy(1:nPart) = std0 * randn(1, nPart);

%Electric potential acceleration
Vxv = dt*C.q_0*V/(C.m_0*nx);
    
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
    title('2-D Particle Projection (1000 iterations)')
    xlim([0 nx])
    ylim([0 ny])
    xlabel('X (m)')
    ylabel('Y (m)')
    col = hsv(nPart);
end

%Main loop
while t < tStop
    %Changing particle positions 
    Vx = Vx + Vxv;
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
ylabel('Temperature (�K)')
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
map = zeros(100,200);
X = ceil(200*x/nx);
Y = ceil(100*y/ny);

for p = 1:nPart
    if Y(p) == 0
        Y(p) = 1;
    end
    if X(p) == 0
        X(p) = 1;
    end
    map(Y(p),X(p)) = map(Y(p),X(p))+1;
end

figure
h1 = mesh(map);
%h1.EdgeColor = 'none';
hold on
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
%h2.EdgeColor = 'none';
hold on
title('Temperature Map')
xlabel('X (nm)')
ylabel('Y (nm)')
zlabel('Temperature (�K)')
hold off

%Output minor calculations
format long
fprintf('Electric field strength is %f N/C\n ',V/nx)
fprintf('Electron Force is %f\n N',C.q_0*V/nx)
fprintf('Electron Acceleration is %f\n m/s^2',C.q_0*V/(C.m_0*nx))
fprintf('Mean free path was calculated to be %f m\n',sum(mfp)/c)
fprintf('Mean collision time is %e s\n',nx/(sum(vth)/c))