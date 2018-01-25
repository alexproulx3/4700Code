function ToRunMC()
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
    dt = 1e-6;
    tStop = 10000 * dt; 

    %Define Box
    nx = 1000;
    ny = nx;
    Box = zeros(nx,ny);

    %Assign position
    x(1:nPart) = nx*rand(1,nPart);
    y(1:nPart) = ny*rand(1,nPart);

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
    Ix = bool(zeros(1,nPart));
    Iy = bool(zeros(1,nPart));
    
    %Main loop
    while t < tStop

        %Changing particle positions
        if Ix == 0
            x = xp + dt * Vx;    
        else
            x = xp - dt * Vx;
        end
        
        if Iy == 0 
            y = yp + dt * Vy;
        else
            y = yp - dt * Vy;
        end
        
        if x > nx
            Vx = -Vx;
        elseif x < 0
            x = xp - 2 * dt * Vx;
        end
        
        if y > ny
            y = yp - 2 * dt * Vy;
        elseif y < 0
            y = yp - 2 * dt * Vy;
        end

        %Iterating loop parameters
        xp = x;
        yp = y;
        c = c + 1;
        t  = t + dt;
        time(c) = t;
     
        %Plot
        plot(x,y,'*')
        
        %Heartbeat
        fprintf('time: %g (%5.3g %%)\n', t, t / tStop * 100);
        pause(0.00001)
    end
end


