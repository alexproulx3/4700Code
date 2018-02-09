function [] = oneDelec(S,P)

% Initialize constants
global C
C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light
C.g = 9.80665;                      % metres (32.1740 ft) per s^2

nSteps = S;
nPart = P;

x(nPart,1:nSteps) = 0;
V(nPart,1:nSteps) = 0;
dF = 1e-30;
dt = 1e-6;

subplot(3,1,1)
hold on
xlabel('t');
ylabel('x');
title('Position')
subplot(3,1,2)
hold on
xlabel('t');
ylabel('V');
title('Velocity')
subplot(3,1,3);
hold on
xlabel('t');
ylabel('F');
title('Force')
grid on

for i = 1:nSteps
    t = i;
    r = rand;
    if t == 1
        V(:,t) = dF/C.m_0;
        dx = V(:,t) * dt;
        x(:,t) = dx;
    elseif r < 0.05
        % Scatter
        V(:,t) = 0;
        dx = V(:,t) * dt;
        x(:,t) = x(:,t-1) + dx;
    else
        V(:,t) = V(:,t-1) + dF/C.m_0;
        dx = V(:,t) * dt;
        x(:,t) = x(:,t-1) + dx;
    end

        subplot(3,1,1)
        plot(t,x(:,t),'-xk');
        hold on
        subplot(3,1,2)
        plot(t,V(:,t),'-or');
        hold on
        subplot(3,1,3);
        plot(t,dF,'og');
        hold on
        grid on

    pause(0.0001)
    
end
hold off
end