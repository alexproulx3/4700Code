%Resistors
R1 = 1;
R2 = 2;
R3 = 10;
R4 = 0.1;
R0 = 1000;
alpha = 100;

G1 = 1/R1;
G2 = 1/R2;
G3 = 1/R3;
G4 = 1/R4;
G0 = 1/R0;

%Capacitors
C1 = 0.25;

%Inductors
L = 0.2;

%Voltages
Vin = 1;
Vinit = [0;
         0;
         0;
         0;
         0;
         0;
         0];
V1 = Vin;
V4 = alpha/R3;

%Transient stuff
dt = 1e-3;
t =  dt:dt:1;

%Set up equation
%V = [V1 V2 V3 V4 V0 IL I3];
G = [1 0 0 0 0 0 0;
    -G2 G1+G2 0 0 0 -1 0;
    0 0 G3 0 0 -1 0;
    0 0 G3 0 0 0 -1;
    0 0 0 -G4 G4+G0 0 0;
    0 1 -1 0 0 0 0;
    0 0 0 1 0 0 -alpha];

C = [0 0 0 0 0 0 0;
    -C1 C1 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 0 0;
    0 0 0 0 0 -L 0;
    0 0 0 0 0 0 0];

%V = [V1; V2; V3; V4; V0; IL; I3];

Foff = [0;
       0;
       0;
       0;
       0;
       0;
       0];
   
Fon = [Vin;
       0;
       0;
       0;
       0;
       0;
       0];

Fstep = ones(1,1000);
Fstep(1:30) = 0;
FTstep = fft(Fstep);
FTstep = fftshift(FTstep);
n = length(Fstep);
freq = (-n/2:n/2-1)*(f/n);
   
%Transient Step Calculation
V = zeros(7,1000);
for i = 1:1:1000
    if i == 1
        V(:,i) = (C./dt+G)\(Foff+C*Vinit/dt);
    elseif i <= 30
        V(:,i) = (C./dt+G)\(Foff+C*Vb4/dt);
    else
        V(:,i) = (C./dt+G)\(Fon+C*Vb4/dt);
    end
    Vb4 = V(:,i);
end

Vft = fft(V(5,:));
Vft = fftshift(Vft);

%Plotting Step Transient
figure
hold on
plot(t,V(5,:))
plot(t,V(1,:))
legend('Output','Input')
grid on
title('Step Transient')
xlabel('Time (s)')
ylabel('Output Voltage (V)')
hold off

%Fourier for Step
figure
hold on
plot(freq,abs(Vft).^2/n)
plot(freq,abs(FTstep).^2/n)
legend('Output','Input')
grid on
title('Step Fourier Transform')
xlabel('Frequency (Hz)')
%xlim([-20 20])
ylabel('Output Voltage (V)')
hold off

%Sine wave
f = 1/0.03;
Vsine = sin(2*pi*f*t);
FTsine = fft(Vsine);
FTsine = fftshift(FTsine);
n = length(Vsine);
freq = (-n/2:n/2-1)*(f/n);

Fsine = zeros(length(G),length(Vsine));
for a = 1:length(Vsine)
    Fsine(:,a) = [Vsine(a);
                    0;
                    0;
                    0;
                    0;
                    0;
                    0];
end

%Transient Sine Calculation
V = zeros(7,1000);
for i = 1:1:1000
    if i == 1
        V(:,i) = (C./dt+G)\(Fsine(:,i)+C*Vinit/dt);
    else
        V(:,i) = (C./dt+G)\(Fsine(:,i)+C*Vb4/dt);
    end
    Vb4 = V(:,i);
end

Vft = fft(V(5,:));
Vft = fftshift(Vft);

%Plotting Sine
figure
hold on
plot(t,V(5,:))
plot(t,V(1,:))
legend('Output','Input')
grid on
title('Sine Transient')
xlabel('Time (s)')
ylabel('Output Voltage (V)')
hold off

%Plotting Fourier for Sine
figure
hold on
plot(freq,abs(Vft).^2/n)
plot(freq,abs(FTsine).^2/n)
legend('Output','Input')
grid on
title('Sine Fourier Transform')
xlabel('Frequency (Hz)')
xlim([-20 20])
ylabel('Output Voltage (V)')
hold off

%Gaussian Function
Vgauss = exp(-(t-0.06).^2/(2*0.03^2));
FTgauss = fft(Vgauss);
FTgauss = fftshift(FTgauss);
n = length(Vgauss);
freq = (-n/2:n/2-1)*(f/n);

Fgauss = zeros(length(G),length(Vgauss));
for a = 1:length(Vgauss)
    Fgauss(:,a) = [Vgauss(a);
                        0;
                        0;
                        0;
                        0;
                        0;
                        0];
end

%Transient Gaussian Calculation
V = zeros(7,1000);
for i = 1:1:1000
    if i == 1
        V(:,i) = (C./dt+G)\(Fgauss(:,i)+C*Vinit/dt);
    else
        V(:,i) = (C./dt+G)\(Fgauss(:,i)+C*Vb4/dt);
    end
    Vb4 = V(:,i);
end

Vft = fft(V(5,:));
Vft = fftshift(Vft);

%Plotting Gaussian
figure
hold on
plot(t,V(5,:))
plot(t,V(1,:))
legend('Output','Input')
grid on
title('Gaussian Transient')
xlabel('Time (s)')
ylabel('Output Voltage (V)')
hold off


%Plotting Fourier for Gaussian
figure
hold on
plot(freq,abs(Vft).^2/n)
plot(freq,abs(FTgauss).^2/n)
legend('Output','Input')
grid on
title('Gaussian Fourier Transform')
xlabel('Frequency (Hz)')
xlim([-20 20])
ylabel('Output Voltage (V)')
hold off