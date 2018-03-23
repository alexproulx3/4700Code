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
Cn = 10e-6;

%Inductors
L = 0.2;

%Voltages
Vin =  1;
Vinit = [0;
         0;
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

%Currents
In = 0.001*randn([1,length(t)]);

%Set up equation
%V = [V1 V2 V3 V4 V0 IL I3 In];
G = [1 0 0 0 0 0 0 0;
    -G2 G1+G2 0 0 0 -1 0 0;
    0 0 G3 0 0 -1 0 1;
    0 0 G3 0 0 0 -1 0;
    0 0 0 -G4 G4+G0 0 0 0;
    0 1 -1 0 0 0 0 0;
    0 0 0 1 0 0 -alpha 0
    0 0 0 0 0 0 0 1];

C = [0 0 0 0 0 0 0 0;
    -C1 C1 0 0 0 0 0 0;
    0 0 Cn 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 0 0 0;
    0 0 0 0 0 -L 0 0;
    0 0 0 0 0 0 0 0
    0 0 0 0 0 0 0 0];

%V = [V1; V2; V3; V4; V0; IL; I3; In];

F = zeros(length(G),length(In));
for a = 1:length(In)
    F(:,a) = [Vin;
        0;
        0;
        0;
        0;
        0;
        0;
        In(a)];
end

%Gaussian Function
Vgauss = exp(-(t-0.06).^2/(2*0.03^2));
FTgauss = fft(Vgauss);
FTgauss = fftshift(FTgauss);
n = length(Vgauss);
freq = (-n/2:n/2-1)*(f/n);

Fgauss = zeros(length(G),length(In));
for a = 1:length(Vgauss)
    Fgauss(:,a) = [Vgauss(a);
                        0;
                        0;
                        0;
                        0;
                        0;
                        0;
                        In(a)];
end

%Transient Gaussian Calculation
V = zeros(length(G),length(In));
for i = 1:1:length(In)
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
