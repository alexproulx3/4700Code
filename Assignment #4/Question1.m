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
Vin = -10:1:10;
V1 = Vin;
V4 = alpha/R3;

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

F = zeros(length(G),length(Vin));
for a = 1:length(Vin)
    F(:,a) = [Vin(a);
        0;
        0;
        0;
        0;
        0;
        0];
end

%DC analysis
Vdc = zeros(length(G),length(Vin));
for a = 1:length(Vin)
    Vdc(:,a) = G\F(:,a);
end

figure
plot(Vin,Vdc(5,:))
hold on
title('DC Sweep of Output Voltage')
xlabel('Input Voltage Vin (V)')
ylabel('Output Voltage Vo (V)')
hold off

figure
plot(Vin,Vdc(3,:))
hold on
title('DC Sweep of Voltage V3')
xlabel('Input Voltage Vin (V)')
ylabel('Output Voltage V3 (V)')
hold off

%AC analysis
%omega = [1 1e1 1e2 1e3 1e4 1e5 1e6 1e7 1e8 1e9 1e10 1e11 1e12];
omega = logspace(0,3,1000);

Vf = zeros(length(G),length(omega));
for b = 1:length(omega)
    Vf(:,b) = (G+1i*omega(b)*C)\F(:,12);
end

%plot Vo
figure
semilogx(omega,abs(Vf(5,:)))
hold on
title('AC Analysis for Output Voltage')
xlabel('Input Frequency (rad/s)')
ylabel('Output Voltage Vo (V)')
hold off

%plot Gain
figure
semilogx(omega,20*log(real(Vf(5,:)./Vf(1,:))))
hold on
title('AC Analysis for Gain')
xlabel('Input Frequency (rad/s)')
ylabel('Gain (dB)')
hold off

%Random Cap redefinition
C2 =  C1 + 0.05*randn([1 100]);
for c = 1:length(C2)
    Crand(:,:,c) = [0 0 0 0 0 0 0;
                    -C2(c) C2(c) 0 0 0 0 0;
                      0 0 -L 0 0 0 0;
                      0 0 0 0 0 0 0;
                      0 0 0 0 0 0 0;
                      0 0 0 0 0 0 0;
                      0 0 0 0 0 0 0];
end

for c = 1:length(C2)
    Vc(:,c) = (G+1i*pi*Crand(:,:,c))\F(:,12);
end

%plot gain histogram
figure
histogram(20*log((real(Vc(5,:)./Vc(1,:)))))
hold on
title('AC Analysis Histogram for Gain with Varying C')
xlabel('Gain (dB)')
ylabel('Bins')
hold off