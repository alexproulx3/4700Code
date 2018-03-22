%% Part 1
Is = 0.01e-12; % Forward Bias Saturation Current
Ib = 0.1e-12;  % Breakdown Bias Saturation Current
Vb = 1.3;      % Breakdown Voltage
Gp = 0.1;      % Parasitic Parallel Conductance
Vt = 1.2/0.025;% Thermal Voltage

% Voltage
V = linspace(-1.95,0.7,200);
% Current
I = Is*(exp(Vt.*V)-1)+Gp*V-Ib*exp(-Vt.*(V+Vb));
Irand = I*0.8+I.*0.2.*rand(size(V));

% Plots
figure(1)
plot(V,I,'.')
title('Diode Voltage vs Current')
xlabel('Current (A)')
ylabel('Voltage (V)')

figure(2)
semilogy(V,abs(I),'.')
title('Diode Voltage vs log(Current)')
xlabel('Current (A)')
ylabel('Voltage (V)')

figure(3)
plot(V,Irand,'.')
title('Diode Voltage vs randCurrent')
xlabel('Current (A)')
ylabel('Voltage (V)')

figure(4)
semilogy(V,abs(Irand),'.')
title('Diode Voltage vs log(randCurrent)')
xlabel('Current (A)')
ylabel('Voltage (V)')

%% Part 2
% Polynomial fit
p4I = polyfit(V,I,4);
Ip4 = polyval(p4I,V);
p8I = polyfit(V,I,8);
Ip8 = polyval(p8I,V);
p4Irand= polyfit(V,I,4);
Irandp4 = polyval(p4Irand,V);
p8Irand= polyfit(V,I,8);
Irandp8 = polyval(p8Irand,V);

% Adding fits to graphs
figure(1)
hold on
plot(V,Ip4)
plot(V,Ip8)
hold off

figure(3)
hold on
plot(V,Irandp4)
plot(V,Irandp8)
hold off

%% Part 3
% Setting non linear fit parameters
fo1 = fittype('A.*(exp(1.2*x/25e-3)-1) + 0.1.*x - C*(exp(1.2*(-(x+1.3))/25e-3))');
ff1 = fit(V',I',fo1);
If1 = ff1(V);

fo2 = fittype('A.*(exp(1.2*x/25e-3)-1) + B.*x - C*(exp(1.2*(-(x+1.3))/25e-3))');
ff2 = fit(V',I',fo2);
If2 = ff2(V);

fo3 = fittype('A.*(exp(1.2*x/25e-3)-1) + B.*x - C*(exp(1.2*(-(x+D))/25e-3))');
ff3 = fit(V',I',fo3);
If3 = ff3(V);

% Adding to new plot
figure(5)
hold on
plot(V,I,'.')
plot(V,If1)
plot(V,If2)
plot(V,If3)
legend('Data','fita','fitb','fitc');
xlim([-2 1])
ylim([-4 4])
title('Non Linear Fit Diode Voltage vs Current')
xlabel('Current (A)')
ylabel('Voltage (V)')
hold off

%% Part 4
% Setting up neural net
inputs = V.';
targets = I.';
hiddenLayerSize = 10;
net = fitnet(hiddenLayerSize);
net.divideParam.trainRatio = 70/100;
net.divideParam.valRatio = 15/100;
net.divideParam.testRatio = 15/100;
[net,tr] = train(net,inputs,targets);
outputs = net(inputs);
errors = gsubtract(outputs,targets);
performance = perform(net,targets,outputs);
%view(net)
Inn = outputs;

% Plotting neural net
figure(6)
hold on
plot(V,I,'.')
plot(V,Inn)
xlim([-2 1])
ylim([-4 4])
title('Neural Net Fit Diode Voltage vs Current')
xlabel('Current (A)')
ylabel('Voltage (V)')
hold off