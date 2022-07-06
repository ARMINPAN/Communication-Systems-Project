%% Communication Systems Project
% implementing a Digital Communication system
%% Armin Panjehpour - 98101288
%% 
clc; clear; close all;
% at first each block will be implemented seperately
%% part 2.1 - Divide/ Combine Blocks for transmitter/ reciever
% functions are implemented at the end of the code

%% part 2.2 - Pulse Shaping - create a pulse from a given sequence
% functions are implemented at the end of the code

%% part 2.3 - Pulse Shaping - Modulation
% functions are implemented at the end of the code

%% part 2.4 - Ideal Channel
% functions are implemented at the end of the code

%% part 2.5 - Demodulation
% functions are implemented at the end of the code

%% part 2.6 - Decode the signal/ Matched Filter
% functions are implemented at the end of the code
%%

%% part 3 - PAM modulation
% initializations 
clc; clear; close all;
Fs = 10^6; % Hz
tPulse = 10^-2; % in s
Fc = 10^4; % in Hz
FcChannel = 10^4; % in Hz
BWchannel = 10^3; % in Hz

% part 3.1
onePulse = ones(1,tPulse*Fs);
zeroPulse = -1*ones(1,tPulse*Fs);

% part 3.1.A - Transmitting the signal
% transmitter
lengthSignal = 10; % arbitrary | should be even
inputSignal = randi([0 1], 1, lengthSignal)
figure;
scatter(1:length(inputSignal),inputSignal,'filled');
grid on; grid minor;
title("Input Digital Signal",'interpreter','latex');
xlabel('numberOfBit','interpreter','latex');

% Devide
[out1, out2] = Devide(inputSignal);
figure;
subplot(1,2,1);
scatter(1:length(out1),out1,'filled');
grid on; grid minor;
title('$b_1[n]$','interpreter','latex');
xlabel('numberOfBit','interpreter','latex');
subplot(1,2,2);
scatter(1:length(out2),out2,'filled');
grid on; grid minor;
title('$b_2[n]$','interpreter','latex');
xlabel('numberOfBit','interpreter','latex');

% Pulse maker
x1 = PulseShaping(out1, zeroPulse, onePulse);
x2 = PulseShaping(out2, zeroPulse, onePulse);
figure;
subplot(1,2,1);
plot((0:1/Fs:length(x1)/Fs-1/Fs)*1000,x1,'LineWidth',2);
grid on; grid minor;
title('$x_1(t)$','interpreter','latex');
xlabel('time(ms)','interpreter','latex');
ylabel('amplitude','interpreter','latex');
xlim([0 50])

subplot(1,2,2);
plot((0:1/10^6:length(x2)/10^6-1/10^6)*1000,x2,'LineWidth',2);
grid on; grid minor;
title('$x_2(t)$','interpreter','latex');
xlabel('time(ms)','interpreter','latex');
ylabel('amplitude','interpreter','latex');
xlim([0 50])

% Modulation
xc = AnalogMod(x1,x2,Fs,Fc);
figure;
plot((0:1/Fs:length(xc)/Fs-1/Fs)*1000,xc)
grid on; grid minor;
title('$x_c(t)$','interpreter','latex');
xlabel('time(ms)','interpreter','latex');
ylabel('amplitude','interpreter','latex');

% receiver
% channel
y = Channel(xc,Fs,FcChannel,BWchannel);
figure
plot((0:1/Fs:length(y)/Fs-1/Fs)*1000,y)
grid on; grid minor;
title('y(t)','interpreter','latex');
xlabel('time(ms)','interpreter','latex');
ylabel('amplitude','interpreter','latex');

% Demodulation
[y1, y2] = AnalogDemod(y,Fs,BWchannel,Fc);
figure;
subplot(1,2,1);
plot((0:1/Fs:length(y1)/Fs-1/Fs)*1000,y1)
grid on; grid minor;
title('$y_1(t)$','interpreter','latex');
xlabel('time(ms)','interpreter','latex');
ylabel('amplitude','interpreter','latex');

subplot(1,2,2);
plot((0:1/Fs:length(y2)/Fs-1/Fs)*1000,y2)
grid on; grid minor;
title('$y_2(t)$','interpreter','latex');
xlabel('time(ms)','interpreter','latex');
ylabel('amplitude','interpreter','latex');
 
% MatchedFilt and threshold
[raw0 raw1 Reconstructed1] = MatchedFilt(y1, zeroPulse, onePulse);
[raw2 raw3 Reconstructed2] = MatchedFilt(y2, zeroPulse, onePulse);

figure;
subplot(1,2,1);
scatter(1:length(Reconstructed1),Reconstructed1,'filled');
grid on; grid minor;
title('$\hat{b_1}[n]$','interpreter','latex');
xlabel('numberOfBit','interpreter','latex');
subplot(1,2,2);
scatter(1:length(Reconstructed2),Reconstructed2,'filled');
grid on; grid minor;
title('$\hat{b_2}[n]$','interpreter','latex');
xlabel('numberOfBit','interpreter','latex');

% Combine
bprime = Combine(Reconstructed1,Reconstructed2)
ReconstructedSignal = bprime;
figure;
subplot(2,1,1)
scatter(1:length(inputSignal),inputSignal,'filled');
grid on; grid minor;
title("Input Digital Signal",'interpreter','latex');
xlabel('numberOfBit','interpreter','latex');
subplot(2,1,2);
scatter(1:length(ReconstructedSignal),ReconstructedSignal,'filled');
grid on; grid minor;
title("Reconstructed Signal",'interpreter','latex');
xlabel('numberOfBit','interpreter','latex');
%% part 3.1.B 
% initializations 
clc; clear; close all;
Fs = 10^6; % Hz
tPulse = 10^-2; % in s
Fc = 10^4; % in Hz
FcChannel = 10^4; % in Hz
BWchannel = 10^3; % in Hz

% part 3.1
onePulse = ones(1,tPulse*Fs);
zeroPulse = -1*ones(1,tPulse*Fs);

% part 3.1.A - Transmitting the signal
% transmitter
lengthSignal = 10; % arbitrary | should be even
inputSignal = randi([0 1], 1, lengthSignal);
%figure;
%scatter(1:length(inputSignal),inputSignal);

% Devide
[out1, out2] = Devide(inputSignal);

% Pulse maker
x1 = PulseShaping(out1, zeroPulse, onePulse);
x2 = PulseShaping(out2, zeroPulse, onePulse);

% Modulation
xc = AnalogMod(x1,x2,Fs,Fc);

% receiver
% channel
y = Channel(xc,Fs,FcChannel,BWchannel);
% noise
error = zeros(1,100);
for i=1:100
    for j=1:5
        
        [i j]
        mu = 0;
        sigma = i;
        noise = normrnd(mu,sigma,1,length(y));
        noisyY = y + noise;


        % Demodulation
        [y1, y2] = AnalogDemod(noisyY,Fs,BWchannel,Fc);

        % MatchedFilt
        [raw0 raw1 Reconstructed1] = MatchedFilt(y1, zeroPulse, onePulse);
        [raw2 raw3 Reconstructed2] = MatchedFilt(y2, zeroPulse, onePulse);

        % Combine
        bprime = Combine(Reconstructed1,Reconstructed2);
        ReconstructedSignal = bprime;
        
       % error calculation
       error(i) = error(i) + sum(abs(bprime-inputSignal));
    end
    error(i) = error(i)/5/10;
end

figure;
stem(1:length(error),error);
grid on; grid minor;
title("Probability of Error",'interpreter','latex');
xlabel('Standard Deviation of noise','interpreter','latex');
ylabel('Probability','interpreter','latex');
%% part 3.1.C
clc; clear; close all;
selectedVars = [10 28 32 40 73 99];
% initializations 
Fs = 10^6; % Hz
tPulse = 10^-2; % in s
Fc = 10^4; % in Hz
FcChannel = 10^4; % in Hz
BWchannel = 10^3; % in Hz

% part 3.1
onePulse = ones(1,tPulse*Fs);
zeroPulse = -1*ones(1,tPulse*Fs);

% part 3.1.A - Transmitting the signal
% transmitter
lengthSignal = 1000; % arbitrary | should be even - big for a beautiful plot
inputSignal = randi([0 1], 1, lengthSignal);


% Devide
[out1, out2] = Devide(inputSignal);

% Pulse maker
x1 = PulseShaping(out1, zeroPulse, onePulse);
x2 = PulseShaping(out2, zeroPulse, onePulse);

% Modulation
xc = AnalogMod(x1,x2,Fs,Fc);

% receiver
% channel
y = Channel(xc,Fs,FcChannel,BWchannel);


figure;
for i = 1:length(selectedVars)
    mu = 0;
    sigma = selectedVars(i);
    noise = normrnd(mu,sigma,1,length(y));
    noisyY = y + noise;


    % Demodulation
    [y1, y2] = AnalogDemod(noisyY,Fs,BWchannel,Fc);

    % MatchedFilt
    [raw0 raw1 Reconstructed1] = MatchedFilt(y1, zeroPulse, onePulse);
    [raw2 raw3 Reconstructed2] = MatchedFilt(y2, zeroPulse, onePulse);

    % Combine
    bprime = Combine(Reconstructed1,Reconstructed2);
    ReconstructedSignal = bprime;  
    
    subplot(2,3,i);
    hold on;
    scatter(raw1,raw3);
    grid on; grid minor;
    title("Signal Constellation for SD = " + selectedVars(i),'interpreter','latex');
    xlabel("$\hat{b_1}[i]$",'interpreter','latex')
    ylabel("$\hat{b_2}[i]$",'interpreter','latex')
end



%% part 3.2 - PSK modulation
% initializations 
clc; clear; close all;
Fs = 10^6; % Hz
tPulse = 10^-2; % in s
Fc = 10^4; % in Hz
FcChannel = 10^4; % in Hz
BWchannel = 10^3; % in Hz

% part 3.1
f = 500; % in Hz
onePulse = cos(2*pi*f*(0:1/Fs:tPulse-1/Fs));
zeroPulse = -cos(2*pi*f*(0:1/Fs:tPulse-1/Fs));

% part 3.1.A - Transmitting the signal
% transmitter
lengthSignal = 10; % arbitrary | should be even
inputSignal = randi([0 1], 1, lengthSignal)
figure;
scatter(1:length(inputSignal),inputSignal,'filled');
grid on; grid minor;
title("Input Digital Signal",'interpreter','latex');
xlabel('numberOfBit','interpreter','latex');

% Devide
[out1, out2] = Devide(inputSignal);
figure;
subplot(1,2,1);
scatter(1:length(out1),out1,'filled');
grid on; grid minor;
title('$b_1[n]$','interpreter','latex');
xlabel('numberOfBit','interpreter','latex');
subplot(1,2,2);
scatter(1:length(out2),out2,'filled');
grid on; grid minor;
title('$b_2[n]$','interpreter','latex');
xlabel('numberOfBit','interpreter','latex');

% Pulse maker
x1 = PulseShaping(out1, zeroPulse, onePulse);
x2 = PulseShaping(out2, zeroPulse, onePulse);
figure;
subplot(1,2,1);
plot((0:1/Fs:length(x1)/Fs-1/Fs)*1000,x1,'LineWidth',2);
grid on; grid minor;
title('$x_1(t)$','interpreter','latex');
xlabel('time(ms)','interpreter','latex');
ylabel('amplitude','interpreter','latex');
xlim([0 50])

subplot(1,2,2);
plot((0:1/10^6:length(x2)/10^6-1/10^6)*1000,x2,'LineWidth',2);
grid on; grid minor;
title('$x_2(t)$','interpreter','latex');
xlabel('time(ms)','interpreter','latex');
ylabel('amplitude','interpreter','latex');
xlim([0 50])

% Modulation
xc = AnalogMod(x1,x2,Fs,Fc);
figure;
plot((0:1/Fs:length(xc)/Fs-1/Fs)*1000,xc)
grid on; grid minor;
title('$x_c(t)$','interpreter','latex');
xlabel('time(ms)','interpreter','latex');
ylabel('amplitude','interpreter','latex');

% receiver
% channel
y = Channel(xc,Fs,FcChannel,BWchannel);
figure
plot((0:1/Fs:length(y)/Fs-1/Fs)*1000,y)
grid on; grid minor;
title('y(t)','interpreter','latex');
xlabel('time(ms)','interpreter','latex');
ylabel('amplitude','interpreter','latex');

% Demodulation
[y1, y2] = AnalogDemod(y,Fs,BWchannel,Fc);
figure;
subplot(1,2,1);
plot((0:1/Fs:length(y1)/Fs-1/Fs)*1000,y1)
grid on; grid minor;
title('$y_1(t)$','interpreter','latex');
xlabel('time(ms)','interpreter','latex');
ylabel('amplitude','interpreter','latex');
xlim([0 50]);

subplot(1,2,2);
plot((0:1/Fs:length(y2)/Fs-1/Fs)*1000,y2)
grid on; grid minor;
title('$y_2(t)$','interpreter','latex');
xlabel('time(ms)','interpreter','latex');
ylabel('amplitude','interpreter','latex');
xlim([0 50]);

% MatchedFilt and threshold
[raw0 raw1 Reconstructed1] = MatchedFilt(y1, zeroPulse, onePulse);
[raw2 raw3 Reconstructed2] = MatchedFilt(y2, zeroPulse, onePulse);

figure;
subplot(1,2,1);
scatter(1:length(Reconstructed1),Reconstructed1,'filled');
grid on; grid minor;
title('$\hat{b_1}[n]$','interpreter','latex');
xlabel('numberOfBit','interpreter','latex');
subplot(1,2,2);
scatter(1:length(Reconstructed2),Reconstructed2,'filled');
grid on; grid minor;
title('$\hat{b_2}[n]$','interpreter','latex');
xlabel('numberOfBit','interpreter','latex');

% Combine
bprime = Combine(Reconstructed1,Reconstructed2)
ReconstructedSignal = bprime;
figure;
subplot(2,1,1)
scatter(1:length(inputSignal),inputSignal,'filled');
grid on; grid minor;
title("Input Digital Signal",'interpreter','latex');
xlabel('numberOfBit','interpreter','latex');
subplot(2,1,2);
scatter(1:length(ReconstructedSignal),ReconstructedSignal,'filled');
grid on; grid minor;
title("Reconstructed Signal",'interpreter','latex');
xlabel('numberOfBit','interpreter','latex');
%% part 3.2.B
% initializations 
clc; clear; close all;
Fs = 10^6; % Hz
tPulse = 10^-2; % in s
Fc = 10^4; % in Hz
FcChannel = 10^4; % in Hz
BWchannel = 10^3; % in Hz

% part 3.1
f = 500; % in Hz
onePulse = cos(2*pi*f*(0:1/Fs:tPulse-1/Fs));
zeroPulse = -cos(2*pi*f*(0:1/Fs:tPulse-1/Fs));

% part 3.1.A - Transmitting the signal
% transmitter
lengthSignal = 10; % arbitrary | should be even
inputSignal = randi([0 1], 1, lengthSignal);
%figure;
%scatter(1:length(inputSignal),inputSignal);

% Devide
[out1, out2] = Devide(inputSignal);

% Pulse maker
x1 = PulseShaping(out1, zeroPulse, onePulse);
x2 = PulseShaping(out2, zeroPulse, onePulse);

% Modulation
xc = AnalogMod(x1,x2,Fs,Fc);

% receiver
% channel
y = Channel(xc,Fs,FcChannel,BWchannel);
% noise
error = zeros(1,100);
for i=1:100
    for j=1:5
        
        [i j]
        mu = 0;
        sigma = i;
        noise = normrnd(mu,sigma,1,length(y));
        noisyY = y + noise;


        % Demodulation
        [y1, y2] = AnalogDemod(noisyY,Fs,BWchannel,Fc);

        % MatchedFilt
        [raw0 raw1 Reconstructed1] = MatchedFilt(y1, zeroPulse, onePulse);
        [raw2 raw3 Reconstructed2] = MatchedFilt(y2, zeroPulse, onePulse);

        % Combine
        bprime = Combine(Reconstructed1,Reconstructed2);
        ReconstructedSignal = bprime;
        
       % error calculation
       error(i) = error(i) + sum(abs(bprime-inputSignal));
    end
    error(i) = error(i)/5/10;
end

figure;
stem(1:length(error),error);
grid on; grid minor;
title("Probability of Error",'interpreter','latex');
xlabel('Standard Deviation of noise','interpreter','latex');
ylabel('Probability','interpreter','latex');

%% part 3.2.C
clc; clear; close all;
selectedVars = [10 28 32 40 73 99];
% initializations 
Fs = 10^6; % Hz
tPulse = 10^-2; % in s
Fc = 10^4; % in Hz
FcChannel = 10^4; % in Hz
BWchannel = 10^3; % in Hz

% part 3.1
f = 500; % in Hz
onePulse = cos(2*pi*f*(0:1/Fs:tPulse-1/Fs));
zeroPulse = -cos(2*pi*f*(0:1/Fs:tPulse-1/Fs));


% part 3.1.A - Transmitting the signal
% transmitter
lengthSignal = 1000; % arbitrary | should be even - big for a beautiful plot
inputSignal = randi([0 1], 1, lengthSignal);


% Devide
[out1, out2] = Devide(inputSignal);

% Pulse maker
x1 = PulseShaping(out1, zeroPulse, onePulse);
x2 = PulseShaping(out2, zeroPulse, onePulse);

% Modulation
xc = AnalogMod(x1,x2,Fs,Fc);

% receiver
% channel
y = Channel(xc,Fs,FcChannel,BWchannel);


figure;
for i = 1:length(selectedVars)
    mu = 0;
    sigma = selectedVars(i);
    noise = normrnd(mu,sigma,1,length(y));
    noisyY = y + noise;


    % Demodulation
    [y1, y2] = AnalogDemod(noisyY,Fs,BWchannel,Fc);

    % MatchedFilt
    [raw0 raw1 Reconstructed1] = MatchedFilt(y1, zeroPulse, onePulse);
    [raw2 raw3 Reconstructed2] = MatchedFilt(y2, zeroPulse, onePulse);

    % Combine
    bprime = Combine(Reconstructed1,Reconstructed2);
    ReconstructedSignal = bprime;  
    
    subplot(2,3,i);
    hold on;
    scatter(raw1,raw3);
    grid on; grid minor;
    title("Signal Constellation for SD = " + selectedVars(i),'interpreter','latex');
    xlabel("$\hat{b_1}[i]$",'interpreter','latex')
    ylabel("$\hat{b_2}[i]$",'interpreter','latex')
end




%% part 3.3 - FSK modulation
% initializations 
clc; clear; close all;
Fs = 10^6; % Hz
tPulse = 10^-2; % in s
Fc = 10^4; % in Hz
FcChannel = 10^4; % in Hz
BWchannel = 10^3; % in Hz

% part 3.1
f1 = 1000; % in Hz
f0 = 1500; % in Hz
onePulse = cos(2*pi*f1*(0:1/Fs:tPulse-1/Fs));
zeroPulse = cos(2*pi*f0*(0:1/Fs:tPulse-1/Fs));

% part 3.1.A - Transmitting the signal
% transmitter
lengthSignal = 10; % arbitrary | should be even
inputSignal = randi([0 1], 1, lengthSignal)
figure;
scatter(1:length(inputSignal),inputSignal,'filled');
grid on; grid minor;
title("Input Digital Signal",'interpreter','latex');
xlabel('numberOfBit','interpreter','latex');

% Devide
[out1, out2] = Devide(inputSignal);
figure;
subplot(1,2,1);
scatter(1:length(out1),out1,'filled');
grid on; grid minor;
title('$b_1[n]$','interpreter','latex');
xlabel('numberOfBit','interpreter','latex');
subplot(1,2,2);
scatter(1:length(out2),out2,'filled');
grid on; grid minor;
title('$b_2[n]$','interpreter','latex');
xlabel('numberOfBit','interpreter','latex');

% Pulse maker
x1 = PulseShaping(out1, zeroPulse, onePulse);
x2 = PulseShaping(out2, zeroPulse, onePulse);
figure;
subplot(1,2,1);
plot((0:1/Fs:length(x1)/Fs-1/Fs)*1000,x1,'LineWidth',2);
grid on; grid minor;
title('$x_1(t)$','interpreter','latex');
xlabel('time(ms)','interpreter','latex');
ylabel('amplitude','interpreter','latex');
xlim([0 50])

subplot(1,2,2);
plot((0:1/10^6:length(x2)/10^6-1/10^6)*1000,x2,'LineWidth',2);
grid on; grid minor;
title('$x_2(t)$','interpreter','latex');
xlabel('time(ms)','interpreter','latex');
ylabel('amplitude','interpreter','latex');
xlim([0 50])

% Modulation
xc = AnalogMod(x1,x2,Fs,Fc);
figure;
plot((0:1/Fs:length(xc)/Fs-1/Fs)*1000,xc)
grid on; grid minor;
title('$x_c(t)$','interpreter','latex');
xlabel('time(ms)','interpreter','latex');
ylabel('amplitude','interpreter','latex');

% receiver
% channel
y = Channel(xc,Fs,FcChannel,BWchannel);
figure
plot((0:1/Fs:length(y)/Fs-1/Fs)*1000,y)
grid on; grid minor;
title('y(t)','interpreter','latex');
xlabel('time(ms)','interpreter','latex');
ylabel('amplitude','interpreter','latex');

% Demodulation
[y1, y2] = AnalogDemod(y,Fs,BWchannel,Fc);
figure;
subplot(1,2,1);
plot((0:1/Fs:length(y1)/Fs-1/Fs)*1000,y1)
grid on; grid minor;
title('$y_1(t)$','interpreter','latex');
xlabel('time(ms)','interpreter','latex');
ylabel('amplitude','interpreter','latex');
xlim([0 50]);

subplot(1,2,2);
plot((0:1/Fs:length(y2)/Fs-1/Fs)*1000,y2)
grid on; grid minor;
title('$y_2(t)$','interpreter','latex');
xlabel('time(ms)','interpreter','latex');
ylabel('amplitude','interpreter','latex');
xlim([0 50]);

% MatchedFilt and threshold
[raw0 raw1 Reconstructed1] = MatchedFilt(y1, zeroPulse, onePulse);
[raw2 raw3 Reconstructed2] = MatchedFilt(y2, zeroPulse, onePulse);

figure;
subplot(1,2,1);
scatter(1:length(Reconstructed1),Reconstructed1,'filled');
grid on; grid minor;
title('$\hat{b_1}[n]$','interpreter','latex');
xlabel('numberOfBit','interpreter','latex');
subplot(1,2,2);
scatter(1:length(Reconstructed2),Reconstructed2,'filled');
grid on; grid minor;
title('$\hat{b_2}[n]$','interpreter','latex');
xlabel('numberOfBit','interpreter','latex');

% Combine
bprime = Combine(Reconstructed1,Reconstructed2)
ReconstructedSignal = bprime;
figure;
subplot(2,1,1)
scatter(1:length(inputSignal),inputSignal,'filled');
grid on; grid minor;
title("Input Digital Signal",'interpreter','latex');
xlabel('numberOfBit','interpreter','latex');
subplot(2,1,2);
scatter(1:length(ReconstructedSignal),ReconstructedSignal,'filled');
grid on; grid minor;
title("Reconstructed Signal",'interpreter','latex');
xlabel('numberOfBit','interpreter','latex');
%% part 3.3.B
% initializations 
clc; clear; close all;
Fs = 10^6; % Hz
tPulse = 10^-2; % in s
Fc = 10^4; % in Hz
FcChannel = 10^4; % in Hz
BWchannel = 10^3; % in Hz

% part 3.1
f1 = 1000; % in Hz
f0 = 1500; % in Hz
onePulse = cos(2*pi*f1*(0:1/Fs:tPulse-1/Fs));
zeroPulse = cos(2*pi*f0*(0:1/Fs:tPulse-1/Fs));

% part 3.1.A - Transmitting the signal
% transmitter
lengthSignal = 10; % arbitrary | should be even
inputSignal = randi([0 1], 1, lengthSignal);
%figure;
%scatter(1:length(inputSignal),inputSignal);

% Devide
[out1, out2] = Devide(inputSignal);

% Pulse maker
x1 = PulseShaping(out1, zeroPulse, onePulse);
x2 = PulseShaping(out2, zeroPulse, onePulse);

% Modulation
xc = AnalogMod(x1,x2,Fs,Fc);

% receiver
% channel
y = Channel(xc,Fs,FcChannel,BWchannel);
% noise
error = zeros(1,100);
for i=1:100
    for j=1:5
        
        [i j]
        mu = 0;
        sigma = i;
        noise = normrnd(mu,sigma,1,length(y));
        noisyY = y + noise;


        % Demodulation
        [y1, y2] = AnalogDemod(noisyY,Fs,BWchannel,Fc);

        % MatchedFilt
        [raw0 raw1 Reconstructed1] = MatchedFilt(y1, zeroPulse, onePulse);
        [raw2 raw3 Reconstructed2] = MatchedFilt(y2, zeroPulse, onePulse);

        % Combine
        bprime = Combine(Reconstructed1,Reconstructed2);
        ReconstructedSignal = bprime;
        
       % error calculation
       error(i) = error(i) + sum(abs(bprime-inputSignal));
    end
    error(i) = error(i)/5/10;
end

figure;
stem(1:length(error),error);
grid on; grid minor;
title("Probability of Error",'interpreter','latex');
xlabel('Standard Deviation of noise','interpreter','latex');
ylabel('Probability','interpreter','latex');

%% part 3.3.C
clc; clear; close all;
selectedVars = [2 10 28 32 40 73];
% initializations 
Fs = 10^6; % Hz
tPulse = 10^-2; % in s
Fc = 10^4; % in Hz
FcChannel = 10^4; % in Hz
BWchannel = 10^3; % in Hz

% part 3.1
f1 = 1000; % in Hz
f0 = 1500; % in Hz
onePulse = cos(2*pi*f1*(0:1/Fs:tPulse-1/Fs));
zeroPulse = cos(2*pi*f0*(0:1/Fs:tPulse-1/Fs));


% part 3.1.A - Transmitting the signal
% transmitter
lengthSignal = 1000; % arbitrary | should be even - big for a beautiful plot
inputSignal = randi([0 1], 1, lengthSignal);


% Devide
[out1, out2] = Devide(inputSignal);

% Pulse maker
x1 = PulseShaping(out1, zeroPulse, onePulse);
x2 = PulseShaping(out2, zeroPulse, onePulse);

% Modulation
xc = AnalogMod(x1,x2,Fs,Fc);

% receiver
% channel
y = Channel(xc,Fs,FcChannel,BWchannel);


figure;
for i = 1:length(selectedVars)
    mu = 0;
    sigma = selectedVars(i);
    noise = normrnd(mu,sigma,1,length(y));
    noisyY = y + noise;


    % Demodulation
    [y1, y2] = AnalogDemod(noisyY,Fs,BWchannel,Fc);

    % MatchedFilt
    [raw0 raw1 Reconstructed1] = MatchedFilt(y1, zeroPulse, onePulse);
    [raw2 raw3 Reconstructed2] = MatchedFilt(y2, zeroPulse, onePulse);

    % Combine
    bprime = Combine(Reconstructed1,Reconstructed2);
    ReconstructedSignal = bprime;  
    
    subplot(2,3,i);
    hold on;
    scatter(raw1,raw3);
    grid on; grid minor;
    title("Signal Constellation for SD = " + selectedVars(i),'interpreter','latex');
    xlabel("$\hat{b_1}[i]$",'interpreter','latex')
    ylabel("$\hat{b_2}[i]$",'interpreter','latex')
end


%% part 4.1
% functions are implemented at the end of the code
%% part 4.2 - it takes a lot to run!
clc; clear; close all;
inputSequence = randi([0 255],1,5);
codedSequence = reshape(SourceGenerator(inputSequence), [1 40]);


Fs = 10^6; % Hz
tPulse = 10^-2; % in s
Fc = 10^4; % in Hz
FcChannel = 10^4; % in Hz
BWchannel = 10^3; % in Hz

% part 3.1
onePulse = ones(1,tPulse*Fs);
zeroPulse = -1*ones(1,tPulse*Fs);

% part 3.1.A - Transmitting the signal
% transmitter
inputSignal = codedSequence;
%figure;
%scatter(1:length(inputSignal),inputSignal);

% Devide
[out1, out2] = Devide(inputSignal);

% Pulse maker
x1 = PulseShaping(out1, zeroPulse, onePulse);
x2 = PulseShaping(out2, zeroPulse, onePulse);

% Modulation
xc = AnalogMod(x1,x2,Fs,Fc);

% receiver
% channel
y = Channel(xc,Fs,FcChannel,BWchannel);
% noise
error = zeros(5,50);
for i=1:20:1000
    for j=1:50
        [i j]
        mu = 0;
        sigma = i;
        noise = normrnd(mu,sigma,1,length(y));
        noisyY = y + noise;


        % Demodulation
        [y1, y2] = AnalogDemod(noisyY,Fs,BWchannel,Fc);

        % MatchedFilt
        [raw0 raw1 Reconstructed1] = MatchedFilt(y1, zeroPulse, onePulse);
        [raw2 raw3 Reconstructed2] = MatchedFilt(y2, zeroPulse, onePulse);

        % Combine
        bprime = Combine(Reconstructed1,Reconstructed2);
        ReconstructedSignal = bprime;
        
       % error calculation
       error(j,(i-1)/20+1) = sum(abs(bprime-inputSignal));
    end
end

reconSeq = OutputDecoder(reshape(ReconstructedSignal,[5 8])); 

figure;
x = var(error,0,1);
plot(1:50,x(1:50));
grid on; grid minor;
title('Var error - SD','interpreter','latex')
xlabel('Standard Deviation of Noise','interpreter','latex');
ylabel('Var of Reconstruction Error','interpreter','latex');
%% part 4.3
clc; clear; close all;
selectedVars = [10 28 32 100 500 700];
inputSequence = randi([0 255],1,5);
codedSequence = reshape(SourceGenerator(inputSequence), [1 40]);


Fs = 10^6; % Hz
tPulse = 10^-2; % in s
Fc = 10^4; % in Hz
FcChannel = 10^4; % in Hz
BWchannel = 10^3; % in Hz

% part 3.1
onePulse = ones(1,tPulse*Fs);
zeroPulse = -1*ones(1,tPulse*Fs);

% part 3.1.A - Transmitting the signal
% transmitter
inputSignal = codedSequence;
%figure;
%scatter(1:length(inputSignal),inputSignal);

% Devide
[out1, out2] = Devide(inputSignal);

% Pulse maker
x1 = PulseShaping(out1, zeroPulse, onePulse);
x2 = PulseShaping(out2, zeroPulse, onePulse);

% Modulation
xc = AnalogMod(x1,x2,Fs,Fc);

% receiver
% channel
y = Channel(xc,Fs,FcChannel,BWchannel);
% noise
error = zeros(5,6);
for i=1:length(selectedVars)
    for j=1:15
        
        [i j]
        mu = 0;
        sigma = selectedVars(i);
        noise = normrnd(mu,sigma,1,length(y));
        noisyY = y + noise;


        % Demodulation
        [y1, y2] = AnalogDemod(noisyY,Fs,BWchannel,Fc);

        % MatchedFilt
        [raw0 raw1 Reconstructed1] = MatchedFilt(y1, zeroPulse, onePulse);
        [raw2 raw3 Reconstructed2] = MatchedFilt(y2, zeroPulse, onePulse);

        % Combine
        bprime = Combine(Reconstructed1,Reconstructed2);
        ReconstructedSignal = bprime;
        
       % error calculation
       error(j,i) = sum(abs(bprime-inputSignal));
    end
    figure;
    histogram(error,100)
    title("Error Distribution - for SD:" + selectedVars(i));
    grid on; grid minor;
end

reconSeq = OutputDecoder(reshape(ReconstructedSignal,[5 8])); 
%% part 5 - Source Coding
clc; clear; close all;
% part 5.5
n = 10;
sequence = InformationSource(n)
% part 5.6
encodedSeq = SourceEncoder(n,sequence)
% part 5.7
decodedSeq = SourceDecoder(encodedSeq)

%% part 5.8
clc; clear; close all;

Hn = zeros(1,5000);
for i = 1:5000
    n = i;
    sequence = InformationSource(n);
    encodedSeq = SourceEncoder(n,sequence);
    Hn(i) = length(encodedSeq)/n;
end
plot(1:5000,Hn);
grid on; grid minor;
title('\(H_n(x)\)','interpreter','latex')
xlabel('n','interpreter','latex');
ylabel('$H_n(X) = \frac{L_B(n)}{n}$','interpreter','latex')
%% ALL FUNCTIONS - functions will be defined here 

% part.2.1 functions
function [b1, b2] = Devide(b)
    b1 = b(1:2:end);
    b2 = b(2:2:end);
end

function b = Combine(b1,b2)
    b = zeros(1,length(b1)+length(b2));
    b(1:2:end) = b1;
    b(2:2:end) = b2;
end

% part.2.2 functions
function pulse = PulseShaping(inputSequence, zeroPulse, onePulse)
    pulse = ((inputSequence.').*onePulse) + ...
        ((1-inputSequence.').*zeroPulse);
    pulseSize = size(pulse);
    pulse = reshape(pulse.',[1 pulseSize(1)*pulseSize(2)]);
end

% part.2.3 functions
function xc = AnalogMod(pulse1,pulse2,fs,fc)
    t = 0:1/fs:(length(pulse1)/fs)-(1/fs);
    xc = pulse1.*cos(2*pi*fc.*t) + pulse2.*sin(2*pi*fc.*t);
end

% part.2.4 functions
function channelOutput = Channel(xc,fs,fcCh,BWch)
    channelOutput = bandpass(xc,[fcCh-BWch/2,fcCh+BWch/2],fs);
end

% part.2.5 functions
function [y1, y2] = AnalogDemod(xc,fs,BWs,fc)
    t = 0:1/fs:(length(xc)/fs)-(1/fs);
    y1 = lowpass(xc.*cos(2*pi.*fc.*t),BWs,fs);
    y2 = lowpass(xc.*sin(2*pi.*fc.*t),BWs,fs);
end

% part.2.6 functions
function [raw0 raw1 Reconstructed] = MatchedFilt(demodulatedSignal, zeroPulse, onePulse)
    Reconstructed = zeros(1,length(demodulatedSignal)/length(zeroPulse));
    raw0 = zeros(1,length(demodulatedSignal)/length(zeroPulse));
    raw1 = zeros(1,length(demodulatedSignal)/length(zeroPulse));
    r0 = conv(demodulatedSignal,zeroPulse);
    r1 = conv(demodulatedSignal,onePulse);
    
    for i=1:length(Reconstructed)
        raw0(i) = r0(i*length(zeroPulse));
        raw1(i) = r1(i*length(zeroPulse));
        if(r1(i*length(zeroPulse)) > r0(i*length(zeroPulse)))
            Reconstructed(i) = 1;
        else
            Reconstructed(i) = 0;
        end
    end
end

% part.4 functions
function output = SourceGenerator(inputSequence)
    output = de2bi(inputSequence,8);
end

function output = OutputDecoder(inputSequence)
    output = bi2de(inputSequence);
end

% part.5 functions
function output = InformationSource(n)
    symbols = ['a','b','c','d','e','f'];
    probabilities = [1/2,1/4,1/8,1/16,1/32,1/64];
    output = randsample(symbols,n, true, probabilities);
end

function encodedSeq = SourceEncoder(n,seq)
    % codeWords after optimization problem
    % 1 - 01 - 001 - 0001 - 00001 - 00000
    % a - b - c - d - e - f
    codeWords = ['1', '01', '001', '0001', '00001', '00000'];
    encodedSeq = [];
    for i=1:n
        if (seq(i) == 'a')
            encodedSeq = [encodedSeq '1'];
        elseif (seq(i) == 'b')
            encodedSeq = [encodedSeq '01'];
        elseif (seq(i) == 'c')
            encodedSeq = [encodedSeq '001'];    
        elseif (seq(i) == 'd')
            encodedSeq = [encodedSeq '0001'];        
        elseif (seq(i) == 'e')
            encodedSeq = [encodedSeq '00001'];
        elseif (seq(i) == 'f')
            encodedSeq = [encodedSeq '00000'];
        end
    end
end

function decodedSeq = SourceDecoder(encodedSeq)
    codeWords = ['1', '01', '001', '0001', '00001', '00000'];
    decodedSeq = [];
    i = 1;
    while i<=length(encodedSeq)
        if (encodedSeq(i) == '1') 
            decodedSeq = [decodedSeq 'a'];
            i = i+1;
        elseif ((encodedSeq(i) == '0')  && (encodedSeq(i+1) == '1'))
            decodedSeq = [decodedSeq 'b'];
            i = i+2;
        elseif ((encodedSeq(i) == '0')  && (encodedSeq(i+1) == '0') && (encodedSeq(i+2) == '1'))
            decodedSeq = [decodedSeq 'c'];
            i = i+3;
        elseif ((encodedSeq(i) == '0')  && (encodedSeq(i+1) == '0') && (encodedSeq(i+2) == '0') && (encodedSeq(i+3) == '1'))
            decodedSeq = [decodedSeq 'd'];
            i = i+4;
        elseif ((encodedSeq(i) == '0')  && (encodedSeq(i+1) == '0') && (encodedSeq(i+2) == '0') && (encodedSeq(i+3) == '0') && (encodedSeq(i+4) == '1'))
            decodedSeq = [decodedSeq 'e'];
            i = i+5;
        elseif ((encodedSeq(i) == '0')  && (encodedSeq(i+1) == '0') && (encodedSeq(i+2) == '0') && (encodedSeq(i+3) == '0') && (encodedSeq(i+4) == '0'))
            decodedSeq = [decodedSeq 'f'];
            i = i+5;
        end
    end   
end