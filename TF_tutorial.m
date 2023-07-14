%% TUTORIAL TIME-FREQUENCY

%this small tutorial will guide us through the analysis and principles of
%time-frequency representation of biological (Electroencephalography) data.
%
%tutorial written by Gianluigi Giannini, MSc. 
%
% have fun!

%% NYQUIST FREQUENCY
%you need to sample at twice the frequency of interest: frequency at half
%of the sampling rate can be recovered.
time = -1:0.001:(1-0.001);
A = 2;
figure; 
plot(time, A*sin(2*pi*2*time)); 
hold on
time = -1:0.01:1; %sampling 100Hz
plot(time,  A*sin(2*pi*2*time), 'r*')
time = -1:0.5:1; %sampling 2Hz
plot(time,  A*sin(2*pi*2*time), 'b*')
legend({'Wave at 2Hz', 'Sampling at 100Hz', 'Sampling at 2Hz'})

%% HOW DOES AN EEG LOOK LIKE
%deep look at the time frequency decomposition and how it relates to the
%signal we measure. 
time = -1:0.001:(1-0.001);
A = 2; 
sin10 = A*sin(2*pi*10*time); 
sin5 = A*sin(2*pi*5*time); 
sin30 = A*sin(2*pi*30*time); 
total_EEG = sin10 + sin5 + sin30; 
figure; 
subplot(4,1,1)
plot(time, sin5); 
title('5 Hz')
subplot(4,1,2)
plot(time, sin10)
title('10 Hz')
subplot(4,1,3)
plot(time, sin30)
title('30 Hz')
subplot(4,1,4)
plot(time, total_EEG)
title('total EEG')

%% CREATE WAVES
% create a simple sine wave over time -1 : 1
time = -1:0.001:(1-0.001);
plot(time, sin(2*pi*time))
%adjust the amplitude (how big is each peak)
A = 2; 
plot(time, A*sin(2*pi*time))
%adjust the frequency (how many peaks in a second)
f = 10; 
plot(time, A*sin(2*pi*f*time))
%add a phase (value of the sin at time 0)
o = pi; 
plot(time, A*sin(2*pi*f*time + o));

figure; 
subplot(4,1,1)
plot(time, sin(2*pi*time))
title('Simple Wave')
subplot(4,1,2)
plot(time, A*sin(2*pi*time))
title('add Amplitude')
subplot(4,1,3)
plot(time, A*sin(2*pi*f*time))
title('add Frequency')
subplot(4,1,4)
plot(time, A*sin(2*pi*f*time + o))
title('add Phase shift')



%% DOT PRODUCT
% another fundamental mathematical concept that we will take advantage of,
% is the dot product. 
%By its algebraical formulation, it is the sum of the multiplications of each element in a
%with each corresponding element in b. 
a = [1 2]; 
b = [2 3]; 
figure; 
plot([0 1], [0 2]);
hold on
plot([0 2], [0 3]); 

sum(a.*b)

%geometrically, the dot product can be interpreted as the length of the
%projection of one vector onto the other
%it's basically a "measure" of how much two vectors are similar

time = -1:0.001:(1-0.001);
c = sin(2*pi*time);
d = sin(2*pi*time*2); 
d = cos(2*pi*time); 
d = cos(2*pi*time +1/2*pi); 

figure; 
plot(time, c)
hold on
plot(time, d)
legend({'1Hz sin wave', '2Hz sin wave'})

sum(c.*d)


%% FOURIER TRANSFORM!
%simply calculate the dot product between our "raw" eeg data and each
%single "simple" wave. 
%in this way, we will estimate "how much" of each single frequency from 0 to 39
%constitutes our data. 
N = numel(time); %basically the number of our datapoint
for i = 1:40
    sine_wave = A*sin(2*pi*(i-1)*time); 
    fourier(i) = sum(sine_wave.*total_EEG); 
end
fourier = fourier / N;
plot((1:40)-1, fourier)
xlabel('Frequences')
ylabel('Amplitude')

%same thing but with fft function
%set up a couple of variables :) 
Fs = 1000; %sampling rate 
T = 1/Fs; %sampling period
L = numel(time); %length of signal
Y = fft(total_EEG);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = Fs*(0:(L/2))/L; 
plot(f, P1)

Y = fft(total_EEG);
plot(abs(Y/L))


for i = 1:40
    sine_wave = A*sin(2*pi*(i-1)*time); 
    sine_wave1 = exp(-1i*2*pi*(i-1).*time); 
    fourier(i) = sum(sine_wave.*total_EEG); 
    fourier1(i) = sum(sine_wave1.*total_EEG);
end
plot((1:40)-1, fourier)
plot(imag(fourier1))

%% CONVOLUTION
% impulse function (all zeros; 1 in the middle)
signal = zeros(1,100);
signal(40:60)=1;

kernel = [1 0.8 0.6 0.4 0.2]; 

% matlab's convolution function
convolved_signal = conv(signal,kernel,'same');

figure; 
subplot(3,1,1)
plot(1:100, signal)
title('Signal')
subplot(3,1,2)
plot(kernel)
xlim([1 100])
title('Kernel')
subplot(3,1,3)
plot(convolved_signal)
title('Convolved signal')

%% MORLET WAVELET
time = -1:0.001:1; 
f = 4; % frequency of sine wave in Hz
% create sine wave 
sine_wave = cos(2*pi*f.*time);
% make a Gaussian
s=4/(2*pi*f); %this is the 
gaussian_win = exp(-time.^2./(2*s^2));

figure;
subplot(3,1,1)
plot(time, sine_wave);
title('Sine wave at 4Hz')
subplot(3,1,2)
plot(time, gaussian_win)
title('Gaussian')
subplot(3,1,3)
plot(time,sine_wave.*gaussian_win)
title('Wavelet')

%create and plot a wavelet of frequency 4 and width 10
f = 4; % frequency of sine wave in Hz
% create sine wave 
a = cos(2*pi*f.*time);
% make a Gaussian
s=10/(2*pi*f); %this is the 
b = exp(-time.^2./(2*s^2));
figure; 
plot(time, a.*b)
title('Wavelet, f: 4, Width: 10')

%create and plot a wavelet of frequency 4 and width 2
f = 4; % frequency of sine wave in Hz
% create sine wave 
a = cos(2*pi*f.*time);
% make a Gaussian
s=2/(2*pi*f); %this is the 
b = exp(-time.^2./(2*s^2));
figure; 
plot(time, a.*b)
title('Wavelet, f: 4, Width: 2')

%also calculate the fft of the wavelets
fft_4_10 = fft(a.*b);
fft_4_10_corr = abs(fft_4_10/numel(time)); 
fft_4_10_corr = fft_4_10_corr(1:numel(time)/2+1);
fft_4_10_corr(2:end-1) = 2*fft_4_10_corr(2:end-1);
f = 1000*(0:(numel(time)/2))/numel(time); 
figure;
plot(f(1:40), fft_4_10_corr(1:40))
title('FFT of a wavelet with f: 4, Width: 10')

%% WAVELET CONVOLUTION
%wavelet convolution as the way to do time-frequency!
%first let's create a new EEG fake signal that has some 4 Hz + 10 Hz + 20
%Hz  
time = 0:0.001:5; 
total_EEG = sin(2*pi*4*time) + sin(2*pi*10*time) + sin(2*pi*20*time);
% sometimes there's also some 15 Hz for 0.5 secs
total_EEG(500:1000) = total_EEG(500:1000) + sin(2*pi*15*[0:0.001:0.5]);
total_EEG(2000:2500) = total_EEG(2000:2500) + sin(2*pi*15*[0:0.001:0.5]);
total_EEG(3000:3500) = total_EEG(3000:3500) + sin(2*pi*15*[0:0.001:0.5]);
%and sometimes there's also 30 Hz for 0.5 secs
total_EEG(1000:1500) = total_EEG(1000:1500) + sin(2*pi*30*[0:0.001:0.5]);
total_EEG(2000:2500) = total_EEG(2000:2500) + sin(2*pi*30*[0:0.001:0.5]);
%plot it
figure; 
plot(time, total_EEG)
title('New EEG signal')

%let's now try to see where the 15 Hz signal is
%first of all, I would be equipped with the right wavelet (or set of
%wavelets)
f = 15; 
centered_time = -2.5:0.001:2.5; 
s=7/(2*pi*f); 
Wavelet_15 = cos(2*pi*f.*centered_time) .* exp(-centered_time.^2./(2*s^2));
%plot(centered_time, Wavelet_15)

%convolve my signal with the wavelet that I just created
TF = conv(total_EEG, Wavelet_15, 'same'); 
%plot(TF)

figure; 
subplot(3,1,1)
plot(total_EEG)
title('Raw EEG')
subplot(3,1,2)
plot(Wavelet_15)
title('Wavelet')
subplot(3,1,3)
plot(TF)
title('Time representation, width 7')

figure
plot(TF)
hold on 
plot(TF.^2)


%% LET'S TRY NOW WITH FIELDTRIP :)
%first put our data in a structure that fieldtrip can recognise
data1 = []; 
data1.fsample = 1000; 
data1.label = {'AUX'};
data1.time = {time}; 
data1.trial = {total_EEG}; 
%then perform wavelet TF calculation
cfg = []; 
cfg.method = 'wavelet';
cfg.foi = 0:1:100; 
cfg.toi = 0:0.01:5; 
TF1 = ft_freqanalysis(cfg, data1); 
%plot it
pcolor(squeeze(TF1.powspctrm))


%% LET'S TRY WITH REAL DATA
%data available at: https://download.fieldtriptoolbox.org/tutorial/timefrequencyanalysis/dataFIC.mat

%load data in matlab (double click) and import
%explore together the data. This is how a data structure in fieldtrip looks
%like
%we see that the sampling rate was 500 Hz and that the timewindow spaces
%from -1 to 1
%TF
cfg = [];
cfg.method = 'wavelet';
cfg.foi = 1:1:30; 
cfg.toi = -1:0.03:2; 
cfg.width = 7;
TF = ft_freqanalysis(cfg, dataFIC); 

%plot them
cfg = []; 
ft_singleplotTFR(cfg, TF)

%calculate spectra in a different way
spectra = mean(TF.powspctrm(:,:,:),[1 3], "omitnan");
plot(1:30, spectra)

%baseline them!
cfg = []; 
cfg.baseline = [-0.5 0];
cfg.baselinetype = 'relative'; 
TF_baselined = ft_freqbaseline(cfg, TF); 

cfg = []; 
ft_singleplotTFR(cfg, TF_baselined)

spectra = mean(TF_baselined.powspctrm(:,:,:),[1 3], "omitnan");
plot(1:30, spectra)


%% A BETTER LOOK INTO BASELINES
%% Types of baseline
%data = whole signal
%meanVals = mean of the signal during baseline

%Absolute
data = data - meanVals;

%Relative
data = data ./ meanVals;

%Relative change
data = (data - meanVals) ./ meanVals;

%Nomative change
data = (data - meanVals) ./ (data + meanVals);

%Decibel conversion (log)
data = 10*log10(data ./ meanVals);

%Z-score
stdVals = repmat(nanstd(data(:,:,baselineTimes),1, 3), [1 1 size(data, 3)]);
data=(data-meanVals)./stdVals;











