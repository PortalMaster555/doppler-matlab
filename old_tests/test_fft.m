clc; clear all; figure(1); clf(1);



fcn = @(t) 5 + 4*cos(10*2*pi*t);
% expected result:
% 5 @ 0Hz
% 4 @ 10Hz
sampleFreq = 48; % Hz
fprintf("Maximum detectable frequency at or just below  %f Hz\n", ...
    sampleFreq/2);
T = 5; % seconds
N = sampleFreq * T; %samples
t = linspace(0,T,N);

y = fft(fcn(t));
absY = abs(y/N); %scaled two-sided amplitudes in f-space
if mod(N,2) == 0
    fourier = absY(1:N/2+1); % take the second half
else
    fourier = absY(1:(N+1)/2); % take the second half
end
fourier(2:end-1) = 2*fourier(2:end-1); % double all but 0 and nyq.


thresholdCond = 0.33*max(fourier);
fourierSmoothIndex= find(fourier < thresholdCond);
fourierSmooth = fourier;
fourierSmooth(fourierSmoothIndex) = 0;

f = sampleFreq/N*(0:(N/2));

figure(1);
splt = subplot(3,1,1);
plot(t,fcn(t),"-k");
title("fcn(t)");
xlabel("t (s)");
ylabel("Amplitude");

splt2 = subplot(3,1,2);
stem(f,fourier,".b");
title("Single-Sided Amplitude Spectrum of fcn(t)")
xlabel("f (Hz)")
ylabel("|fourier(f)|")
axis([-1, 1.1*max(f), 0, 1.05*max(fourier)]);

splt3 = subplot(3,1,3);
stem(f,fourierSmooth, ".r");
title("Threshold Single-Sided Amplitude Spectrum of fcn(t)")
xlabel("f (Hz)")
ylabel("|fourier(f)|")
axis([-1, 1.1*max(f), 0, 1.05*max(fourier)]);





















