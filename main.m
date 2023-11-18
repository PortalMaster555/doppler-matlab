clc; printSplash(); clear all;

audioFile = "test_440_48khz.wav";

[Amps, Fs] = audioread(audioFile);
N = size(Amps,1); % samples
fprintf("Audio file loaded: %s (length %d samples)\n", audioFile, N);

fprintf("Playing sample of audio if available.\n")
audio = audioplayer(Amps,Fs); play(audio); pause(1.5); stop(audio);

fprintf("Maximum detectable frequency at or just below  %f Hz\n", Fs/2);

figure(1);
t = linspace(0,N/Fs,N);
stem(Fs*t/25, Amps, ".k");
xlabel("Time (not to scale)"); ylabel("Amplitude");
axis([0, max(t), -1, 1]);

y = fft(Amps);
absY = abs(y/N); %scaled two-sided amplitudes in f-space
if mod(N,2) == 0
    fourier = absY(1:N/2+1); % take the second half
else
    fourier = absY(1:(N+1)/2); % take the second half
end
fourier(2:end-1) = 2*fourier(2:end-1); % double all but 0 and nyq.

f = linspace(0,Fs/2,size(fourier,1));
figure(2); clf(2); plot(f,fourier,".b");



