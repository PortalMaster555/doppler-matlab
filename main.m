clc; printSplash(); clear all;

audioFile = "test_440_48khz.wav";

[Amps, Fs] = audioread(audioFile);
N = size(Amps,1); % number of samples
fprintf("Audio file loaded: %s (length %d samples)\n", audioFile, N);

fprintf("Playing sample of audio if available.\n")
audio = audioplayer(Amps,Fs); play(audio); pause(1.5); stop(audio);

fprintf("Maximum detectable frequency at or just below  %f Hz\n", Fs/2);

% Plot scaled input audio file vs. time
scaleConst = Fs/25;
figure(1);
t = linspace(0,N/Fs,N);
fprintf("Length of audio file is %f seconds\n", max(t));
stem(t*scaleConst, Amps, ".k");
xlabel("Time (not to scale)"); ylabel("Amplitude");
axis([0, max(t), -1, 1]);

% Compute full FFT for the input audio.
fourier = pos_fft(Amps);
f = linspace(0,Fs/2,size(fourier,1));
figure(2); clf(2); plot(f,fourier,".b");



