clc; printSplash(); clear all;

audioFile = "50mphobserver.wav";

[Amps, Fs] = audioread(audioFile);
N = size(Amps,1); % number of samples
fprintf("Audio file loaded: %s (length %d samples)\n", audioFile, N);

fprintf("Playing sample of audio:\n")
audio = audioplayer(Amps,Fs); play(audio); 
% pause(1.5); stop(audio);
fprintf("Playback complete.\n")

fprintf("Maximum detectable frequency at or just below  %f Hz\n", Fs/2);

% Plot scaled input audio file vs. time
figure(1);
t = linspace(0,N/Fs,N);
fprintf("Length of audio file is %f seconds\n", max(t));
stem(t, Amps, ".k");
xlabel("Time (s)"); ylabel("Amplitude");
axis([0, max(t), -1, 1]);

% Compute and plot full FFT for the input audio
fourier = pos_fft(Amps);
f = linspace(0,Fs/2,size(fourier,1));
figure(2); clf(2); stem(f,fourier,".b");
xlabel("Frequency (Hz)"); ylabel("Amplitude");



