clc; printSplash(); clear all;

audioFile = "test_440_48khz.wav";

[Amps, Fs] = audioread(audioFile);
T = size(Amps,1) / Fs; % seconds
fprintf("Audio file loaded: %s (length %f s)\n", audioFile, T);

fprintf("Playing sample of audio if available.\n")
audio = audioplayer(Amps,Fs);
play(audio); pause(1); stop(audio);

fprintf("Maximum detectable frequency at or just below  %f Hz\n", Fs/2);

figure(1);

