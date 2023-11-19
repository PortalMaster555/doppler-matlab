clc; figure(1); clf(1); clear all;

audioFile = "50mphobserver.wav";
[Amps, Fs] = audioread(audioFile);
%audio = audioplayer(Amps,Fs); play(audio); 
N = size(Amps,1); % number of samples

WINDOW = hann(2048);

stfourier = stft(Amps, Fs, FrequencyRange="onesided", ...
    Window=WINDOW);
ampStft = abs(stfourier);

imagesc(ampStft);
set(gca,'YDir','normal');
xlabel("Time (nts)");
ylabel("Frequency (nts)");
axis([0, size(stfourier,2), 0, floor(0.25*size(stfourier,1))]);