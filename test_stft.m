clc; 
%figure(2); clf(2); 
figure(1); clf(1); clear all;

audioFile = "test_10mps.wav";
[Amps, Fs] = audioread(audioFile);
audio = audioplayer(Amps,Fs); play(audio); 
N = size(Amps,1); % number of samples

WINDOW = hann(2048);
[stfourier, f, t] = stft(Amps, Fs, FrequencyRange="onesided", ...
    Window=WINDOW);
ampStft = abs(stfourier);

deltaT = t(end)-t(end-1);
deltaF = f(end)-f(end-1);
%f,t,ampStft

c = "default"; colormap(c);
imagesc(t, f, ampStft);
hold on;
set(gca,'YDir','normal');
xlabel("Time (s)");
ylabel("Frequency (Hz)");

changeIndices = findchangepts(ampStft,MaxNumChanges=2,Statistic="rms");
beginT = changeIndices(1)*deltaT;
endT = changeIndices(2)*deltaT;

%plot vertical lines at largest changes
xline(beginT, "r", LineWidth=1);
xline(endT, "r", LineWidth=1);
axis([0, max(t), 0, floor(max(f)*0.15)])
