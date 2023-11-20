clc; printSplash(); figure(1); clf(1); clear all;

audioFile = "50mphobserver.wav";

[Amps, Fs] = audioread(audioFile);
N = size(Amps,1); % number of samples
fprintf("Audio file loaded: %s (length %d samples)\n", audioFile, N);

%fprintf("Playing sample of audio...\n")
% audio = audioplayer(Amps,Fs); play(audio); 
% pause(1); stop(audio);
%fprintf("Playback complete.\n")

fprintf("Maximum detectable frequency at or just below  %f Hz\n", Fs/2);

% Compute the short-time Fourier transform
% Windowing function: length 2048 Hann
% One-sided frequency range (absolute amplitude not preserved)
WINDOW = hann(2048);
[stfourier, f, t] = stft(Amps, Fs, FrequencyRange="onesided", ...
    Window=WINDOW);
ampStft = abs(stfourier);

deltaT = t(end)-t(end-1);
deltaF = f(end)-f(end-1);
changeIndices = findchangepts(ampStft,MaxNumChanges=2,Statistic="rms");
beginT = changeIndices(1)*deltaT;
endT = changeIndices(2)*deltaT;

% Plot spectrogram
fprintf("Plotting spectrogram...\n")
imagesc(t, f, ampStft); hold on;
set(gca,'YDir','normal');
xlabel("Time (s)"); ylabel("Frequency (Hz)");
axis([0, max(t), 0, floor(max(f)*0.15)])
cb = colorbar; cb.TicksMode = "manual";
ylabel(cb,'Amplitude', Rotation=270);

% Plot vertical lines at largest changes
fprintf("Plotting time position of change points...\n")
xline(beginT, "r", LineWidth=1); xline(endT, "r", LineWidth=1);

% Get mean frequency up to before the first changepoint


