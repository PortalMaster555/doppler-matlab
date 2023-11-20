clc; printSplash(); figure(1); clf(1); clear all;

% CONSTANTS (tweak as needed)
smoothingConstant = 20;
num_to_take = 10;
audioFile = "horn.ogg";

[Amps, Fs] = audioread(audioFile);
N = size(Amps,1); % number of samples
fprintf("Audio file loaded: %s (length %d samples)\n", audioFile, N);

fprintf("Playing sample of audio...\n")
audio = audioplayer(Amps,Fs); play(audio); 
pause(1); stop(audio);
fprintf("Playback complete.\n")

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

% Plot spectrogram
fprintf("Plotting spectrogram...\n")
imagesc(t, f, ampStft); hold on;
set(gca,'YDir','normal');
xlabel("Time (s)"); ylabel("Frequency (Hz)");

cb = colorbar; cb.TicksMode = "manual";
ylabel(cb,'Amplitude', Rotation=270);

% Plot vertical lines at largest changes
changeIndices = findchangepts(ampStft,MaxNumChanges=2,Statistic="rms");
if size(changeIndices, 2) < 2
    %fixes a bug where it doesn't find a second change
    changeIndices(2) = changeIndices(1);
end
beginT = changeIndices(1)*deltaT;
endT = changeIndices(2)*deltaT;
fprintf("Plotting time position of change points...\n")
xline(beginT, "r", LineWidth=1); xline(endT, "r", LineWidth=1);

% Draw a colormap of the unsmoothed audio spectrogram
fprintf("Drawing unsmoothed colormap...\n")
figure(2); clf(2);
c = "default"; colormap(c);
imagesc(t, f, ampStft);
hold on;
set(gca,'YDir','normal');
xlabel("Time (s)");
ylabel("Frequency (Hz)");
axis([0, max(t), 0, 5000]);
colorbar;

% Perform very aggressive smoothing on the image across time
fprintf("Smoothing over time...\nAdjust smoothing constant if" + ...
    " results are not valid.\n")

smoothAmps = movmean(ampStft', smoothingConstant);
smoothAmps = smoothAmps';

figure(3); clf(3);
fprintf("Plotting smoothed spectrogram...\n")
c = "default"; colormap(c);
imagesc(t, f, smoothAmps);
hold on;
set(gca,'YDir','normal');
xlabel("Time (s)");
ylabel("Frequency (Hz)");
axis([0, max(t), 0, 5000]);
colorbar;
xline(beginT, "r", LineWidth=1);
xline(endT, "r", LineWidth=1);

fprintf("Detecting edges...\n")
smoothAmps = edge(smoothAmps, 'Canny');

% Isolate ends (outside of the middle transition)
begRange = 1:changeIndices(1);
endRange = changeIndices(2):size(smoothAmps,2);
begAmps = smoothAmps(:, begRange);
endAmps = smoothAmps(:, endRange);
tbeg = t(begRange); tend = t(endRange);
%imagesc(tbeg, f, begAmps);
%imagesc(tend, f, endAmps);

% Calculate mean of each frequency band across time
% outside of the middle region to find most constant frequencies
% frequencies across time (index 2)
fprintf("Calculating means of all frequencies...\n")
fBegSums = mean(begAmps, 2);
fEndSums = mean(endAmps, 2);
figure(4); clf(4);
stem(f, fBegSums, ".r");
hold on;
stem(f, fEndSums, ".b");
axis([0, 4000, 0, max(fBegSums)]);

% Take highest-amp frequencies (constant in beginning of file)

fprintf("Choosing and plotting highest %d means...\n", num_to_take);

[maxAmpsB, begMaxFIndices] = maxk(fBegSums, num_to_take);
highestFBeg = f(begMaxFIndices);
plot(highestFBeg, maxAmpsB, "om", LineStyle="none");
highestFBeg = sort(highestFBeg);

[maxAmpsE, endMaxFIndices] = maxk(fEndSums, num_to_take);
highestFEnd = f(endMaxFIndices);
plot(highestFEnd, maxAmpsE, "om",LineStyle="none");
highestFEnd = sort(highestFEnd);

c = vSound();
fprintf("Speed of sound taken as %f m/s.\n", c);

fprintf("Analyzing velocity based on highest %d means...\n", num_to_take);
options = optimset('Display','off');
for i = 1:size(highestFBeg, 1)
    app(i) = highestFBeg(i);
    rec(i) = highestFEnd(i);
    velFcn = @(v) (c-v)*app(i) - (c+v)*rec(i);
    sourceV(i) = fsolve(velFcn, 0, options);
    %fprintf("Approach: %f Hz\nRecede: %f Hz\nVel Est.: %f m/s\n\n", ...
    %app(i), rec(i), sourceV(i));
end

fprintf("Finding closest two velocities from the %d and averaging...\n", ...
    num_to_take);
minDiffInd = find(abs(diff(sourceV))==min(abs(diff(sourceV))));
%extract this index, and it's neighbor index from A
v1 = sourceV(minDiffInd);
v2 = sourceV(minDiffInd+1);

avgClosestV = mean([v1, v2]);
fprintf("Estimated velocity: %f m/s.\n", avgClosestV);

% Convert to MPH from m/s
avgClosestV = 2.237 * avgClosestV;
fprintf("Estimated velocity: %f mph.\n", avgClosestV);