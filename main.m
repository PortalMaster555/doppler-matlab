clc; clear all;
% Conversion functions
mps2mph = @(v) v * 2.237;
mph2mps = @(v) v / 2.237;

% CONSTANTS (tweak as needed)

WINDOW_SIZE = 2048; %larger seems to work better, but the changepoint 
                    % calculation hangs sometimes

% Used to find the speed of sound, leave as NaN if unknown.                    
TEMPERATURE_F = NaN; %F

% Used to make comparisons, leave as NaN if unknown.   
TRUE_V_MPH = NaN; % mph

% audioFile = "./audio/50mphobserver.wav"; TRUE_V_MPH = 50; 
% audioFile = "horn.ogg"; trueV = NaN;
% audioFile = "test_10mps.wav";
% audioFile = "test_440_48khz.wav";
% audioFile = "./audio/56mph44F.wav"; TRUE_V_MPH = 56; TEMPERATURE_F = 44; %F
audioFile = "./audio/appx40to50.wav";
% audioFile = "./audio/onecar/30kph.wav"; TRUE_V_MPH = 18.64;
% audioFile = "./audio/onecar/50kph.wav"; TRUE_V_MPH = 31.07;
% audioFile = "./audio/onecar/70kph.wav"; TRUE_V_MPH = 43.50;

% Minimum mean frequency amplitude prominence to be considered a peak
% E.g. 0.5
MIN_FREQ_PROMINENCE = 0.5;

% Minimum jump in velocity estimate curve to be considered a peak
% E.g. 5 (may need to be tweaked based on available sample data)
MIN_DIFF_PROMINENCE = 1; %mph

% Calculate speed of sound based on temperature
if ~isnan(TEMPERATURE_F)
     C = vSound(TEMPERATURE_F);
else
     C = vSound();
end

% Doppler shift not detectable in supersonic objects
CUTOFF_VELOCITY_MPH = mps2mph(C); 

% Additional cutoff if results make no sense
% 100 mph reasonable for automobiles
CUTOFF_VELOCITY_MPH = 100;


%


%   %   %   %   %
% BEGIN PROGRAM %
%   %   %   %   %

printSplash();

[amps, Fs] = audioread(audioFile);
N = size(amps,1); % number of samples
fprintf("Audio file loaded: %s (length %d samples)\n", audioFile, N);
fprintf("Maximum detectable frequency at or just below  %f Hz\n", Fs/2);

% Comment if playback not needed
playAudio(amps, Fs, 2);

% Compute the short-time Fourier transform
% Windowing function: Hann
% One-sided frequency range (absolute amplitude not preserved)

[stfourier, f, t] = stft(amps, Fs, FrequencyRange="onesided", ...
    Window=hann(WINDOW_SIZE));
ampsStft = abs(stfourier);

% Smallest resolvable timesteps and frequency steps
deltaT = t(end)-t(end-1);
deltaF = f(end)-f(end-1);

% Plot vertical lines at largest changes
changeIndices = findchangepts(ampsStft,MaxNumChanges=2,Statistic="rms");
if size(changeIndices, 2) < 2
    % Fixes a bug where it doesn't find a second change for some files
    changeIndices(2) = changeIndices(1);
end
% Beginning and ending times of the transition interval (used for plots)
beginT = changeIndices(1)*deltaT;
endT = changeIndices(2)*deltaT;

% Plot spectrogram
fprintf("Plotting spectrogram...\n")
figure(1); clf(1);
imagesc(t, f, ampsStft); hold on;
set(gca,'YDir','normal');
xlabel("Time (s)"); ylabel("Frequency (Hz)");
cb = colorbar; cb.TicksMode = "manual";
ylabel(cb,'Amplitude', Rotation=270);
xline(beginT, "r", LineWidth=1);
xline(endT, "r", LineWidth=1);

% % plot frequency breakdown in 3d
% figure(2); clf(2);
% surf(t, f, ampsStft);
% shading interp

% figure(3); clf(3);
% fprintf("Plotting zoomed spectrogram...\n")
% cmap = "default"; colormap(cmap);
% imagesc(t, f, ampsStft);
% hold on;
% set(gca,'YDir','normal');
% xlabel("Time (s)");
% ylabel("Frequency (Hz)");
% axis([0, max(t), 0, 5000]);
% colorbar;
% xline(beginT, "r", LineWidth=1);
% xline(endT, "r", LineWidth=1);

fprintf("Detecting edges...\n");
edges = edge(ampsStft, 'Canny');

% Isolate ends (outside of the middle transition)
begRange = 1:changeIndices(1);
endRange = changeIndices(2):size(edges,2);

figure(4); clf(4);
subplot(2,1,1);
imagesc(t(begRange), f, edges(:, begRange));
set(gca,'YDir','normal');
subplot(2,1,2);
imagesc(t(endRange), f, edges(:, endRange));
set(gca,'YDir','normal');

fprintf("Plotting beginning and end spectrograms...\n");
% Plot beginning frequencies
figure(11); clf(11);
begFiltAmps = ampsStft(:,begRange);
endFiltAmps = ampsStft(:,endRange);
imagesc(t(begRange), f, begFiltAmps);
set(gca,'YDir','normal');

% Plot end frequencies
figure(12); clf(12);
imagesc(t(endRange), f, endFiltAmps);
set(gca,'YDir','normal');

fprintf("Computing averages...\n");
% Find the averages for each frequency section
begFreqAvgs = mean(begFiltAmps, 2);
endFreqAvgs = mean(endFiltAmps, 2);

% Plot the averages and the peaks
figure(13); clf(13);
subplot(2,1,1);
stem(f, begFreqAvgs);
xlim([0 0.5e4]);
[pksBeg, locsBeg] = findpeaks(begFreqAvgs, MinPeakProminence=MIN_FREQ_PROMINENCE);


subplot(2,1,2);
stem(f, endFreqAvgs);
xlim([0 0.5e4]);
[pksEnd, locsEnd] = findpeaks(endFreqAvgs, MinPeakProminence=MIN_FREQ_PROMINENCE);

% Plot the highest prominence peaks on charts for reference

figure(14); clf(14);
pltB = subplot(2,1,1);
Ampa = begFreqAvgs(locsBeg);
stem(f(locsBeg), Ampa)

pltE = subplot(2,1,2);
Ampr = endFreqAvgs(locsEnd);
stem(f(locsEnd), Ampr)

fprintf("Constructing initial/final frequency pairs...\n");
% Construct pairs for each combination of FREQUENCIES.
Fa = f(locsBeg);
Fr = f(locsEnd);
% columns are F_a
[FaM, FrM] = ndgrid(Fa, Fr);
% straighten out into columns and merge into an Nx2 matrix
FaM = FaM(:);
FrM = FrM(:);
inputMatrix = [FaM FrM];


% Find velocities for all combinations
fprintf("Computing all combinations of velocities...\n");
for i = 1:size(inputMatrix, 1)
    options = optimset('Display','off');
    velFcn = @(v) (C-v)*inputMatrix(i,1) - (C+v)*inputMatrix(i,2);
    sourceV(i) = fsolve(velFcn, 0, options);
end
fprintf("Computed %d combinations of velocities.\n", size(sourceV,2));

% Plot all (including negative) velocities
figure(15); clf(15);
plot(mps2mph(sourceV), "b");
hold on; yline(50, "r");
ylabel("Velocity (MPH)");
title(audioFile);

% Remove negative velocities (non-physical)
possibleI = find(sourceV>=0);
sourceV = sourceV(possibleI);

% Plot all possible velocities
sortMphVels = mps2mph(sort(sourceV));
figure(16); clf(16);
% stem(mps2mph(sourceV), "b");
plot(sortMphVels, ".b");
ylabel("Velocity (MPH)");
title(audioFile);
% saveas(gcf, "50mph.png");


% EXTRACTING VELOCITY FROM THE CURVE
fprintf("Finding differences in possible velocities...\n");
% Find differences and peaks
figure(17); clf(17);
diffVector = diff(sortMphVels);
[peakDiffs,peakDiffLocs] = findpeaks(diffVector, ...
    MinPeakProminence=MIN_DIFF_PROMINENCE);
%plot
findpeaks(diffVector, ...
    MinPeakProminence=MIN_DIFF_PROMINENCE);

fprintf("Estimating mean velocities...\n");
% two point mean estimate
format short g
beforeVJumps = sortMphVels(peakDiffLocs);
afterVJumps = sortMphVels(peakDiffLocs+1);

dif = beforeVJumps(1);
meanV = 0 + dif/2;
for i = 1:size(beforeVJumps,2)-1
    dif(i+1) = beforeVJumps(i+1)-afterVJumps(i);
    meanV(i+1) = afterVJumps(i) + dif(i+1)/2;
end

%find horizontal span of the nonpeak regions
differenceArray = [0 peakDiffLocs];
for i = 1:size(differenceArray,2)-1
    horizDiffs(i) = differenceArray(i+1)-differenceArray(i);
end

pairs = [meanV; horizDiffs]';

fprintf("Performing cutoff at %f mph...\n", CUTOFF_VELOCITY_MPH);
% Perform cutoff checks
indVec = pairs(:,1) < CUTOFF_VELOCITY_MPH;
pairs = [pairs(indVec, 1) pairs(indVec, 2)];

% Place additional arbitrary limitations here.

pairsSort = flip(sortrows(pairs, 2));
figure(16); hold on; yline(pairsSort(1,1), "m"); yline(TRUE_V_MPH, "r");
legend("","Estimated V", "True V", Location="northwest");

fprintf("The most likely velocity of this object is %f mph.\n" + ...
    "\t(%d contiguous entries, cutoff of %f mph.\n", ...
    pairsSort(1,1), pairsSort(1,2),  CUTOFF_VELOCITY_MPH);
if ~isnan(TRUE_V_MPH)
    fprintf("This is %f mph ", abs(TRUE_V_MPH - pairsSort(1,1)))
    if TRUE_V_MPH - pairsSort(1,1) <= 0
        fprintf("above the true velocity, %f mph.\n", TRUE_V_MPH);
    else
        fprintf("below the true velocity, %f mph.\n", TRUE_V_MPH);
    end
    fprintf("Possible velocities, contiguous entries, and distance from trueV:\n");
    disp([pairsSort TRUE_V_MPH-pairsSort(:,1)]);
else
    fprintf("Possible velocities and contiguous entries:\n");
    disp(pairsSort);
end