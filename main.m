clc; clear all;
% Conversion functions
mps2mph = @(v) v * 2.237;
mph2mps = @(v) v / 2.237;
mph2kph = @(v) v * 1.609;
kph2mph = @(v) v / 1.609;

% CONSTANTS (tweak as needed)

% Short-time Fourier transform window size: 2048 works well
WINDOW_SIZE = 2048; 

% Used to find the speed of sound, leave as NaN if unknown.                    
TEMPERATURE_F = NaN; %F

% Used to make comparisons, leave as NaN if unknown.   
TRUE_V_MPH = NaN; % mph

% These work remarkably well!
% audioFile = "./audio/50mphobserver.wav"; TRUE_V_MPH = 50; 
% audioFile = "./audio/appx40to50.wav"; TRUE_V_MPH = 45;
    % Source estimate was approximately 40 to 50 mph, so within tolerance.
    % Line drawn in the middle of the range
% audioFile = "./audio/onecar/70kph.wav"; TRUE_V_MPH = kph2mph(70);
% audioFile = "./audio/onecar/50kph.wav"; TRUE_V_MPH = kph2mph(50);
    % Above seems to imply that the 50 is the wrong value (63kph instead???)

% TO-DO: Make difference code find continuous regions AND lowest slope.
audioFile = "./audio/onecar/30kph.wav"; TRUE_V_MPH = kph2mph(30); 

% This one is okay, but without true value it is uncertain (curved "flat" area).
% audioFile = "./audio/horn.ogg";

% Incredibly low number of velocity estimates (quiet, almost no overtones!)
% audioFile = "./audio/56mph44F.wav"; TRUE_V_MPH = 56; TEMPERATURE_F = 44; %F

% These break everything (TO-DO: Patch the difference code to skip when only one vel.)
% audioFile = "./audio/test_10mps.wav"; TRUE_V_MPH = mps2mph(10);
% audioFile = "./audio/test_440_48khz.wav"; TRUE_V_MPH = 0;


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
fprintf("Audio file loaded: %s\n", audioFile);
fprintf("Sample rate: %d Hz\nLength: %d samples\n", Fs, N);
fprintf("Maximum detectable frequency at or just below %f Hz.\n\n", Fs/2);

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
title("Beginning Edge Detection");
set(gca,'YDir','normal');
subplot(2,1,2);
imagesc(t(endRange), f, edges(:, endRange));
title("End Edge Detection");
set(gca,'YDir','normal');

fprintf("Plotting beginning and end spectrograms...\n");
% Plot beginning frequencies
figure(11); clf(11);
subplot(2,2,1);
begFiltAmps = ampsStft(:,begRange);
imagesc(t(begRange), f, begFiltAmps);
title("Beginning Frequencies vs. Time");
set(gca,'YDir','normal');

% Plot end frequencies
subplot(2,2,3);
endFiltAmps = ampsStft(:,endRange);
imagesc(t(endRange), f, endFiltAmps);
title("End Frequencies vs. Time");
set(gca,'YDir','normal');

fprintf("Computing averages...\n");
% Find the averages for each frequency section
begFreqAvgs = mean(begFiltAmps, 2);
endFreqAvgs = mean(endFiltAmps, 2);

% Plot the averages and the peaks
subplot(2,2,2);
stem(f, begFreqAvgs);
title("Beginning Frequency Averages");
xlim([0 0.5e4]);
[pksBeg, locsBeg] = findpeaks(begFreqAvgs, MinPeakProminence=MIN_FREQ_PROMINENCE);


subplot(2,2,4);
stem(f, endFreqAvgs);
title("End Frequency Averages");
xlim([0 0.5e4]);
[pksEnd, locsEnd] = findpeaks(endFreqAvgs, MinPeakProminence=MIN_FREQ_PROMINENCE);

% Plot the highest prominence peaks on charts for reference

figure(14); clf(14);
pltB = subplot(2,1,1);
Ampa = begFreqAvgs(locsBeg);
stem(f(locsBeg), Ampa);
title("Highest Prominence Beginning Peaks");

pltE = subplot(2,1,2);
Ampr = endFreqAvgs(locsEnd);
stem(f(locsEnd), Ampr);
title("Highest Prominence End Peaks");

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
plot(mps2mph(sourceV), ".b");
ylabel("Velocity (MPH)");
title("All Estimated Velocities");

% Remove negative velocities (non-physical)
possibleI = find(sourceV>=0);
sourceV = sourceV(possibleI);

% Plot all possible velocities
sortMphVels = mps2mph(sort(sourceV));
figure(16); clf(16);
% stem(mps2mph(sourceV), "b");
plot(sortMphVels, ".b");
ylabel("Velocity (MPH)");
title("All Positive Velocity Estimates");
% saveas(gcf, "50mph.png");

save("yetanothertest.mat");

% EXTRACTING VELOCITY FROM THE CURVE
% fprintf("Finding differences in possible velocities...\n");
% % Find differences and peaks
% figure(17); clf(17);
% diffVector = diff(sortMphVels);
% [peakDiffs,peakDiffLocs] = findpeaks(diffVector, ...
%     MinPeakProminence=MIN_DIFF_PROMINENCE);
% %plot
% findpeaks(diffVector, ...
%     MinPeakProminence=MIN_DIFF_PROMINENCE);
% title("Velocity Estimate Difference Plot")
% 
% fprintf("Estimating mean velocities...\n");
% % two point mean estimate
% format short g
% beforeVJumps = sortMphVels(peakDiffLocs);
% afterVJumps = sortMphVels(peakDiffLocs+1);
% 
% dif = beforeVJumps(1);
% meanV = 0 + dif/2;
% for i = 1:size(beforeVJumps,2)-1
%     dif(i+1) = beforeVJumps(i+1)-afterVJumps(i);
%     meanV(i+1) = afterVJumps(i) + dif(i+1)/2;
% end
% 
% %find horizontal span of the nonpeak regions
% differenceArray = [0 peakDiffLocs];
% for i = 1:size(differenceArray,2)-1
%     horizDiffs(i) = differenceArray(i+1)-differenceArray(i);
% end
% 
% pairs = [meanV; horizDiffs]';

% fprintf("Performing cutoff at %f mph...\n\n", CUTOFF_VELOCITY_MPH);
% % Perform cutoff checks
% indVec = pairs(:,1) < CUTOFF_VELOCITY_MPH;
% pairs = [pairs(indVec, 1) pairs(indVec, 2)];

% Place additional arbitrary limitations here.

% pairsSort = flip(sortrows(pairs, 2));
% figure(16); hold on; yline(pairsSort(1,1), "m");
% if ~isnan(TRUE_V_MPH)
%      yline(TRUE_V_MPH, "r");
%      legend("","Estimated V", "True V", Location="northwest");
% else
%     legend("","Estimated V", Location="northwest");
% end


% fprintf("The most likely velocity of this object is %f mph.\n" + ...
%     "\t(%d contiguous estimates, cutoff of %f mph.\n", ...
%     pairsSort(1,1), pairsSort(1,2),  CUTOFF_VELOCITY_MPH);
% if ~isnan(TRUE_V_MPH)
%     fprintf("This is %f mph ", abs(TRUE_V_MPH - pairsSort(1,1)))
%     if TRUE_V_MPH - pairsSort(1,1) <= 0
%         fprintf("above the true velocity, %f mph.\n", TRUE_V_MPH);
%     else
%         fprintf("below the true velocity, %f mph.\n", TRUE_V_MPH);
%     end
%     fprintf("Possible velocities, contiguous estimates, and distance" + ...
%         " from the true velocity:\n");
%     disp([pairsSort TRUE_V_MPH-pairsSort(:,1)]);
% else
%     fprintf("Possible velocities and contiguous estimates:\n");
%     disp(pairsSort);
% end