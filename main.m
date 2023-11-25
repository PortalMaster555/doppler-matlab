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

% Minimum mean frequency amplitude prominence to be considered a peak
% E.g. 0.5
MIN_FREQ_PROMINENCE = 0.5;

% Minimum jump in velocity estimate curve to be considered a peak
% E.g. 5 (may need to be tweaked based on available sample data)
MIN_DIFF_PROMINENCE = 1; %mph

% Set to 1 to interpolate, 0 to use original velocity estimate curve.
interp_flag = 1;
% Multiplies number of points. Be careful for files with many tones!
% Minimum of 1 should be used.
INTERPOLATION_VALUE = 1;

% Set to 1 to smooth, 0 to use original gradient curve.
smooth_flag = 0;


% audioFile = "./audio/onecar/70kph.wav"; TRUE_V_MPH = kph2mph(70);
% audioFile = "./audio/onecar/50kph.wav"; TRUE_V_MPH = kph2mph(50);
    % Above seems to imply that the 50 is the wrong value (64kph instead???)
audioFile = "./audio/onecar/30kph.wav"; TRUE_V_MPH = kph2mph(30);

% audioFile = "./audio/appx40to50.wav"; TRUE_V_MPH = 45;
    % Source estimate was approximately 40 to 50 mph, so acceptable.
    % Line drawn in the middle of the range

% Requires cutoff and interpolation
% audioFile = "./audio/50mphobserver.wav"; TRUE_V_MPH = 50; 
% INTERPOLATION_VALUE = 10;

% Slightly inaccurate; overtones actually help the program
% audioFile = "./audio/test_10mps.wav"; TRUE_V_MPH = mps2mph(10);
% interp_flag = 0;

% I would have been surprised if this were anything but 0.
% audioFile = "./audio/test_440_48khz.wav"; TRUE_V_MPH = 0;


% BROKEN
% Incredibly low number of velocity estimates;
% audioFile = "./audio/56mph44F.wav"; 
% TRUE_V_MPH = 56; TEMPERATURE_F = 44; INTERPOLATION_VALUE = 10;

% Likely requires a lower limit (true velocity not 0, may be ~20 mph)
% audioFile = "./audio/horn.ogg";


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
CUTOFF_VELOCITY_MPH = 150;


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

% sort list of velocities
sortMphVels = mps2mph(sort(sourceV));
figure(16); clf(16);
% stem(mps2mph(sourceV), "b");
plot(sortMphVels, ".b");
ylabel("Velocity (MPH)");
title("All Positive Velocity Estimates");

if size(sortMphVels, 2) == 1
    finalVelEst = sortMphVels(1,1);
else
    if interp_flag
        % defines additional precision
        step = 1/INTERPOLATION_VALUE;
        xq = 0:step:size(sortMphVels, 2);
        sortMphVels = interp1(sortMphVels, xq);
        plot(sortMphVels, ".r");
    end
    
    % saveas(gcf, "50mph.png");
    
    % Apply cutoff (does this need to happen later?)
    sortMphVels = sortMphVels(sortMphVels < CUTOFF_VELOCITY_MPH);
    hold on;
    
    % EXTRACTING VELOCITY FROM THE CURVE
    
    fprintf("Calculating gradient of velocity estimate curve...\n");
    figure(21); clf(21);
    gradVels = gradient(sortMphVels);
    
    plot(gradVels, "--r");
    hold on;
    title("Gradient of the possible velocity curve");
    ylabel("Point-to-Point Change in Velocity");
    
    % if smooth_flag
    %     gradVels = smoothdata(gradVels);
    %     plot(gradVels, "-b");
    % end
    
    fprintf("Calculating peaks of gradient of velocity estimate curve...\n");
    figure(22); clf(22);
    [gradPeaks, gradPeakLocs] = findpeaks(gradVels, ...
        MinPeakProminence=MIN_DIFF_PROMINENCE);
    % Plot
    findpeaks(gradVels, MinPeakProminence=MIN_DIFF_PROMINENCE);
    title("Peaks of the gradient of the possible velocity curve");
    
    % Idea: Sum of the points in the regions between peaks
    %       should be close to 0 for the region to be FLAT.
    
    
    
    figure(16);
    hold on;
    yline(sortMphVels(gradPeakLocs), "--k");
    
    fprintf("Calculating flattest region on curve...\n");
    % The most overcomplicated method for what is probably just a mean ever
    regionCounter = 1;
    regionSum = 0;
    regionLength = 1;
    for i = 2:size(sortMphVels, 2)
        % If still part of the current region
        if ~any(i == gradPeakLocs)
            regionSum(regionCounter) = regionSum(regionCounter) + gradVels(i);
            regionLength(regionCounter) = regionLength(regionCounter) + 1;
        else %if peak (new region)
            regionCounter = regionCounter + 1;
            % for some reason it doesn't want to expand the array itself
            regionSum = [regionSum 0];
            regionLength(regionCounter) = 1;
            % add first peak to the new region
            regionSum(regionCounter) = regionSum(regionCounter) + gradVels(i);
        end
    end
    
    regionPairs = [regionSum; regionLength]';
    % Effectively a measure of flatness per length?
    regionMeans = regionSum ./ regionLength;
    [lowestMean, lowestRegionIndex] = min(regionMeans);
    textLocs = [0 gradPeakLocs];
    % text(textLocs, regionSum(textLocs));
    
    fprintf("Finding all velocity estimates in region...\n");
    % First point in the region
    beginningIndex = sum(regionLength(1:lowestRegionIndex-1)) + 1;
    % Last point in the region
    endingIndex = sum(regionLength(1:lowestRegionIndex));
    xline(beginningIndex, "--r", DisplayName="");
    xline(endingIndex, "--r", DisplayName="");
    
    % Find final estimate
    velEstimates = sortMphVels(beginningIndex:endingIndex);
    finalVelEst = mean(velEstimates);
end
% Draw some lines
estLine = yline(finalVelEst, "b");
cutLine = yline(CUTOFF_VELOCITY_MPH, "m");
l = legend([estLine, cutLine], "Estimate", "Cutoff");
l.Direction = "reverse";
l.Location = "northwest";

fprintf("\nFinal velocity estimate is %f mph. (%f km/h, %f m/s)\n" ...
    , finalVelEst, mph2kph(finalVelEst), mph2mps(finalVelEst));

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

if ~isnan(TRUE_V_MPH)
    fprintf("This is %f mph ", abs(TRUE_V_MPH - finalVelEst))
    if TRUE_V_MPH - finalVelEst <= 0
        fprintf("above the true velocity, %f mph.\n", TRUE_V_MPH);
    else
        fprintf("below the true velocity, %f mph.\n", TRUE_V_MPH);
    end
    fprintf("Percentage error: %f%%\n", 100*(finalVelEst-TRUE_V_MPH)/TRUE_V_MPH)
end