clc; clear all;
trueV = NaN;
% matstring = "alternateanalysisroute.mat"; trueV = 50;
matstring = "alternateanalysisroute56mph44F.mat"; trueV = 56;
load(matstring);

% conversion functions
mps2mph = @(v) v * 2.237;
mph2mps = @(v) v / 2.237;

% CONSTANTS
% Calculate speed of sound based on temperature
if isnumeric(TEMPERATURE)
     c = vSound(TEMPERATURE);
else
     c = vSound();
end

% Doppler shift not detectable in supersonic objects
CUTOFF_VELOCITY_MPH = mps2mph(c); 

% Additional cutoff if results make no sense
CUTOFF_VELOCITY_MPH = 150; % reasonable for automobiles

% Minimum mean frequency amplitude prominence to be considered a peak
% E.g. 0.5
MIN_FREQ_PROMINENCE = 0.5;

% Minimum jump in velocity estimate curve to be considered a peak
% E.g. 5
MIN_DIFF_PROMINENCE = 10; %mph


% BEGIN 

% Plot beginning frequencies
figure(1); clf(1);
begFiltAmps = filteredAmps(:,begRange);
endFiltAmps = filteredAmps(:,endRange);
imagesc(tbeg, f, begFiltAmps);
ylim([0 4000]);
set(gca,'YDir','normal');

% Plot end frequencies
figure(2); clf(2);
imagesc(tbeg, f, endFiltAmps);
ylim([0 4000]);
set(gca,'YDir','normal');

% Find the averages for each frequency section
begFreqAvgs = mean(begFiltAmps, 2);
endFreqAvgs = mean(endFiltAmps, 2);

% Plot the averages and the peaks

figure(3); clf(3);
pltB = subplot(2,1,1);
stem(f, begFreqAvgs);
xlim([0 0.5e4]);
[pksBeg, locsBeg] = findpeaks(begFreqAvgs, MinPeakProminence=MIN_FREQ_PROMINENCE);


pltE = subplot(2,1,2);
stem(f, endFreqAvgs);
xlim([0 0.5e4]);
[pksEnd, locsEnd] = findpeaks(endFreqAvgs, MinPeakProminence=MIN_FREQ_PROMINENCE);

% Plot the highest prominence peaks on charts for reference

figure(4); clf(4);
pltB = subplot(2,1,1);
Ampa = begFreqAvgs(locsBeg);
stem(f(locsBeg), Ampa)

pltE = subplot(2,1,2);
Ampr = endFreqAvgs(locsEnd);
stem(f(locsEnd), Ampr)

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

for i = 1:size(inputMatrix, 1)
    options = optimset('Display','off');
    velFcn = @(v) (c-v)*inputMatrix(i,1) - (c+v)*inputMatrix(i,2);
    sourceV(i) = fsolve(velFcn, 0, options);
end


% Plot all (including negative) velocities
figure(5); clf(5);
stem(mps2mph(sourceV), "b");
hold on; yline(50, "r");
ylabel("Velocity (MPH)");
title(matstring);

% Remove negative velocities (non-physical)
possibleI = find(sourceV>=0);
sourceV = sourceV(possibleI);

% Plot all possible velocities
sortMphVels = mps2mph(sort(sourceV));
figure(6); clf(6);
% stem(mps2mph(sourceV), "b");
plot(sortMphVels, ".b");
ylabel("Velocity (MPH)");
title(matstring);
% saveas(gcf, "50mph.png");


% EXTRACTING VELOCITY FROM THE CURVE

% Find differences and peaks
figure(7); clf(7);
diffVector = diff(sortMphVels);
[peakDiffs,peakDiffLocs] = findpeaks(diffVector, ...
    MinPeakProminence=MIN_DIFF_PROMINENCE);
%plot
findpeaks(diffVector, ...
    MinPeakProminence=MIN_DIFF_PROMINENCE);

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

% Perform cutoff checks
indVec = pairs(:,1) < CUTOFF_VELOCITY_MPH;
pairs = [pairs(indVec, 1) pairs(indVec, 2)];

% Place additional arbitrary limitations here.

pairsSort = flip(sortrows(pairs, 2));
figure(6); hold on; yline(pairsSort(1,1), "m"); yline(trueV, "r");
legend("","Estimated V", "True V", Location="northwest");

fprintf("The most likely velocity of this object is %f mph.\n" + ...
    "\t(%d contiguous entries, cutoff of %f mph.\n", ...
    pairsSort(1,1), pairsSort(1,2),  CUTOFF_VELOCITY_MPH);
if ~isnan(trueV)
    fprintf("This is %f mph ", abs(trueV - pairsSort(1,1)))
    if trueV - pairsSort(1,1) <= 0
        fprintf("above the true velocity, %f mph.\n", trueV);
    else
        fprintf("below the true velocity, %f mph.\n", trueV);
    end
    fprintf("Possible velocities, contiguous entries, and distance from trueV:\n");
    disp([pairsSort trueV-pairsSort(:,1)]);
else
    fprintf("Possible velocities and contiguous entries:\n");
    disp([pairsSort]);
end