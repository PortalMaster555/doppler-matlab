clc; clear all;
load("alternateanalysisroute.mat");

begFiltAmps = filteredAmps(:,begRange);
endFiltAmps = filteredAmps(:,endRange);

% Plot beginning frequencies
figure(1); clf(1);

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

MIN_PROMINENCE = 1;

figure(3); clf(3);
pltB = subplot(2,1,1);
stem(f, begFreqAvgs);
xlim([0 0.5e4]);
[pksBeg, locsBeg] = findpeaks(begFreqAvgs, MinPeakProminence=MIN_PROMINENCE);
Fa = begFreqAvgs(locsBeg);

pltE = subplot(2,1,2);
stem(f, endFreqAvgs);
xlim([0 0.5e4]);
[pksEnd, locsEnd] = findpeaks(endFreqAvgs, MinPeakProminence=MIN_PROMINENCE);
Fr = endFreqAvgs(locsEnd);

% Plot the highest prominence peaks on charts for reference

figure(4); clf(4);
pltB = subplot(2,1,1);
findpeaks(begFreqAvgs, MinPeakProminence=MIN_PROMINENCE);

pltE = subplot(2,1,2);
findpeaks(endFreqAvgs, MinPeakProminence=MIN_PROMINENCE);

figure(5); clf(5);
pltB = subplot(2,1,1);
stem(f(locsBeg), Fa);
pltE = subplot(2,1,2);
stem(f(locsBeg), Fr);

% Construct pairs for each combination.

% columns are F_a
[FaM, FrM] = ndgrid(Fa, Fr);
% straighten out into columns and merge into an Nx2 matrix
FaM = FaM(:);
FrM = FrM(:);
inputMatrix = [FaM FrM];

% Calculate speed of sound based on temperature
if isnumeric(TEMPERATURE)
     c = vSound(TEMPERATURE);
else
     c = vSound();
end

% Find velocities for all combinations

for i = 1:size(inputMatrix, 1)
    options = optimset('Display','off');
    velFcn = @(v) (c-v)*inputMatrix(i,1) - (c+v)*inputMatrix(i,2);
    sourceV(i) = fsolve(velFcn, 0, options);
end

% conversion functions
mps2mph = @(v) v * 2.237;
mph2mps = @(v) v / 2.237;

% Remove negative velocities (non-physical)
possibleI = find(sourceV>=0);
sourceV = sourceV(possibleI);

% Plot all possible velocities
figure(6);
bar(mps2mph(sort(sourceV)), "b");
ylabel("Velocity (MPH)");

% debug 

hold on; yline(50, "r");


