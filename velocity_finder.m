clc; clear all;
matstring = "alternateanalysisroute.mat";
% matstring = "alternateanalysisroute56mph44F.mat";

load(matstring);

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

MIN_PROMINENCE = 0.5;

figure(3); clf(3);
pltB = subplot(2,1,1);
stem(f, begFreqAvgs);
xlim([0 0.5e4]);
[pksBeg, locsBeg] = findpeaks(begFreqAvgs, MinPeakProminence=MIN_PROMINENCE);


pltE = subplot(2,1,2);
stem(f, endFreqAvgs);
xlim([0 0.5e4]);
[pksEnd, locsEnd] = findpeaks(endFreqAvgs, MinPeakProminence=MIN_PROMINENCE);

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

% Plot all (including negative) velocities
figure(6); clf(6);
stem(mps2mph(sourceV), "b");
hold on; yline(50, "r");
ylabel("Velocity (MPH)");
title(matstring);

% Remove negative velocities (non-physical)
possibleI = find(sourceV>=0);
sourceV = sourceV(possibleI);

% Plot all possible velocities
figure(7); clf(7);
stem(mps2mph(sourceV), "b");
stem(mps2mph(sort(sourceV)), "b");
hold on; yline(56, "r");
ylabel("Velocity (MPH)");
title(matstring);
% saveas(gcf, "50mph.png");
