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
stem(f(locsEnd), Fr);

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

% Plot all (including negative) velocities
figure(6); clf(6);
subplot(2,1,1);
stem(mps2mph(sourceV), "b");
hold on; yline(50, "r");
ylabel("Velocity (MPH)");
title(matstring);

% Remove negative velocities (non-physical)
possibleI = find(sourceV>=0);
sourceV = sourceV(possibleI);

% Plot all possible velocities
subplot(2,1,2);
stem(mps2mph(sourceV), "b");
stem(mps2mph(sort(sourceV)), "b");
hold on; yline(50, "r");
ylabel("Velocity (MPH)");
title(matstring);
% figure(7); clf(7);
% plot(diff(mps2mph(sort(sourceV))));

%
%
%

% Pair each beginning frequency with the closest end frequency below it
% Repeats are allowed

for i = 1:size(highestFBeg, 1)
    ind = find(highestFEnd <= highestFBeg(i,1));
    % fprintf("highestVal: %f\n", highestFBeg(i));
    % fprintf("correspondingEnd: %f\n\n", highestFEnd(ind(end)));
    if isempty(ind)
        %Fixes cases where there is no lower end value than a beginning
        %value (sets beginning and end equal, so v=0)
        app(i) = highestFBeg(i);
        rec(i) = highestFBeg(i);
    else
        app(i) = highestFBeg(i);
        rec(i) = highestFEnd(ind(end));
    end
end



