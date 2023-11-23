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
% Plot the averages

figure(3); clf(3);
pltB = subplot(2,1,1);
stem(f, begFreqAvgs);
xlim([0 0.5e4]);
findpeaks(begFreqAvgs);

pltE = subplot(2,1,2);
stem(f, endFreqAvgs);
xlim([0 0.5e4]);
findpeaks(endFreqAvgs);