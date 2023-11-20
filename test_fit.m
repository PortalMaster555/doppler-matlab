clc; clear all;  figure(2); clf(2); figure(1); clf(1);
load("testvars50.mat");

% Draw a colormap of the unsmoothed audio spectrogram
c = "default"; colormap(c);
imagesc(t, f, ampStft);
hold on;
set(gca,'YDir','normal');
xlabel("Time (s)");
ylabel("Frequency (Hz)");
axis([0, max(t), 0, 5000]);
colorbar;

% Perform very aggressive smoothing on the image across time
smoothAmps = movmean(ampStft', 50);
smoothAmps = smoothAmps';

figure(2);
c = "default"; colormap(c);
%imagesc(t, f, smoothAmps);
hold on;
set(gca,'YDir','normal');
xlabel("Time (s)");
ylabel("Frequency (Hz)");
axis([0, max(t), 0, 5000]);
colorbar;
xline(beginT, "r", LineWidth=1);
xline(endT, "r", LineWidth=1);

% Calculate sum of each frequency band across time
% outside of the middle region
begRange = 1:changeIndices(1);
endRange = changeIndices(2):size(smoothAmps,2);
begAmps = smoothAmps(:, begRange);
endAmps = smoothAmps(:, endRange);
tbeg = t(begRange); tend = t(endRange);
imagesc(tbeg, f, begAmps);
imagesc(tend, f, endAmps);





















