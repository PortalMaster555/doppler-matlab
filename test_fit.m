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
smoothAmps = movmean(ampStft', 100);
smoothAmps = smoothAmps';

figure(2); clf(2);
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
fBegSums = mean(begAmps, 2);
fEndSums = mean(endAmps, 2);
figure(3); clf(3);
stem(f, fBegSums, ".r");
hold on;
stem(f, fEndSums, ".b");
axis([0, 4000, 0, max(fBegSums)]);

% Take highest-amp frequencies
num_to_take = 3;

[maxAmpsB, begMaxFIndices] = maxk(fBegSums, num_to_take);
highestFBeg = f(begMaxFIndices);
plot(highestFBeg, maxAmpsB, "om", LineStyle="none");
highestFBeg = sort(highestFBeg);

[maxAmpsE, endMaxFIndices] = maxk(fEndSums, num_to_take);
highestFEnd = f(endMaxFIndices);
plot(highestFEnd, maxAmpsE, "om",LineStyle="none");
highestFEnd = sort(highestFEnd);

save("test_highestf");
















