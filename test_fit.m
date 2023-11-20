clc; clear all;  figure(2); clf(2); figure(1); clf(1);
load("testvars.mat");

c = "default"; colormap(c);
imagesc(t, f, ampStft);
hold on;
set(gca,'YDir','normal');
xlabel("Time (s)");
ylabel("Frequency (Hz)");
axis([0, max(t), 0, 1000]);

BW = imbinarize(ampStft, "adaptive");

figure(2);
c = "default"; colormap(c);
imagesc(t, f, BW);
hold on;
set(gca,'YDir','normal');
xlabel("Time (s)");
ylabel("Frequency (Hz)");
axis([0, max(t), 0, 1000]);

xline(beginT, "r", LineWidth=1);
xline(endT, "r", LineWidth=1);