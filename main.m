clear all;
audioFile = "./audio/50mphobserver.wav";
% audioFile = "horn.ogg";
% audioFile = "test_10mps.wav";
% audioFile = "test_440_48khz.wav";
% audioFile = "56mph44F.wav";

clc; printSplash();

% CONSTANTS (tweak as needed)
% NUM_TO_TAKE = 1; % decrease if the audio has fewer overtones or if the 
                    % window size is very large

WINDOW_SIZE = 2048; %larger seems to work better, but the changepoint 
                    % calculation hangs sometimes
TEMPERATURE = "UNKNOWN";
% TEMPERATURE = 44; %F
% smoothingConstant = 10;



[Amps, Fs] = audioread(audioFile);
N = size(Amps,1); % number of samples
fprintf("Audio file loaded: %s (length %d samples)\n", audioFile, N);
% 
% fprintf("Playing sample of audio...\n")
% audio = audioplayer(Amps, Fs); 
% play(audio); 
% % pause(2); stop(audio);
% fprintf("Playback complete.\n")

fprintf("Maximum detectable frequency at or just below  %f Hz\n", Fs/2);

% Compute the short-time Fourier transform
% Windowing function: length 2048 Hann
% One-sided frequency range (absolute amplitude not preserved)
[stfourier, f, t] = stft(Amps, Fs, FrequencyRange="onesided", ...
    Window=hann(WINDOW_SIZE));
ampStft = abs(stfourier);

deltaT = t(end)-t(end-1);
deltaF = f(end)-f(end-1);

% Plot vertical lines at largest changes
changeIndices = findchangepts(ampStft,MaxNumChanges=2,Statistic="rms");
if size(changeIndices, 2) < 2
    %fixes a bug where it doesn't find a second change
    changeIndices(2) = changeIndices(1);
end
beginT = changeIndices(1)*deltaT;
endT = changeIndices(2)*deltaT;

% Plot spectrogram
fprintf("Plotting spectrogram...\n")
figure(1); clf(1);
imagesc(t, f, ampStft); hold on;
set(gca,'YDir','normal');
xlabel("Time (s)"); ylabel("Frequency (Hz)");
cb = colorbar; cb.TicksMode = "manual";
ylabel(cb,'Amplitude', Rotation=270);
xline(beginT, "r", LineWidth=1);
xline(endT, "r", LineWidth=1);

figure(2); clf(2);
surf(t, f, ampStft); %plot frequency breakdown in 3d
shading interp

freqHigh = 2000;
% filteredAmps = highpass(ampStft,freqHigh,Fs);
filteredAmps = ampStft;

figure(3); clf(3);
% fprintf("Plotting smoothed spectrogram...\n")
fprintf("Plotting spectrogram...\n")
c = "default"; colormap(c);
imagesc(t, f, filteredAmps);
hold on;
set(gca,'YDir','normal');
xlabel("Time (s)");
ylabel("Frequency (Hz)");
axis([0, max(t), 0, 5000]);
colorbar;
xline(beginT, "r", LineWidth=1);
xline(endT, "r", LineWidth=1);

fprintf("Detecting edges...\n")
edges = edge(filteredAmps, 'Canny');

% Filter
% filteredAmps = filterAudio(ampStft);

% Isolate ends (outside of the middle transition)
begRange = 1:changeIndices(1);
endRange = changeIndices(2):size(edges,2);
begAmps = edges(:, begRange);
endAmps = edges(:, endRange);
tbeg = t(begRange); tend = t(endRange);
%imagesc(tbeg, f, begAmps);
%imagesc(tend, f, endAmps);

save("alternateanalysisroute");

% Weirdest part of the data analysis

% Calculate mean of each frequency band across time
% outside of the middle region to find most constant frequencies
% frequencies across time (index 2)
fprintf("Calculating means of all frequencies...\n")
fBegSums = mean(begAmps, 2);
fEndSums = mean(endAmps, 2);
figure(4); clf(4);
stem(f, fBegSums, ".r");
hold on;
stem(f, fEndSums, ".b");
axis([0, 4000, 0, max(fBegSums)]);


for i = 1:20
    avgClosestV(i) = mainfcn(i, fBegSums, fEndSums, f, TEMPERATURE);
end
figure(5); clf(5);
stem(avgClosestV, ".k");
title(audioFile);
xlabel("Points chosen"); ylabel("Speed (MPH)");
text(1:length(avgClosestV),avgClosestV,num2str(avgClosestV'), ...
    'vert','bottom','horiz','center'); 
%savefig(strcat(audioFile,'.fig'));


function avgClosestV = mainfcn(NUM_TO_TAKE, fBegSums, fEndSums, f, ...
    TEMPERATURE)
% Take highest-amp frequencies (constant in beginning of file)

fprintf("Choosing and plotting highest %d means...\n", NUM_TO_TAKE);

[maxAmpsB, begMaxFIndices] = maxk(fBegSums, NUM_TO_TAKE);
highestFBeg = f(begMaxFIndices);
plot(highestFBeg, maxAmpsB, "om", LineStyle="none");
highestFBeg = sort(highestFBeg);

[maxAmpsE, endMaxFIndices] = maxk(fEndSums, NUM_TO_TAKE);
highestFEnd = f(endMaxFIndices);
plot(highestFEnd, maxAmpsE, "om",LineStyle="none");
highestFEnd = sort(highestFEnd);

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

% Calculate speed of sound based on temperature
if isnumeric(TEMPERATURE)
     c = vSound(TEMPERATURE);
else
     c = vSound();
end

fprintf("Speed of sound taken as %f m/s.\n", c);
fprintf("Analyzing velocity based on highest %d means...\n", NUM_TO_TAKE);
options = optimset('Display','off');
for i = 1:size(highestFBeg, 1)
    velFcn = @(v) (c-v)*app(i) - (c+v)*rec(i);
    sourceV(i) = fsolve(velFcn, 0, options);
    %fprintf("Approach: %f Hz\nRecede: %f Hz\nVel Est.: %f m/s\n\n", ...
    %app(i), rec(i), sourceV(i));
end

if size(sourceV,2) >= 2
    fprintf("Finding closest two velocities from the %d and averaging...\n", ...
        NUM_TO_TAKE);
    minDiffInd = find(abs(diff(sourceV))==min(abs(diff(sourceV))));
    % extract this index and its neighbor index
    v1 = sourceV(minDiffInd);
    v2 = sourceV(minDiffInd+1);
    avgClosestV = mean([v1, v2]);
else
    avgClosestV = sourceV;
end

fprintf("Estimated velocity: %f m/s.\n", avgClosestV);

% Convert to MPH from m/s
avgClosestV = 2.237 * avgClosestV;
fprintf("Estimated velocity: %f mph.\n", avgClosestV);

end