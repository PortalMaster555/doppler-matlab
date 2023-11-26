clear all;
load("yetanothertest.mat");
MIN_DIFF_PROMINENCE = 1;
figure(16); clf(16);
% stem(mps2mph(sourceV), "b");
plot(sortMphVels, ".b");
ylabel("Velocity (MPH)");
title("All Positive Velocity Estimates");

figure(21); clf(21);
gradVels = gradient(sortMphVels);
plot(gradVels, "-r");
title("Gradient of the possible velocity curve");
ylabel("Point-to-Point Change in Velocity");

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
yline(sortMphVels(gradPeakLocs));

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

% First point in the region
beginningIndex = sum(regionLength(1:lowestRegionIndex-1)) + 1;
% Last point in the region
endingIndex = sum(regionLength(1:lowestRegionIndex));
xline(beginningIndex, "r");
xline(endingIndex, "r");

velEstimates = sortMphVels(beginningIndex:endingIndex);

finalVelocityEstimate = mean(velEstimates);
disp(finalVelocityEstimate);













