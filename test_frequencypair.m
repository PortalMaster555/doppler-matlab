load("test_freqpairbeg")

% Pair each beginning frequency with the closest end frequency below it
% Repeats are allowed

fullFArray = [highestFBeg highestFEnd];
% col 1, fBeg   col 2, fEnd

for i = 1:size(highestFBeg, 1)
    ind = find(highestFEnd < highestFBeg(i,1));
    fprintf("highestVal: %f\n", highestFBeg(i));
    fprintf("correspondingEnd: %f\n\n", highestFEnd(ind(end)));
    app(i) = highestFBeg(i);
    rec(i) = highestFEnd(ind(end));
end

c = vSound();
fprintf("Speed of sound taken as %f m/s.\n", c);
fprintf("Analyzing velocity based on highest %d means...\n", num_to_take);
options = optimset('Display','off');
for i = 1:size(highestFBeg, 1)
    velFcn = @(v) (c-v)*app(i) - (c+v)*rec(i);
    sourceV(i) = fsolve(velFcn, 0, options);
    %fprintf("Approach: %f Hz\nRecede: %f Hz\nVel Est.: %f m/s\n\n", ...
    %app(i), rec(i), sourceV(i));
end

fprintf("Finding closest two velocities from the %d and averaging...\n", ...
    num_to_take);
minDiffInd = find(abs(diff(sourceV))==min(abs(diff(sourceV))));
% extract this index and its neighbor index
v1 = sourceV(minDiffInd);
v2 = sourceV(minDiffInd+1);

avgClosestV = mean([v1, v2]);
fprintf("Estimated velocity: %f m/s.\n", avgClosestV);

% Convert to MPH from m/s
avgClosestV = 2.237 * avgClosestV;
fprintf("Estimated velocity: %f mph.\n", avgClosestV);