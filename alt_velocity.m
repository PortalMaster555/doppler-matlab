load("yetanothertest.mat");

figure(16); clf(16);
% stem(mps2mph(sourceV), "b");
plot(sortMphVels, ".b");
ylabel("Velocity (MPH)");
title("All Positive Velocity Estimates");