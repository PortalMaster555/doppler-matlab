clc; clear all;
load("test_highestf.mat");
% Generate sample signal


% Approaching = negative vSource; Receding = positive vSource
fObserver = @(vSource, vObserver, vSound, fSource) ...
    fSource*(vSound+vObserver)/(vSound+vSource);

fSource = 440; %Hz
trueV = 20; %m/s when receding
c = vSound();
app = fObserver(-trueV, 0, c, fSource);
rec = fObserver(trueV, 0, c, fSource);
% SUBSTITUTE MEASURED app, rec VALUES HERE


% only works for stationary observer in line with the vehicle
% technically there should be a "radial" term for v

options = optimset('Display','off');
for i = 1:size(highestFBeg, 1)
    app(i) = highestFBeg(i);
    rec(i) = highestFEnd(i);
    velFcn = @(v) (c-v)*app(i) - (c+v)*rec(i);
    sourceV(i) = fsolve(velFcn, 0, options);
    fprintf("Approach: %f Hz\nRecede: %f Hz\nVel Est.: %f m/s\n\n", ...
    app(i), rec(i), sourceV(i));
end

%finds the index with the minimal difference in A
minDiffInd = find(abs(diff(sourceV))==min(abs(diff(sourceV))));
%extract this index, and it's neighbor index from A
v1 = sourceV(minDiffInd);
v2 = sourceV(minDiffInd+1);

avgClosestV = mean([v1, v2]);
% Convert to MPH from m/s
avgClosestV = 2.237 * avgClosestV;

fprintf("Estimated velocity: %f mph\n", avgClosestV);

%vSList();

% FUNCTIONS
 
function vS = vSound(TF)
% Find speed of sound in air based on temperature if provided
% http://hyperphysics.phy-astr.gsu.edu/hbase/Sound/souspe3.html
  if nargin < 1
    vS = 343; % m/s
  else
    TK = @(TF) ((TF-32) * 5/9) + 273.15;
    vS = 20.05 * sqrt(TK(TF));
  end
end

function vSList()
    for tempF = -459.67:20:100
        fprintf("%f m/s at %f degrees Fahrenheit\n", vSound(tempF), tempF);
    end
    fprintf("Default vSound() output: %f\n", vSound());
end

