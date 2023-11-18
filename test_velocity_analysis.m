clc; clear all;

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
velFcn = @(v) (c-v)*app - (c+v)*rec;
sourceV = fsolve(velFcn, 0);

fprintf("Approach: %f Hz\nRecede: %f Hz\nTrue: %f Hz\n", ...
    app, rec, fSource);
fprintf("True velocity: %f m/s\nEstimated velocity: %f m/s\n", ...
    trueV, sourceV);

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

