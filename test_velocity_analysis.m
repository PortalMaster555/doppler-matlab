clc; clear all;

% Generate sample signal


% Approaching = positive vSource; Receding = positive vSource
fObserver = @(vSource, vObserver, vSound, fSource) ...
    fSource*(vSound+vObserver)/(vSound+vSource);

fSource = 440; %Hz
app = fObserver(-10, 0, 343, fSource);
rec = fObserver(10, 0, 343, fSource);
avg = mean([app, rec]);

fprintf("Approach: %f\nRecede: %f\nAvg: %f\nTrue: %f\n", ...
    app, rec, avg, fSource);
vSList();

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