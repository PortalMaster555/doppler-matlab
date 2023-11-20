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
