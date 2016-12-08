function meanAmplitude = MeanAmp(proteinTrajectory)
%finds mean amplitude of the trajectory of the protein

%get extrema from trajectory. Since this gives both max and min points, can average maxdat to find mean amplitude.
[maxdat, ~] = extrema(smooth(proteinTrajectory, 750));
meanAmplitude = mean(maxdat);