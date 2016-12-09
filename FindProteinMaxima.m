function [maxIndices, maxProtein] = FindProteinMaxima(generation, proteinTrajectory, SpecifyDataOrModel)

%find maxima of peaks of protein levels for a given trajectory
if strcmp(SpecifyDataOrModel,'model') == 1
    [maxdat, maxidx] = extrema(smooth(proteinTrajectory, 750));
elseif strcmp(SpecifyDataOrModel,'data') == 1
    [maxdat, maxidx] = extrema(proteinTrajectory);
end
maxima = maxdat(maxdat > mean(maxdat));
maxima_idx = generation(maxidx(maxdat > mean(maxdat)));

%Sort the maxima by the order of when they occur, not how high they are.
%Thus, sort by the maxima indices.
[~, j] = sort(maxima_idx);
maxIndices = maxima_idx(j);
maxProtein = maxima(j);
