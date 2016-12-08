function [maxIndices, maxProtein] = FindProteinMaxima(proteinTrajectory)
%find maxima of peaks of protein levels for a given trajectory
[maxdat, maxidx] = extrema(smooth(proteinTrajectory, 750));
maxima = maxdat_tetR(maxdat > mean(maxdat));
maxima_idx = generation(maxidx(maxdat > mean(maxdat)));

%Sort the maxima by the order of when they occur, not how high they are.
%Thus, sort by the maxima indices.
[~, j] = sort(maxima_idx);
maxIndices = maxima_idx(j);
maxProtein = maxima(j);