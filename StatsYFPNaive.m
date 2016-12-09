%function [meanPeriodYFPNaive, stdPeriodYFPNaive, meanAmpYFPNaive] = StatsYFPNaive()

dataYFPNaive = dlmread('yfpwithoutsponge.txt');

[maxdat, maxidx] = extrema(smooth(dataYFPNaive(:, 2),50));
maxima = maxdat(maxdat > 1);
maxima_idx = generation(maxidx(maxdat > 1));

%Sort the maxima by the order of when they occur, not how high they are.
%Thus, sort by the maxima indices.
[~, j] = sort(maxima_idx);
maxIndices = maxima_idx(j);
maxProtein = maxima(j);

period = [];
for mm = 1 : (length(maxIndices) - 1)
    period(mm) = maxIndices(mm + 1) - maxIndices(mm);
end

%discard transients
period = period(period > 5);

meanPeriodYFPNaive = mean(period);
stdPeriodYFPNaive = std(period);
meanAmpYFPNaive = mean(maxdat);
