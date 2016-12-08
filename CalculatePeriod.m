function period = CalculatePeriod(generation, proteinTrajectory)
%returns stats analyses

[maxIndices, maxProtein] = FindProteinMaxima(generation, proteinTrajectory);

period = zeros(1 : (length(maxIndices) - 1));
for mm = 1 : (length(maxIndices) - 1)
    period(mm) = maxIndices(mm + 1) - maxIndices(mm);
end

%discard transients
period = period(period > 5);
