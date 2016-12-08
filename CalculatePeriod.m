function period = CalculatePeriod(generation, proteinTrajectory)

%returns periods to models (NaiveRepressilator or ImprovedRepressilator)
[maxIndices, ~] = FindProteinMaxima(generation, proteinTrajectory);

period = [];
for mm = 1 : (length(maxIndices) - 1)
    period(mm) = maxIndices(mm + 1) - maxIndices(mm);
end

%discard transients
period = period(period > 5);
