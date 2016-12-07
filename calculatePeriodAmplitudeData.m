function [mean_period, std_period, mean_amp] = calculatePeriodAmplitudeData(maxIndices, maxProtein)
%returns stats analyses
period = zeros(1 : (length(maxIndices) - 1));
for mm = 1 : length(maxIndices) - 1
    period(mm) = maxIndices(mm - 1) - maxIndices(mm);
end

%discard transients
period = period(period > 5);

%mean and standard deviation for comparison
mean_period = mean(period);
std_period = std(period);
%mean amplitude of oscillations is just:
mean_amp = mean(maxProtein);