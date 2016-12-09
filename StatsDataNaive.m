function [meanPeriod, stdPeriod, meanAmp] = StatsDataNaive(Filename)

if strcmp(Filename,'yfpwithoutsponge.txt') == 1
    dataYFPNaive = dlmread('yfpwithoutsponge.txt');
    generation = dataYFPNaive(:,1);

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

    meanPeriod = mean(period);
    stdPeriod = std(period);
    meanAmp = mean(maxdat);
    
elseif strcmp(Filename,'cfpwithoutsponge.txt') == 1
    dataCFPNaive = dlmread('cfpwithoutsponge.txt');
    generation = dataCFPNaive(:,1);

    [maxdat, maxidx] = extrema(smooth(dataCFPNaive(:, 2),50));
    maxima = maxdat(maxdat > 0.5);
    maxima_idx = generation(maxidx(maxdat > 0.5));

    %Sort the maxima by the order of when they occur, not how high they are.
    %Thus, sort by the maxima indices.
    [~, j] = sort(maxima_idx);
    maxIndices = maxima_idx(j);
    maxIndices = maxIndices(2:end);
    maxProtein = maxima(j);

    period = [];
    for mm = 1 : (length(maxIndices) - 1)
        period(mm) = maxIndices(mm + 1) - maxIndices(mm);
    end

    %discard transients
    period = period(period > 5);

    meanPeriod = mean(period);
    stdPeriod = std(period);
    meanAmp = mean(maxdat);
    
elseif strcmp(Filename,'rfpwithoutsponge.txt') == 1
    dataRFPNaive = dlmread('rfpwithoutsponge.txt');
    PeriodRFPNaive = CalculatePeriod(dataRFPNaive(:,1), dataRFPNaive(:, 2),'data');
    meanPeriod = mean(PeriodRFPNaive);
    stdPeriod = std(PeriodRFPNaive);
    meanAmp = MeanAmp(dataRFPNaive(:, 2));
end
