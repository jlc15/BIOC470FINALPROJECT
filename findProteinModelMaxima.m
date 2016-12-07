function [genAligned, proteinModelAligned] = findProteinModelMaxima(proteinData, generation, proteinModel)
%returns maxima for all timepts for a given protein model
%maybe discard transients from model here??
%shift data
shiftedProteinData(:, 2) = proteinData(:, 2) - mean(proteinData(:, 2));

%Note: didn't smooth since reverse-engineered data are single pts
[~, iDataMax] = extrema(shiftedProteinData);
proteinFirstMax = min(iDataMax);

%find maxima of peaks of protein levels in model
[maxDataModel, maxIdxModel] = extrema(smooth(proteinModel, 750));

%remove mins
maxImaModel = maxDataModel(maxDataModel > mean(maxDataModel));
%return indices of max points in model
maxImaIdxModel = generation(maxIdxModel(maxDataModel > mean(maxDataModel)));

%Sort the maxima by index not amplitude
[~, j] = sort(maxImaIdxModel);
maxProteinIdx = maxImaIdxModel(j);

modelFirstMaxAfterFirstDataMax = maxProteinIdx(maxProteinIdx > proteinFirstMax);
firstMaxGreaterThanDataMax = min(modelFirstMaxAfterFirstDataMax);

%line up maxima
ndat = length(proteinData(:, 1));

firstpoint = firstMaxGreaterThanDataMax - proteinFirstMax + 1;

lastpoint = firstMaxGreaterThanDataMax - proteinFirstMax + ndat;

%get model points corresponding to data points to line up
proteinModelAligned = proteinModel(firstpoint : lastpoint);
modelGenAligned = generation(firstpoint : lastpoint);
genAligned = modelGenAligned - modelGenAligned(1);