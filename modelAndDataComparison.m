function modelAndDataComparison(proteinData, generationNaive, naiveProteinModel, generationImproved, improvedProteinModel)
%Function returns models of paperData, naiveProteinModel, and
%improvedProtein Model superimposed

readIn = dlmread(proteinData);

%for a single protein align each model type and the data
[genNaive, proteinNaive] = findProteinModelMaxima(readIn, generationNaive, naiveProteinModel);
[genImproved, proteinImproved] = findProteinModelMaxima(readIn, generationImproved, improvedProteinModel);

%ttest each model using code that I submitted (sponge v nonsponge)

%[is_sig, pval] = kstest2(xx, yy)

 
%plot(CFPSponge(:, 1), CFPSponge(:, 2), 'r.'); hold on
%plot MODELS HERE

%legend('CFP with Sponge', 'Model');

figure; 
plot(readIn(:, 1), readIn(:, 2), 'r.-'); hold on;
plot(genNaive, proteinNaive,'b.-');
plot(genImproved, proteinImproved, 'm.-'); hold off;
