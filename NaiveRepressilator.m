%protein degradation rate = 1 (as noted in 4.1 Supplementary Info from 
%Synchronous long-term oscillations in a synthetic gene circuit paper)

%initialize p(i) as they are constant for both models
dataRFPNaive = dlmread('rfpwithoutsponge.txt');
dataYFPNaive = dlmread('yfpwithoutsponge.txt');
dataCFPNaive = dlmread('cfpwithoutsponge.txt');
%for plotting
periodTetR = [];
periodLambdaCl = [];
periodLacL = [];
ampStorageTetR = [];
ampStorageLambdaCl = [];
ampStorageLacL = [];
generationOutput = [];
tetROutput = [];
lambdaClOutput = [];
lacLOutput = [];

for ii = 1:100
    %parameters are empirically determined in paper
    %production rate of protein
    lambda = 2000; 
    %at half Vmax
    kTetR = 100; 
    kLambdaCl = 100; 
    kLacL = 100;
    hillCoeff = 4; 
    initTetR = 20;
    initLambdaCl = 20;
    initLacL = 20;
    
    tetR = 0;
    lambdaCl = 0;
    lacL = 0;
    generation = 0;

    generationMax = 1000;
    generation(1) = 0;

    %protein amounts are stored in tetR, lambdaCl, lacL
    tetR(1) = initTetR;
    lambdaCl(1) = initLambdaCl; 
    lacL(1) = initLacL;
    
    %counter
    currentTime = 2;

    while generation(currentTime - 1) < generationMax
        %In the form -> dp(tetR | lambdaCl | lacL)/dt = lambda times kTetR to the 
        %hillCoeff divided by kTetR to the hillCoeff + amount of protein that represses 
        %specific the protein to the hillCoeff all 
        rateup_TetR = ((lambda * (kTetR ^ hillCoeff)) / ( (kTetR ^ hillCoeff) + (lacL(currentTime - 1) ^ hillCoeff))); 

        ratedown_TetR = tetR(currentTime - 1);
    
        rateup_LambdaCl = ((lambda * (kLambdaCl ^ hillCoeff)) / ( (kLambdaCl ^ hillCoeff) + (tetR(currentTime - 1) ^ hillCoeff)));

        ratedown_LambdaCl = lambdaCl(currentTime - 1);
    
        rateup_LacL = ((lambda * (kLacL ^ hillCoeff)) / ( (kLacL ^ hillCoeff) + (lambdaCl(currentTime - 1) ^ hillCoeff))); 

        ratedown_LacL = lacL(currentTime - 1);
    
        total_rxn_rates = rateup_TetR + ratedown_TetR + rateup_LambdaCl + ratedown_LambdaCl + rateup_LacL + ratedown_LacL;
    
        rand1 = rand();
    
        if rand1 < rateup_TetR / total_rxn_rates
            %production rxn for TetR was chosen
            tetR(currentTime) = tetR(currentTime - 1) + 1;
            lambdaCl(currentTime) = lambdaCl(currentTime - 1);
            lacL(currentTime) = lacL(currentTime - 1);
        
        elseif rand1 < ((rateup_TetR + ratedown_TetR) / total_rxn_rates)
            %degradation rxn for TetR was chosen
            tetR(currentTime) = tetR(currentTime - 1) - 1;
            lambdaCl(currentTime) = lambdaCl(currentTime - 1);
            lacL(currentTime) = lacL(currentTime - 1);
           
        elseif rand1 < ((rateup_TetR + ratedown_TetR + rateup_LambdaCl) / total_rxn_rates)
            %production rxn for LambdaCl was chosen
            tetR(currentTime) = tetR(currentTime - 1);
            lambdaCl(currentTime) = lambdaCl(currentTime - 1) + 1;
            lacL(currentTime) = lacL(currentTime - 1);
             
        elseif rand1 < ((rateup_TetR + ratedown_TetR + rateup_LambdaCl + ratedown_LambdaCl) / total_rxn_rates)
            %degradation rxn for LambdaCl was chosen
            tetR(currentTime) = tetR(currentTime - 1);
            lambdaCl(currentTime) = lambdaCl(currentTime - 1) - 1;
            lacL(currentTime) = lacL(currentTime - 1);
             
        elseif rand1 < ((rateup_TetR + ratedown_TetR + rateup_LambdaCl + ratedown_LambdaCl + rateup_LacL) / total_rxn_rates)
            %production rxn for LacL was chosen
            tetR(currentTime) = tetR(currentTime - 1);
            lambdaCl(currentTime) = lambdaCl(currentTime - 1);
            lacL(currentTime) = lacL(currentTime - 1) + 1;
        else  
            %degradation rxn for LacL was chosen if nothing else is true
            %rateup_TetR + ratedown_TetR + rateup_LambdaCl + ratedown_LambdaCl

            tetR(currentTime) = tetR(currentTime - 1);
            lambdaCl(currentTime) = lambdaCl(currentTime - 1);
            lacL(currentTime) = lacL(currentTime - 1) - 1;
        end
        %randomly assign amount of time spent in each reaction (corresponds to
        %generation/tau)
        rand2 = rand();
        generation(currentTime) = generation(currentTime - 1) + 1/((total_rxn_rates) * log(1 / rand2));
    
        %advance the time
        currentTime = currentTime + 1;
    end

    %store the traces from each iteration of gillespie for averaging later and then plotting
    subsetOfGeneration = floor(linspace(10, currentTime - 1, 2000));
    for count = 1 : length(subsetOfGeneration)
        generation_trace(count) = generation(subsetOfGeneration(count));
        tetR_trace(count) = tetR(subsetOfGeneration(count));
        lambdaCl_trace(count) = lambdaCl(subsetOfGeneration(count));
        lacL_trace(count) = lacL(subsetOfGeneration(count));
    end
    generationOutput = [generationOutput; generation_trace];
    tetROutput = [tetROutput; tetR_trace];
    lambdaClOutput = [lambdaClOutput; lambdaCl_trace];
    lacLOutput = [lacLOutput; lacL_trace];

    %clear traces that were stored in Output matrices
    tetR_trace = 0;
    lambdaCl_trace = 0;
    lacL_trace = 0;
 
    %input to histograms
    periodTetR = [periodTetR, CalculatePeriod(tetR)];
    ampStorageTetR = [ampStorageTetR, MeanAmp(tetR)];
    periodLambdaCl = [periodLambdaCl, CalculatePeriod(lambdaCl)];
    ampStorageLambdaCl = [ampStorageLambdaCl, MeanAmp(lambdaCl)];
    periodLacL = [periodLacL, CalculatePeriod(lacL)];
    ampStorageLacL = [ampStorageLacL, MeanAmp(lacL)];
end
%mean and standard deviation for comparison
meanPeriodTetR = mean(periodTetR);
stdPeriodTetR = std(periodTetR);
meanPeriodLambdaCl = mean(periodLambdaCl);
stdPeriodLambdaCl = std(periodLambdaCl);
meanPeriodLacL = mean(periodLacL);
stdPeriodLacL = std(periodLacL);

hist(periodTetR);
hist(periodLambdaCl);
hist(periodLacL);

%stats for mean amplitude
meanAmpTetR = mean(ampStorageTetR);
stdAmpTetR = std(ampStorageTetR);
meanAmpLambdaCl = mean(ampStorageLambdaCl);
stdAmpLambdaCl = std(ampStorageLambdaCl);
meanAmpLacL = mean(ampStorageLacL);
stdAmpLacL = std(ampStorageLacL);

   
figure;
title('Repressilator Protein Concentration as a Fxn of Generation Number');
hold on;
plot(generationOutput, tetROutput, 'r');
plot(generationOutput, lambdaClOutput, 'y');
plot(generationOutput, lacLOutput, 'b');
xlabel('Generation Number'); ylabel('Protein Molecules');
legend('tetR', 'lambdaCl', 'lacL');
set(gca, 'FontSize', 24);
hold off

