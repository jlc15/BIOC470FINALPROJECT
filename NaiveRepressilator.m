%protein degradation rate = 1 (as noted in 4.1 Supplementary Info from 
%Synchronous long-term oscillations in a synthetic gene circuit paper)

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
    subsetOfGeneration = floor(linspace(10, currentTime - 1, 2000)); %get 2000 evenly spaced points that can be used to index generation for an axis to plot
    for count = 1 : length(subsetOfGeneration) %use this to index generation and model traces' values at 2000 evenly spaced intervals
        generation_trace(count) = generation(subsetOfGeneration(count));
        tetR_trace(count) = tetR(subsetOfGeneration(count));
        lambdaCl_trace(count) = lambdaCl(subsetOfGeneration(count));
        lacL_trace(count) = lacL(subsetOfGeneration(count));
    end
    generationOutput = [generationOutput; generation_trace]; %stores generation (x-axis) values for each iteration of Gillespie
    tetROutput = [tetROutput; tetR_trace]; %stores trace of protein levels
    lambdaClOutput = [lambdaClOutput; lambdaCl_trace];
    lacLOutput = [lacLOutput; lacL_trace];

    %clear traces that were stored in Output matrices
    tetR_trace = 0;
    lambdaCl_trace = 0;
    lacL_trace = 0;
 
    %input to histograms of mean period and amplitude of protein level
    %oscillations for all trajectories
    periodTetR = [periodTetR, CalculatePeriod(generation, tetR,'model')];
    ampStorageTetR = [ampStorageTetR, MeanAmp(tetR)];
    periodLambdaCl = [periodLambdaCl, CalculatePeriod(generation, lambdaCl,'model')];
    ampStorageLambdaCl = [ampStorageLambdaCl, MeanAmp(lambdaCl)];
    periodLacL = [periodLacL, CalculatePeriod(generation, lacL,'model')];
    ampStorageLacL = [ampStorageLacL, MeanAmp(lacL)];
end

%mean and standard deviation for comparison
meanPeriodTetR = mean(periodTetR);
stdPeriodTetR = std(periodTetR);
meanPeriodLambdaCl = mean(periodLambdaCl);
stdPeriodLambdaCl = std(periodLambdaCl);
meanPeriodLacL = mean(periodLacL);
stdPeriodLacL = std(periodLacL);

%plot histograms
figure;
title('Distributions of Proteins'' Periods of Oscillation');
subplot(1,3,1); hist(periodTetR,10);
xlabel('TetR period'); ylabel('frequency');
subplot(1,3,2); hist(periodLambdaCl,10);
xlabel('LambdaCl period'); ylabel('frequency');
subplot(1,3,3); hist(periodLacL,10);
xlabel('LacL period'); ylabel('frequency');

%stats for mean amplitude
meanAmpTetR = mean(ampStorageTetR);
stdAmpTetR = std(ampStorageTetR);
meanAmpLambdaCl = mean(ampStorageLambdaCl);
stdAmpLambdaCl = std(ampStorageLambdaCl);
meanAmpLacL = mean(ampStorageLacL);
stdAmpLacL = std(ampStorageLacL);

%read in the data from the paper's simulation of original model to analyze
%stats
dataRFPNaive = dlmread('rfpwithoutsponge.txt');
dataCFPNaive = dlmread('cfpwithoutsponge.txt'); 

%find mean period and std of data from paper
PeriodRFPNaive = CalculatePeriod(dataRFPNaive(:,1), dataRFPNaive(:, 2),'data');
PeriodCFPNaive = CalculatePeriod(dataCFPNaive(:, 1), dataCFPNaive(:, 2),'data');

meanPeriodRFPNaive = mean(PeriodRFPNaive);
stdPeriodRFPNaive = std(PeriodRFPNaive);
meanPeriodCFPNaive = mean(PeriodCFPNaive);
stdPeriodCFPNaive = std(PeriodCFPNaive);

%find mean amplitude from data
meanAmpRFPNaive = MeanAmp(dataRFPNaive(:, 2));
meanAmpCFPNaive = MeanAmp(dataCFPNaive(:, 2));

%YFP graph from paper needed specific processing to obtain stats
[meanPeriodYFPNaive, stdPeriodYFPNaive, meanAmpYFPNaive] = StatsYFPNaive();

%Plot model and paper data side by side
figure;
title('Protein as a Fxn of Generation Number');
subplot(1,2,1);
hold on;
title('Data');
plot(dataCFPNaive(:, 1), dataCFPNaive(:, 2), 'b.');
plot(dataYFPNaive(:, 1), dataYFPNaive(:, 2), 'y.');
plot(dataRFPNaive(:, 1), dataRFPNaive(:, 2), 'r.');
legend('tetR', 'lambdaCl', 'lacL');
xlabel('Generation Number'); ylabel('Protein Molecules');
hold off;
subplot(1,2,2);
hold on;
title('Model');
plot(mean(generationOutput), mean(lacLOutput), 'b');
plot(mean(generationOutput), mean(lambdaClOutput), 'y');
plot(mean(generationOutput), mean(tetROutput), 'r');
legend('tetR', 'lambdaCl', 'lacL');
xlabel('Generation Number'); ylabel('Protein Molecules');
hold off;
