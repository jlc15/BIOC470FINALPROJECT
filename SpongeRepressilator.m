%Updated Model described in "Synchronous long-term oscillations in a
%synthetic gene circuit paper

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

for ii=1:10
    
    %parameters from paper (all based on Supplementary Info 4.3):

    %production rate of protein
    lambda = 60; 
    %at half Vmax
    kTetR = 5; 
    kLambdaCl = 10; 
    kLacL = 10;

    n = 1.5; 

    %initialize p(i)
    initTetR = 20;
    initLambdaCl = 20;
    initLacL = 20;

    tetR=0;
    lambdaCl=0;
    lacL=0;
    generation = 0;

    generationMax = 200;
    generation(1) = 0;

    %total protein amounts are stored in tetR, lambdaCl, lacL
    tetR(1) = initTetR;
    lambdaCl(1) = initLambdaCl; 
    lacL(1) = initLacL;

    %sponge parameters
    N0(1)=20; %number of repressilator plasmids
    Nt(1)=20; %number of sponge plasmids (GFP promoters)

    %counter
    currentTime = 2;

    while generation(currentTime - 1) < generationMax

        %In the model that incorporates the sponge, number of plasmids fluctuates,
        %so to account for this, run a gillespie that is independent of the protein
        %quantification gillespie, but runs parallel in time.

        rateupN0 = 10;
        ratedownN0 = N0(currentTime-1);
        rateupNt = 40;
        ratedownNt = Nt(currentTime-1);

        rateNx_total = rateupN0 + ratedownN0 + rateupNt + ratedownNt;
        rand0 = rand();

        if rand0 < rateupN0/rateNx_total
            N0(currentTime) = N0(currentTime - 1) + 1;
            Nt(currentTime) = Nt(currentTime - 1);
        elseif rand0 < (rateupN0 + ratedownN0)/rateNx_total
            N0(currentTime) = N0(currentTime - 1) - 1;
            Nt(currentTime) = Nt(currentTime - 1);
        elseif rand0 < (rateupN0 + ratedownN0 + rateupNt)/rateNx_total
            Nt(currentTime) = Nt(currentTime - 1) + 1;
            N0(currentTime) = N0(currentTime - 1);
        elseif rand0 < (rateupN0 + ratedownN0 + rateupNt + ratedownNt)/rateNx_total
            Nt(currentTime) = Nt(currentTime - 1) - 1;
            N0(currentTime) = N0(currentTime - 1);
        end

        %Protein production is repressed specifically by the FREE repressor
        %proteins. The formula relating total (T) to free (F) repressor protein is:
        %T = F + 2N0*( (F^n)/(F^n + K^n) ) ==> F^(n+1) - (T+2*N0)*F^n + F*K^n -
        %T*K^n = 0. (For tetR, N0 becomes N0+Nt since the sponge plasmids
        %contribute to the number of promoter regions tetR can bind to.

        %functions for free protein in terms of total protein, hill coeff, N0, Nt,
        %and K values.
        findFreeTetR = @(x) (x^(n+1) + (2*(N0(currentTime-1)+Nt(currentTime-1))-tetR(currentTime-1))*x^n + (kTetR^n)*x - (kTetR^n)*tetR(currentTime-1));
        findFreeLambdaCl = @(x) (x^(n+1) + (2*N0(currentTime-1) - lambdaCl(currentTime-1))*x^n + (kLambdaCl^n)*x - (kLambdaCl^n)*lambdaCl(currentTime-1));
        findFreeLacL = @(x) (x^(n+1) + (2*N0(currentTime-1) - lacL(currentTime-1))*x^n + (kLacL^n)*x - (kLacL^n)*lacL(currentTime-1));

        %calculate free protein
        tetR_freeTemp = fsolve(findFreeTetR,tetR(currentTime-1),optimset('Display','off'));
        tetR_free(currentTime-1) = tetR_freeTemp(tetR_freeTemp == real(tetR_freeTemp));

        lambdaCl_freeTemp = fsolve(findFreeLambdaCl,lambdaCl(currentTime-1),optimset('Display','off'));
        lambdaCl_free(currentTime-1) = lambdaCl_freeTemp(lambdaCl_freeTemp == real(lambdaCl_freeTemp));

        lacL_freeTemp = fsolve(findFreeLacL,lacL(currentTime-1),optimset('Display','off'));
        lacL_free(currentTime-1) = lacL_freeTemp(lacL_freeTemp == real(lacL_freeTemp));

        %Now calculate rates of producing a protein (p-->p+b, p-->p-1) where b is a
        %geometrically distributed variable with mean 10, according to the paper.
        %In the form -> dp(tetR | lambdaCl | lacL)/dt = lambda times kTetR to the 
        %hillCoeff divided by kTetR to the hillCoeff + amount of protein that represses 
        %specific the protein to the hillCoeff all times amount of protein at the last recorded timepoint

        rateup_TetR = (lambda*(kTetR ^ n)) / ( (kTetR ^ n) + (lacL_free(currentTime - 1) ^ n) ); 

        ratedown_TetR = tetR(currentTime - 1);

        rateup_LambdaCl = (lambda*(kLambdaCl ^ n)) / ( (kLambdaCl ^ n) + (tetR_free(currentTime - 1) ^ n));

        ratedown_LambdaCl = lambdaCl(currentTime - 1);

        rateup_LacL = (lambda * (kLacL ^ n)) / ( (kLacL ^ n) + (lambdaCl_free(currentTime - 1) ^ n)); 

        ratedown_LacL = lacL(currentTime - 1);

        total_rxn_rates = rateup_TetR + ratedown_TetR + rateup_LambdaCl + ratedown_LambdaCl + rateup_LacL + ratedown_LacL;
        rand1 = rand();

        if rand1 < rateup_TetR / total_rxn_rates
            %production rxn for TetR was chosen
            tetR(currentTime) = tetR(currentTime - 1) + geornd(0.1);
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
            lambdaCl(currentTime) = lambdaCl(currentTime - 1) + geornd(0.1);
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
            lacL(currentTime) = lacL(currentTime - 1) + geornd(0.1);
        else  
            %degradation rxn for LacL was chosen if nothing else is true
            %rateup_TetR + ratedown_TetR + rateup_LambdaCl + ratedown_LambdaCl

            tetR(currentTime) = tetR(currentTime - 1);
            lambdaCl(currentTime) = lambdaCl(currentTime - 1);
            lacL(currentTime) = lacL(currentTime - 1) - 1;
        end

        rand2 = rand();
        generation(currentTime) = generation(currentTime - 1) + 1/((total_rxn_rates) * log(1/rand2));
        %clear vars
        tetR_freeTemp=0;
        lambdaCl_freeTemp=0;
        lacCl_freeTemp=0;
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
    periodTetR = [periodTetR, CalculatePeriod(generation, tetR)];
    ampStorageTetR = [ampStorageTetR, MeanAmp(tetR)];
    periodLambdaCl = [periodLambdaCl, CalculatePeriod(generation, lambdaCl)];
    ampStorageLambdaCl = [ampStorageLambdaCl, MeanAmp(lambdaCl)];
    periodLacL = [periodLacL, CalculatePeriod(generation, lacL)];
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

%read in the data from the paper's simulation of sponge model
dataRFPSponge = dlmread('rfpsponge.txt');
dataYFPSponge = dlmread('yfpsponge.txt');
dataCFPSponge = dlmread('cfpsponge.txt');   

figure;
title('Protein as a Fxn of Generation Number');
subplot(1,2,1);
hold on;
title('Data');
plot(dataRFPSponge(:, 1), dataRFPSponge(:, 2), 'r.');
plot(dataYFPSponge(:, 1), dataYFPSponge(:, 2), 'y.');
plot(dataCFPSponge(:, 1), dataCFPSponge(:, 2), 'b.');
legend('tetR', 'lambdaCl', 'lacL');
xlabel('Generation Number'); ylabel('Protein Molecules');
hold off;
subplot(1,2,2);
hold on;
title('Model');
plot(mean(generationOutput), mean(tetROutput), 'r');
plot(mean(generationOutput), mean(lambdaClOutput), 'y');
plot(mean(generationOutput), mean(lacLOutput), 'b');
legend('tetR', 'lambdaCl', 'lacL');
xlabel('Generation Number'); ylabel('Protein Molecules');
hold off;
