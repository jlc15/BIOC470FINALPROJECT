function [tetROutput, lambdaClOutput, laclOutput] = gillepsieRepressilator(initTetR, initLambdaCl, initLacL)
%all parameters are from the paper "Synchronous long-term oscillations in a synthetic gene circuit"
%protein degradation rate = 1 (as noted in 4.1 Supplementary Info)

%production rate of protein
lambda = 2000; 
%at half Vmax
kTetR = 100; 
kLambdaCl = 100; 
kLacL = 100;
hillCoeff = 4; 

%will contain average trajectories for each protein
generationOutput = cell(1: 100);
tetROutput = cell(1 : 100);
lambdaClOutput = cell(1 : 100);
lacLOutput = cell(1 : 100);

for ii = 1 : 100
    %clear prev protein amounts
    tetR = 0;
    lambdaCl = 0;
    lacL = 0;
    generation = 0;
    generationMax = 45;
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
            %production rxn for TetR was randomly chosen
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
        %corresponds to tau (amount of time spent in particular reaction)
        rand2 = rand();
        generation(currentTime) = generation(currentTime - 1) + 1/((total_rxn_rates) * log(1 / rand2));
        %advance the time
        currentTime = currentTime + 1;
        generationOutput(ii) = [generationOutput(currentTime - 1) generation(currentTime)];
        tetROutput(ii) = [tetROutput(currentTime - 1) tetR(currentTime)];
        lambdaClOutput(ii) = [lambdaClOutput(currentTime - 1) lambdaClOutput(currentTime)];
        lacLOutput(ii) =  [lacLOutput(currentTime - 1) lacLOutput(currentTime)]
    end     
end



