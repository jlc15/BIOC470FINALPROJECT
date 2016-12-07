%function gillepsieRepressilator()
%protein degradation rate = 1 (as noted in 4.1 Supplementary Info from 
%Synchronous long-term oscillations in a synthetic gene circuit paper)

for ii = 1:100
%parameters from paper:

%production rate of protein
lambda = 2000; 
%at half Vmax
kTetR = 100; 
kLambdaCl = 100; 
kLacL = 100;

hillCoeff = 4; 

%initialize p(i)
initTetR = 20;
initLambdaCl = 20;
initLacL = 20;

tetR=0;
lambdaCl=0;
lacL=0;
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
%specific the protein to the hillCoeff all times amount of protein at the last recorded timepoint
    rateup_TetR = ((lambda * (kTetR ^ hillCoeff)) / ( (kTetR ^ hillCoeff) + (lacL(currentTime - 1) ^ hillCoeff))); 
%     * tetR(currentTime - 1);
    ratedown_TetR = tetR(currentTime - 1);
    
    rateup_LambdaCl = ((lambda * (kLambdaCl ^ hillCoeff)) / ( (kLambdaCl ^ hillCoeff) + (tetR(currentTime - 1) ^ hillCoeff)));
%     * lambdaCl(currentTime - 1); 
    ratedown_LambdaCl = lambdaCl(currentTime - 1);
    
    rateup_LacL = ((lambda * (kLacL ^ hillCoeff)) / ( (kLacL ^ hillCoeff) + (lambdaCl(currentTime - 1) ^ hillCoeff))); 
%     * lacL(currentTime - 1)
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
    
    rand2 = rand();
    generation(currentTime) = generation(currentTime - 1) + 1/((total_rxn_rates) * log(1 / rand2));
    %advance the time
    currentTime = currentTime + 1;
end

% figure;
% title('Repressilator Protein Concentration as a Fxn of Generation Number');
% hold on
% plot(generation, tetR, 'r');
% plot(generation, lambdaCl, 'g');
% plot(generation, lacL, 'b');
% xlabel('Generation Number'); ylabel('Protein Concentration');
% legend('tetR Concentration', 'lambdaCl Concentration', 'lacL Concentration');
% set(gca, 'FontSize', 24);
% hold off

%find maxima of peaks of protein levels
[maxdat_tetR,maxidx_tetR] = extrema(smooth(tetR,750));
maxima_tetR = maxdat_tetR(maxdat_tetR>mean(maxdat_tetR));
maxima_idx_tetR = generation(maxidx_tetR(maxdat_tetR>mean(maxdat_tetR)));

%Sort the maxima by the order of when they occur, not how high they are.
%Thus, sort by the maxima indices.
[~, j] = sort(maxima_idx_tetR);
maxima_idx_tetR = maxima_idx_tetR(j);
maxima_tetR = maxima_tetR(j);

% figure; 
% plot(generation, tetR,'.-r',maxima_idx_tetR,maxima_tetR,'kx');

for mm=1:length(maxima_idx_tetR)-1
    period(mm) = maxima_idx_tetR(mm+1)-maxima_idx_tetR(mm);
end
%discard transients
period = period(period>5);

%mean and standard deviation for comparison
mean_period(ii) = mean(period);
std_period(ii) = std(period);
%mean amplitude of oscillations is just:
mean_amp(ii) = mean(maxdat_tetR);

%clear vars after use
period = 0;
maxima_idx_tetR=0;
maxima_tetR=0;
end

mean_period_overall = mean(~isnan(mean_period));
mean_std_overall = mean(~isnan(std_period));
mean_amp_overall = mean(~isnan(mean_amp));
