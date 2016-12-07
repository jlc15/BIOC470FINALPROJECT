%function gillepsieRepressilator()
%protein degradation rate = 1 (as noted in 4.1 Supplementary Info from 
%Synchronous long-term oscillations in a synthetic gene circuit paper)

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

generationMax = 10;
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
rateupNt = 10;
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
tetR_freeTemp = fsolve(findFreeTetR,tetR(currentTime-1));
tetR_free(currentTime-1) = tetR_freeTemp(tetR_freeTemp == real(tetR_freeTemp));

lambdaCl_freeTemp = fsolve(findFreeLambdaCl,lambdaCl(currentTime-1));
lambdaCl_free(currentTime-1) = lambdaCl_freeTemp(lambdaCl_freeTemp == real(lambdaCl_freeTemp));

lacL_freeTemp = fsolve(findFreeLacL,lacL(currentTime-1));
lacL_free(currentTime-1) = lacL_freeTemp(lacL_freeTemp == real(lacL_freeTemp));

%Now calculate rates of producing a protein (p-->p+b, p-->p-1) where b is a
%geometrically distributed variable with mean 10, according to the paper.
%In the form -> dp(tetR | lambdaCl | lacL)/dt = lambda times kTetR to the 
%hillCoeff divided by kTetR to the hillCoeff + amount of protein that represses 
%specific the protein to the hillCoeff all times amount of protein at the last recorded timepoint

    rateup_TetR = ((lambda*(kTetR ^ n)) / ( (kTetR ^ n) + (lacL_free(currentTime - 1) ^ n))); 

    ratedown_TetR = tetR(currentTime - 1);
    
    rateup_LambdaCl = ((lambda*(kLambdaCl ^ n)) / ( (kLambdaCl ^ n) + (tetR_free(currentTime - 1) ^ n)));
 
    ratedown_LambdaCl = lambdaCl(currentTime - 1);
    
    rateup_LacL = ((lambda * (kLacL ^ n)) / ( (kLacL ^ n) + (lambdaCl_free(currentTime - 1) ^ n))); 

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

figure;
title('Repressilator Protein Concentration as a Fxn of Generation Number');
hold on
plot(generation, tetR, 'r');
plot(generation, lambdaCl, 'g');
plot(generation, lacL, 'b');
xlabel('Generation Number'); ylabel('Protein Concentration');
%legend('tetR Concentration', 'lambdaCl Concentration', 'lacL Concentration');
set(gca, 'FontSize', 24);
hold off
