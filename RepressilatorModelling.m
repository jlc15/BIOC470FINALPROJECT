function RepressilatorModelling()
%initialize p(i) as they are constant for both models

modelAndDataComparison('rfpwithoutsponge.txt', NaiveRepressilator(20, 20, 20));
% modelAndDataComparison('rfpwithoutsponge.txt', gillespieRepressilatorOutput, SpongeOutput);
% 
% modelAndDataComparison('yfpwithoutsponge.txt', gillespieRepressilatorOutput, SpongeOutput);
% modelAndDataComparison('yfpsponge.txt', gillespieRepressilatorOutput, SpongeOutput);
% 
% modelAndDataComparison('cfpwithoutsponge.txt', gillespieRepressilatorOutput, SpongeOutput);
% modelAndDataComparison('cfpsponge.txt', gillespieRepressilatorOutput, SpongeOutput);