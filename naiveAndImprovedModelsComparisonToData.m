function naiveAndImprovedModelsComparisonToData()
modelAndDataComparison('rfpsponge.txt', gillespieRepressilatorOutput, SpongeOutput);
modelAndDataComparison('rfpwithoutsponge.txt', gillespieRepressilatorOutput, SpongeOutput);

modelAndDataComparison('yfpwithoutsponge.txt', gillespieRepressilatorOutput, SpongeOutput);
modelAndDataComparison('yfpsponge.txt', gillespieRepressilatorOutput, SpongeOutput);

modelAndDataComparison('cfpwithoutsponge.txt', gillespieRepressilatorOutput, SpongeOutput);
modelAndDataComparison('cfpsponge.txt', gillespieRepressilatorOutput, SpongeOutput);

