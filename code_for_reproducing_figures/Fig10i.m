
% Load preprocessed data
evr = importdata('./data/explainedVarianceReactivationAnalysis.mat');

% Do statistical tests
[~, p_interaction, p_medium, p_high] = evr.mixedEffectsNestedAnova();

% Plot
evr.plotPeakPreRewardCoactivity(p_medium, p_high, p_interaction);