
% Load preprocessed data
evr = importdata('./data/explainedVarianceReactivationAnalysis.mat');

% Do statistical tests
[~, p_interaction, p_medium, p_high] = evr.mixedEffectsNestedAnova('cellPairs', 'vStr');

% Plot
evr.plotPeakPeriRewardCoactivity(p_medium, p_high, p_interaction, 'cellPairs', 'vStr');
ylim([-10, 120])