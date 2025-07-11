
% Load preprocessed data
evr = importdata('./data/explainedVarianceReactivationAnalysis.mat');

% Plot
evr = evr.getReactivatedCellPairActivity('rewardArrival', 'trialSubset', 'rewarded_only');