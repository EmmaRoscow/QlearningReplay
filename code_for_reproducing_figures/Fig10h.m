
% Load preprocessed data
evr = importdata('./data/explainedVarianceReactivationAnalysis.mat');

% Plot
evr.plotReactivatedCellPairActivity('rewardArrival');