
% Load preprocessed data
evr = importdata('./data/explainedVarianceReactivationAnalysis.mat');

% Plot (these take a long time; progress bars indicate how long)
evr.plotRunningSpeedAndPopulationActivity('centralPlatform', 'all', 'plotProgress', true);
evr.plotRunningSpeedAndPopulationActivity('rewardArrival', 'all', 'plotProgress', true);