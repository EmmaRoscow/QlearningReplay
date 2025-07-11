
% Load stored config parameters
config;

% Set up analysis
evr = ExplainedVarianceReactivation(inclusionCriteria, analysisConfig);
evr.selectSessions();

% Analyse population activity during the sessions
evr.plotRunningSpeedAndPopulationActivity('centralPlatform', 'all', 'plotProgress', true);
evr.plotRunningSpeedAndPopulationActivity('rewardArrival', 'all', 'plotProgress', true);
evr.countRewardResponsiveCells();

% Calculate reactivation and each inter-brain-area contribution to reactivation
evr = evr.calculateCA1EVREV();
evr = evr.calculateVStrEVREV();
evr = evr.calculateInterAreaEVREV();
evr.plotEVREV();

% Analyse the tuning of reactivated CA1-vStr pairs
evr = evr.getSignificantlyReactivatedCellPairs('cellPairs', 'interarea');
evr.countRewardResponsivenessInCellPairs();
evr = evr.getReactivatedCellPairActivity('centralPlatform');
evr.plotReactivatedCellPairActivity('centralPlatform');
evr = evr.getReactivatedCellPairActivity('rewardArrival', 'trialSubset', 'rewarded_only');
evr.plotReactivatedCellPairActivity('rewardArrival');
evr = evr.getReactivatedCellPairActivity('rewardArrival', 'trialSubset', 'all');
evr.plotReactivatedCellPairActivity('rewardArrival');
evr = evr.getEventTriggeredActivity('cellPairs', 'interarea');
[~, p_interaction, p_medium, p_high] = evr.mixedEffectsNestedAnova();
evr.plotPeakPeriRewardCoactivity(p_medium, p_high, p_interaction);

% Analyse the tuning of reactivated vStr-vStr pairs
evr = evr.getSignificantlyReactivatedCellPairs('cellPairs', 'vStr');
evr = evr.getReactivatedCellPairActivity('rewardArrival', 'cellPairs', 'vStr', 'trialSubset', 'rewarded_only');
evr.plotReactivatedCellPairActivity('rewardArrival', 'cellPairs', 'vStr');
evr = evr.getEventTriggeredActivity('cellPairs', 'vStr');
[~, p_interaction, p_medium, p_high] = evr.mixedEffectsNestedAnova('cellPairs', 'vStr');
evr.plotPeakPeriRewardCoactivity(p_medium, p_high, p_interaction, 'cellPairs', 'vStr');
