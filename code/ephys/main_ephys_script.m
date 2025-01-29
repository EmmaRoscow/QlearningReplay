
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

% Analyse the tuning of reactivated cell pairs
evr.getSignificantlyReactivatedCellPairs();
evr = evr.getReactivatedCellPairActivity('centralPlatform');
evr.plotReactivatedStrActivity('centralPlatform');
evr = evr.getReactivatedCellPairActivity('rewardArrival');
evr.plotReactivatedCellPairActivity('rewardArrival');
evr = evr.getEventTriggeredActivity();
[table, p_interaction, p_medium, p_high] = evr.mixedEffectsNestedAnova();
evr.plotPeakPreRewardCoactivity(p_medium, p_high, p_interaction);
