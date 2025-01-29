
classdef RewardRelatedFiring < EphysAnalysis
    
    properties (Access = public)
        initialLearningSessions
        sessionsForSpecificAnalysis
        rewardFiring
        figures
    end
    
    % Helper methods (should not be inherited)
    methods(Access = private)
        
        function rewardFiring = initialiseRewardFiring(obj)
            
            % Create empty structure with all necessary subfields
            rewardFiring = struct('highRewardExpectation', cell2struct(cell(1, length(obj.inclusionCriteria.rats)), cellfun(@(x) x, obj.inclusionCriteria.rats, 'UniformOutput', false), 2),...
                'mediumRewardExpectation', cell2struct(cell(1, length(obj.inclusionCriteria.rats)), cellfun(@(x) x, obj.inclusionCriteria.rats, 'UniformOutput', false), 2),...
                'perCellHighRewardExpectation', cell2struct(cell(1, length(obj.inclusionCriteria.rats)), cellfun(@(x) x, obj.inclusionCriteria.rats, 'UniformOutput', false), 2),...
                'perCellMediumRewardExpectation', cell2struct(cell(1, length(obj.inclusionCriteria.rats)), cellfun(@(x) x, obj.inclusionCriteria.rats, 'UniformOutput', false), 2),...
                'rewarded', cell2struct(cell(1, length(obj.inclusionCriteria.rats)), cellfun(@(x) x, obj.inclusionCriteria.rats, 'UniformOutput', false), 2),...
                'unrewarded', cell2struct(cell(1, length(obj.inclusionCriteria.rats)), cellfun(@(x) x, obj.inclusionCriteria.rats, 'UniformOutput', false), 2),...
                'perCellRewarded', cell2struct(cell(1, length(obj.inclusionCriteria.rats)), cellfun(@(x) x, obj.inclusionCriteria.rats, 'UniformOutput', false), 2),...
                'perCellUnrewarded', cell2struct(cell(1, length(obj.inclusionCriteria.rats)), cellfun(@(x) x, obj.inclusionCriteria.rats, 'UniformOutput', false), 2));
            
        end
        
        function [armChoice, armLegitimacy, rewardOutcome, spiketimes, rewardArrivalTimes, CPEntryTimes, binnedfr] = loadSessionData(obj, filepath, behavData, session)
            
            % Behavioural data
            armChoice = behavData.reward_probs{session}(behavData.actions{session});
            armLegitimacy = cellfun(@(x, y) x(y), behavData.arm_values{session}, mat2cell(behavData.actions{session}, 1, ones(behavData.n_trials(session), 1)));
            armLegitimacy = cellfun(@(x) ~strcmp(x, 'Illegitimate'), armLegitimacy);
            rewardOutcome = behavData.rewarded{session};
            
            % Other data
            spiketimes = obj.loadVentralStriatalSpikes(filepath);
            [rewardArrivalTimes, CPEntryTimes, ~] = obj.loadRewardCPTimes(filepath);
            data = obj.loadData(filepath);
            binnedfr = data.binnedfr;
            
        end
        
        function [mediumExpTrials, highExpTrials, rewardedTrials, unrewardedTrials] = processTrials(~, armChoice, armLegitimacy, rewardOutcome)
                    
            % Divide trials into mid- and high- probability arms
            mediumExpTrials = strcmp(armChoice, 'medium')  & armLegitimacy;
            highExpTrials = strcmp(armChoice, 'high') & armLegitimacy;

            % Divide trials into rewarded and unrewarded trials
            rewardedTrials = boolean(rewardOutcome) & armLegitimacy;
            unrewardedTrials = ~boolean(rewardOutcome) & armLegitimacy;
            
        end
        
        function [firingRate, ts] = computeAndZScoreFiringRates(obj, includeCells, mediumExpTrials, highExpTrials, rewardedTrials, unrewardedTrials, rewardArrivalTimes, spiketimes, binnedfr)
                    
            % Ensure booleans
            mediumExpTrials = boolean(mediumExpTrials);
            highExpTrials = boolean(highExpTrials);
            includeCells = boolean(includeCells);
            
            timeBefore = 2;
            timeAfter = 0;
            ts = (timeBefore - obj.BinSize/2) : obj.BinSize : (timeAfter - obj.BinSize/2);
            
            % Get firing rates for each trial type
            firingRate = struct();
                firingRate.mediumExpTrials = obj.getEventTriggeredFiring(rewardArrivalTimes(mediumExpTrials)',...
                    spiketimes(includeCells),...
                    timeBefore, timeAfter);
                firingRate.highExpTrials = obj.getEventTriggeredFiring(rewardArrivalTimes(highExpTrials)',...
                    spiketimes(includeCells),...
                    timeBefore, timeAfter);
                firingRate.rewardedTrials = obj.getEventTriggeredFiring(rewardArrivalTimes(rewardedTrials)',...
                    spiketimes(includeCells),...
                    timeBefore, timeAfter);
                firingRate.unrewardedTrials = obj.getEventTriggeredFiring(rewardArrivalTimes(unrewardedTrials)',...
                    spiketimes(includeCells),...
                    timeBefore, timeAfter);
                    
                % Z-score the firing rates
                firingRate.mediumExpTrials = obj.zScore(firingRate.mediumExpTrials, binnedfr(includeCells, :));
                firingRate.highExpTrials = obj.zScore(firingRate.highExpTrials, binnedfr(includeCells, :));
                firingRate.rewardedTrials = obj.zScore(firingRate.rewardedTrials, binnedfr(includeCells, :));
                firingRate.unrewardedTrials = obj.zScore(firingRate.unrewardedTrials, binnedfr(includeCells, :));
        
        end
        
        function rewardFiring = aggregateFiringRates(~, rewardFiring, firingRate, ratName, mediumExpTrials, highExpTrials, rewardedTrials, unrewardedTrials)
                                        
            % Stack per-cell average firing rates together
            nTS = size(firingRate.mediumExpTrials, 3);
            if (sum(mediumExpTrials) > 1) && (sum(highExpTrials) > 1)
                rewardFiring.perCellHighRewardExpectation.(ratName) = [rewardFiring.perCellHighRewardExpectation.(ratName);...
                    reshape(mean(firingRate.highExpTrials, 1), [], nTS)];
                rewardFiring.perCellMediumRewardExpectation.(ratName) = [rewardFiring.perCellMediumRewardExpectation.(ratName);...
                    reshape(mean(firingRate.mediumExpTrials, 1), [], nTS)];
            end
            if (sum(rewardedTrials) > 1) && (sum(unrewardedTrials) > 1)
                rewardFiring.perCellRewarded.(ratName) = [rewardFiring.perCellRewarded.(ratName);...
                    reshape(mean(firingRate.rewardedTrials, 1), [], nTS)];
                rewardFiring.perCellUnrewarded.(ratName) = [rewardFiring.perCellUnrewarded.(ratName);...
                    reshape(mean(firingRate.unrewardedTrials, 1), [], nTS)];
            end

            % Stack the firing rate arrays together
            nTS = size(firingRate.mediumExpTrials, 3);
            firingRate.mediumExpTrials = reshape(firingRate.mediumExpTrials, [], nTS);
            firingRate.highExpTrials = reshape(firingRate.highExpTrials, [], nTS);
            rewardFiring.highRewardExpectation.(ratName) = [rewardFiring.highRewardExpectation.(ratName);...
                firingRate.highExpTrials];
            rewardFiring.mediumRewardExpectation.(ratName) = [rewardFiring.mediumRewardExpectation.(ratName);...
                firingRate.mediumExpTrials];
            firingRate.rewardedTrials = reshape(firingRate.rewardedTrials, [], nTS);
            firingRate.unrewardedTrials = reshape(firingRate.unrewardedTrials, [], nTS);
            rewardFiring.rewarded.(ratName) = [rewardFiring.rewarded.(ratName);...
                firingRate.rewardedTrials];
            rewardFiring.unrewarded.(ratName) = [rewardFiring.unrewarded.(ratName);...
                firingRate.unrewardedTrials];
                    
        end
        
        function plotMeanSEMFiringRate(obj, arm, ratName, colourscheme)
            
            % Calculate mean and standard error from data to plot
            arm(1) = upper(arm(1));
            firing = obj.rewardFiring.(sprintf('perCell%sRewardExpectation', arm)).(ratName);
            firingMean = nanmean(firing);
            firingStd = nanstd(firing);
            nObs = size(firing, 1);
            firingSem = firingStd / sqrt(nObs);
            
            % Get other parameters
            ts = (-5 + obj.BinSize/2) : obj.BinSize : (5 - obj.BinSize/2);
            FaceColor = colourscheme.(sprintf('%s_exp_reward_primary', lower(arm)));
            LineColor = colourscheme.(sprintf('%s_exp_reward_secondary', lower(arm)));
            
            % Plot shaded area for standard error
            patch([ts flip(ts)], [firingMean + firingSem flip(firingMean - firingSem)], 'k', 'FaceColor', FaceColor', 'EdgeColor', 'none')
            
            % Plot line for mean
            plot(ts, firingMean, 'Color', LineColor)
            
        end
        
    end
        
    % Main methods
    methods (Access = public)
        
        % Create inclusion criteria and analysis configurations if not already specified
        function obj = RewardRelatedFiring(inclusionCriteria, analysisConfig)
            
            % Assigns parameter values to the superclass analysis object, applying
            % default values wherever they are not specified by the user
            % 
            % Inputs
            %   - inclusionCriteria:        optional; structure with the following fields
            %     - rats:                       cell array of names of rats to include in analysis,
            %                                    which should correspond to file paths and names
            %     - sessions:                structure wiith the following fields
            %         - initialLearningOnly:  boolean; if true, includes only sessions before reversal learning
            %         - skipEarliest:        boolean; if true, sessions with poor behavioural performance are excluded from analysis
            %         - significantPerformanceThreshold:    float; only sessions with behavioural performance significantly above 
            %                                    this threshold are included in analysis
            %         - cells:                 structure with the following fields
            %             - rewardResponsive:    boolean; if true, only cells that are statistically significantly responsive to reward times
            %                                    are included in analysis
            %             - minimumFiringRate:  float; only cells whose mean firing rate over the whole recording session
            %                                   meets or exceeds this threshold are included in analysis
            %   - analysisConfig:       optional; structure with the following fields
            %     - gaussianWindowWidth:    float; used for smoothing firing rates                
            
            % Specify defaults, to use if not user-specifiied
            defaultInclusionCriteria = struct(...
                    'rats', {{'Quirinius', 'Severus', 'Trevor'}},...
                    'sessions', struct(...
                        'initialLearningOnly', true,...
                        'skipEarliest', true,...
                        'significantPerformanceThreshold', 1/3),...       % Chance: one-third
                    'cells', struct(...
                        'rewardResponsive', false,...
                        'minimumFiringRate', 0));
                    
            defaultAnalysisConfig = struct(...
                    'gaussianWindowWidth', 0.25);
                                        
            % If no inclusion criteria given, use default
            if (nargin < 1) || (isempty(inclusionCriteria)) || (~isstruct(inclusionCriteria))
                inclusionCriteria = defaultInclusionCriteria;
            else
                inclusionCriteria = obj.mergeStructs(defaultInclusionCriteria, inclusionCriteria);
            end
                            
            % If no analysis config given, use default
            if (nargin < 2) || (isempty(analysisConfig)) || (~isstruct(analysisConfig))
                analysisConfig = defaultAnalysisConfig;
            else
                analysisConfig = obj.mergeStructs(defaultAnalysisConfig, analysisConfig);                
            end
            
            % Add both as a properties of the class object
            obj.inclusionCriteria = inclusionCriteria;
            obj.analysisConfig = analysisConfig;
            
        end
        
        
        function merged = mergeStructs(~, defaults, user)
                
            % Merges two structs together, iterating over fields and subfields 
            % and giving precedence to user-defined fields
            % 
            % Inputs
            %   - defaults:      a struct containing a default value for all required fields and subfields
            %                          
            %   - user:           a struct containing at least some of the same fields and subfields as 'default'
            % 
            % Outputs
            %  - merged:        a struct containing all of the fields and subfields as 'default', with values from the 'user' structure where given

            requiredFieldNames = fieldnames(defaults);
            userFieldNames = fieldnames(user);
            merged = user;
            
            % Check if user has provided any field names that are not in
            % the default (if so, it is likely to be a typo)
            for iField = 1:length(userFieldNames)
                field = userFieldNames{iField};
                if ~ismember(field, requiredFieldNames)
                    warning('"%s" is not an expected field name; check that it is spelled correctly.', field);
                end
            end

            for iField = 1:length(requiredFieldNames)
                field = requiredFieldNames{iField};

                % If missing from user input, use default
                if ~isfield(merged, field)
                    merged.(field) = defaults.(field);

                % If not missing, check for subfields
                elseif isstruct(defaults.(field))
                    disp(['Using user-defined value for ' field])
                    requiredSubFieldNames = fieldnames(defaults.(field));
                    for iSubfield = 1:length(requiredSubFieldNames)
                        subfield = requiredSubFieldNames{iSubfield};

                        % If subfield is missing from user input, use default
                        if ~isfield(merged.(field), subfield)
                            merged.(field).(subfield) = defaults.(field).(subfield);
                        else
                            disp(['Using user-defined value for ' field ' - ' subfield])
                        end
                    end
                end
            end
        end
        

        % Select sessions for analysis
        function obj = selectSessions(obj, varargin)
            
            % Selects sessions to include for analysis based on various
            % criteria specified by the superclass property inclusionCriteria, 
            % and assigns them to the property significantSessions
            % 
            % Inputs
            %   - dataPath:     optionally, 'dataPath' along with the directory where behavioural data is stored (otherwise uses default path)
            
            function sessions = getAllAvailableSessions(obj)
                
                % All sessions for which data is available
                sessions.Quirinius = 1:17;
                sessions.Severus = [5:7 9:13 15:17];
                sessions.Trevor = [1:2 4:10 13:20];
                
                % Check that there are no rats unaccounted for
                for ratName = obj.inclusionCriteria.rats
                    if ~ismember(ratName{1}, fieldnames(sessions))
                        warning(sprintf('No learning sessions known for rat "%s"', ratName{1}))
                    end
                end
                
            end
            
            function sessions = skipEarliestSessions(sessions)
                % Excluding first 2-5 sessions before consistent behaviour was noted as established by experimenter
                sessions.Quirinius = sessions.Quirinius(sessions.Quirinius >= 3);
                sessions.Severus = sessions.Severus(sessions.Severus >= 5);
                sessions.Trevor = sessions.Trevor(sessions.Trevor >= 6);
            end
            
            function sessions = initialLearningOnly(sessions)
                % Excluding sessions from the reversal learning point onwards
                sessions.Quirinius = sessions.Quirinius(sessions.Quirinius <= 12);
                sessions.Severus = sessions.Severus(sessions.Severus <= 12);
                sessions.Trevor = sessions.Trevor(sessions.Trevor <= 15);
            end
            
            function significantSessions = significantPerformanceOnly(obj, sessions, dataPath)
                
                % Skip if minimum performance requirement is 0
                if obj.inclusionCriteria.sessions.significantPerformanceThreshold == 0
                    significantSessions = sessions;
                    return;
                end
                
                % Create a structure for storing significant sessions
                nRats = length(obj.inclusionCriteria.rats);
                pChance = obj.inclusionCriteria.sessions.significantPerformanceThreshold;
                significantSessions = cell2struct(cell(1, nRats), obj.inclusionCriteria.rats, 2);
                
                for iRat = 1:nRats
                    ratName = obj.inclusionCriteria.rats{iRat};
                    filepath = fullfile(dataPath, sprintf('behaviour_%s.mat', lower(ratName)));
                    behavData = importdata(filepath);
                    for iSess = sessions.(ratName)
                        % Perform binomial test to check if rate of optimal arm choice is significantly above user-defined threshold
                        % at alpha level 0.05
                        n = behavData.n_trials(iSess);
                        x = behavData.n_trials(iSess) * behavData.proportion_optimal(iSess);
                        pval = 1 - binocdf(x - 1, n, pChance);
                        if pval < 0.05
                            significantSessions.(ratName) = horzcat( significantSessions.(ratName), iSess );
                        end
                    end
                end
                
            end
            
            % Process path where data is stored
            p = inputParser;
            addParameter(p, 'dataPath', fullfile('.', 'data', 'behavioural_data'), @(x) ischar(x) || isstring(x));
            parse(p, varargin{:});
            dataPath = p.Results.dataPath;
                            
            % Sessions to evaluate: all available sessions, then exclude as needed
            allSessions = getAllAvailableSessions(obj);
            initialLearningSessions = allSessions;
            sessionsForSpecificAnalysis = allSessions;
            if obj.inclusionCriteria.sessions.initialLearningOnly
                % Skip later sessions with reversal learning
                initialLearningSessions = initialLearningOnly(allSessions);
                sessionsForSpecificAnalysis = initialLearningSessions;
            end
            
            if obj.inclusionCriteria.sessions.skipEarliest
                % Also skip the first few sessions before good behaviour established
                sessionsWithConsistentBehaviour = skipEarliestSessions(sessionsForSpecificAnalysis);
                initialLearningSessions = sessionsWithConsistentBehaviour;
                sessionsForSpecificAnalysis = sessionsWithConsistentBehaviour;
            end
            
            % Restrict to sessions with significant behavioural performance
            significantSessions = significantPerformanceOnly(obj, sessionsForSpecificAnalysis, dataPath);
            
            obj.sessionsForSpecificAnalysis = significantSessions;            
            obj.initialLearningSessions = initialLearningSessions;
            
        end

        
        % Select cells with high firing rates
        function includeCells = selectHighFiringRateCells(obj, binnedfr)
            
            % Identifies which rows of a binned firing rate matrix
            % represent an average firing rate greater than or equal to the
            % minimum set by the inclusionCriteria property
            % 
            % Inputs
            %   - binnedfr:      N x T array of N neurons' firing rates in T time bins
            % 
            % Outputs
            %   - includeCells: boolean array of length N, true wherever the mean firing rate
            %                        is at least the minimum set by the inclusionCriteria 
            
            includeCells = boolean(ones(1, size(binnedfr, 1)));
            
            % Select cells with high enough firing rates overall
            meanfr = mean(binnedfr, 2);
            includeCells(meanfr < obj.inclusionCriteria.cells.minimumFiringRate) = false;
        
        end
                
        
        % Get reward-related firing
        function obj = getRewardRelatedFiring(obj, varargin)
            
            % Calculates the z-scored firing rate of ventral striatal
            % neurons in the two seconds before arrival at the reward
            % location on trials where rats made a "legitimate" choice
            % (i.e. did not return the arm visited on the previous trial).
            % Splits this data by action choice (high- and
            % medium-probability arms) and reward outcome, and assigns it
            % to the property rewardFiring
            % 
            % Inputs
            %   - behaviouralDataPath:      optional; specifies the path where .mat files of behavioural data are stored
            %   - ephysDataPath:             optional; specifies the path where .mat files of electrophysiological and related data are stored
            % 
            
            % Process paths where data is stored
            p = inputParser;
            addParameter(p, 'behaviouralDataPath', fullfile('.', 'data', 'behavioural_data'), @(x) ischar(x) || isstring(x));
            addParameter(p, 'ephysDataPath', fullfile('.', 'data', 'ephys_data'), @(x) ischar(x) || isstring(x));
            parse(p, varargin{:});
            behaviouralDataPath = p.Results.behaviouralDataPath;
            ephysDataPath = p.Results.ephysDataPath;
            
            % Create a structure for storing the results
            rewardFiring = obj.initialiseRewardFiring();
            
            for iRat = 1:length(obj.inclusionCriteria.rats)
                
                ratName = obj.inclusionCriteria.rats{iRat};
                                
                % Load behavioural data
                behavData = importdata(fullfile(behaviouralDataPath, ['behaviour_' lower(ratName) '.mat']));
                
                ratSessions = obj.sessionsForSpecificAnalysis.(ratName);
                nSessions = length(ratSessions);
                
                for iSess = 1:nSessions
                    
                    session = ratSessions(iSess);
                    filepath = fullfile(ephysDataPath, ['Rat_' ratName(1)], ['Session' num2str(session)]);
                    
                    % Load session-specific data
                    [armChoice, armLegitimacy, rewardOutcome, spiketimes, rewardArrivalTimes, CPEntryTimes, binnedfr] = obj.loadSessionData(filepath, behavData, session);
                    
                    % Select cells according to inclusion criteria
                    includeCells = obj.selectRewardResponsiveCells(spiketimes, rewardArrivalTimes, CPEntryTimes);
                    binnedfr = binnedfr(find(includeCells), :);
                    includeCells = obj.selectHighFiringRateCells(binnedfr);
                    
                    % Process trials
                    [mediumExpTrials, highExpTrials, rewardedTrials, unrewardedTrials] = obj.processTrials(armChoice, armLegitimacy, rewardOutcome);
                                        
                    % Get firing of ventral striatal cells 2 seconds around arrival time at reward location
                    [firingRate, ts] = obj.computeAndZScoreFiringRates(includeCells, mediumExpTrials, highExpTrials, rewardedTrials, unrewardedTrials, rewardArrivalTimes, spiketimes, binnedfr);
                    
                    % Aggregate data
                    rewardFiring = obj.aggregateFiringRates(rewardFiring, firingRate, ratName, mediumExpTrials, highExpTrials, rewardedTrials, unrewardedTrials);
                    
                end
                
            end
            
            rewardFiring.ts = ts;
            obj.rewardFiring = rewardFiring;
            
        end
        
        
        % Plot reward-related firing
        function obj = plotRewardRelatedFiring(obj)
            
            % Plots line graphs of trial-averaged +- sem firing rates around the point of arrival at the reward location, 
            % for high- and medium-probability arms separately, and on a separate figure for each rat
            
            % Check that there is data to plot
            assert(isstruct(obj.rewardFiring) && ~isempty(obj.rewardFiring), 'RewardRelatedFiring:dataNotAssigned', 'No data to plot. Please run the getRewardRelatedFiring() method to generate data')
            for ratName = obj.inclusionCriteria.rats
                assert(isfield(obj.rewardFiring.perCellMediumRewardExpectation, ratName{1}), 'RewardRelatedFiring:dataNotAssigned', sprintf('Rat "%s" is missing from rewardFiring data. Please re-run the getRewardFiring() method', ratName{1}))
                assert(isfield(obj.rewardFiring.perCellHighRewardExpectation, ratName{1}), 'RewardRelatedFiring:dataNotAssigned', sprintf('Rat "%s" is missing from rewardFiring data. Please re-run the getRewardFiring() method', ratName{1}))
            end
            
            % Set design properties
            colour_scheme;
            axes_properties;
            
            figures = {};
            
            for iRat = 1:length(obj.inclusionCriteria.rats)
                
                ratName = obj.inclusionCriteria.rats{iRat};
            
                % Create figure of desired size
                figures{iRat} = figure;

                % Plot firing related to medium reward expectation (left subplot)
                subplot(1, 2, 1); hold on
                obj.plotMeanSEMFiringRate('medium', ratName, colourscheme)
                ylimits = [ylim];
                
                 % Plot firing related to high reward expectation (right subplot)
                subplot(1, 2, 2); hold on
                obj.plotMeanSEMFiringRate('high', ratName, colourscheme)
                ylimits = [ylim];                
                
                % Axis formatting
                for i = 1:2
                    subplot(1, 2, i)
                    xlabel('time (s)')
                    set(gca, 'XTick', [-5 0 5])
                    ylim([min(ylimits(:, 1)) max(ylimits(:, 2))])
                    line([0 0], [0 1]*max(ylimits(:, 2)), 'Color', colourscheme.zero, 'LineStyle', '--', 'LineWidth', 2)
                end
                subplot(1, 2, 1)
                ylabel(sprintf('Rat %s\nz-scored FR', ratName(1)))
                
                
            end
            
            obj.figures = figures;
            
        end
        
        % Plot bar graph
        function [h] = plotBar(~, data, colour_name, x, h)
            
            % Plots a single bar with a horizontal line indicating the median,
            % a box for the interquartile range, and whiskers extending from the quartiles to the full range
            % 
            % Inputs
            %   - data:             vector of numeric data whose distribution to plot
            %   - colour_name: corresponding to a primary and secondary colour specified in the colourscheme
            %   - x:                 float; where on the x-axis to centre the bar (allows calling the function multiple times to plot multiple bars)
            %   - h:                 figure handle
            % 
            % Outputs
            %   - h:                 figure handle with the bar added

            axes_properties;
            colour_scheme;

            % Calculate quartiles
            assert(min(size(data))==1, 'ExplainedVarianceReactivation:invalidData', 'Data must be a vector.')
            quartiles = quantile(data, [0.25 0.5 0.75]);

            % Get primary and secondary colours from colours scheme
            evalc(['colour_primary = colourscheme.'  colour_name '_primary']);
            evalc(['colour_secondary = colourscheme.'  colour_name '_secondary']);

            % Plot whiskers showing full range
            line([0 0]+x, [min(data) max(data)], 'Color', 'k', 'LineWidth', 2)

            % Plot bar showing interquartile range
            patch([-0.25 0.25 0.25 -0.25]+x, [quartiles(1) quartiles(1) quartiles(3) quartiles(3)], colour_primary, 'EdgeColor', colour_secondary)

            % Plot line showing median
            line([-0.25 0.25]+x, [quartiles(2) quartiles(2)], 'Color', colour_secondary, 'LineWidth', 2)

        end
        
        % Plot indication of significance if significant
        function [h, p] = plotTTestSignificance(~, group1, group2, xData, h)
            
            % Performs a paired-sample t-test and plots lines and a significant marker: * if significant at alpha=0.05, and n.s. otherwise
            % 
            % Inputs
            %   - group1:          vector of numeric data on which to perform the t-test
            %   - group2:          second vector of numeric data on which to perform the t-test
            %   - xData:           vector of length 2, corresponding to where the two groups might be plotted; this is used to position the significance marker
            %   - h:                 figure handle
            % 
            % Outputs
            %   - h:                 figure handle with the bar added
            %   - p:                 p-value of the paired t-test

            % Perform t-test
            [~, p] = ttest(group1, group2);

            % Plot lines to show comparison
            y = ylim;
            assert(numel(xData) <= 2, 'ExplainedVarianceReactivation:invalidData', 'xData must not have more than 2 elements; a vector of length 2 is recommended.')
            line([xData(1) xData(1)], [[0.9 0.95] * y(2)], 'Color', 'k', 'LineWidth', 1.5)
            line([xData(2) xData(2)], [[0.9 0.95] * y(2)], 'Color', 'k', 'LineWidth', 1.5)
            line(xData, [[0.95 0.95] * y(2)], 'Color', 'k', 'LineWidth', 1.5)
            
            if p < 0.05
                text(mean(xData), 0.96 * y(2), '*', 'FontSize', 24, 'HorizontalAlignment', 'center')
            else
                text(mean(xData), 0.98 * y(2), 'n.s.', 'FontSize', 13, 'HorizontalAlignment', 'center')
            end

        end

        function peakFiring = getPeakFiringRates(obj, arm, trialType)
            
            % Get firing rates for the 2 seconds before arrival at reward location
            % 
            % Inputs
            %   - arm:           'high' or 'medium'; must be a field of the perCellRewardExpectation property
            %   - trialType:    'reactivated' 'control' or NaN, if not Nan, must be a subfield of the arm field
            % 
            % Outputs
            %   - peakFiring:   peak firing rate in the 2 seconds prior to arrival at the reward location
            
            % Check that there is data to plot
            arm(1) = upper(arm(1));
            requiredField = sprintf('perCell%sRewardExpectation', arm);
            assert(isstruct(obj.rewardFiring) && ~isempty(obj.rewardFiring), 'RewardRelatedFiring:dataNotAssigned', 'No data to plot. Please run the getRewardRelatedFiring() method to generate data.')
            assert(isfield(obj.rewardFiring, requiredField) && ~isempty(obj.rewardFiring.(requiredField)), 'RewardRelatedFiring:dataNotAssigned', sprintf('Poperty "rewardFiring" is missing data for the "%s" field. Please run the getRewardRelatedFiring() method to generate data.', requiredField))
            if ~isnan(trialType)
                assert(isfield(obj.rewardFiring.(requiredField), trialType) && ~isempty(obj.rewardFiring.(requiredField).(trialType)), 'RewardRelatedFiring:dataNotAssigned', sprintf('Poperty "rewardFiring.%s" is missing data for the "%s" field. Please run the getRewardRelatedFiring() method to generate data.', requiredField, trialType))
            end
            assert(isfield(obj.rewardFiring, 'ts') && ~isempty(obj.rewardFiring.ts), 'RewardRelatedFiring:dataNotAssigned', 'Poperty "rewardFiring" is missing timestamp data for the "ts" field. Please run the getRewardRelatedFiring() method to generate data.')

            ts = obj.rewardFiring.ts;
            if isnan(trialType)
                firing = obj.rewardFiring.(sprintf('perCell%sRewardExpectation', arm));
            else
                firing = obj.rewardFiring.(sprintf('perCell%sRewardExpectation', arm)).(trialType);
            end
            firing = cell2mat(struct2cell(firing));
            assert(size(firing, 2) == numel(ts), 'RewardRelatedFiring:invalidData', 'Reward firing data does not match timestamps.')
            preArrivalFiring = firing(:, (ts > -2) & (ts < 0));
            
            % Calculate peak
            peakFiring = max(preArrivalFiring, [], 2);
                        
        end

        
        % Plot peak pre-reward firing rate
        function obj = plotPeakPreRewardFiring(obj)
            
            % Plots a box-and-whisker plot of peak per-trial firing rates in the 2 seconds before arrival at reward locations, 
            % showing median and quartiles, for high- and medium-probability arms
            
            % Check that there is data to plot
            assert(isstruct(obj.rewardFiring) && ~isempty(obj.rewardFiring), 'RewardRelatedFiring:dataNotAssigned', 'No data to plot. Please run the getRewardRelatedFiring() method to generate data')
            assert(~isempty(obj.inclusionCriteria.rats), 'RewardRelatedFiring:dataNotAssigned', 'Rat names are missing from inclusionCriteria property')
            for ratName = obj.inclusionCriteria.rats
                assert(isfield(obj.rewardFiring.perCellMediumRewardExpectation, ratName{1}), 'RewardRelatedFiring:dataNotAssigned', sprintf('Rat "%s" is missing from rewardFiring data. Please re-run the getRewardFiring() method', ratName{1}))
                assert(isfield(obj.rewardFiring.perCellHighRewardExpectation, ratName{1}), 'RewardRelatedFiring:dataNotAssigned', sprintf('Rat "%s" is missing from rewardFiring data. Please re-run the getRewardFiring() method', ratName{1}))
            end
            
            % Set design properties
            colour_scheme;
            axes_properties;
            
            fig = figure('Position', [680 558 230 420]); hold on
                            
            % Calculate peak firing rate in the 2 seconds prior to arrival at medium-expectation reward location
            mediumRewardExpectationPeaks = obj.getPeakFiringRates('medium', NaN);

            % Calculate peak firing rate in the 2 seconds prior to arrival at high-expectation reward location
            highRewardExpectationPeaks = obj.getPeakFiringRates('high', NaN);
                
            % Plot bars
            fig = obj.plotBar(mediumRewardExpectationPeaks, 'medium_exp_reward', 1, fig);
            fig = obj.plotBar(highRewardExpectationPeaks, 'high_exp_reward', 2, fig);

            % T-test
            [fig, p] = obj.plotTTestSignificance(highRewardExpectationPeaks, mediumRewardExpectationPeaks, [1, 2], fig);
            
            % Axis formatting
            ylabel('z-scored FR')
            set(gca, 'XTick', [1 2])
            xlim([0.5 2.5])
            set(gca, 'XTickLabel', {'medium', 'high'})
            xlabel('expected reward')
            
            obj.figures{end+1} = fig;
            
        end

    end
    
end

