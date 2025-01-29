classdef ExplainedVarianceReactivation < RewardRelatedFiring
    
    properties (Access = public)
        CA1_vStr
        perCellContributions
        vStr
        CA1
        data
        meanEventTriggeredActivity
        taskCoactivity
    end
    
    methods(Access = public)
        
        % Create inclusion criteria if not already specified
        function obj = ExplainedVarianceReactivation(inclusionCriteria, analysisConfig)
            
            % Assigns parameter values to the superclass analysis object, applying
            % default values wherever they are not specified by the user
            % 
            % Inputs
            %   - inclusionCriteria:                  optional; structure with the following fields
            %     - rats:                                cell array of names of rats to include in analysis,
            %                                             which should correspond to file paths and names
            %     - sessions:                         structure wiith the following fields
            %         - initialLearningOnly:        boolean; if true, includes only sessions before reversal learning
            %         - skipEarliest:                  boolean; if true, sessions with poor behavioural performance are excluded from analysis
            %         - significantPerformanceThreshold:    float; only sessions with behavioural performance significantly above 
            %                                              this threshold are included in analysis
            %         - minimumNCA1Cells:      integer; sessions below this number of CA1 neurons are included for broad analyses
            %                                              but excluded from specific EV-REV analyses
            %         - minimumNVStrCells:      integer; sessions below this number of ventral striatal neurons are included for broad analyses
            %                                              but excluded from specific EV-REV analyses
            %         - cells:                           structure with the following fields
            %             - rewardResponsive:    boolean; if true, only cells that are statistically significantly responsive to reward times
            %                                              are included in analysis
            %             - minimumFiringRate:   float; only cells whose mean firing rate over the whole recording session
            %                                              meets or exceeds this threshold are included in analysis
            %   - analysisConfig:                   optional; structure with the following fields
            %     - postDuration:                    float; total duration of activity during POST to analyse EV-REV for, in seconds
            %     - nRipples:                          integer, or array of two integers; the total number of ripples to analyse EV-REV for
            %                                              If length 1, applies to POST only; if length 2, the first applies to PRE ripples and the second to POST.
            %                                              Applies only when 'rippleActivityMode' is not 'none'
            %     - rippleActivityMode:            'none', 'variable', or 'equalise'; if 'none', all spiking activity during PRE and POST is used for EV-REV analysis
            %                                              (within and between sharp-wave ripples); if 'variable', spikes fired between the start and end of
            %                                              every ripple are used; if 'equalise', spikes fired during a window of 200ms from the ripple peak is used
            %     - controlCriterion:                 'minimalContributions' or 'negativeContributions'; if 'negativeContributions', control (non-reactivated)
            %                                              cell-pairs are defined as those that have the lowest, most negative contributions to overall EV-REV; 
            %                                              if 'minimalContributions', cell-pairs are defined as those that have the smallest, closest-to-zero contributions to
            %                                              overall EV-REV
            %     - reactivationQuantileThreshold:  float, between 0 and 1 exclusive; reactivated and control cell pairs are defined by their contributions to EV-REV
            %                                              being in the top reactivationQuantileThreshold or bottom reactivationQuantileThreshold percentile, respectively
            %     - gaussianWindowWidth:       float; used for smoothing firing rates                
            
            % Specify defaults, to use if not user-specifiied
            defaultInclusionCriteria = struct(...
                    'rats', {{'Quirinius', 'Severus', 'Trevor'}},...
                    'sessions', struct(...
                        'initialLearningOnly', false,...
                        'skipEarliest', false,...
                        'significantPerformanceThreshold', 1/3,...
                        'minimumNCA1Cells', 1,...
                        'minimumNVStrCells', 1),...
                    'cells', struct(...
                        'rewardResponsive', false,...
                        'minimumFiringRate', 0));
            
            % Specify defaults, to use if not user-specifiied
            default_analysisConfig = struct(...
                'postDuration', NaN,...
                'nRipples', NaN,...
                'rippleActivityMode', 'variable',...
                'controlCriterion', 'minimalContributions',...
                'reactivationQuantileThreshold', 0.1,...
                'gaussianWindowWidth', 0.25);
                
            % If no inclusion criteria given, use default
            if (nargin < 1) || (isempty(inclusionCriteria)) || (~isstruct(inclusionCriteria))
                inclusionCriteria = defaultInclusionCriteria;
            else
                inclusionCriteria = obj.mergeStructs(defaultInclusionCriteria, inclusionCriteria);
            end
                            
            % If no analysis config given, use default
            if (nargin < 2) || (isempty(analysisConfig)) || (~isstruct(analysisConfig))
                analysisConfig = default_analysisConfig;
            else
                analysisConfig = obj.mergeStructs(default_analysisConfig, analysisConfig);
            end
                            
            % Add both as a properties of the class object
            obj.inclusionCriteria = inclusionCriteria;
            obj.analysisConfig = analysisConfig;
            
        end
        
        
        % Select sessions
        function obj = selectSessions(obj, varargin)
            
            % Selects sessions to include for analysis based on various
            % criteria specified by the superclass property inclusionCriteria, 
            % and assigns them to the property significantSessions
            % 
            % Inputs
            %   - behaviouralDataPath:     optional; specifies the path where .mat files of behavioural data are stored
            %   - ephysDataPath:            optional; specifies the path where .mat files of electrophysiological and related data are stored
            
            % Process paths where data is stored
            p = inputParser;
            addParameter(p, 'behaviouralDataPath', fullfile('.', 'data', 'behavioural_data'), @(x) ischar(x) || isstring(x));
            addParameter(p, 'ephysDataPath', fullfile('.', 'data', 'ephys_data'), @(x) ischar(x) || isstring(x));
            parse(p, varargin{:});
            behaviouralDataPath = p.Results.behaviouralDataPath;
            ephysDataPath = p.Results.ephysDataPath;
            
            % Identify sessions based on behavioural criteria (learning stage and performance)
            objWithInitialSessionSelection = selectSessions@RewardRelatedFiring(obj, 'dataPath', behaviouralDataPath);
            sessionsForBroadAnalysis = objWithInitialSessionSelection.initialLearningSessions;
            significantSessions = objWithInitialSessionSelection.sessionsForSpecificAnalysis;
            
            % Remove sessions with too few CA1 and/or ventral striatal units
            sessionsForSpecificAnalysis = significantSessions;
            for rat = obj.inclusionCriteria.rats
                ratName = rat{1};
                for session = significantSessions.(ratName)
                    filepath = fullfile(ephysDataPath, sprintf('Rat_%s', ratName(1)), sprintf('Session%s', num2str(session)));
                    sessionData = obj.loadData(filepath);
                    if sessionData.nCA1Units < obj.inclusionCriteria.sessions.minimumNCA1Cells
                        sessionsForSpecificAnalysis.(ratName) = sessionsForSpecificAnalysis.(ratName)(sessionsForSpecificAnalysis.(ratName) ~= session);
                    end
                    if sessionData.nNAcUnits < obj.inclusionCriteria.sessions.minimumNVStrCells
                        sessionsForSpecificAnalysis.(ratName) = sessionsForSpecificAnalysis.(ratName)(sessionsForSpecificAnalysis.(ratName) ~= session);
                    end
                end
            end
                    
            obj.initialLearningSessions = objWithInitialSessionSelection.initialLearningSessions;
            obj.sessionsForBroadAnalysis = sessionsForBroadAnalysis;
            obj.sessionsForSpecificAnalysis = sessionsForSpecificAnalysis;

        end
        
        % Filter ripples according to analysis configurations
        function [rippleStart, rippleStop, ripplePeak] = filterRipples(obj, data)
            
            % Filters ripples in PRE and POST-task rest for one session according to the number specified by inclusion criteria
            % 
            % Inputs
            %   - data:                            structure containing the following fields
            %       - rippleStart:                vector of timestamps marking the start of all identified sharp-wave ripples in the session
            %       - rippleStop:                vector of timestamps marking the end of all identified sharp-wave ripples in the session;
            %                                        must be the same length as 'rippleStart'
            %       - ripplePeak:                vector of timestamps marking the time of peak amplitude of all identified sharp-wave ripples in the session;
            %                                        must be the same length as 'rippleStart'
            %       - startEndTimes:          3 x 2 array of timestamps in seconds,
            %                                        marking the start and end of the pre-task rest period (first row), 
            %                                        start and end of the task period (second row), 
            %                                        and start and end of the post-task rest period (third row)
            % 
            % Outputs
            %   - rippleStart:                    vector of timestamps marking the start of sharp-wave ripples that meet the criteria
            %   - rippleStop:                    vector of timestamps marking the end of sharp-wave ripples that meet the criteria
            %   - ripplePeak:                    vector of timestamps marking the time of peak amplitude of sharp-wave ripples that meet the criteria
            
            rippleStart = data.rippleStart;
            rippleStop = data.rippleStop;
            ripplePeak = data.ripplePeak;
            
            % Skip if nRipples config parameter is NaN
            if isnan(obj.analysisConfig.nRipples)
                return
            end
            
            % Ensure ripple data is vectors
            assert(min(size(rippleStart)) == 1, 'ExplainedVarianceReactivation:invalidData', 'Ripple start times must be a vector')
            assert(min(size(rippleStop)) == 1, 'ExplainedVarianceReactivation:invalidData', 'Ripple end times must be a vector')
            assert(min(size(ripplePeak)) == 1, 'ExplainedVarianceReactivation:invalidData', 'Ripple peak times must be a vector')
            
            % Ensure ripple data vectors are the same length
            assert(numel(rippleStart)==numel(rippleStop), 'ExplainedVarianceReactivation:invalidData', 'Vectors for the start and end times of ripples must be the same length')
            assert(numel(rippleStart)==numel(ripplePeak), 'ExplainedVarianceReactivation:invalidData', 'Vectors for the start and peak times of ripples must be the same length')
            
            % Ensure nRipples parameter is valid
            assert(isnumeric(obj.analysisConfig.nRipples), 'ExplainedVarianceReactivation:invalidConfigParameter', 'nRipples parameter must be a single integer or a matrix of two integers')
            assert(all(obj.analysisConfig.nRipples == ceil(obj.analysisConfig.nRipples)), 'ExplainedVarianceReactivation:invalidConfigParameter', 'nRipples parameter must be a single integer or a matrix of two integers')
            assert(numel(obj.analysisConfig.nRipples)==1 || numel(obj.analysisConfig.nRipples)==2, 'ExplainedVarianceReactivation:invalidConfigParameter', 'nRipples parameter must be a single integer or a matrix of two integers')
            
            if length(obj.analysisConfig.nRipples)==1
                % Restrict to all the ripples of PRE and the first n ripples of POST
                keep = false(1, length(ripplePeak));
                keep(find(rippleStop < data.startEndTimes(1, 2), obj.analysisConfig.nRipples, 'last')) = true;
                assert(sum(keep) > 0, 'ExplainedVarianceReactivation:invalidData', 'There are no ripples during PRE.')
                keep(find(rippleStart > data.startEndTimes(3, 1), obj.analysisConfig.nRipples, 'first')) = true;
                assert(sum(keep) > 0, 'ExplainedVarianceReactivation:invalidData', 'There are no ripples during POST.')
                rippleStart = rippleStart(keep);
                rippleStop = rippleStop(keep);
                ripplePeak = ripplePeak(keep);
                if length(rippleStart) < obj.analysisConfig.nRipples * 2
                    warning('There are not enough ripples in PRE and POST to satisfy analysis of %i ripples', obj.analysisConfig.nRipples)
                end
                assert(length(rippleStart) <= obj.analysisConfig.nRipples * 2, '%i ripples included but analysis configuration specifies %i', length(rippleStart), obj.analysisConfig.nRipples * 2)

            elseif length(obj.analysisConfig.nRipples)==2
                % Restrict to last n(1) ripples of PRE and first n(2) ripples of POST
                keep_pre = false(1, length(ripplePeak));
                keep_post = false(1, length(ripplePeak));
                keep_pre(find(rippleStop < data.startEndTimes(1, 2), obj.analysisConfig.nRipples(1), 'last')) = true;
                keep_post(find(rippleStart > data.startEndTimes(3, 1), obj.analysisConfig.nRipples(2), 'first')) = true;
                fprintf('%i ripples from PRE and %i ripples from POST\n', sum(keep_pre), sum(keep_post))
                assert(sum(keep_pre) > 0, 'ExplainedVarianceReactivation:invalidData', 'There are no ripples during PRE.')
                assert(sum(keep_post) > 0, 'ExplainedVarianceReactivation:invalidData', 'There are no ripples during POST.')
                keep = find(keep_pre | keep_post);
                rippleStart = rippleStart(keep);
                rippleStop = rippleStop(keep);
                ripplePeak = ripplePeak(keep);
                if length(rippleStart) < sum(obj.analysisConfig.nRipples)
                    warning('There are not enough ripples in PRE and POST to satisfy analysis of %i PRE and %i POST ripples', obj.analysisConfig.nRipples(1), obj.analysisConfig.nRipples(2))
                end
                assert(length(rippleStart) <= sum(obj.analysisConfig.nRipples), '%i ripples included but analysis configuration specifies %i', length(rippleStart), sum(obj.analysisConfig.nRipples))

            end

        end
        
        % Bin session-wide spiketimes
        function [binnedfrPRE, binnedfrPOST, binnedfrTASK] = binActivity(obj, data, timestamps, analysisConfig)
            
            % Filters ripples in PRE and POST-task rest according to the number specified by inclusion criteria
            % 
            % Inputs
            %   - data:                            structure containing the following fields
            %       - spiketimes:               cell array of spiketimes in seconds; one cell per neuron
            %       - binnedfr:                   pre-computed firing rates in 50ms bins,
            %                                        calculated from the spiketimes cell array; one row per unit
            %   - timestamps:                  structure containing the following fields
            %       - startEndTimes:          3 x 2 array of timestamps in seconds,
            %                                        marking the start and end of the pre-task rest period (first row), 
            %                                        start and end of the task period (second row), 
            %                                        and start and end of the post-task rest period (third row)
            %       - rippleStart:                timestamps in seconds marking the start of each identified sharp-wave ripple
            %       - rippleStop:                timestamps in seconds marking the end of each identified sharp-wave ripple
            %       - ripplePeak:                timestamps in seconds marking the peak amplitude of each identified sharp-wave ripple
            %   - analysisConfig:             structure containing the following fields
            %       - postDuration:            total duration of activity during POST to analyse EV-REV for, in seconds
            %       - rippleActivityMode:    'none', 'variable', or 'equalise'; if 'none', all spiking activity during PRE and POST is used for EV-REV analysis
            %                                        (within and between sharp-wave ripples); if 'variable', spikes fired between the start and end of
            %                                        every ripple are used; if 'equalise', spikes fired during a window of 200ms from the ripple peak is used
            % 
            % Outputs
            %   - binnedfrPRE:                firing rates in 50ms bins for the PRE-task epoch, computed from the spiketimes
            %   - binnedfrPOST:              firing rates in 50ms bins for the TASK epoch
            %   - binnedfrTASK:              firing rates in 50ms bins for the POST-task epoch, computed from the spiketimes
            
            % Process inputs
            dataFields = {'spiketimes', 'binnedfr'};
            timestampFields = {'startEndTimes', 'rippleStart', 'rippleStop', 'ripplePeak'};
            analysisconfigFields = {'postDuration', 'rippleActivityMode'};
            for field = dataFields
                fieldName = field{1};
                if isfield(data, fieldName)
                    eval([fieldName ' = data.' fieldName ';']);
                else
                    eval([fieldName ' = [];']);
                end
            end
            for field = timestampFields
                fieldName = field{1};
                if isfield(timestamps, fieldName)
                    eval([fieldName ' = timestamps.' fieldName ';']);
                else
                    eval([fieldName ' = [];']);
                end
            end
            for field = analysisconfigFields
                fieldName = field{1};
                if isfield(analysisConfig, fieldName)
                    eval([fieldName ' = analysisConfig.' fieldName ';']);
                else
                    eval([fieldName ' = [];']);
                end
            end
    
            % Binsize in seconds
            binsizePrePost = obj.BinSize;
            binsizeTask = obj.BinSize;
            
            % Convolve the spiketimes. Three options: equalise the time around each
            % ripple, use the calculated ripple length, or ignore ripples
            % altogether and use the full duration
            switch rippleActivityMode

                % Don't restrict to ripples at all
                case 'none'
                    if isnan(postDuration)
                        binnedfrPRE = squeeze(obj.convolveAndBin(startEndTimes(1, 1), diff(startEndTimes(1, :)), binsizePrePost, spiketimes));
                        binnedfrPOST = squeeze(obj.convolveAndBin(startEndTimes(3, 1), diff(startEndTimes(3, :)), binsizePrePost, spiketimes));
                    else
                        binnedfrPRE = squeeze(obj.convolveAndBin(startEndTimes(1, 2), diff(startEndTimes(1, :)), binsizePrePost, spiketimes));
                        binnedfrPOST = squeeze(obj.convolveAndBin(startEndTimes(3, 1), postDuration, binsizePrePost, spiketimes));
                    end
                    clear ripplePeakPRE rippleStartPRE rippleStopPRE ripplePeakPOST rippleStartPOST rippleStopPOST
            
                % Restrict the ripples to the last X minutes of PRE and the first X minutes of POST
                otherwise
                    if isnan(postDuration)
                        includePRE = (rippleStart >= (startEndTimes(1, 1))) & (rippleStop <= startEndTimes(1, 2));
                        includePOST = (rippleStart >= startEndTimes(3, 1)) & (rippleStop <= (startEndTimes(3, 2)));
                    else
                        includePRE = (rippleStart >= (startEndTimes(1, 1))) & (rippleStop <= startEndTimes(1, 2));
                        includePOST = (rippleStart >= startEndTimes(3, 1)) & (rippleStop <= (startEndTimes(3, 1) + postDuration));
                    end

                    switch rippleActivityMode

                        % Take 200ms of activity around each ripple peak
                        case 'equalise'            
                            
                            ripplePeakPRE = ripplePeak(includePRE);
                            ripplePeakPOST = ripplePeak(includePOST);
                            rippleStopPRE = rippleStop(includePRE);
                            rippleStopPOST  = rippleStop(includePOST);
                            rippleStartPRE = rippleStart(includePRE);
                            rippleStartPOST = rippleStart(includePOST);
                            
                            timeBefore = 0.00; timeAfter = 0.20;
                            if length(rippleStartPRE) == 1 || length(spiketimes) == 1
                                dontPermute = true;
                            else
                                dontPermute = false;
                            end
                    
                            % PRE ripples
                            nRipples = length(rippleStartPRE);
                            binnedfrPRE = cell(1, nRipples);
                            iRip = 1; 
                            while iRip <= nRipples

                                windowStart = rippleStartPRE(iRip) - timeBefore;
                                windowStop = rippleStartPRE(iRip) + timeAfter;

                                % Merge any other ripples within timeAfter to the same window
                                while iRip < nRipples && (rippleStartPRE(iRip + 1) - timeBefore) <= windowStop
                                    iRip = iRip + 1;
                                    windowStop = max(windowStop, rippleStartPRE(iRip) + timeAfter);
                                end

                                % Convolve spiketimes for the window of one or more ripples
                                windowDuration = windowStop - windowStart;
                                popV = obj.convolveAndBin(windowStart, windowDuration, binsizePrePost, spiketimes);
                                if ndims(popV) == 3 && ~dontPermute
                                    binnedfrPRE{iRip} = reshape(permute(popV, [2, 3, 1]), length(spiketimes), []);
                                else
                                    binnedfrPRE{iRip} = squeeze(popV)';
                                end

                                iRip = iRip + 1;
                            end
                            binnedfrPRE = cat(2, binnedfrPRE{1:iRip-1});
                    
                            % POST ripples
                            nRipples = length(rippleStartPOST);
                            binnedfrPOST = cell(1, nRipples);
                            iRip = 1; 
                            while iRip <= nRipples

                                windowStart = rippleStartPOST(iRip) - timeBefore;
                                windowStop = rippleStartPOST(iRip) + timeAfter;

                                % Merge any other ripples within timeAfter to the same window
                                while iRip < nRipples && (rippleStartPOST(iRip + 1) - timeBefore) <= windowStop
                                    iRip = iRip + 1;
                                    windowStop = max(windowStop, rippleStartPOST(iRip) + timeAfter);
                                end

                                % Convolve spiketimes for the window of one or more ripples
                                windowDuration = windowStop - windowStart;
                                popV = obj.convolveAndBin(windowStart, windowDuration, binsizePrePost, spiketimes);
                                if ndims(popV) == 3 && ~dontPermute
                                    binnedfrPOST{iRip} = reshape(permute(popV, [2, 3, 1]), length(spiketimes), []);
                                else
                                    binnedfrPOST{iRip} = squeeze(popV)';
                                end

                                iRip = iRip + 1;
                            end
                            binnedfrPOST = cat(2, binnedfrPOST{1:iRip-1});

                            clear popV ripplePeakPRE rippleStartPRE rippleStopPRE ripplePeakPOST rippleStartPOST rippleStopPOST

                        % Use the true ripple duration 
                        case 'variable'

                            rippleStartPRE = rippleStart(includePRE);
                            rippleStopPRE = rippleStop(includePRE);
                            rippleStartPOST = rippleStart(includePOST);
                            rippleStopPOST = rippleStop(includePOST);

                            if length(rippleStartPRE) > length(rippleStopPRE)
                                rippleStopPRE = [rippleStopPRE; startEndTimes(1,2)];
                            end
                            if length(rippleStartPOST) > length(rippleStopPOST)
                                rippleStopPRE = [rippleStopPOST; startEndTimes(1,2)];
                            end

                            binnedfrPRE = [];
                            for i = 1:length(rippleStartPRE)
                                popV = obj.convolveAndBin(rippleStartPRE(i), rippleStopPRE(i) - rippleStartPRE(i), binsizePrePost, spiketimes);
                                binnedfrPRE = [binnedfrPRE reshape(permute(popV, [2, 3, 1]), length(spiketimes), [])];
                            end
                            binnedfrPOST = [];
                            for i = 1:length(rippleStartPOST)
                                popV = obj.convolveAndBin(rippleStartPOST(i), rippleStopPOST(i) - rippleStartPOST(i), binsizePrePost, spiketimes);
                                binnedfrPOST = [binnedfrPOST reshape(permute(popV, [2, 3, 1]), length(spiketimes), [])];
                            end
                            clear popV ripplePeakPRE rippleStartPRE rippleStopPRE ripplePeakPOST rippleStartPOST rippleStopPOST

                    end

            end

            % Convolve the spiketimes during TASK
            if numel(binnedfr) > 0
                binnedfrTASK = binnedfr(:,ceil(startEndTimes(2,1)/binsizeTask)+1:ceil(startEndTimes(2,2)/binsizeTask));
            else
                binnedfrTASK = squeeze(obj.convolveAndBin(startEndTimes(2, 1), diff(startEndTimes(2, :)), binsizeTask, spiketimes));
            end
            clear binnedFR
            
        end
        
        
        % Calculate explained variance
        function [EV, REV, binnedfrPRE, binnedfrPOST, binnedfrTASK, preCC, taskCC, postCC] = explainedVariance(obj, data, timestamps, analysisConfig)
            
            % Calculates the explained variance and reverse explained variance for one recording session
            % 
            % Inputs
            %   - data:                            structure containing the following fields
            %       - spiketimes:               optional; cell array of spiketimes in seconds for one session; one cell per neuron. If calculating 
            %                                        explained variance within one brain area, must contain spiketimes only from that area.
            %                                        Must be provided if any of binnedfrPRE, binnedfrPOST and binnedfrTASK are not
            %       - binnedfr:                   optional; pre-computed firing rates in 50ms bins, calculated from the
            %                                        spiketimes cell array; one row per unit. If calculating  explained variance within one brain area,
            %                                        must contain data only from that area. Makes computation of binnedfrTASK faster if provided
            %       - binnedfrPRE:            optional; pre-computed firing rates in 50ms bins, calculated from the spiketimes cell
            %                                       array; one row per unit. If calculating explained variance within one brain area, must contain data only
            %                                       from that area. Can be calculated from the entire PRE epoch or a subset of it (e.g. concatenated ripple times)
            %       - binnedfrTASK:          optional; pre-computed firing rates in 50ms bins, calculated from the spiketimes cell
            %                                       array; one row per unit. If calculating explained variance within one brain area, must contain data only
            %                                       from that area. Can be calculated from the entire TASK epoch or a subset of it
            %       - binnedfrPOST           optional; pre-computed firing rates in 50ms bins, calculated from the spiketimes cell
            %                                       array; one row per unit. If calculating explained variance within one brain area, must contain data only
            %                                       from that area. Can be calculated from the entire POST epoch or a subset of it (e.g. concatenated ripple times)
            %       - preCC:                     optional; matrix of correlation coefficients between neurons during PRE epoch,
            %                                        of size NxN or NxM or MxM where N is number of ventral striatal neurons and M is number 
            %                                        of CA1 neurons. Makes computation faster if provided
            %       - taskCC:                   optional; matrix of correlation coefficients between neurons during TASK epoch,  
            %                                        of same size as preCC. Makes computation faster if provided
            %       - postCC:                   optional; matrix of correlation coefficients between neurons during POST epoch; 
            %                                        of same size as preCC. Makes computation faster if provided
            %       - nNAcUnits:               optional; number of ventral striatal neurons in the spiketimes and binnedfr data 
            %                                        if calculating explained variance between brain areas, or empty if within brain area.
            %                                        Required if preCC, taskCC and postCC are not provided
            %   - timestamps:                  structure containing the following fields
            %       - startEndTimes:          optional; 3 x 2 array of timestamps in seconds,
            %                                        marking the start and end of the pre-task rest period (first row), 
            %                                        start and end of the task period (second row), 
            %                                        and start and end of the post-task rest period (third row). Required if any of
            %                                        binnedfrPRE, binnedfrPOST and binnedfrTASK are not provided 
            %       - rippleStart:                optional; timestamps in seconds marking the start of each identified sharp-wave ripple.
            %                                        Required if any ofbinnedfrPRE, binnedfrPOST and binnedfrTASK are not provided 
            %       - rippleStop:                optional; timestamps in seconds marking the end of each identified sharp-wave ripple. 
            %                                        Required if any of binnedfrPRE, binnedfrPOST and binnedfrTASK are not provided 
            %       - ripplePeak:                optional; timestamps in seconds marking the peak amplitude of each identified sharp-wave ripple. 
            %                                        Required if any of binnedfrPRE, binnedfrPOST and binnedfrTASK are not provided 
            %   - analysisConfig:             structure containing the following fields
            %       - postDuration:            optional; total duration of activity during POST to analyse EV-REV for, in seconds.
            %                                        Required if any of binnedfrPRE, binnedfrPOST and binnedfrTASK are not provided 
            %       - rippleActivityMode:    optional; 'none', 'variable', or 'equalise'; if 'none', all spiking activity during PRE and POST is used for EV-REV analysis
            %                                        (within and between sharp-wave ripples); if 'variable', spikes fired between the start and end of
            %                                        every ripple are used; if 'equalise', spikes fired during a window of 200ms from the ripple peak is used.
            %                                        Required if any of binnedfrPRE, binnedfrPOST and binnedfrTASK are not provided 
            %      - exclude:                    optional; Nx2 matrix of indices of neurons to exclude from the explained
            %                                        variance calculation
            % 
            % Outputs
            %   - EV:                              explained variance for the session
            %   - REV:                           reverse explained variance for the session
            %   - binnedfrPRE:                firing rates in 50ms bins for the PRE-task epoch, computed from the spiketimes
            %   - binnedfrPOST:              firing rates in 50ms bins for the TASK epoch
            %   - binnedfrTASK:              firing rates in 50ms bins for the POST-task epoch, computed from the spiketimes
            %   - preCC:                        matrix of correlation coefficients between neurons during PRE epoch
            %   - taskCC:                       matrix of correlation coefficients between neurons during TASK epoch
            %   - postCC:                       matrix of correlation coefficients between neurons during POST epoch

            function corrCoeff = nancorr2(A, B)
                % Ensure A and B are the same size
                assert(all(size(A) == size(B)), 'ExplainedVarianceReactivation:invalidData', 'Matrices A and B must be the same size.');

                % Create mask for non-NaN values in both A and B
                validMask = ~isnan(A) & ~isnan(B);

                % Extract values from A and B that are not NaN in either A or B
                validA = A(validMask);
                validB = B(validMask);

                % Calculate means
                meanA = mean(validA);
                meanB = mean(validB);

                % Subtract means from validA and validB
                devA = validA - meanA;
                devB = validB - meanB;

                % Calculate standard deviations
                stdA = std(validA);
                stdB = std(validB);

                % Compute correlation
                corrCoeff = sum(devA .* devB) / (length(validA) - 1) / (stdA * stdB);
            end
            
            % If binned firing rates are not provided, create them
            if (~isfield(data, 'binnedfrPRE')) || (~isfield(data, 'binnedfrPOST')) || (~isfield(data, 'binnedfrTASK'))

                % Check inputs
                if (~isfield(data, 'spiketimes')) || (~isfield(analysisConfig, 'postDuration'))
                    error('Must provide either binned firing rates for each task epoch, or spiketimes, ripple times and duration to calculate them')
                elseif ((~isfield(timestamps, 'rippleStart')) ||  (~isfield(timestamps, 'rippleStop'))) && strcmp(analysisConfig.rippleActivityMode, 'variable')
                    error('To extract activity during ripples with variable duration, you need to provide ripple start and stop times')
                elseif (~isfield(timestamps, 'ripplePeak')) && strcmp(analysisConfig.rippleActivityMode, 'equalise')
                    error('To extract activity around ripple peaks with equal duration, you need to provide ripple peak times')
                else
                    [data.binnedfrPRE, data.binnedfrPOST, data.binnedfrTASK] = obj.binActivity(data, timestamps, analysisConfig);
                end

            end
            
            % Check that binned firing rates contain the same number of
            % neurons
            nNeuronsPre = size(data.binnedfrPRE, 1);
            nNeuronsTask = size(data.binnedfrTASK, 1);
            nNeuronsPost = size(data.binnedfrPOST, 1);
            assert(nNeuronsPre==nNeuronsTask, 'ExplainedVarianceReactivation:invalidData', 'Provided matrices for PRE and TASK binned firing rates have incompatible sizes: %i and %i rows', nNeuronsPre, nNeuronsTask)
            assert(nNeuronsPre==nNeuronsPost, 'ExplainedVarianceReactivation:invalidData', 'Provided matrices for PRE and POST binned firing rates have incompatible sizes: %i and %i rows', nNeuronsPre, nNeuronsPost)
            assert(nNeuronsTask==nNeuronsPost, 'ExplainedVarianceReactivation:invalidData', 'Provided matrices for TASK and POST binned firing rates have incompatible sizes: %i and %i rows', nNeuronsTask, nNeuronsPost)

            % If no data for any of the splits, return NaN
            if numel(data.binnedfrPRE)==0 || numel(data.binnedfrTASK)==0 || numel(data.binnedfrPOST)==0
                EV = NaN;
                REV = NaN;
                binnedfrPRE = data.binnedfrPRE;
                binnedfrTASK = data.binnedfrTASK;
                binnedfrPOST = data.binnedfrPOST;
                preCC = NaN;
                taskCC = NaN;
                postCC = NaN;
                return
            end

            % Check if analysis should be done between brain areas or within
            if isempty(data.nNAcUnits)
                intraregion = 1;
            else
                intraregion = 0;
            end

            % Intra-region correlation matrices
            if isfield(data, 'preCC') && isfield(data, 'taskCC') && isfield(data, 'postCC')
                preCC = data.preCC;
                taskCC = data.taskCC;
                postCC = data.postCC;

            else
                if intraregion
                    preCC = corrcoef(data.binnedfrPRE'); preCC(isnan(preCC)) = 0;
                    taskCC = corrcoef(data.binnedfrTASK'); taskCC(isnan(taskCC)) = 0;
                    postCC = corrcoef(data.binnedfrPOST'); postCC(isnan(postCC)) = 0;

                elseif numel(data.binnedfrPRE) == 0 || numel(data.binnedfrPOST) == 0 || numel(data.binnedfrTASK) == 0
                    preCC = NaN;
                    taskCC = NaN;
                    postCC = NaN;

                % Inter-region correlation matrices
                else
                    nNAcUnits = data.nNAcUnits;
                    nHPCUnits = size(data.binnedfrPRE, 1) - nNAcUnits;
                    preCC = zeros(nNAcUnits, nHPCUnits);
                    taskCC = zeros(nNAcUnits, nHPCUnits);
                    postCC = zeros(nNAcUnits, nHPCUnits);

                    for iNAc = 1:nNAcUnits
                        for iHPC = 1:nHPCUnits
                            indexHPC = iHPC + nNAcUnits; 
                            c = corrcoef(data.binnedfrPRE(iNAc,:), data.binnedfrPRE(indexHPC,:));
                            preCC(iNAc, iHPC) = c(2);
                            c = corrcoef(data.binnedfrTASK(iNAc,:), data.binnedfrTASK(indexHPC,:));
                            taskCC(iNAc, iHPC) = c(2);
                            c = corrcoef(data.binnedfrPOST(iNAc,:), data.binnedfrPOST(indexHPC,:));
                            postCC(iNAc, iHPC) = c(2);
                        end
                    end

                end

                preCC(isnan(preCC)) = 0;
                taskCC(isnan(taskCC)) = 0;
                postCC(isnan(postCC)) = 0;

            end

            % Exclude pairs if needed
            if isfield(analysisConfig, 'exclude')
                for i = 1:size(analysisConfig.exclude,1)
                    preCC(analysisConfig.exclude(i,1),analysisConfig.exclude(i,2)) = NaN;
                    postCC(analysisConfig.exclude(i,1),analysisConfig.exclude(i,2)) = NaN;
                    taskCC(analysisConfig.exclude(i,1),analysisConfig.exclude(i,2)) = NaN;
                end
            end

            % Get correlation between PRE window, TASK and POST window (only for units that fired in both windows)
            taskPreCC = nancorr2(taskCC,preCC);
            taskPostCC = nancorr2(taskCC,postCC);
            prePostCC = nancorr2(preCC,postCC);

            % Compute EV
            EV = ( (taskPostCC - taskPreCC * prePostCC) / sqrt( (1 - taskPreCC.^2) * (1 - prePostCC.^2) ) ).^2;

            % Compute REV
            REV = ( (taskPreCC - taskPostCC * prePostCC) / sqrt( (1 - taskPostCC.^2) * (1 - prePostCC.^2) ) ).^2;

            if isinf(EV) || isinf(REV)
                EV = 0;
                REV = 0;
            end

            binnedfrPRE = data.binnedfrPRE;
            binnedfrPOST = data.binnedfrPOST;
            binnedfrTASK = data.binnedfrTASK;

        end
                
        
        % Calculate contribution of each cell-pair to EV-REV
        function [perPairContribution] = calculatePerCellPairContributions(obj, data, includeCells, intraarea, EV_REV0)
            
            % Removes each cell pair in turn, recalculates explained variance and reverse explained variance, and uses the
            % difference to define the cell pair's contribution to EV-REV reactivation
            % 
            % Inputs
            %   - data:                            structure containing the following fields
            %       - spiketimes:               cell array of spiketimes in seconds for one session; one cell per neuron
            %       - binnedfrPRE:            matrix of pre-computed firing rates in 50ms bins, calculated from the spiketimes cell
            %                                        array; one row per unit. If calculating explained variance within one brain area, must contain data only
            %                                        from that area. Can be calculated from the entire PRE epoch or a subset of it (e.g. concatenated ripple times)
            %       - binnedfrTASK:           matrix of pre-computed firing rates in 50ms bins, calculated from the spiketimes cell
            %                                        array; one row per unit. If calculating explained variance within one brain area, must contain data only
            %                                        from that area. Can be calculated from the entire TASK epoch or a subset of it
            %       - binnedfrPOST           matrix of pre-computed firing rates in 50ms bins, calculated from the spiketimes cell
            %                                        array; one row per unit. If calculating explained variance within one brain area, must contain data only
            %                                        from that area. Can be calculated from the entire POST epoch or a subset of it (e.g. concatenated ripple times)
            %       - nNAcUnits:               empty matrix for calculations of reactivation within a brain area, or
            %                                        number of ventral striatal units, corresponding to the
            %                                        first cells of spiketimes and the first rows of binnedfr, for reactivation between brain areas
            %       - nUnits                      total number of neurons in the session (must be at least nNAcUnits and equal to
            %                                        the length of spiketimes and the first dimension of binnedfr)
            %   - includeCells:                 vector of integers of length nUnits, containing the indices of ventral striatal
            %                                        neurons to include in the analysis
            %   - intraarea:                      boolean; true for calculating explained variance within one brain area, false
            %                                        for between brain areas
            %   - EV_REV0:                    EV-REV for the whole session including all neurons; this is used to compute the
            %                                       contribution of each cell pair
            % 
            % Outputs
            %   - perPairContribution:       1xN struct for N cell pairs with the following fields
            %       - cells:                       pairs of indices of cells to which the per-cell contribution corresponds
            %       - contribution:             the per-cell contribution, relative to EV_REV0
            %       - correlations:             3x1 matrix of the correlation coefficient between the two cells' binned firing rate vectors for PRE, TASK and POST respectively

            
            % Check data inputs
            assert(isstruct(data), 'ExplainedVarianceReactivation:invalidData', 'Argument  "data" must be a struct');
            requiredDataFields = {'nUnits', 'nNAcUnits', 'spiketimes', 'binnedfrPRE', 'binnedfrTASK', 'binnedfrPOST'};
            for field = requiredDataFields
                assert(isfield(data, field{1}), 'ExplainedVarianceReactivation:dataNotAssigned', 'Argument "data" must contain a field for %s', field{1});
            end
            assert(numel(data.spiketimes) == data.nUnits, 'ExplainedVarianceReactivation:invalidData', 'Spiketimes data must contain one cell per neuron');
            assert(size(data.binnedfrPRE, 1) == data.nUnits, 'ExplainedVarianceReactivation:invalidData', 'Binned firing rate data for PRE must contain one cell per neuron');
            assert(size(data.binnedfrTASK, 1) == data.nUnits, 'ExplainedVarianceReactivation:invalidData', 'Binned firing rate data for TASK must contain one cell per neuron');
            assert(size(data.binnedfrPOST, 1) == data.nUnits, 'ExplainedVarianceReactivation:invalidData', 'Binned firing rate data for POST must contain one cell per neuron');
            assert(length(includeCells) <= data.nNAcUnits, 'ExplainedVarianceReactivation:invalidData', 'Argument "includeCells" must contain no more than the number of accumbens units');
            assert(all(max(includeCells) <= data.nNAcUnits), 'ExplainedVarianceReactivation:invalidData', 'Argument "includeCells" must contain no more than the number of accumbens units');
            assert(length(unique(includeCells)) == length(includeCells), 'ExplainedVarianceReactivation:invalidData', 'Argument "includeCells" must not have repeating values');
            assert(islogical(intraarea))
            data.nCA1Units = data.nUnits - data.nNAcUnits;
            if intraarea
                assert(all(size(data.preCC)==[data.nNAcUnits, data.nNAcUnits]), 'ExplainedVarianceReactivation:invalidData', 'If analysing intra-area cell pairs, correlation coefficient matrix for PRE must be a square matrix of order nNAcUnits');
                assert(all(size(data.taskCC)==[data.nNAcUnits, data.nNAcUnits]), 'ExplainedVarianceReactivation:invalidData', 'If analysing intra-area cell pairs, correlation coefficient matrix for TASK must be a square matrix of order nNAcUnits');
                assert(all(size(data.postCC)==[data.nNAcUnits, data.nNAcUnits]), 'ExplainedVarianceReactivation:invalidData', 'If analysing intra-area cell pairs, correlation coefficient matrix for POST must be a square matrix of order nNAcUnits');
            else
                assert(all(size(data.preCC)==[data.nNAcUnits, data.nCA1Units]), 'ExplainedVarianceReactivation:invalidData', 'If analysing intra-area cell pairs, correlation coefficient matrix for PRE must be a matrix of size nNAcUnits x nCA1Units');
                assert(all(size(data.taskCC)==[data.nNAcUnits, data.nCA1Units]), 'ExplainedVarianceReactivation:invalidData', 'If analysing intra-area cell pairs, correlation coefficient matrix for TASK must be a matrix of size nNAcUnits x nCA1Units');
                assert(all(size(data.postCC)==[data.nNAcUnits, data.nCA1Units]), 'ExplainedVarianceReactivation:invalidData', 'If analysing intra-area cell pairs, correlation coefficient matrix for POST must be a matrix of size nNAcUnits x nCA1Units');
            end
            
            perPairContribution = struct('cells', {}, 'contribution', {}, 'correlations', {});
            
            % Get correlation matrices between all pairs of units
            if intraarea
                nNAcUnits = [];
                spiketimes = data.spiketimes(1:data.nNAcUnits);
                binnedfrPRE = data.binnedfrPRE(1:data.nNAcUnits, :);
                binnedfrPOST = data.binnedfrPOST(1:data.nNAcUnits, :);
                binnedfrTASK = data.binnedfrTASK(1:data.nNAcUnits, :);
            else
                nNAcUnits = data.nNAcUnits;
                spiketimes = data.spiketimes;
                binnedfrPRE = data.binnedfrPRE;
                binnedfrPOST = data.binnedfrPOST;
                binnedfrTASK = data.binnedfrTASK;
            end

           EVData = struct('spiketimes', {spiketimes},...
                'nNAcUnits', nNAcUnits,...
                'binnedfrPRE', binnedfrPRE,...
                'binnedfrPOST', binnedfrPOST,...
                'binnedfrTASK', binnedfrTASK);
            EVTimestamps = struct();
            EVAnalysisConfig = struct('exclude', [],...
                'duration', duration,...
                'rippleActivityMode', obj.analysisConfig.rippleActivityMode);
            preCC = data.preCC;
            taskCC = data.taskCC;
            postCC = data.postCC;
                           
            % Iterate through cell pairs
            if intraarea
               [Y, X] = ndgrid(includeCells, includeCells);
               exclude_cellPairs = [X(:), Y(:)];

           else
                [Y, X] = ndgrid(1:data.nCA1Units, includeCells);
                exclude_cellPairs = [X(:), Y(:)];

           end
                
            n_exclude = size(exclude_cellPairs, 1);
            for i = 1:n_exclude
                pair = exclude_cellPairs(i, :);

                % Remove pair-to-be-excluded from correlation matrices
                EVData.preCC = preCC;
                EVData.taskCC = taskCC;
                EVData.postCC = postCC;
                EVData.preCC(pair(1), pair(2)) = NaN;
                EVData.taskCC(pair(1), pair(2)) = NaN;
                EVData.postCC(pair(1), pair(2)) = NaN;

                % Get difference in EV-REV without this pair
                [ev, rev] = obj.explainedVariance(EVData, EVTimestamps, EVAnalysisConfig);
                perPairContribution(i).cells = pair;
                perPairContribution(i).contribution = (EV_REV0) - (ev-rev);
                perPairContribution(i).correlations = [preCC(pair(1), pair(2)); taskCC(pair(1), pair(2)); postCC(pair(1), pair(2))];

           end
           
        end
        
        % Restrict to immobile periods
        function [restStartPRE, restStopPRE, restStartPOST, restStopPOST] = getFixedDurationImmobility(obj, data)
            
            % Returns timestamps for the start and end of bouts of immobility (rest) that do not exceed a specificed duration
            % limit when concatenated
            % 
            % Inputs
            %   - data:                            structure containing the following fields
            %       - startEndTimes:          3 x 2 array of timestamps in seconds,
            %                                        marking the start and end of the pre-task rest period (first row), 
            %                                        start and end of the task period (second row), 
            %                                        and start and end of the post-task rest period (third row)
            %       - restStart:                  vector of timestamps in seconds marking the start of each identified bout of rest
            %       - restStop:                  vector of timestamps in seconds marking the end of each identified bout of rest
            % 
            % Outputs
            %   - restStartPRE:                vector of timestamps in seconds marking the start of each bout of rest in the
            %                                        PRE-task epoch that combined do not exceed the specified duration
            %   - restStopPRE:                vector of timestamps in seconds marking the end of each bout of rest in the
            %                                        PRE-task epoch that combined do not exceed the specified duration
            %   - restStartPOST:              vector of timestamps in seconds marking the start of each bout of rest in the
            %                                        POST-task epoch that combined do not exceed the specified duration
            %   - restStopPOST:              vector of timestamps in seconds marking the start of each bout of rest in the
            %                                        POST-task epoch that combined do not exceed the specified duration

            
            if isnan(obj.analysisConfig.postDuration)
                obj.analysisConfig.postDuration = inf;
            end
        
            % Find rest periods during desired duration: PRE-task rest
            restPRE = find((data.restStart < data.startEndTimes(1, 2)) & (data.restStop > data.startEndTimes(1, 1)));
            restDuration = 0;
            restStart = [];
            restStop = [];
            for iRest = 0:length(restPRE)-1
                immobilePeriodDuration = data.restStop(restPRE(end-iRest)) - data.restStart(restPRE(end-iRest));
                restDuration = restDuration + immobilePeriodDuration;
                restStart = [data.restStart(restPRE(end-iRest)) restStart];
                restStop = [data.restStop(restPRE(end-iRest)) restStop];
                if restDuration > obj.analysisConfig.postDuration
                    break
                end
            end

            restStartPRE = restStart;
            restStopPRE = restStop;

            % Find rest periods during desired duration: POST-task rest
            restPOST = find((data.restStart < data.startEndTimes(3, 2)) & (data.restStop > data.startEndTimes(3, 1)));
            restDuration = 0;
            restStart = [];
            restStop = [];
            for iRest = 1:length(restPOST)
                immobilePeriodDuration = data.restStop(restPOST(iRest)) - data.restStart(restPOST(iRest));
                restDuration = restDuration + immobilePeriodDuration;
                restStart = [restStart data.restStart(restPOST(iRest))];
                restStop = [restStop data.restStop(restPOST(iRest))];
                if restDuration > obj.analysisConfig.postDuration
                    break
                end
            end
            if restDuration < obj.analysisConfig.postDuration && ~isinf(obj.analysisConfig.postDuration)
                warning(['total pre-task immobility of ' num2str(restDuration) 's is less than desired duration of ' num2str(obj.analysisConfig.postDuration) 's'])
            end
            
            % Clip end of immobility period
            if ~isinf(obj.analysisConfig.postDuration) && restDuration > obj.analysisConfig.postDuration
                restStop(end) = restStop(end) - (restDuration - obj.analysisConfig.postDuration);
                if restStop(end) <= restStart(end)
                    restStop = restStop(1:end-1);
                    restStart = restStart(1:end-1);
                end
                assert (round(sum(restStop - restStart), 2) == round(obj.analysisConfig.postDuration, 2))
            end
            restStartPOST = restStart;
            restStopPOST = restStop;
            
        end

        
        
        % Calculate EV-REV for inter-area (CA1-vStr) pairs
        function obj = calculateInterAreaEVREV(obj, varargin)
            
            % Calculates the explained variance and reverse explained variance between CA1 and ventral striatum across recording sessions
            % 
            % Inputs
            %   - dataPath:            optional; specifies the path where .mat files of electrophysiological and related data are stored

            
            % Process paths where data is stored
            p = inputParser;
            addParameter(p, 'dataPath', fullfile('.', 'data', 'ephys'), @(x) ischar(x) || isstring(x));
            parse(p, varargin{:});
            dataPath = p.Results.dataPath;
            
            % Check there are sessions to analyse
            assert(~isempty(obj.sessionsForBroadAnalysis), 'ExplainedVarianceReactivation:sessions:noSessions', 'No sessions selected for analysis. Please revise the inclusion criteria')
            
             % Struct for storing EV and REV scores
            EV = cell2struct(cell(1, length(obj.inclusionCriteria.rats)), cellfun(@(x) x, obj.inclusionCriteria.rats, 'UniformOutput', false), 2);
            REV = cell2struct(cell(1, length(obj.inclusionCriteria.rats)), cellfun(@(x) x, obj.inclusionCriteria.rats, 'UniformOutput', false), 2);
            perCellContributions = struct();
            
            for iRat = 1:length(obj.inclusionCriteria.rats)
                
                ratName = obj.inclusionCriteria.rats{iRat};
                
                for iSess = 1:length(obj.sessionsForBroadAnalysis.(ratName))
                    
                    session = obj.sessionsForBroadAnalysis.(ratName)(iSess);
                    filepath = fullfile(dataPath, sprintf('Rat_%s', ratName(1)), sprintf('Session%i', session));
                    
                    % Load data
                    data = obj.loadData(filepath);
                    [rewardArrivalTimes, CPEntryTimes, CPExitTimes] = obj.loadRewardCPTimes(filepath);
                    
                    % Select cells according to inclusion criteria
                    includeCellsNAc = true(1, data.nNAcUnits);
                    includeCellsCA1 = true(1, data.nCA1Units);
                    meanfrNAc = mean(data.binnedfr(1:data.nNAcUnits, :), 2);
                    meanfrCA1 = mean(data.binnedfr(data.nNAcUnits+1:end, :), 2);
                    includeCellsNAc(meanfrNAc < obj.inclusionCriteria.cells.minimumFiringRate) = false;
                    includeCellsCA1(meanfrCA1 < obj.inclusionCriteria.cells.minimumFiringRate) = false;
                    includeCells = [includeCellsNAc includeCellsCA1];
                    spiketimes = data.spiketimes(includeCells);
                    binnedfr = data.binnedfr(includeCells, :);
                    nNAcUnits = sum(includeCellsNAc);
                    nUnits = sum(includeCells);
                    
                    % Filter ripples according to analysis config
                    [rippleStart, rippleStop, ripplePeak] = obj.filterRipples(data);
                                        
                    % Calculate EV-REV
                    exclude = [];
                    EVData = struct('spiketimes', {spiketimes},...
                        'nNAcUnits', nNAcUnits,...
                        'binnedfr', binnedfr);
                    EVTimestamps = struct('startEndTimes', data.startEndTimes,...
                        'rippleStart', rippleStart,...
                        'rippleStop', rippleStop,...
                        'ripplePeak', ripplePeak,...
                        'CPEntryTimes', CPEntryTimes,...
                        'rewardArrivalTimes', rewardArrivalTimes');
                    EVAnalysisConfig = struct('exclude', exclude,...
                        'postDuration', obj.analysisConfig.postDuration,...
                        'rippleActivityMode', obj.analysisConfig.rippleActivityMode);
                    [EV0, REV0, data.binnedfrPRE, data.binnedfrPOST, data.binnedfrTASK, data.preCC, data.taskCC, data.postCC] = obj.explainedVariance(EVData, EVTimestamps, EVAnalysisConfig);
                    
                    % Store result
                    EV_REV0 = EV0-REV0;
                    EV.(ratName) = [EV.(ratName) EV0];
                    REV.(ratName) = [REV.(ratName) REV0];
                    
                    % Calculate per-cell ventral striatal contribution to EV and REV (only if this session has positive reactivation)
                    if EV_REV0 <= 0
                        warning('EV-REV negative')
                    end
                    
                    % Calculate contribution of each cell pair to the session's EV-REV
                    if EV_REV0 > 0 && ismember(session, obj.sessionsForSpecificAnalysis.(ratName))
                                includeStrCells = obj.selectRewardResponsiveCells(data.spiketimes(1:data.nNAcUnits), rewardArrivalTimes, CPEntryTimes);
                                includeCells = find(boolean(includeStrCells));
                                perCellContributions.(ratName){session} = obj.calculatePerCellPairContributions(data, includeCells, false, EV_REV0);                        
                   else
                     perCellContributions.(ratName){session} = NaN;
                                
                   end
                    
                end
                
            end
            
            obj.CA1_vStr.EV = EV;
            obj.CA1_vStr.REV = REV;
            obj.perCellContributions= perCellContributions;
            
         end
        
        
        % Calculate EV-REV for vStr  pairs
        function obj = calculateVStrEVREV(obj, varargin)
            
            % Calculates the explained variance and reverse explained variance within ventral striatum across recording sessions
            % 
            % Inputs
            %   - dataPath:            optional; specifies the path where .mat files of electrophysiological and related data are stored

            % Process paths where data is stored
            p = inputParser;
            addParameter(p, 'dataPath', fullfile('.', 'data', 'ephys'), @(x) ischar(x) || isstring(x));
            parse(p, varargin{:});
            dataPath = p.Results.dataPath;
            
            % Check there are sessions to analyse
            assert(~isempty(obj.sessionsForBroadAnalysis), 'ExplainedVarianceReactivation:sessions:noSessions', 'No sessions selected for analysis. Please revise the inclusion criteria')
            
            % Struct for storing EV and REV scores
            EV = cell2struct(cell(1, length(obj.inclusionCriteria.rats)), cellfun(@(x) x, obj.inclusionCriteria.rats, 'UniformOutput', false), 2);
            REV = cell2struct(cell(1, length(obj.inclusionCriteria.rats)), cellfun(@(x) x, obj.inclusionCriteria.rats, 'UniformOutput', false), 2);
            
            for iRat = 1:length(obj.inclusionCriteria.rats)
                
                ratName = obj.inclusionCriteria.rats{iRat};
                
                for iSess = 1:length(obj.sessionsForBroadAnalysis.(ratName))
                    
                    session = obj.sessionsForBroadAnalysis.(ratName)(iSess);
                    filepath = fullfile(dataPath, sprintf('Rat_%s', ratName(1)), sprintf('Session%i', session));
                    
                    % Load data
                    data = obj.loadData(filepath);
                    [rewardArrivalTimes, CPEntryTimes, CPExitTimes] = obj.loadRewardCPTimes(filepath);
                    
                    % Select cells according to inclusion criteria
                    includeCellsNAc = true(1, data.nNAcUnits);
                    includeCellsCA1 = false(1, data.nCA1Units);
                    meanfrNAc = mean(data.binnedfr(1:data.nNAcUnits, :), 2);
                    includeCellsNAc(meanfrNAc < obj.inclusionCriteria.cells.minimumFiringRate) = false;
                    includeCells = [includeCellsNAc includeCellsCA1];
                    spiketimes = data.spiketimes(includeCells);
                    binnedfr = data.binnedfr(includeCells, :);
                    nNAcUnits = nNAcUnits;
                                        
                    % Filter ripples according to analysis config
                    [rippleStart, rippleStop, ripplePeak] = obj.filterRipples(data);

                    % Calculate EV-REV
                    exclude = [];
                    EVData = struct('spiketimes', {spiketimes},...
                        'nNAcUnits', [],...
                        'binnedfr', binnedfr);
                    EVTimestamps = struct('startEndTimes', data.startEndTimes,...
                        'rippleStart', rippleStart,...
                        'rippleStop', rippleStop,...
                        'ripplePeak', ripplePeak,...
                        'CPEntryTimes', CPEntryTimes,...
                        'rewardArrivalTimes', rewardArrivalTimes');
                    EVAnalysisConfig = struct('exclude', exclude,...
                        'postDuration', obj.analysisConfig.postDuration,...
                        'rippleActivityMode', obj.analysisConfig.rippleActivityMode);
                    [EV0, REV0, data.binnedfrPRE, data.binnedfrPOST, data.binnedfrTASK] = obj.explainedVariance(EVData, EVTimestamps, EVAnalysisConfig);
                    
                    % Store result
                    EV.(ratName) = [EV.(ratName) EV0];
                    REV.(ratName) = [REV.(ratName) REV0];
                                                            
                end
                
            end
            
            obj.vStr.EV = EV;
            obj.vStr.REV = REV;
            
        end
        
        
        % Calculate EV-REV for CA1  pairs
        function obj = calculateCA1EVREV(obj, varargin)
            
            % Calculates the explained variance and reverse explained variance between CA1 and ventral striatum across recording sessions
            % 
            % Inputs
            %   - dataPath:            optional; specifies the path where .mat files of electrophysiological and related data are stored
            
            % Process paths where data is stored
            p = inputParser;
            addParameter(p, 'dataPath', fullfile('.', 'data', 'ephys'), @(x) ischar(x) || isstring(x));
            parse(p, varargin{:});
            dataPath = p.Results.dataPath;
            
            % Check there are sessions to analyse
            assert(~isempty(obj.sessionsForBroadAnalysis), 'ExplainedVarianceReactivation:sessions:noSessions', 'No sessions selected for analysis. Please revise the inclusion criteria')
            
             % Struct for storing EV and REV scores
            EV = cell2struct(cell(1, length(obj.inclusionCriteria.rats)), cellfun(@(x) x, obj.inclusionCriteria.rats, 'UniformOutput', false), 2);
            REV = cell2struct(cell(1, length(obj.inclusionCriteria.rats)), cellfun(@(x) x, obj.inclusionCriteria.rats, 'UniformOutput', false), 2);
            
            for iRat = 1:length(obj.inclusionCriteria.rats)
                
                ratName = obj.inclusionCriteria.rats{iRat};
                
                for iSess = 1:length(obj.sessionsForBroadAnalysis.(ratName))
                    
                    session = obj.sessionsForBroadAnalysis.(ratName)(iSess);
                    filepath = fullfile(dataPath, sprintf('Rat_%s', ratName(1)), sprintf('Session%i', session));
                    
                    % Load data
                    data = obj.loadData(filepath);
                    [rewardArrivalTimes, CPEntryTimes, CPExitTimes] = obj.loadRewardCPTimes(filepath);
                    
                    % Select cells according to inclusion criteria
                    includeCellsNAc = false(1, data.nNAcUnits);
                    includeCellsCA1 = true(1, data.nCA1Units);
                    meanfrCA1 = mean(data.binnedfr(data.nNAcUnits+1:end, :), 2);
                    includeCellsCA1(meanfrCA1 < obj.inclusionCriteria.cells.minimumFiringRate) = false;
                    includeCells = [includeCellsNAc includeCellsCA1];
                    spiketimes = data.spiketimes(includeCells);
                    binnedfr = data.binnedfr(includeCells, :);
                    nNAcUnits = sum(includeCellsNAc);
                    nUnits = sum(includeCells);
                    
                    % Filter ripples according to analysis config
                    [rippleStart, rippleStop, ripplePeak] = obj.filterRipples(data);
                    
                    % Calculate EV-REV
                    exclude = [];
                    EVData = struct('spiketimes', {spiketimes},...
                        'nNAcUnits', [],...
                        'binnedfr', binnedfr);
                    EVTimestamps = struct('startEndTimes', data.startEndTimes,...
                        'rippleStart', rippleStart,...
                        'rippleStop', rippleStop,...
                        'ripplePeak', ripplePeak,...
                        'CPEntryTimes', CPEntryTimes,...
                        'rewardArrivalTimes', rewardArrivalTimes');
                    EVAnalysisConfig = struct('exclude', exclude,...
                        'postDuration', obj.analysisConfig.postDuration,...
                        'rippleActivityMode', obj.analysisConfig.rippleActivityMode);
                    [EV0, REV0, data.binnedfrPRE, data.binnedfrPOST, data.binnedfrTASK] = obj.explainedVariance(EVData, EVTimestamps, EVAnalysisConfig);
                    
                    % Store result
                    EV.(ratName) = [EV.(ratName) EV0];
                    REV.(ratName) = [REV.(ratName) REV0];
                    
                end
                
            end
            
            obj.CA1.EV = EV;
            obj.CA1.REV = REV;
            
        end
                
        % Plot EV-REV
        function obj = plotEVREV(obj)
            
            % Plots three pairs of bars for EV and REV within and between brain areas, summarising EV and REV across rats and sessions.
            % Also plots statistical significance or lack of significance of each pair
            
            % Check there is data to plot
            areas = {'CA1', 'vStr', 'CA1_vStr'};
            errorMessages = {...
                'No data to plot for CA1 reactivation. Please run the calculateCA1EVREV() method to generate data',...
                'No data to plot for vStr reactivation. Please run the calculateVStrEVREV() method to generate data',...
                'No data to plot for CA1-vStr reactivation. Please run the calculateInterAreaEVREV() method to generate data'};
            for iArea = 1:numel(areas)
                area = areas{iArea};
                assert(isstruct(obj.(area)), 'ExplainedVarianceReactivation:dataNotAssigned', errorMessages{iArea})
                assert(~isempty(obj.(area)), 'ExplainedVarianceReactivation:dataNotAssigned', errorMessages{iArea})
                for explainedVar = {'EV', 'REV'}
                    assert(~isempty(obj.(area).(explainedVar{1})), 'ExplainedVarianceReactivation:dataNotAssigned', errorMessages{iArea})
                end
            end
            
            colour_scheme;
            axes_properties();
            
            % Plot boxes
            h = figure('Position', [680 558 454 420]); hold on
            x = 1;
            for area = {'CA1', 'vStr', 'CA1_vStr'}
                for explainedVar = {'EV', 'REV'}
                    data = cell2mat(struct2cell(obj.(area{1}).(explainedVar{1})')');
                    obj.plotBar(data, explainedVar{1}, x, h);
                    x = x + 1;
                end
            end
            
            % Re-colour
            ax = gca();
            ax.Children(end).Color = colourscheme.CA1; ax.Children(end-1).FaceColor = colourscheme.CA1;     % Box and whiskers for CA1 EV: green
            ax.Children(end-3).Color = colourscheme.CA1; ax.Children(end-4).EdgeColor = colourscheme.CA1; ax.Children(end-4).FaceColor = 'w'; ax.Children(end-5).Color = colourscheme.CA1;     % Box for CA1 REV: white with green outline; whiskers green
            ax.Children(end-6).Color = colourscheme.vStr; ax.Children(end-7).FaceColor = colourscheme.vStr;  % Box and whiskers for vStr EV: purple
             ax.Children(end-9).Color = colourscheme.vStr; ax.Children(end-10).EdgeColor = colourscheme.vStr; ax.Children(end-10).FaceColor = 'w'; ax.Children(end-11).Color = colourscheme.vStr;  % Box for vStr REV: white with purple outline; whiskers purple
            
            % T-tests
            ylim([0 1.1])
            [h, ~] = obj.plotTTestSignificance(cell2mat(struct2cell(obj.CA1.EV)'), cell2mat(struct2cell(obj.CA1.REV)'), [1 2], h);
            [h, ~] = obj.plotTTestSignificance(cell2mat(struct2cell(obj.vStr.EV)'), cell2mat(struct2cell(obj.vStr.REV)'), [3 4], h);
            [h, ~] = obj.plotTTestSignificance(cell2mat(struct2cell(obj.CA1_vStr.EV)'), cell2mat(struct2cell(obj.CA1_vStr.REV)'), [5 6], h);
            
            % Format axes
            xlim([0 7])
            set(gca, 'YTick', [0 1])
            set(gca, 'XTick', [1.5 3.5 5.5])
            set(gca, 'XTickLabels', {sprintf('CA1-CA1\\newline    pairs'), sprintf('vStr-vStr\\newline    pairs'), sprintf('CA1-vStr\\newline    pairs')})
            ylabel('(R)EV')

        end
        
        % Use drop-one-cell-pair-out to identify cell pairs that are
        % significantly reactivated and controls with which to compare them
        function obj = getSignificantlyReactivatedCellPairs(obj)
            
            % Uses the previously calculated contribution of each cell pair to the session's EV-REV reactivation metric to identify cell
            % pairs that are significantly coreactivated during POST, and matching control cell pairs that are not coreactivated.
            % Returns a list of cell pairs that belong to the two groups
                                    
            % Find upper and lower quartiles of contributions to EV-REV, per session
            if ~isstruct(obj.perCellContributions)
                error('Run EV-REV analysis first to get per-cell or per-cell-pair contributions to EV-REV')
            end
            fields = fieldnames(obj.perCellContributions);
            rats = fields(ismember(fields, obj.inclusionCriteria.rats));
            lowContributingCellPairs = cell2struct(cell(1, length(rats)), cellfun(@(x) x, rats, 'UniformOutput', false), 2);
            highContributingCellPairs = cell2struct(cell(1, length(rats)), cellfun(@(x) x, rats, 'UniformOutput', false), 2);
            
            % Quantile method, per session
            for iRat = 1:length(rats)
                ratName = rats{iRat};
                for session = 1:length(obj.perCellContributions.(ratName))
                    if isstruct(obj.perCellContributions.(ratName){session})
                        
                        % Get a cell array of all the correlation coefficients
                        corrs = mat2cell([obj.perCellContributions.(ratName){session}.correlations], [3], ones(size(obj.perCellContributions.(ratName){session}, 2), 1));
                        % Get a vector of the differences
                        corrIncrease = cellfun(@(x) x(3) - x(1), corrs);
                        
                        allCellPairs = reshape([obj.perCellContributions.(ratName){session}.cells], 2, [])';
                        allCellPairContributions = [obj.perCellContributions.(ratName){session}.contribution];
 
                        q = obj.analysisConfig.reactivationQuantileThreshold;
                        switch obj.analysisConfig.controlCriterion

                            case 'negativeContributions'
                                lq = quantile(allCellPairContributions, q);
                                lowContributingCellPairs.(ratName){session} = allCellPairs(allCellPairContributions < lq & corrIncrease > 0, :);
                                uq = quantile(allCellPairContributions, 1-q);
                                highContributingCellPairs.(ratName){session} = allCellPairs(find(allCellPairContributions > uq & corrIncrease > 0), :);

                            case 'minimalContributions'
                                lq = quantile(abs(allCellPairContributions), q);
                                lowContributingCellPairs.(ratName){session} = allCellPairs(abs(allCellPairContributions) < lq, :);
                                uq = quantile(allCellPairContributions, 1-q);
                                highContributingCellPairs.(ratName){session} = allCellPairs(find(allCellPairContributions > uq & corrIncrease > 0), :);
                        end

                        
                    end
                end
            end
            
            obj.perCellContributions.controlCellPairs = lowContributingCellPairs;
            obj.perCellContributions.reactivatedCellPairs = highContributingCellPairs;
                                
        end

                                                
        function obj = getReactivatedCellPairActivity(obj, event, varargin)
            
            % Calculates the trial-averaged, event-triggered coactivity of pairs of neurons across sessions, divided into trials of
            % different types.
            % 
            % Inputs
            %   - event:                        'centralPlatform' or 'rewardArrival'; trial event around which to calculate the coactivity
            %   - behaviouralDataPath:   optional; specifies the path where .mat files of behavioural data are stored
            %   - ephysDataPath:          optional; specifies the path where .mat files of electrophysiological and related data are stored
            
            
            % Check the event input
            assert(nargin >= 2, 'Not enough input arguments.');
            acceptableEvents = {'centralPlatform', 'rewardArrival'};
            assert(ismember(event, acceptableEvents), sprintf('Value of input argument event "%s" is not valid. Must be "centralPlatform" or "rewardArrival".', event));
            
            % Check that there are reactivated cell pairs
            assert(~isempty(obj.perCellContributions), 'Reactivted cell pairs have not been identified. Please run getSignificantlyReactivatedCellPairs() first.')
            
            % Process paths where data is stored
            p = inputParser;
            addParameter(p, 'behaviouralDataPath', fullfile('.', 'data', 'behavioural_data'), @(x) ischar(x) || isstring(x));
            addParameter(p, 'ephysDataPath', fullfile('.', 'data', 'ephys_data'), @(x) ischar(x) || isstring(x));
            parse(p, varargin{:});
            behaviouralDataPath = p.Results.behaviouralDataPath;
            ephysDataPath = p.Results.ephysDataPath;
            
            fields = fieldnames(obj.perCellContributions.reactivatedCellPairs);
            rats = fields(ismember(fields, obj.inclusionCriteria.rats));
            
            % Create a structure for storing the results
            taskCoactivity = struct(...
                'highRewardExpectation', struct(...
                    'reactivated', cell2struct(cell(1, length(rats)), cellfun(@(x) x, rats, 'UniformOutput', false), 2),...
                    'control', cell2struct(cell(1, length(rats)), cellfun(@(x) x, rats, 'UniformOutput', false), 2)),...
                'mediumRewardExpectation', struct(...
                    'reactivated', cell2struct(cell(1, length(rats)), cellfun(@(x) x, rats, 'UniformOutput', false), 2),...
                    'control', cell2struct(cell(1, length(rats)), cellfun(@(x) x, rats, 'UniformOutput', false), 2)),...
                'perCellHighRewardExpectation', struct(...
                    'reactivated', cell2struct(cell(1, length(rats)), cellfun(@(x) x, rats, 'UniformOutput', false), 2),...
                    'control', cell2struct(cell(1, length(rats)), cellfun(@(x) x, rats, 'UniformOutput', false), 2)),...
                'perCellMediumRewardExpectation', struct(...
                    'reactivated', cell2struct(cell(1, length(rats)), cellfun(@(x) x, rats, 'UniformOutput', false), 2),...
                    'control', cell2struct(cell(1, length(rats)), cellfun(@(x) x, rats, 'UniformOutput', false), 2)));
            
            for iRat = 1:length(rats)
                ratName = rats{iRat};
                
                % Load behavioural data
                behav_data = importdata(fullfile(behaviouralDataPath, sprintf('behaviour_%s.mat', lower(ratName))));
                
                mediumExpTaskCoactivityHigh = [];
                mediumExpTaskCoactivityLow = [];
                highExpTaskCoactivityHigh = [];
                highExpTaskCoactivityLow = [];
                
                for session = 1:length(obj.perCellContributions.reactivatedCellPairs.(ratName))
                    
                    if isempty(obj.perCellContributions.reactivatedCellPairs.(ratName){session})
                        continue
                    end
                    
                    filepath = fullfile(ephysDataPath, sprintf('Rat_%s', ratName(1)), sprintf('Session%i', session));
                    data = obj.loadData(filepath);
                    
                    reactivatedCellPairs = obj.perCellContributions.reactivatedCellPairs.(ratName){session};
                    controlCellPairs = obj.perCellContributions.controlCellPairs.(ratName){session};
                    reactivatedCellPairs = mat2cell(reactivatedCellPairs, ones(size(reactivatedCellPairs, 1), 1), 2);
                    controlCellPairs = mat2cell(controlCellPairs, ones(size(controlCellPairs, 1), 1), 2);
                    reactivatedCellPairs = cellfun(@(x) [x(1) x(2)+data.nNAcUnits], reactivatedCellPairs, 'UniformOutput', false);
                    controlCellPairs = cellfun(@(x) [x(1) x(2)+data.nNAcUnits], controlCellPairs, 'UniformOutput', false);
                    
                    % Load behavioural info
                    armChoice = behav_data.reward_probs{session}(behav_data.actions{session});
                    armLegitimacy = cellfun(@(x, y) x(y), behav_data.arm_values{session}, mat2cell(behav_data.actions{session}, 1, ones(behav_data.n_trials(session), 1)));
                    armLegitimacy = cellfun(@(x) ~strcmp(x, 'Illegitimate'), armLegitimacy);
                    
                    % Load spiking data
                    [rewardArrivalTimes, CPEntryTimes, ~] = obj.loadRewardCPTimes(filepath);
                    includeStrCells = find(obj.selectRewardResponsiveCells(data.spiketimes(1:data.nNAcUnits), rewardArrivalTimes, CPEntryTimes));
                    includeCells = [includeStrCells, ones(1, data.nCA1Units)];
                    spiketimes = data.spiketimes(boolean(includeCells));
                    binnedfr = data.binnedfr(boolean(includeCells), :);
                    
                    % Get activity of all members of cell pairs during TASK
                    spiketimes1Reactivated = spiketimes(cellfun(@(x) x(1), reactivatedCellPairs));
                    spiketimes2Reactivated = spiketimes(cellfun(@(x) x(2), reactivatedCellPairs));
                    spiketimes1Control = spiketimes(cellfun(@(x) x(1), controlCellPairs));
                    spiketimes2Control = spiketimes(cellfun(@(x) x(2), controlCellPairs));
         
                    binnedfr1Reactivated = binnedfr(cellfun(@(x) x(1), reactivatedCellPairs), :);
                    binnedfr2Reactivated = binnedfr(cellfun(@(x) x(2), reactivatedCellPairs), :);
                    binnedfr1Control = binnedfr(cellfun(@(x) x(1), controlCellPairs), :);
                    binnedfr2Control = binnedfr(cellfun(@(x) x(2), controlCellPairs), :);
                                
                    % Divide trials into mid- and high- probability arms
                    mediumExpTrials = strcmp(armChoice, 'medium') & armLegitimacy;
                    highExpTrials = strcmp(armChoice, 'high') & armLegitimacy;
                    
                    % Divide trials into mid- and high- probability arms
                    switch event
                        case 'centralPlatform'
                            armChoice = armChoice(2:end);
                            armLegitimacy = armLegitimacy(2:end);
                    end
                    mediumExpTrials = strcmp(armChoice, 'medium')  & armLegitimacy;
                    highExpTrials = strcmp(armChoice, 'high') & armLegitimacy;

                    % Get coactivity rates -5 to +5 seconds around arrival, divided into
                    % mid- and high-probability arms, for each rat
                    mediumExpTrials = boolean(mediumExpTrials);
                    highExpTrials = boolean(highExpTrials);
                    coactivity = struct();
                    switch event
                        case 'rewardArrival'
                            timestamps = rewardArrivalTimes;
                        case 'centralPlatform'
                            timestamps = CPEntryTimes(2:end);
                    end
                    coactivity.mediumExpTrials.reactivated = obj.getEventTriggeredCoactivity(timestamps(mediumExpTrials)', spiketimes1Reactivated, spiketimes2Reactivated, 5, 5, obj.analysisConfig.gaussianWindowWidth);
                    coactivity.mediumExpTrials.control = obj.getEventTriggeredCoactivity(timestamps(mediumExpTrials)', spiketimes1Control, spiketimes2Control, 5, 5, obj.analysisConfig.gaussianWindowWidth);
                    coactivity.highExpTrials.reactivated = obj.getEventTriggeredCoactivity(timestamps(highExpTrials)', spiketimes1Reactivated, spiketimes2Reactivated, 5, 5, obj.analysisConfig.gaussianWindowWidth);
                    coactivity.highExpTrials.control = obj.getEventTriggeredCoactivity(timestamps(highExpTrials)', spiketimes1Control, spiketimes2Control, 5, 5, obj.analysisConfig.gaussianWindowWidth);
                    
                    % Z-score the coactivity
                    coactivity.mediumExpTrials.reactivated = obj.zScore(coactivity.mediumExpTrials.reactivated, min(binnedfr1Reactivated, binnedfr2Reactivated));
                    coactivity.mediumExpTrials.control = obj.zScore(coactivity.mediumExpTrials.control, min(binnedfr1Control, binnedfr2Control));
                    coactivity.highExpTrials.reactivated = obj.zScore(coactivity.highExpTrials.reactivated, min(binnedfr1Reactivated, binnedfr2Reactivated));
                    coactivity.highExpTrials.control = obj.zScore(coactivity.highExpTrials.control, min(binnedfr1Control, binnedfr2Control));
                                                            
                    % Stack per-cell average firing rates together
                    if (sum(mediumExpTrials) > 1) && (sum(highExpTrials) > 1)
                        for condition = {'reactivated', 'control'}
                            condition = condition{1};
                            taskCoactivity.perCellHighRewardExpectation.(condition).(ratName) = [taskCoactivity.perCellHighRewardExpectation.(condition).(ratName);...
                                permute(mean(coactivity.highExpTrials.(condition), 1), [2, 3, 1])];
                            taskCoactivity.perCellMediumRewardExpectation.(condition).(ratName) = [taskCoactivity.perCellMediumRewardExpectation.(condition).(ratName);...
                                permute(mean(coactivity.mediumExpTrials.(condition), 1), [2, 3, 1])];
                        end
                    end
                    
                    % Stack the firing rate arrays together
                    for condition = {'reactivated', 'control'}
                        condition = condition{1};
                        n_ts = size(coactivity.mediumExpTrials.(condition), 3);
                        coactivity.mediumExpTrials.(condition) = reshape(coactivity.mediumExpTrials.(condition), [], n_ts);
                        coactivity.highExpTrials.(condition) = reshape(coactivity.highExpTrials.(condition), [], n_ts);   
                        taskCoactivity.highRewardExpectation.(condition).(ratName) = [taskCoactivity.highRewardExpectation.(condition).(ratName);...
                            coactivity.highExpTrials.(condition)];
                        taskCoactivity.mediumRewardExpectation.(condition).(ratName) = [taskCoactivity.mediumRewardExpectation.(condition).(ratName);...
                            coactivity.mediumExpTrials.(condition)];
                    end
                    
                end
                
            end
            
            % Add the time bin centres relative to events as timestamps for further analysis
            taskCoactivity.ts = (-5 + obj.BinSize/2) : obj.BinSize : (5 - obj.BinSize/2);
            if ~isempty(taskCoactivity.highRewardExpectation.reactivated.(ratName))
                assert(size(taskCoactivity.highRewardExpectation.reactivated.(ratName), 2) == size(taskCoactivity.ts, 2), 'Dimensions of cell pair coactivity and timestamps are not consistent. Please revise the window length and bin size.')
            end
            
            obj.taskCoactivity = taskCoactivity;
                        
        end

                
        function obj = plotReactivatedCellPairActivity(obj, event, varargin)
            
            % Plots line graphs of the mean and standard error of event-triggered coactivity of 
            % two types of cell pairs (reactivated and control) on two types of trial
            % (high-reward-expectation and medium-reward-expectation)
            % 
            % Inputs
            %   - event:            'rewardArrival' or 'centralPlatform', used for labelling the x-axis
            %   - perRat:           optional; logical which plots a separate figure for each rat if true, or pools cell pairs
            %                          from all rats together if false
            
            % Process input
            p = inputParser;
            addParameter(p, 'perRat', false, @(x) islogical(x));
            parse(p, varargin{:});
            perRat = p.Results.perRat;
            
            if ~isempty(obj.taskCoactivity)
                activity = obj.taskCoactivity;
                ts = obj.taskCoactivity.ts;
            else
                error('No activity has been calculated, so cannot be plotted. Please run getReactivatedCellPairActivity().')
            end
                        
            % Plot mean and SEM
            
            fields = fieldnames(obj.perCellContributions.reactivatedCellPairs);
            rats = fields(ismember(fields, obj.inclusionCriteria.rats));
            figures = {};
            
            if perRat
                
                for iRat = 1:length(rats)
                    ratName = rats{iRat};

                    % Skip if no control / reactivated firing computed
                    if isempty(activity.perCellHighRewardExpectation.reactivated.(ratName))...
                        || isempty(activity.perCellHighRewardExpectation.control.(ratName))
                        continue
                    end

                    mediumControl = activity.perCellMediumRewardExpectation.control.(ratName);
                    mediumReactivated = activity.perCellMediumRewardExpectation.reactivated.(ratName);
                    highControl = activity.perCellHighRewardExpectation.control.(ratName);
                    highReactivated = activity.perCellHighRewardExpectation.reactivated.(ratName);
                    fig = obj.plotActivity(mediumControl, mediumReactivated, highControl, highReactivated, ts);
                    suptitle(ratName)
                    figures{iRat} = fig;

                end

                obj.figures = figures;
                
            else
                mediumControl = cell2mat(struct2cell(activity.perCellMediumRewardExpectation.control));
                mediumReactivated = cell2mat(struct2cell(activity.perCellMediumRewardExpectation.reactivated));
                highControl = cell2mat(struct2cell(activity.perCellHighRewardExpectation.control));
                highReactivated = cell2mat(struct2cell(activity.perCellHighRewardExpectation.reactivated));
                fig = obj.plotActivity(mediumControl, mediumReactivated, highControl, highReactivated, ts, event);
                fig.Units = 'centimeters';
                fig.Position = [14.8, 16.6, 17.9, 8.7];
                figures{1} = fig;
                
            end

        end
        
        % Plot 
        function fig = plotActivity(~, mediumControl, mediumReactivated, highControl, highReactivated, ts, event)
            
            % Plots line graphs of the mean and standard error of event-triggered coactivity of 
            % two types of cell pairs (reactivated and control) on two types of trial
            % (high-reward-expectation and medium-reward-expectation)
            % 
            % Inputs
            %   - mediumControl:        vector of binned coactivity values of control cell pairs on medium-reward-expectation trials
            %   - mediumReactivated: vector of binned coactivity values of reactivated cell pairs on medium-reward-expectation trials
            %   - highControl:             vector of binned coactivity values of control cell pairs on high-reward-expectation trials
            %   - highReactivated:      vector of binned coactivity values of reactivated cell pairs on high-reward-expectation trials
            %   - ts:                          vector of timestamps used for labelling the x-axis; must be the same length as the medium and high control and
            %                                   reactivated vectors
            %   - event:                     'centralPlatform' or 'rewardArrival'; used for labelling the x-axis
            % 
            % Outputs
            %   - fig:                         handle for the resulting figure
            
            colour_scheme;
            axes_properties;

            fig = figure;
            subplot(1, 2, 1); hold on

            % Plot control cell-pair activity on medium-reward-expectation trials
            FaceColor = colourscheme.insignificant_secondary*1.5;
            LineColor = colourscheme.insignificant_primary;
            firing = mediumControl;
            if size(firing, 2)==1
                firing = firing';
            end
            firing_mean = nanmean(firing, 1);
            firing_std = nanstd(firing, [], 1);
            n_obs = size(firing, 1);
            firing_sem = firing_std / sqrt(n_obs);
            patch([ts flip(ts)], [firing_mean + firing_sem flip(firing_mean - firing_sem)], 'k', 'FaceColor', FaceColor', 'EdgeColor', 'none')
            plot(ts, firing_mean, 'Color', LineColor)

            % Plot reactivated cell-pair activity on medium-reward-expectation trials
            FaceColor = colourscheme.medium_exp_reward_primary;
            LineColor = colourscheme.medium_exp_reward_secondary;
            firing = mediumReactivated;
            if size(firing, 2)==1
                firing = firing';
            end
            firing_mean = nanmean(firing, 1);
            firing_std = nanstd(firing, [], 1);
            n_obs = size(firing, 1);
            firing_sem = firing_std / sqrt(n_obs);
            patch([ts flip(ts)], [firing_mean + firing_sem flip(firing_mean - firing_sem)], 'k', 'FaceColor', FaceColor', 'EdgeColor', 'none', 'FaceAlpha', 0.75)
            plot(ts, firing_mean, 'Color', LineColor)

            title('Medium reward expectation')
            ylabel('Coactivity (z-score)')

            subplot(1, 2, 2); hold on

            % Plot control cell-pair activity on high-reward-expectation trials
            FaceColor = colourscheme.insignificant_primary*1.5;
            LineColor = colourscheme.insignificant_secondary;
            firing = highControl;
            if size(firing, 2)==1
                firing = firing';
            end
            firing_mean = nanmean(firing, 1);
            firing_std = nanstd(firing, [], 1);
            n_obs = size(firing, 1);
            firing_sem = firing_std / sqrt(n_obs);
            patch([ts flip(ts)], [firing_mean + firing_sem flip(firing_mean - firing_sem)], 'k', 'FaceColor', FaceColor', 'EdgeColor', 'none')
            plot(ts, firing_mean, 'Color', LineColor)

            % Plot reactivated cell-pair activity on high-reward-expectation trials
            FaceColor = colourscheme.high_exp_reward_primary;
            LineColor = colourscheme.high_exp_reward_secondary;
            firing = highReactivated;
            if size(firing, 2)==1
                firing = firing';
            end
            firing_mean = nanmean(firing, 1);
            firing_std = nanstd(firing, [], 1);
            n_obs = size(firing, 1);
            firing_sem = firing_std / sqrt(n_obs);
            patch([ts flip(ts)], [firing_mean + firing_sem flip(firing_mean - firing_sem)], 'k', 'FaceColor', FaceColor', 'EdgeColor', 'none', 'FaceAlpha', 0.75)
            plot(ts, firing_mean, 'Color', LineColor)

            title('High reward expectation')

            % Equalise y-limits
            subplot(1, 2, 1); y1 = ylim;
            subplot(1, 2, 2); y2 = ylim;
            for i = 1:2
                subplot(1, 2, i)
                ylim([min([y1(1) y2(1)]) max([y1(2) y2(2)])])
                xlim([-5 5])
                plot([0 0], [min([y1(1) y2(1)]) max([y1(2) y2(2)])], 'k', 'LineStyle', '--')
                switch event
                    case 'rewardArrival'
                        xlabel({'Time from arrival', 'at reward location (s)'})
                    case 'centralPlatform'
                        xlabel({'Time from arrival', 'at central platform (s)'})
                end
                ax = gca(); pos = ax.Position; pos(2) = pos(2) + 0.15; pos(4) = pos(4) - 0.2; ax.Position = pos;
            end

        end

                
            
        % Get event-triggered activity
        function obj = getEventTriggeredActivity(obj)

            % Calculates the mean coactivity in the 2 seconds prior to event, averaged by rat, reward expectation, and
            % reactivated/control cell type
            
            if ~isempty(obj.taskCoactivity)
                ts = obj.taskCoactivity.ts;
            else
                error('No activity has been calculated, so cannot be averaged. Please run getReactivatedCellPairActivity() first.')
            end
            meanEventTriggeredActivity = struct('high', struct('reactivated', struct(), 'control', struct()),...
                'medium', struct('reactivated', struct(), 'control', struct()));

            for ratName = obj.inclusionCriteria.rats
                for trialType = {'High', 'Medium'}
                    f = cell2mat(['perCell' trialType 'RewardExpectation']);
                    for cellType = {'reactivated', 'control'}
                        if numel(obj.taskCoactivity.(f).(cellType{1}).(ratName{1})) > 0
                            firing = obj.taskCoactivity.(f).(cellType{1}).(ratName{1})(:, (ts > -2) & (ts < 0));
                            meanEventTriggeredActivity.(lower(trialType{1})).(cellType{1}).(ratName{1}) = mean(firing, 2);
                        else
                            meanEventTriggeredActivity.(lower(trialType{1})).(cellType{1}).(ratName{1}) = [];
                        end
                    end
                end
            end

            obj.meanEventTriggeredActivity = meanEventTriggeredActivity;

        end
        
        
        function [table, p_interaction, p_medium, p_high] = mixedEffectsNestedAnova(obj)
            
            % Performs a mixed-effects nested ANOVA with cell-pair type (reactivated or control) and trial type (high- or
            % medium-reward-expectation) as fixed effects, cell-pair identity nested within rat identity as random effects, and
            % mean event-triggered coactivity as the dependent variable.
            % Also performs post-hoc tests of main effects.
            % 
            % Outputs
            %   - table:               full ANOVA table with all results
            %   - p_interaction:    p-value associated with the interaction effect of cell-pair type and trial type
            %   - p_medium:       p-value associated with post-hoc t-test on medium-reward-expectation trials
            %   - p_high:       p-value associated with post-hoc t-test on high-reward-expectation trials
            
            % Concatenate data from all conditions
            reactivatedHigh = cell2mat(struct2cell(obj.meanEventTriggeredActivity.high.reactivated))';
            reactivatedMid = cell2mat(struct2cell(obj.meanEventTriggeredActivity.medium.reactivated))';
            controlHigh = cell2mat(struct2cell(obj.meanEventTriggeredActivity.high.control))';
            controlMid = cell2mat(struct2cell(obj.meanEventTriggeredActivity.medium.control))';
            Y =  [reactivatedHigh reactivatedMid controlHigh controlMid];

            % Specify reactivated or not
            react = [ones(1, length(reactivatedHigh)) * 1, ...
                        ones(1, length(reactivatedMid)) * 1, ...
                        ones(1, length(controlHigh)) * 2, ...
                        ones(1, length(controlMid)) * 2];

            % Specify high or medium expected reward
            exp_reward = [ones(1, length(reactivatedHigh)) * 1, ...
                                ones(1, length(reactivatedMid)) * 2, ...
                                ones(1, length(controlHigh)) * 1, ...
                                ones(1, length(controlMid)) * 2];

            % Specify unit number
            nReactivated = length(reactivatedHigh);
            nControl = length(controlHigh);
            nTotal = nReactivated + nControl;
            neuron = [1:nReactivated, ...
                          1:nReactivated, ...
                          nReactivated + 1 : nTotal, ...
                          nReactivated + 1 : nTotal];

            % Specify rat number
            rat = [];
            for cellType = {'reactivated', 'control'}
                for trialType = {'high', 'medium'} 
                    for iRat = 1:length(obj.inclusionCriteria.rats)
                        ratName = obj.inclusionCriteria.rats{iRat};
                        nCells = length(obj.meanEventTriggeredActivity.(trialType{1}).(cellType{1}).(ratName));
                        rat = [rat ones(1, nCells) * iRat];
                    end
                end
            end

            nesting = [0, 0, 0, 0;
                           0, 0, 0, 0;
                           1, 0, 0, 1;
                           0, 0, 0, 0];

            % Compute ANOVA
            [p, table, stats] = anovan(Y, {react, exp_reward, neuron, rat}, ...
                                                'model', 'interaction', ...
                                                'random', [3, 4], ...
                                                'nested', nesting, ...
                                                'varnames', {'reactivated', 'expected_reward', 'neuron', 'rat'});
            p_interaction = table(find(strcmp(table(:, 1), 'reactivated*expected_reward')), find(strcmp(table(1, :), 'Prob>F')));
            p_interaction = p_interaction{1};
                                            
            % Compute post-hoc tests
            highReactivated = Y(:, exp_reward==1 & react==1);
            highControl = Y(:, exp_reward==1 & react==2);
            mediumReactivated = Y(:, exp_reward==2 & react==1);
            mediumControl = Y(:, exp_reward==2 & react==2);
            [~, p_high] = ttest2(highReactivated, highControl);
            [~, p_medium] = ttest2(mediumReactivated, mediumControl);
                                            
        end
        
        
        function obj = plotCellPairCoactivityExamples(~, ratName, session, data, varargin)
            
            % Plots three illustrative examples of coactivity of high-contributing cell pairs.
            % 
            % Inputs
            %   - ratName:               string, corresponding to the rats listed in the inclusion criteria for which data has been analysed
            %   - session:                integer, corresponding to a session for which data has been analysed
            %   - data:                     struct containing the spiketimes and number of ventral striatal units for the session
            
            all_contributions = [obj.perCellContributions.(ratName){session}(:).contribution];
            [~, order] = sort(all_contributions, 'descend');

            for iPair = 1:3

                % Get coactivity for both arms
                cellPair = obj.perCellContributions.(ratName){session}(order(iPair)).cells;
                spiketimes1 = data.spiketimes(cellPair(1));
                spiketimes2 = data.spiketimes(cellPair(2) + data.nNAcUnits);
                mediumCoactivity = obj.getEventTriggeredCoactivity(rewardArrivalTimes(mediumExpTrials)', spiketimes1, spiketimes2, 5, 5);
                highCoactivity = obj.getEventTriggeredCoactivity(rewardArrivalTimes(highExpTrials)', spiketimes1, spiketimes2, 5, 5);

                % Overall minimum and maximum coactivity to scale c-limits
                pairMaxFR = max([mediumCoactivity(:); highCoactivity(:)]);
                pairMinFR = min([mediumCoactivity(:); highCoactivity(:)]);

                % Plot coactivity on medium-probability arm
                figure('Position', [530, 600, 1100, 260]);
                binsize = 0.05;
                subplot(1, 2, 1);
                [nTrials, ~, nBins] = size(mediumCoactivity);
                xData = linspace(-5+binsize/2, 5-binsize/2, nBins);
                yData = find(mediumExpTrials);
                imagesc(xData, 1:nTrials, squeeze(mediumCoactivity), [pairMinFR, pairMaxFR])
                y = ylim; line([0, 0], y, 'Color', 'k');
                set(gca, 'YTick', 1:length(yData));
                set(gca, 'YTickLabel', yData)
                xlabel('Time from arrival at reward location (s)')
                ylabel('Trial #')
                title('Medium-probability arm')

                % Plot coactivity on high-probability arm
                subplot(1, 2, 2);
                [nTrials, ~, nBins] = size(highCoactivity);
                xData = linspace(-5+binsize/2, 5-binsize/2, nBins);
                yData = find(highExpTrials);
                imagesc(xData, 1:nTrials, squeeze(highCoactivity), [pairMinFR, pairMaxFR])
                y = ylim; line([0, 0], y, 'Color', 'k');
                set(gca, 'YTick', 1:length(yData));
                set(gca, 'YTickLabel', yData)
                xlabel('Time from arrival at reward location (s)')
                ylabel('Trial #')
                title('High-probability arm')

            end
            
        end
        
        
        function obj = plotCellPairRippleCoactivity(obj, varargin)
            
            % Plots ripple-triggered average coactivity of reactivated and control cell pairs.
            % 
            % Inputs
            %   - ephysDataPath:     optional; specifies the path where .mat files of electrophysiological and related data are stored

            % Process paths where data is stored, and progress tracking
            p = inputParser;
            addParameter(p, 'ephysDataPath', fullfile('.', 'data', 'ephys_data'), @(x) ischar(x) || isstring(x));
            parse(p, varargin{:});
            ephysDataPath = p.Results.ephysDataPath;
            
            fields = fieldnames(obj.perCellContributions.reactivatedCellPairs);
            rats = fields(ismember(fields, obj.inclusionCriteria.rats));

            % Create a structure for storing the results
            rippleCoactivity = struct(...
                'PRE', struct(...
                    'reactivated', cell2struct(cell(1, length(rats)), cellfun(@(x) x, rats, 'UniformOutput', false), 2),...
                    'control', cell2struct(cell(1, length(rats)), cellfun(@(x) x, rats, 'UniformOutput', false), 2)),...
                'POST', struct(...
                    'reactivated', cell2struct(cell(1, length(rats)), cellfun(@(x) x, rats, 'UniformOutput', false), 2),...
                    'control', cell2struct(cell(1, length(rats)), cellfun(@(x) x, rats, 'UniformOutput', false), 2)),...
                'perCellPairPre', struct(...
                    'reactivated', cell2struct(cell(1, length(rats)), cellfun(@(x) x, rats, 'UniformOutput', false), 2),...
                    'control', cell2struct(cell(1, length(rats)), cellfun(@(x) x, rats, 'UniformOutput', false), 2)),...
                'perCellPairPost', struct(...
                    'reactivated', cell2struct(cell(1, length(rats)), cellfun(@(x) x, rats, 'UniformOutput', false), 2),...
                    'control', cell2struct(cell(1, length(rats)), cellfun(@(x) x, rats, 'UniformOutput', false), 2)));

            for iRat = 1:length(rats)
                ratName = rats{iRat};

                for session = 1:length(obj.perCellContributions.reactivatedCellPairs.(ratName))

                    if isempty(obj.perCellContributions.reactivatedCellPairs.(ratName){session})
                        continue
                    end

                    % Load spiketimes
                    filepath = fullfile(ephysDataPath, sprintf('Rat_%s', ratName(1)), sprintf('Session%i', session));
                    data = obj.loadData(filepath);

                    % Load ripple times during PRE and POST
                    ripplePeakPre = data.ripplePeak(data.ripplePeak > data.startEndTimes(1, 1) & data.ripplePeak < data.startEndTimes(1, 2));
                    ripplePeakPost = data.ripplePeak(data.ripplePeak > data.startEndTimes(3, 1) & data.ripplePeak < data.startEndTimes(3, 2));

                    % Get spiketimes1 and spiketimes 2 for reactivated and control cell pairs
                    reactivatedCellPairs = obj.perCellContributions.reactivatedCellPairs.(ratName){session};
                    controlCellPairs = obj.perCellContributions.controlCellPairs.(ratName){session};
                    reactivatedCellPairs = mat2cell(reactivatedCellPairs, ones(size(reactivatedCellPairs, 1), 1), 2);
                    controlCellPairs = mat2cell(controlCellPairs, ones(size(controlCellPairs, 1), 1), 2);
                    reactivatedCellPairs = cellfun(@(x) [x(1) x(2)+data.nNAcUnits], reactivatedCellPairs, 'UniformOutput', false);
                    controlCellPairs = cellfun(@(x) [x(1) x(2)+data.nNAcUnits], controlCellPairs, 'UniformOutput', false);
                    spiketimes1Reactivated = data.spiketimes(cellfun(@(x) x(1), reactivatedCellPairs));
                    spiketimes2Reactivated = data.spiketimes(cellfun(@(x) x(2), reactivatedCellPairs));
                    spiketimes1Control = data.spiketimes(cellfun(@(x) x(1), controlCellPairs));
                    spiketimes2Control = data.spiketimes(cellfun(@(x) x(2), controlCellPairs));
                    binnedfr1Reactivated = data.binnedfr(cellfun(@(x) x(1), reactivatedCellPairs), :);
                    binnedfr2Reactivated = data.binnedfr(cellfun(@(x) x(2), reactivatedCellPairs), :);
                    binnedfr1Control = data.binnedfr(cellfun(@(x) x(1), controlCellPairs), :);
                    binnedfr2Control = data.binnedfr(cellfun(@(x) x(2), controlCellPairs), :);

                    % Get coactivity of spiketimes for 500ms around each of a sample of 100 ripples
                    preIndex = round(linspace(1, length(ripplePeakPre), 100));
                    if length(ripplePeakPost) < 100
                        warning(sprintf('Only %i ripples in POST', length(ripplePeakPost)))
                        postIndex = 1:length(ripplePeakPost);
                    else
                        postIndex = round(linspace(1, length(ripplePeakPost), 100));
                    end
                    rippleCoactivity.PRE.reactivated = obj.getEventTriggeredCoactivity(ripplePeakPre(preIndex)', spiketimes1Reactivated, spiketimes2Reactivated, 0.5, 0.5, 0.05);
                    rippleCoactivity.PRE.control = obj.getEventTriggeredCoactivity(ripplePeakPre(preIndex)', spiketimes1Control, spiketimes2Control, 0.5, 0.5, 0.05);
                    rippleCoactivity.POST.reactivated = obj.getEventTriggeredCoactivity(ripplePeakPost(postIndex)', spiketimes1Reactivated, spiketimes2Reactivated, 0.5, 0.5, 0.05);
                    rippleCoactivity.POST.control = obj.getEventTriggeredCoactivity(ripplePeakPost(postIndex)', spiketimes1Control, spiketimes2Control, 0.5, 0.5, 0.05);

                    % Z-score coactivity
                    rippleCoactivity.PRE.reactivated = obj.zScore(rippleCoactivity.PRE.reactivated, min(binnedfr1Reactivated, binnedfr2Reactivated));
                    rippleCoactivity.PRE.control = obj.zScore(rippleCoactivity.PRE.control, min(binnedfr1Control, binnedfr2Control));
                    rippleCoactivity.POST.reactivated = obj.zScore(rippleCoactivity.POST.reactivated, min(binnedfr1Reactivated, binnedfr2Reactivated));
                    rippleCoactivity.POST.control = obj.zScore(rippleCoactivity.POST.control, min(binnedfr1Control, binnedfr2Control));

                    % Store
                    for condition = {'reactivated', 'control'}
                        condition = condition{1};
                        rippleCoactivity.perCellPairPre.(condition).(ratName) = [rippleCoactivity.perCellPairPre.(condition).(ratName);...
                            reshape(rippleCoactivity.PRE.(condition), [], 20)];
                        rippleCoactivity.perCellPairPost.(condition).(ratName) = [rippleCoactivity.perCellPairPost.(condition).(ratName);...
                            reshape(rippleCoactivity.POST.(condition), [], 20)];
                    end

                end

            end

            % Plot
            controlPre = cell2mat(struct2cell(rippleCoactivity.perCellPairPre.control));
            controlPost = cell2mat(struct2cell(rippleCoactivity.perCellPairPost.control));
            reactivatedPre = cell2mat(struct2cell(rippleCoactivity.perCellPairPre.reactivated));
            reactivatedPost = cell2mat(struct2cell(rippleCoactivity.perCellPairPost.reactivated));
            ts = -475:50:475;
            fig = obj.plotActivity(reactivatedPre, reactivatedPost, reactivatedPre, reactivatedPost, ts);
            
        end
        
        function obj = plotRunningSpeedAndPopulationActivity(obj, event, scope, varargin)
            
            % Plots a figure of three subplots showing trial-averaged running speed, CA1 population firing rate,
            % and ventral striatum population firing rate, respectively
            % 
            % Inputs
            %   - event:                        'centralPlatform' or 'rewardArrival'; trial event around which to calculate the coactivity
            %   - scope:                        'all', or 'rewardExpectancy'; if 'all' then all trials are pooled together; if 'rewardExpectancy' 
            %                                      then activity on high- and medium-reward-expectation trials are plotted separately
            %   - behaviouralDataPath:   optional; specifies the path where .mat files of behavioural data are stored
            %   - ephysDataPath:          optional; specifies the path where .mat files of electrophysiological and related data are stored
            %   - plotProgress:              optional; plot a progress bar (this function takes a long time to run)
            
             % Check the event input
            assert(nargin >= 2, 'Not enough input arguments.');
            acceptableEvents = {'centralPlatform', 'rewardArrival'};
            assert(ismember(event, acceptableEvents), sprintf('Value of input argument event "%s" is not valid. Must be "centralPlatform" or "rewardArrival".', event));
            
            % Check the scope input
            assert(nargin >= 3, 'Not enough input arguments.');
            acceptableScopes = {'rewardExpectancy', 'all'};
            assert(ismember(scope, acceptableScopes), sprintf('Value of input argument scope "%s" is not valid. Must be "rewardExpectancy" or "all".', scope));
                        
            % Process paths where data is stored, and progress tracking
            p = inputParser;
            addParameter(p, 'behaviouralDataPath', fullfile('.', 'data', 'behavioural_data'), @(x) ischar(x) || isstring(x));
            addParameter(p, 'ephysDataPath', fullfile('.', 'data', 'ephys_data'), @(x) ischar(x) || isstring(x));
            addParameter(p, 'plotProgress', false, @(x) islogical(x));
            parse(p, varargin{:});
            behaviouralDataPath = p.Results.behaviouralDataPath;
            ephysDataPath = p.Results.ephysDataPath;
            plotProgress = p.Results.plotProgress;
            
            % Create a struct for storing the data
            switch scope
                case 'rewardExpectancy'
                    runningSpeed = struct(...
                        'medium', cell2struct(repmat({{}}, 1, length(obj.inclusionCriteria.rats)), obj.inclusionCriteria.rats, 2),...
                        'high', cell2struct(repmat({{}}, 1, length(obj.inclusionCriteria.rats)), obj.inclusionCriteria.rats, 2));
                    vStrFiring = struct(...
                        'medium', cell2struct(repmat({{}}, 1, length(obj.inclusionCriteria.rats)), obj.inclusionCriteria.rats, 2),...
                        'high', cell2struct(repmat({{}}, 1, length(obj.inclusionCriteria.rats)), obj.inclusionCriteria.rats, 2));
                    CA1Firing = struct(...
                        'medium', cell2struct(repmat({{}}, 1, length(obj.inclusionCriteria.rats)), obj.inclusionCriteria.rats, 2),...
                        'high', cell2struct(repmat({{}}, 1, length(obj.inclusionCriteria.rats)), obj.inclusionCriteria.rats, 2));
                case 'all'
                    runningSpeed = struct('all', cell2struct(repmat({{}}, 1, length(obj.inclusionCriteria.rats)), obj.inclusionCriteria.rats, 2));
                    vStrFiring = struct('all', cell2struct(repmat({{}}, 1, length(obj.inclusionCriteria.rats)), obj.inclusionCriteria.rats, 2));
                    CA1Firing = struct('all', cell2struct(repmat({{}}, 1, length(obj.inclusionCriteria.rats)), obj.inclusionCriteria.rats, 2));
            end
            
            if plotProgress
                switch scope
                    case 'rewardExpectancy'
                        nSessions = length(cell2mat(struct2cell(obj.sessionsForSpecificAnalysis)'));
                    case 'all'
                        nSessions = length(cell2mat(struct2cell(obj.sessionsForBroadAnalysis)'));
                end
                progressBar = ['|' repmat('-', 1, nSessions+4) '|'];
                disp(progressBar)
            end

            % For all the selected sessions...
            for iRat = 1:length(obj.inclusionCriteria.rats)

                ratName = obj.inclusionCriteria.rats{iRat};

                % Load behavioural data
                behav_data = importdata(fullfile(behaviouralDataPath, sprintf('behaviour_%s.mat', lower(ratName))));

                switch scope
                    case 'rewardExpectancy'
                        sessions = obj.sessionsForSpecificAnalysis.(ratName);
                    case 'all'
                        sessions = obj.sessionsForBroadAnalysis.(ratName);
                end

                for iSess = 1:length(sessions)

                    session = sessions(iSess);
                    filepath = fullfile(ephysDataPath, sprintf('Rat_%s', ratName(1)), sprintf('Session%i', session));

                    % Load spiking data
                    data = obj.loadData(filepath);

                    % Load speed
                    speed = importdata(fullfile(filepath, 'runningSpeed.mat'));

                    % Load timestamps of trial events
                    [rewardArrivalTimes, CPEntryTimes, ~] = obj.loadRewardCPTimes(filepath);
                    switch event
                        case 'rewardArrival'
                            eventTimes = rewardArrivalTimes;
                        case 'centralPlatform'
                            eventTimes = CPEntryTimes(2:end);
                    end


                    % Divide trials into mid- and high- probability arms
                    if strcmp(scope, 'rewardExpectancy')
                        armChoice = behav_data.reward_probs{session}(behav_data.actions{session});
                        armLegitimacy = cellfun(@(x, y) x(y), behav_data.arm_values{session}, mat2cell(behav_data.actions{session}, 1, ones(behav_data.n_trials(session), 1)));
                        armLegitimacy = cellfun(@(x) ~strcmp(x, 'Illegitimate'), armLegitimacy);
                        if strcmp(event, 'centralPlatform')
                            armChoice = armChoice(2:end);
                            armLegitimacy = armLegitimacy(2:end);
                        end
                        mediumExpTrials = strcmp(armChoice, 'medium')  & armLegitimacy;
                        highExpTrials = strcmp(armChoice, 'high') & armLegitimacy;
                    end

                    % Get speed around the reward arrival times (use reward arrival as an input
                    % argument, so we can do the same for central platform entry if necessary)
                    switch scope
                        case 'rewardExpectancy'
                            [runningSpeed.medium.(ratName){session}, ~] = obj.getEventTriggeredRunningSpeed(eventTimes(boolean(mediumExpTrials)), speed.speed, speed.timestamps, 5, 5);
                            [runningSpeed.high.(ratName){session}, ts] = obj.getEventTriggeredRunningSpeed(eventTimes(boolean(highExpTrials)), speed.speed, speed.timestamps, 5, 5);
                        case 'all'
                            [runningSpeed.all.(ratName){session}, ts] = obj.getEventTriggeredRunningSpeed(eventTimes, speed.speed, speed.timestamps, 5, 5);
                    end

                    % Also get NAc and CA1 population firing
                    switch scope
                        case 'rewardExpectancy'
                            vStrFiring.medium.(ratName){session} = obj.zScore(obj.getEventTriggeredFiring(eventTimes(boolean(mediumExpTrials)), data.spiketimes(1:data.nNAcUnits), 5, 5),...
                                data.binnedfr(1:data.nNAcUnits, :));
                            CA1Firing.medium.(ratName){session} = obj.zScore(obj.getEventTriggeredFiring(eventTimes(boolean(mediumExpTrials)), data.spiketimes(data.nNAcUnits+1:end), 5, 5),...
                                data.binnedfr(data.nNAcUnits+1:end, :));
                            vStrFiring.high.(ratName){session} = obj.zScore(obj.getEventTriggeredFiring(eventTimes(boolean(highExpTrials)), data.spiketimes(1:data.nNAcUnits), 5, 5),...
                                data.binnedfr(1:data.nNAcUnits, :));
                            CA1Firing.high.(ratName){session} = obj.zScore(obj.getEventTriggeredFiring(eventTimes(boolean(highExpTrials)), data.spiketimes(data.nNAcUnits+1:end), 5, 5),...
                                data.binnedfr(data.nNAcUnits+1:end, :));

                        case 'all'
                            vStrFiring.all.(ratName){session} = obj.zScore(obj.getEventTriggeredFiring(eventTimes, data.spiketimes(1:data.nNAcUnits), 5, 5),...
                                data.binnedfr(1:data.nNAcUnits, :));
                            CA1Firing.all.(ratName){session} = obj.zScore(obj.getEventTriggeredFiring(eventTimes, data.spiketimes(data.nNAcUnits+1:end), 5, 5),...
                                data.binnedfr(data.nNAcUnits+1:end, :));
                            
                    end
                    conditions = fieldnames(vStrFiring);
                    for i = 1:length(conditions)
                        condition = conditions{i};
                        CA1Firing.(condition).(ratName){session} = reshape(CA1Firing.(condition).(ratName){session}, [], 200);
                        vStrFiring.(condition).(ratName){session} = reshape(vStrFiring.(condition).(ratName){session}, [], 200);
                    end
                    
                    if plotProgress
                        i = strfind(progressBar, '-'); progressBar(i(1)) = '+'; disp(progressBar)
                    end


                end

            end


            % Create a figure with three subplots
            colour_scheme;
            axes_properties;
            f = figure('Units', 'centimeters',...
                'Position', [9.7, 7.9, 9, 17.9]);

            switch scope
                case 'rewardExpectancy'
                    nCols = 2;
                case 'all'
                    nCols = 1;
            end

            % Top subplot: running speed
            conditions = fieldnames(vStrFiring);
            plots = [];
            for i = 1:nCols
                condition = conditions{i};
                subplot(3, nCols, i); hold on
                data = struct2cell(runningSpeed.(condition)); data = cell2mat(vertcat(transpose([data{:}])));
                if isempty(data)
                    continue
                    end
                speedMean = nanmean(data, 1);
                speedSEM = nanstd(data, [], 1) / sqrt(size(data, 1));
                patch([ts flip(ts)], [speedMean + speedSEM flip(speedMean - speedSEM)], [0.5, 0.5, 0.5], 'FaceAlpha', 0.2, 'EdgeColor', 'None')
                plot(ts, speedMean, 'k', 'LineWidth', 2.5)
                switch scope
                    case 'rewardExpectancy'
                        title([upper(condition(1)), condition(2:end)  ' exp. reward'])
                end
                set(gca, 'XTick', -5:5)
                set(gca, 'XTickLabel', [])
                xlim([-5 5])
            end
            subplot(3, nCols, 1); ylabel({'\bf{Running}', '\bf{speed (cm/s)}'})
            axis manual;
            if plotProgress
                i = strfind(progressBar, '-'); progressBar(i(1)) = '+'; disp(progressBar)
            end

            % Middle: CA1 firing
            plots = [];
            for i = 1:nCols
                condition = conditions{i};
                subplot(3, nCols, nCols+i); hold on
                data = struct2cell(CA1Firing.(condition)); data = cell2mat(vertcat(transpose([data{:}])));
                ca1Mean = mean(data);
                ca1SEM = std(data) / sqrt(size(data, 1));
                ts2 = -4.975:0.05:4.975;
                patch([ts2 flip(ts2)], [ca1Mean + ca1SEM flip(ca1Mean - ca1SEM)], 'k', 'FaceAlpha', 0.2, 'FaceColor', colourscheme.CA1,  'EdgeColor', 'None')
                plot(ts2, ca1Mean, 'Color', colourscheme.CA1, 'LineWidth', 2.5)
                set(gca, 'XTick', -5:5)
                set(gca, 'XTickLabel', [])
            end
            subplot(3, nCols, nCols+1); ylabel({'\bf{CA1}', 'Firing rate (z-score)'}, 'FontSize', 14)
            axis manual;
            if plotProgress
                i = strfind(progressBar, '-'); progressBar(i(1)) = '+'; disp(progressBar)
            end

            % Bottom: NAc firing
            plots = [];
            for i = 1:nCols
                condition = conditions{i};
                subplot(3, nCols, (nCols*2)+i); hold on
                data = struct2cell(vStrFiring.(condition)); data = cell2mat(vertcat(transpose([data{:}])));
                vstrMean = mean(data);
                vstrSEM = std(data) / sqrt(size(data, 1));
                ts2 = -4.975:0.05:4.975;
                patch([ts2 flip(ts2)], [vstrMean + vstrSEM flip(vstrMean - vstrSEM)], 'k', 'FaceAlpha', 0.2, 'FaceColor', colourscheme.vStr, 'EdgeColor', 'None')
                plot(ts2, vstrMean, 'Color', colourscheme.vStr, 'LineWidth', 2.5)
                set(gca, 'XTick', -5:5)
                switch event
                    case 'rewardArrival'
                        xlabel({'Time from', 'arrival at reward location (s)'})
                    case 'centralPlatform'
                        xlabel({'Time from', 'central platform entry (s)'})
                end
            end
            subplot(3, nCols, (2*nCols)+1); ylabel({'\bf{vStr}', 'Firing rate (z-score)'}, 'FontSize', 14)
            axis manual;
            for i = 1:nCols*3
                subplot(3, nCols, i)
                y = ylim;
                plot([0 0], y, 'k', 'LineStyle', '--')
            end
            if plotProgress
                i = strfind(progressBar, '-'); progressBar(i(1)) = '+'; disp(progressBar)
            end

            switch scope
                case 'rewardExpectancy'
                    linkaxes(f.Children(1:2))
                    linkaxes(f.Children(3:4))
                    linkaxes(f.Children(5:6))
            end
            if plotProgress
                i = strfind(progressBar, '-'); progressBar(i(1)) = '+'; disp(progressBar)
            end
            
        end
        
        function peakCoactivity = getPeakCoactivity(obj, arm, trialType)
            
            % Get coactivity for the 2 seconds before an event
            % 
            % Inputs
            %   - arm:                   'high' or 'medium'; must be a field of the perCellRewardExpectation property
            %   - trialType:            'reactivated' 'control' or NaN, if not Nan, must be a subfield of the arm field
            % 
            % Outputs
            %   - peakCoactivity:   peak coactivity in the 2 seconds prior to event

            obj.rewardFiring = obj.taskCoactivity;
            peakCoactivity = getPeakFiringRates(obj, arm, trialType);
            
        end

        function obj = plotPeakPreRewardCoactivity(obj, p_medium, p_high, p_interaction)
            
            % Check that there is data to plot
            assert(isstruct(obj.taskCoactivity) && ~isempty(obj.taskCoactivity), 'ExplainedVarianceReactivation:dataNotAssigned', 'No data to plot. Please run the getReactivatedCellPairActivity() method to generate data.')
            obj.rewardFiring = obj.taskCoactivity;
            
            % Set design properties
            colour_scheme;
            axes_properties
            
            fig = figure('Units', 'centimeters',...
                'Position', [18.0, 17.2, 12.2, 8.7]); hold on
            
             % Calculate peak coactivity in the 2 seconds prior to event on medium-expectation trials
             % for reactivated and control cell pairs
            mediumReactivatedRewardExpectationPeaks = obj.getPeakCoactivity('medium', 'reactivated');
            mediumControlRewardExpectationPeaks = obj.getPeakCoactivity('medium', 'control');

            % Calculate peak coactivity in the 2 seconds prior to event on high-expectation trials
             % for reactivated and control cell pairs
            highReactivatedRewardExpectationPeaks = obj.getPeakCoactivity('high', 'reactivated');
            highControlRewardExpectationPeaks = obj.getPeakCoactivity('high', 'control');
                
            % Plot bars
            fig = obj.plotBar(mediumControlRewardExpectationPeaks, 'insignificant', 1, fig);
            fig = obj.plotBar(mediumReactivatedRewardExpectationPeaks, 'medium_exp_reward', 2, fig);
            fig = obj.plotBar(highControlRewardExpectationPeaks, 'insignificant', 4, fig);
            fig = obj.plotBar(highReactivatedRewardExpectationPeaks, 'high_exp_reward', 5, fig);

            % Significance markers
            y = ylim;
             y(2) = y(2)*1.2; ylim(y)
            p_values = [p_medium, p_high, p_interaction];
            xData = {[1, 2], [4, 5], [1.5, 4.5]};
            yData = {[0.85 0.9, 0.91, 0.935], [0.85, 0.9, 0.91, 0.935], [0.96, 1, 1.01, 1.035]};
            for i = 1:3
                line([xData{i}(1) xData{i}(1)], [yData{i}(1:2) * y(2)], 'Color', 'k', 'LineWidth', 1.5)
                line([xData{i}(2) xData{i}(2)], [yData{i}(1:2) * y(2)], 'Color', 'k', 'LineWidth', 1.5)
                line(xData{i}, [[yData{i}(2) yData{i}(2)] * y(2)], 'Color', 'k', 'LineWidth', 1.5)
                if p_values(i) < 0.05
                    text(mean(xData{i}), yData{i}(3) * y(2), '*', 'FontSize', 24, 'HorizontalAlignment', 'center')
                else
                    text(mean(xData{i}), yData{i}(3) * y(2), 'n.s.', 'FontSize', 13, 'HorizontalAlignment', 'center')
                end
            end
            
            % Axis formatting
            ylabel('Coactivity (z-score)')
            set(gca, 'XTick', [1.5 4.5])
            xlim([0.5 5.5])
            set(gca, 'XTickLabel', {'Medium', 'High'})
            xlabel('Reward expectation')

            obj.figures{end+1} = fig;
            
        end
        
        function obj = plotCoactivityIllustration(obj, varargin)
            
            % Plots an example of a significantly reactivated CA1-vStr cell
            % pair, their event-triggered coactivity on medium- and
            % high-reward-expectation trials, and the peri-spike histogram
            % of ventral striatal firing rates relative to CA1 spikes in
            % PRE, TASK and POST.
            
            % Process inputs
            p = inputParser;
            addParameter(p, 'ratName', 'Quirinius', @(x) ischar(x) || isstring(x));
            addParameter(p, 'session', 8, @(x) isnumeric(x) && floor(x)==x);
            addParameter(p, 'cellPairIndex', 7, @(x) isnumeric(x) && floor(x)==x);
            addParameter(p, 'ephysDataPath', './data/ephys_data', @(x) ischar(x) || isstring(x));
            addParameter(p, 'behaviouralDataPath', './data/behavioural_data', @(x) ischar(x) || isstring(x));
            parse(p, varargin{:});
            ratName = p.Results.ratName;
            session = p.Results.session;
            cellPairIndex = p.Results.cellPairIndex;
            ephysDataPath = p.Results.ephysDataPath;
            behaviouralDataPath = p.Results.behaviouralDataPath;
            
            % Select pair to visualise
            reactivatedCells = obj.perCellContributions.reactivatedCellPairs.(ratName){session};
            pair = reactivatedCells(cellPairIndex, :);

            % Get spiketimes
            filepath = fullfile(ephysDataPath, sprintf('Rat_%s', ratName(1)), sprintf('Session%i', session));
            data = obj.loadData(filepath);
            pair(2) = pair(2) + data.nNAcUnits;

            % Get behavioural timestamps
            [rewardArrivalTimes, ~, ~] = obj.loadRewardCPTimes(filepath);
            behav_data = importdata(fullfile(behaviouralDataPath, sprintf('behaviour_%s.mat', lower(ratName))));
            arm_choice = behav_data.reward_probs{session}(behav_data.actions{session});
            arm_legitimacy = cellfun(@(x, y) x(y), behav_data.arm_values{session}, mat2cell(behav_data.actions{session}, 1, ones(behav_data.n_trials(session), 1)));
            arm_legitimacy = cellfun(@(x) ~strcmp(x, 'Illegitimate'), arm_legitimacy);
            medium_exp_trials = strcmp(arm_choice, 'medium')  & arm_legitimacy;
            high_exp_trials = strcmp(arm_choice, 'high') & arm_legitimacy;
            highExpArrivalTimes = rewardArrivalTimes(high_exp_trials);
            mediumExpArrivalTimes = rewardArrivalTimes(medium_exp_trials);

            % Plot high
            axes_properties;
            colour_scheme;
            f_high = figure('Units', 'centimeters',...
                'Position', [18 14.7 14.8 11.1]);
            subplot(4, 1, 1:3); hold on; set(gca, 'YDir', 'reverse'); title('High reward expectation trials')
            tickHeight = 0.5; coactivityScalingFactor = 5;
            vStr_spiketimes = data.spiketimes{pair(1)};
            CA1_spiketimes = data.spiketimes{pair(2)};
            for iTrial = 1:length(highExpArrivalTimes)

                trialSpikes_vStr = vStr_spiketimes(vStr_spiketimes > (highExpArrivalTimes(iTrial) - 5) & vStr_spiketimes < (highExpArrivalTimes(iTrial) + 5)) - highExpArrivalTimes(iTrial);
                trialSpikes_CA1 = CA1_spiketimes(CA1_spiketimes > (highExpArrivalTimes(iTrial) - 5) & CA1_spiketimes < (highExpArrivalTimes(iTrial) + 5)) - highExpArrivalTimes(iTrial);
                yValue = iTrial*3-3;
                for iSpike = 1:length(trialSpikes_vStr)
                    line([trialSpikes_vStr(iSpike), trialSpikes_vStr(iSpike)], [0,  tickHeight] + yValue, 'Color', [0.9 0.5 0.9])
                end
                yValue = iTrial*3-2.5;
                for iSpike = 1:length(trialSpikes_CA1)
                    line([trialSpikes_CA1(iSpike), trialSpikes_CA1(iSpike)], [0,  tickHeight] + yValue, 'Color', [0.5 0.9 0.5])
                end

                coactivity = squeeze(obj.getEventTriggeredCoactivity(highExpArrivalTimes(iTrial), {vStr_spiketimes}, {CA1_spiketimes}, 5, 5, obj.analysisConfig.gaussianWindowWidth));
                coactivityTransformed = 1 - (coactivity / max(coactivity));
                yData = [coactivityTransformed' ones(1, length(coactivityTransformed))] + (iTrial*3-2);
                xData = [-4.975:0.05:4.975 flip(-4.975:0.05:4.975)];
                patch(xData, yData, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.25)

            end

            % Axis formatting
            xlim([-5 5])
            ylim([-0.5 iTrial*3-1])
            set(gca, 'XTick', [-5, 0, 5])
            set(gca, 'XTickLabel', [])
            set(gca, 'XColor', [1 1 1])
            set(gca, 'YTick', [0 : iTrial-1]*3 + 0.75)
            set(gca, 'YTickLabel', 1:iTrial)
            ylabel('Trial')
            text(4.5, 0, 'vStr', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.9, 0.5, 0.9], 'HorizontalAlignment', 'right')
            text(4.5, 1.25, 'CA1', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0.9, 0.5], 'HorizontalAlignment', 'right')
            text(4.75, 2.25, 'coactivity', 'FontSize', 16, 'Color', [0.6, 0.6, 0.6], 'HorizontalAlignment', 'right')

            % Plot mean z-scored coactivity
            meanCoactivity = squeeze(mean(obj.getEventTriggeredCoactivity(highExpArrivalTimes', {vStr_spiketimes}, {CA1_spiketimes}, 5, 5, obj.analysisConfig.gaussianWindowWidth), 1));
            zScoredMeanCoactivity = obj.zScore(meanCoactivity, min(data.binnedfr(pair(1), :), data.binnedfr(pair(2), :)));
            z_high = subplot(4, 1, 4);  plot(-4.975:0.05:4.975, zScoredMeanCoactivity, 'k', 'LineWidth', 3);

            axis tight; xlim([-5 5]); box off
            xlabel({'Time from arrival', 'at reward location (s)'})
            ylabel({'Mean', 'coactivity', '(z-score)'})
            drawnow;

            % Align axes
            pos1 = max([f_high.Children(1).Position(1) f_high.Children(2).Position(1)]);
            pos3 = min([f_high.Children(1).Position(3) f_high.Children(2).Position(3)]);
            for i = 1:2
                f_high.Children(i).Position(1) = pos1;
                f_high.Children(i).Position(3) = pos3;
            end

            % Plot medium
            f_medium = figure('Units', 'centimeters',...
                'Position', [18 14.7 14.8 11.1]);
            subplot(4, 1, 1:3); hold on; set(gca, 'YDir', 'reverse'); title('Medium reward expectation trials')
            tickHeight = 0.5; coactivityScalingFactor = 5;
            vStr_spiketimes = data.spiketimes{pair(1)};
            CA1_spiketimes = data.spiketimes{pair(2)};
            for iTrial = 1:length(mediumExpArrivalTimes)

                trialSpikes_vStr = vStr_spiketimes(vStr_spiketimes > (mediumExpArrivalTimes(iTrial) - 5) & vStr_spiketimes < (mediumExpArrivalTimes(iTrial) + 5)) - mediumExpArrivalTimes(iTrial);
                trialSpikes_CA1 = CA1_spiketimes(CA1_spiketimes > (mediumExpArrivalTimes(iTrial) - 5) & CA1_spiketimes < (mediumExpArrivalTimes(iTrial) + 5)) - mediumExpArrivalTimes(iTrial);
                yValue = iTrial*3-3;
                for iSpike = 1:length(trialSpikes_vStr)
                    line([trialSpikes_vStr(iSpike), trialSpikes_vStr(iSpike)], [0,  tickHeight] + yValue, 'Color', [0.9 0.5 0.9])
                end
                yValue = iTrial*3-2.5;
                for iSpike = 1:length(trialSpikes_CA1)
                    line([trialSpikes_CA1(iSpike), trialSpikes_CA1(iSpike)], [0,  tickHeight] + yValue, 'Color', [0.5 0.9 0.5])
                end

                coactivity = squeeze(obj.getEventTriggeredCoactivity(mediumExpArrivalTimes(iTrial), {vStr_spiketimes}, {CA1_spiketimes}, 5, 5, obj.analysisConfig.gaussianWindowWidth));
                coactivityTransformed = 1 - (coactivity / max(coactivity));
                yData = [coactivityTransformed' ones(1, length(coactivityTransformed))] + (iTrial*3-2);
                xData = [-4.975:0.05:4.975 flip(-4.975:0.05:4.975)];
                patch(xData, yData, 'k', 'EdgeColor', 'none', 'FaceAlpha', 0.25)

            end

            xlim([-5 5])
            ylim([-0.5 iTrial*3-1])
            set(gca, 'XTick', [-5, 0, 5])
            set(gca, 'XTickLabel', [])
            set(gca, 'XColor', [1 1 1])
            set(gca, 'YTick', [0 : iTrial-1]*3 + 0.75)
            set(gca, 'YTickLabel', 1:iTrial)
            ylabel('Trial')

            text(4.5, 0, 'vStr', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.9, 0.5, 0.9], 'HorizontalAlignment', 'right')
            text(4.5, 1.25, 'CA1', 'FontSize', 16, 'FontWeight', 'bold', 'Color', [0.5, 0.9, 0.5], 'HorizontalAlignment', 'right')
            text(4.75, 2.25, 'coactivity', 'FontSize', 16, 'Color', [0.6, 0.6, 0.6], 'HorizontalAlignment', 'right')

            % Plot mean z-scored coactivity
            meanCoactivity = squeeze(mean(obj.getEventTriggeredCoactivity(mediumExpArrivalTimes', {vStr_spiketimes}, {CA1_spiketimes}, 5, 5, obj.analysisConfig.gaussianWindowWidth), 1));
            zScoredMeanCoactivity = obj.zScore(meanCoactivity, min(data.binnedfr(pair(1), :), data.binnedfr(pair(2), :)));
            z_medium = subplot(4, 1, 4); plot(-4.975:0.05:4.975, zScoredMeanCoactivity, 'k', 'LineWidth', 3);

            axis tight; xlim([-5 5]); box off
            xlabel({'Time from arrival', 'at reward location (s)'})
            ylabel({'Mean', 'coactivity', '(z-score)'})
            drawnow;

            % Align axes
            pos1 = max([f_medium.Children(1).Position(1) f_medium.Children(2).Position(1)]);
            pos3 = min([f_medium.Children(1).Position(3) f_medium.Children(2).Position(3)]);
            for i = 1:2
                f_medium.Children(i).Position(1) = pos1;
                f_medium.Children(i).Position(3) = pos3;
            end
            if z_high.YLim(2) > z_medium.YLim(2)
                linkaxes([z_high z_medium])
            else
                linkaxes([z_medium z_high])
            end
            z_high.YTick = z_high.YTick([1, end]);
            z_medium.YTick = z_medium.YTick([1, end]);
            z_medium.Position(4) = z_medium.Position(4) + 0.05;
            z_high.Position(4) = z_high.Position(4) + 0.05;

            % Take a random selection of 1000 ripples in PRE, and find all  CA1 & vStr spikes wtihin 200ms around riple peak
            rng(0),
            ripplesPRE = data.ripplePeak(data.ripplePeak > data.startEndTimes(1, 1) & data.ripplePeak < data.startEndTimes(1, 2));
            CA1Spikes = cell2mat(arrayfun(@(x) CA1_spiketimes(CA1_spiketimes > (x - 0.1) & CA1_spiketimes < (x + 0.1)), ripplesPRE, 'UniformOutput', false)');
            vStrSpikes = arrayfun(@(x) vStr_spiketimes(vStr_spiketimes > (x - 0.1) & vStr_spiketimes < (x + 0.1)) - x, CA1Spikes, 'UniformOutput', false);
            vStrFiringRatePRE = cell2mat(cellfun(@(spikes) squeeze(obj.convolveAndBin(-0.1, 0.2, 0.01, {spikes})), vStrSpikes, 'UniformOutput', false)');

            % Repeat for POST
            rng(0),
            ripplesPOST = data.ripplePeak(data.ripplePeak > data.startEndTimes(3, 1) & data.ripplePeak < data.startEndTimes(3, 2));
            CA1Spikes = cell2mat(arrayfun(@(x) CA1_spiketimes(CA1_spiketimes > (x - 0.1) & CA1_spiketimes < (x + 0.1)), ripplesPOST, 'UniformOutput', false)');
            vStrSpikes = arrayfun(@(x) vStr_spiketimes(vStr_spiketimes > (x - 0.1) & vStr_spiketimes < (x + 0.1)) - x, CA1Spikes, 'UniformOutput', false);
            vStrFiringRatePOST = cell2mat(cellfun(@(spikes) squeeze(obj.convolveAndBin(-0.1, 0.2, 0.01, {spikes})), vStrSpikes, 'UniformOutput', false)');

            % Repeat for all activity in TASK
            rng(0);
            CA1Spikes = CA1_spiketimes(CA1_spiketimes > data.startEndTimes(2, 1) & CA1_spiketimes < data.startEndTimes(2, 2));
            vStrSpikes = arrayfun(@(x) vStr_spiketimes(vStr_spiketimes > (x - 0.1) & vStr_spiketimes < (x + 0.1)) - x, CA1Spikes, 'UniformOutput', false);
            vStrFiringRateTASK = cell2mat(cellfun(@(spikes) squeeze(obj.convolveAndBin(-0.1, 0.2, 0.01, {spikes})), vStrSpikes, 'UniformOutput', false)');

            % Plot
            f = figure;
            meanPRE = mean(vStrFiringRatePRE, 2);
            semPRE = std(vStrFiringRatePRE, [], 2) / sqrt(size(vStrFiringRatePRE, 2));
            meanTASK = mean(vStrFiringRateTASK, 2);
            semTASK = std(vStrFiringRateTASK, [], 2) / sqrt(size(vStrFiringRateTASK, 2));
            meanPOST = mean(vStrFiringRatePOST, 2);
            semPOST = std(vStrFiringRatePOST, [], 2) / sqrt(size(vStrFiringRatePOST, 2));
            subplot(1, 3, 1); errorbar(-0.095:0.01:0.095, meanPRE, semPRE, 'Color', colourscheme.vStr, 'LineWidth', 3); title('PRE'); xlabel({'Time from', 'CA1 spike (s)'}); ylabel({'vStr firing', 'rate (Hz)'})
            subplot(1, 3, 2); errorbar(-0.095:0.01:0.095, meanTASK, semTASK, 'Color', colourscheme.vStr, 'LineWidth', 3); title('TASK'); xlabel({'Time from', 'CA1 spike (s)'})
            subplot(1, 3, 3); errorbar(-0.095:0.01:0.095, meanPOST, semPOST, 'Color', colourscheme.vStr, 'LineWidth', 3); title('POST'); xlabel({'Time from', 'CA1 spike (s)'})
            [~, order] = sort([f.Children(1).YLim(2) f.Children(2).YLim(2) f.Children(3).YLim(2)], 'descend');
            linkaxes(f.Children(order))

        end
        
    end
    
end
