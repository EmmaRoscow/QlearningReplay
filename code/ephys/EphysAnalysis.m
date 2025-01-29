

classdef EphysAnalysis
    
    properties (Access = public)
        inclusionCriteria
        analysisConfig
        sessionsForBroadAnalysis
    end
    
    
    properties (Constant)
        BinSize = 0.05;
    end

    
    methods (Access = public)
                
        % Load data
        function data = loadData(~, filepath)
            
            % Loads .mat files containing spiketimes and other data for electrophysiological analysis
            % 
            %   Inputs
            %       - filepath: path of the directory where .mat files are stored
            % 
            %   Outputs
            %       -  data:     structure containing data with the following fields
            %           - spiketimes:      cell array of spiketimes in seconds; one cell per neuron
            %           - nNAcUnits:      number of nucleus accumbens (ventral
            %                                   striatal) neurons in the data. Cells 1:nNAcUnits of
            %                                   spiketimes correspond to accumbens units; the rest are
            %                                   hippocampal units
            %           - nCA1Units:      number of hippocampal (CA1) neurons in the spiketimes data
            %           - nUnits:            total number of neurons in the spiketimes data
            %           - binnedfr:          pre-computed firing rates in 50ms bins,
            %                                   calculated from the spiketimes cell array; one row per unit
            %           - startEndTimes: 3 x 2 array of timestamps in seconds,
            %                                   marking the start and end of the pre-task rest period (first row), 
            %                                   start and end of the task period (second row), 
            %                                   and start and end of the post-task rest period (third row)
            %           - rippleStart:       timestamps in seconds marking the start of each identified sharp-wave ripple
            %           - rippleStop:       timestamps in seconds marking the end of each identified sharp-wave ripple
            %           - ripplePeak:      timestamps in seconds marking the peak amplitude of each identified sharp-wave ripple
            
            % Check that the filepath exists
            if exist(filepath, 'dir') ~= 7
                error('EphysAnalysis:load:couldNotFindPath', 'Specified filepath "%s" does not exist.', filepath)
            end
            
            % Check that all required files are present
            requiredFiles = {'spiketimes.mat', 'nNAcUnits.mat', 'startEndTimes.mat', 'ripples3std.mat'};
            for k = 1:length(requiredFiles)
                if exist(fullfile(filepath, requiredFiles{k}), 'file') ~= 2
                    error('EphysAnalysis:load:couldNotReadFile', 'Required file "%s" not found in path "%s".', requiredFiles{k}, filepath);
                end
            end

            % Load files
            load([filepath, filesep, 'spiketimes.mat'])
            load([filepath, filesep, 'nNAcUnits.mat'])
            try
                binnedfr = importdata([filepath, filesep, 'binnedfr.mat']);
            catch
                load([filepath, filesep, 'binlessfr50ms.mat'])
                binnedfr = binlessfr;
                clear binlessfr
            end
            load([filepath, filesep, 'startEndTimes.mat'])
            load([filepath, filesep, 'ripples3std.mat'])
            nUnits = size(spiketimes, 2);
            nCA1Units = nUnits - nNAcUnits;
            
            % Validate loaded variables
            requiredVars = {'spiketimes', 'nNAcUnits', 'startEndTimes', 'binnedfr', ...
                            'rippleStart', 'rippleStop', 'ripplePeak'};
            for var = requiredVars
                if exist(var{1}) ~= 1
                    error('EphysAnalysis:load:missingData', 'Variable "%s" is missing in the loaded data.', var{1});
                end
            end
            
            data.spiketimes = spiketimes;
            data.nNAcUnits = nNAcUnits;
            data.binnedfr = binnedfr;
            data.startEndTimes = startEndTimes;
            data.rippleStart = rippleStart;
            data.rippleStop = rippleStop;
            data.ripplePeak = ripplePeak;
            data.nUnits = nUnits;
            data.nCA1Units = nCA1Units;
            
        end

        
        % Load ventral striatal spiking data
        function spiketimes = loadVentralStriatalSpikes(~, filepath)
            
            % Check that the filepath exists
            if exist(filepath, 'dir') ~= 7
                error('EphysAnalysis:load:couldNotFindPath', 'Specified filepath "%s" does not exist.', filepath)
            end
            
            % Check that all required files are present
            requiredFiles = {'spiketimes.mat', 'nNAcUnits.mat'};
            for k = 1:length(requiredFiles)
                if exist(fullfile(filepath, requiredFiles{k}), 'file') ~= 2
                    error('EphysAnalysis:load:couldNotReadFile', 'Required file "%s" not found in path "%s".', requiredFiles{k}, filepath);
                end
            end

            spiketimes = importdata(fullfile(filepath, 'spiketimes.mat'));
            nNAcUnits = importdata(fullfile(filepath, 'nNAcUnits.mat'));
            
            % Validate loaded variables
            requiredVars = {'spiketimes', 'nNAcUnits'};
            for var = requiredVars
                if exist(var{1}) ~= 1
                    error('EphysAnalysis:load:missingData', 'Variable "%s" is missing in the loaded data.', var{1});
                end
            end
            isNonNegativeInteger = @(x)numel(x)==1 && isfinite(x) && x==floor(x) && x>=0;
            assert(isNonNegativeInteger(nNAcUnits), 'EphysAnalysis:load:invalidData', 'nNAcUnits must be a non-negative integer')
            assert(iscell(spiketimes), 'EphysAnalysis:load:invalidData', 'Spiketimes must be a cell array')
            assert(length(spiketimes) >= nNAcUnits, 'EphysAnalysis:load:invalidData', sprintf('There are %i ventral striatal units but only %i cells of spiketimes', nNAcUnits, length(spiketimes)))
            spiketimes = spiketimes(1:nNAcUnits);
        end
        
        
        % Load behavioural timestamps
        function [rewardArrivalTimes, CPEntryTimes, CPExitTimes] = loadRewardCPTimes(~, filepath)
            
            % Loads .mat files containing timestamps for rat-initiatiated trial events
            % 
            %   Inputs
            %       - filepath:                     path of the directory where .mat files are stored
            % 
            %   Outputs
            %       -  rewardArrivalTimes:    timestamps in seconds of arrival at the reward locations
            %       - CPEntryTimes:           timestamps in seconds of entry to the central platform
            %       - CPExitTimes:             timestamps in seconds of exit from the central platform
            
            rewardArrivalTimes = importdata([filepath filesep 'RarrivalTimes.mat']);
            CPEntryTimes = importdata([filepath filesep 'CPentryTimes.mat']);
            CPExitTimes = importdata([filepath filesep 'CPexitTimes.mat']);
            
            % Check that there is roughly one of each per trial. Note that
            % rat may enter the central platform but not leave it at the
            % end of the session, or leave it but not reach the reward
            % location, so there can be a legitimate difference of 1
            assert(ismember(length(CPExitTimes) - length(CPEntryTimes), [0, 1]), 'EphysAnalysis:load:invalidData', 'Number of timestamps for entry to and exit from the central platform do not match')
            assert(ismember(length(CPExitTimes) - length(rewardArrivalTimes), [0, 1]), 'EphysAnalysis:load:invalidData', 'Number of timestamps for arrival at reward location does not match number of exits from central platform')
        end
        
        
        % Select reward-responsive cells for analysis
        function rewardResponsive = selectRewardResponsiveCells(obj, spiketimes, rewardArrivalTimes, CPEntryTimes)
            
            % Computes a binary vector indicating which cells are
            % statistically responsive to reward-related events, if used as
            % an inclusion criterion for analysis; otherwise returns a
            % vector of ones
            % 
            %   Inputs
            %       - spiketimes:               cell array of spiketimes, one cell per neuron
            %       - rewardArrivalTimes:   timestamps of reward-related events
            %       - CPEntryTimes:         timestamps of control events, specifically entry to the central platform
            % 
            %   Outputs
            %       -  rewardResponsive:   binary vector with the same length as spiketimes, indicating the cells to be
            %                                        included according to the inclusion criteria for reward-responsivity
            
            nCells = length(spiketimes);
            
            if obj.inclusionCriteria.cells.rewardResponsive
                
                rewardResponsive = obj.identifyRewardResponsiveCells(spiketimes, rewardArrivalTimes, CPEntryTimes);
                
            else
                
                rewardResponsive = boolean(ones(1, nCells));
                
            end
            
        end
            
        
        % Identify reward-responsive cells for analysis
        function rewardResponsive = identifyRewardResponsiveCells(~, spiketimes, rewardArrivalTimes, CPEntryTimes)
            
            % Computes a binary vector indicating which cells are statistically responsive to reward-related events, 
            % by counting spikes in each of 250ms bins 1s before to 1s after arrival time,
            % as per https://doi.org/10.1371/journal.pbio.1000173,
            % compared to three control bins before entry to central platform
            % 
            %   Inputs
            %       - spiketimes:               cell array of spiketimes, one cell per neuron
            %       - rewardArrivalTimes:   timestamps of reward-related events
            %       - CPEntryTimes:         timestamps of control events, specifically entry to the central platform
            % 
            %   Outputs
            %       -  rewardResponsive:   binary vector with the same length as spiketimes, indicating the cells to be
            %                                        included according to the inclusion criteria for reward-responsivity
            
            function bins = getBins(timestamps, offsets)
                % Computes bin edges as offsets relative to trial timestamps
                bins = arrayfun(@(x) x + offsets, timestamps, 'UniformOutput', false);
                bins = cell2mat(bins);
            end
            
            function spikeCount = getBinnedSpikeCounts(spiketimes, binStarts, binEnds)
                % Computes the number of spikes in each cell of spiketimes,
                % in each bin defined by binStarts and binEnds
                spikeCount = arrayfun(@(x, y) sum(spiketimes > x & spiketimes < y), binStarts, binEnds);
            end
                        
            nCells = length(spiketimes);
            nTrials = length(rewardArrivalTimes);
            
                % Find reward-responsive cells
                rewardResponsive = zeros(1, nCells);
                for iCell = 1:nCells
                    
                    % Pre-allocate spike counts for all relevant bins
                    rewardSpikeCount = zeros(nTrials, 8);
                    controlSpikeCount = zeros(nTrials, 3);
                                        
                    % Define reward bins (1s before to 1s after reward arrival time) 
                    % and control bins (0.75s to 0s before  central
                    % platform arrival time) for all trials
                    RBinStarts = getBins(rewardArrivalTimes(:), -1:0.25:0.75);
                    RBinEnds = getBins(rewardArrivalTimes(:), -0.75:0.25:1);
                    CPBinStarts = getBins(CPEntryTimes(:), -0.75:0.25:-0.25);
                    CPBinEnds = getBins(CPEntryTimes(:), -0.5:0.25:0);
                    
                    % Count spikes in each bin
                    for iTrial = 1:nTrials
                        rewardSpikeCount(iTrial, :) = getBinnedSpikeCounts(spiketimes{iCell}, RBinStarts(iTrial, :), RBinEnds(iTrial, :));
                        controlSpikeCount(iTrial, :) = getBinnedSpikeCounts(spiketimes{iCell}, CPBinStarts(iTrial, :), CPBinEnds(iTrial, :));
                    end
                    
                    % At least one reward-related bin must be significantly
                    % different from all three control bins for cell to be
                    % labelled as reward-responsive, as per above paper
                    for iBin = 1:size(RBinStarts, 2)
                        pvals = cellfun(@(x) ranksum(rewardSpikeCount(:, iBin), x), mat2cell(controlSpikeCount, nTrials, ones(3, 1)));
                        alpha = 0.05;
                        if sum(pvals < alpha)==3
                            rewardResponsive(iCell) = 1;
                        end
                    end
                    
                end
                    
            fprintf('%i cells out of %i are reward-responsive\n', sum(rewardResponsive), length(rewardResponsive))
            rewardResponsive = boolean(rewardResponsive);
        end
        
        
        % Calculate the proportion of reward-responsive ventral striatal cells
        function perRatRewardResponsiveRate = countRewardResponsiveCells(obj, varargin)
            
            % Identifies how many ventral striatal cells in each session
            % are significantly responsive to reward-related events, 
            % prints the mean and standard error overall and for each rat,
            % and returns the proportion per rat
            %
            % Inputs
            %   - varargin:                                  optionally, 'DataPath' along with the directory where data 
            %                                                   for spike times and other data are
            %                                                   stored (otherwise uses default path)
            %   
            % Outputs
            %   - perRatRewardResponsiveRate:  structure with one field for each rat name per inclusion criteria, corresponding to each
            %                                                   rat's total number of reward-responsive striatal cells as a proportion of 
            %                                                   rat's total number of striatal cells across sessions per inclusion criteria
            
            % Process path where data is stored
            p = inputParser;
            addParameter(p, 'DataPath', fullfile('.', 'data', 'ephys_data'), @(x) ischar(x) || isstring(x));
            parse(p, varargin{:});
            datapath = p.Results.DataPath;
            
            % Create variables for computing rate of reward-responsiveness
            nRats = length(obj.inclusionCriteria.rats);
            perRatRewardResponsiveRate = cell2struct(mat2cell(zeros(1, nRats), 1, ones(1, nRats)), obj.inclusionCriteria.rats, 2);
            overallRewardResponsiveTotal = 0;
            overallTotal = 0;

            for iRat = 1:length(obj.inclusionCriteria.rats)

                ratName = obj.inclusionCriteria.rats{iRat};
                ratSessions = obj.sessionsForBroadAnalysis.(ratName);

                perRatRewardResponsive = 0;
                perRatTotal = 0;

                for iSess = 1:numel(ratSessions)

                    session = ratSessions(iSess);

                    % Load session data
                    filepath = fullfile(datapath, ['Rat_' ratName(1)], ['Session' num2str(session)]);
                    data = obj.loadData(filepath);
                    vStrSpiketimes = data.spiketimes(1:data.nNAcUnits);
                    [rewardArrivalTimes, CPEntryTimes, ~] = obj.loadRewardCPTimes(filepath);
                    
                    % Identify reward-responsive cells for the session
                    rewardResponsive = obj.identifyRewardResponsiveCells(vStrSpiketimes, rewardArrivalTimes, CPEntryTimes);
                    
                    % Update per-rat counts
                    perRatRewardResponsive = perRatRewardResponsive + sum(rewardResponsive);
                    perRatTotal = perRatTotal + length(rewardResponsive);

                end
                
                % Update overall counts
                overallRewardResponsiveTotal = overallRewardResponsiveTotal + perRatRewardResponsive;
                overallTotal = overallTotal + perRatTotal;
                perRatRewardResponsiveRate.(ratName) = perRatRewardResponsive / perRatTotal;

            end

            % Compute statistics
            meanPerRatRewardRepsonsiveRate = mean(cell2mat(struct2cell(perRatRewardResponsiveRate)));
            semPerRateRewardResponsiveRate = std(cell2mat(struct2cell(perRatRewardResponsiveRate))) / sqrt(nRats);

            % Print statistics
            fprintf('%i out of %i (%.0f%%) of ventral striatal cells show significant reward-related firing\n',...
                overallRewardResponsiveTotal,...
                overallTotal,...
                overallRewardResponsiveTotal / overallTotal * 100)

            fprintf('On average %.2f%% ± %.1f%% ventral striatal units per rat show significant reward-related firing\n',...
                meanPerRatRewardRepsonsiveRate * 100,...
                semPerRateRewardResponsiveRate * 100)
                
        end
        
        
        % Convolve spike trains using binless method
        function [popV] = convolveAndBin(~, startTimes, duration, binSize, spiketimes, varargin)
            
            % Convolves and "bins" spike trains in a binless way, based on Kruskal et al. (2007) https://doi.org/10.1002/sim.2946
            % 
            % Inputs
            %   - startTimes:   vector of timestamps from which to calculate firing rates (e.g. start of trials)
            %   - duration:       desired duration of firing rate calculations, such that the returned firing rates
            %                        correspond to the period (startTimes) to (startTimes + duration)
            %   - binSize:        compute firing rates equivalent to binning with this window; returned data corresponds to
            %                        timestamps with intervals of binSize
            %   - spiketimes:   cell array of spike times, one cell per neuron
            %   - varargin:       optionally, a parameter for calculating the Gaussian kernel width
            %
            % Outputs
            %   - popV:           population vector(s) of firing rates of size T x N x B where T is length of startTimes, 
            %                        N is length of spiketimes, and B is duration / binSize

            if isempty(varargin)
                % Calculate sigma for Gaussian from binsize
                sigma = binSize/sqrt(12);
            else
                sigma = varargin{1};
            end
    
            tolerance = 1e-3;
    
            for iTrial = 1:size(startTimes,1)

                % Get desired timestamps (i.e. bins) starting from desired start times
                 t = [(binSize/2):binSize:(duration-binSize/2)] + startTimes(iTrial);
                 finalBinSize = rem(duration, binSize);
                 if finalBinSize > tolerance
                    if isempty(t)
                        finalT = startTimes(iTrial) + finalBinSize/2;
                    else
                         finalT = t(end) + finalBinSize/2;
                    end
                    finalBinSigma = sigma * finalBinSize/binSize;
                 end

                 for iUnit = 1:size(spiketimes,2)
            
                    % Find spike times that are within 5 seconds of these timestamps
                    st = arrayfun(@(x,y) x{1}(x{1} >= (y-5) & x{1} < (y+5)), repmat({spiketimes{iUnit}},1,length(t)),t,'UniformOutput',false);

                    % Convolve
                    conv = arrayfun(@(x,y) (1/sqrt(2*pi*sigma.^2))*exp(-((y-x{1}).^2)/(2*sigma.^2)), st, t, 'UniformOutput', false);
            
                % Convolve for final, smaller bin
                if finalBinSize > tolerance
                    st = arrayfun(@(x,y) x{1}(x{1} >= (y-5) & x{1} < (y+5)), spiketimes(iUnit), finalT, 'UniformOutput',false);
                    finalConv = arrayfun(@(x,y) (1/sqrt(2*pi*finalBinSigma.^2))*exp(-((y-x{1}).^2)/(2*finalBinSigma.^2)), st, finalT, 'UniformOutput', false);
                    conv = [conv, finalConv];
                end
            
                % Add up influence of all spikes for each timestamp
                popV(iTrial, iUnit,:) = cellfun(@(x) sum(x), conv);
                 end
            
            end
        
        end

        
        % Get event-triggered activity
        function firingRate = getEventTriggeredFiring(obj, eventTimes, spiketimes, timeBefore, timeAfter)
            
            % Calculates peri-event firing rates equivalent to 50ms bins
            % 
            % Inputs
            %   - eventTimes:   vector of timestamps around which to calculate firing rates (e.g. reward delivery), in seconds
            %   - spiketimes:    cell array of spike times, one cell per neuron
            %   - timeBefore:   double indicating offset in seconds before eventtimes from which to calculate firing rates, 
            %                         e.g. 1 to get firing rates starting 1 second before
            %   - timeAfter:     double indicating offset in seconds after eventtimes from which to calculate firing rates, 
            %                         e.g. 1 to get firing rates until 1 second after
            % 
            % Outputs
            %   - firingRate:    array of convolved firing rates of size T x N x B where T is length of eventtimes, 
            %                        N is length of spiketimes, and B is number of 50 ms bins between timeBefore and timeAfter
            
            start = eventTimes(:) - timeBefore;
            duration = timeBefore + timeAfter;
            nEvents = length(eventTimes);
            nCells = length(spiketimes);
            nBins = duration / obj.BinSize;
            if length(eventTimes) == 0
                firingRate = double.empty(0, nCells, nBins);
            elseif length(spiketimes) == 0
                firingRate = double.empty(0, nCells, nBins);
            else
                firingRate = obj.convolveAndBin(start, duration, obj.BinSize, spiketimes, obj.analysisConfig.gaussianWindowWidth);
                assert(size(firingRate, 1)==nEvents)
                assert(size(firingRate, 2)==nCells)
                assert(size(firingRate, 3)==nBins)
            end
        end
        
        
        % Get event-triggered running speed
        function [speed, ts] = getEventTriggeredRunningSpeed(~, eventTimes, speedData, timestamps, timeBefore, timeAfter)
            
            % Gets peri-event running speed
            % 
            % Inputs
            %   - eventTimes:  vector of timestamps around which to calculate running speed (e.g. reward delivery), in seconds
            %   - speedData:   vector of running speed values
            %   - timestamps:  vector of timestamps in seconds, corresponding to speedData, with same length
            %   - timeBefore:   double indicating offset in seconds before eventtimes from which to get running speed, 
            %                         e.g. 1 to get firing rates starting 1 second before
            %   - timeAfter:     double indicating offset in seconds after eventtimes from which to get running speed, 
            %                         e.g. 1 to get firing rates until 1 second after
            % 
            % Outputs
            %   - speed:          vector of running speed values
            %   - ts:                vector of timestamps corresponding to the speed values
            
            % Check inputs
            assert(length(speedData)==length(timestamps), 'Number of timestamps does not match length of running speed data');
            
            % Convert timestamps to miliseconds for better precision
            eventTimes = int64(eventTimes*1000);
            timestamps = int64(timestamps*1000);
            timeBefore = int64(timeBefore*1000);
            timeAfter = int64(timeAfter*1000);
            
            speed = [];
            for iEvent = 1:length(eventTimes)
                inWindow = (timestamps >= (eventTimes(iEvent) - timeBefore)) & (timestamps < (eventTimes(iEvent) + timeAfter));
                if iEvent > 1 && size(speed, 2) ~= sum(inWindow)
                    error('Dimensions of timestamps do not match between event times. Try splining to get evenly sampled timestamps, or padding with NaN values before the first even or after the last event.')
                else
                    speed = [speed; speedData(inWindow)];
                end
            end
            
            % Convert timestamps to seconds
            if ~isempty(eventTimes)
                ts = double(timestamps(inWindow) - eventTimes(end)) / 1000;
            else
                ts = [];
            end
            
        end

        
        % Get event-triggered coactivity of two spike trains
        function coactivity = getEventTriggeredCoactivity(obj, eventTimes, spiketimes1, spiketimes2, timeBefore, timeAfter, gaussianWindowWidth)
            
            % Calculates peri-event coactivity between firing rates of two neurons equivalent to 50ms bins
            % 
            % Inputs
            %   - eventTimes:                    vector of timestamps around which to calculate firing rates (e.g. reward delivery), in seconds
            %   - spiketimes1:                    cell of spike times for one neuron
            %   - spiketimes2:                    cell of spike times for a second neuron
            %   - timeBefore:                     double indicating offset in seconds before eventtimes from which to calculate firing rates, 
            %                                           e.g. 1 to get firing rates starting 1 second before
            %   - timeAfter:                       double indicating offset in seconds after eventtimes from which to calculate firing rates, 
            %                                          e.g. 1 to get firing rates until 1 second after
            %   - gaussianWindowWidth:    used for smoothing firing rates
            % 
            % Outputs
            %   - coactivity:                      array of coactivity values, calculated as the minimum firing rate of the two
            %                                          spike trains, of size T x 1 x B where T is length of eventtimes, 
            %                                          and B is number of 50 ms bins between timeBefore and timeAfter

            % Check inputs: cell arrays of spike times must have the same numbers of cells
            assert (size(spiketimes1, 2) == size(spiketimes2, 2))
            
            start = eventTimes - timeBefore;
            duration = timeBefore + timeAfter;
            nBins = duration / obj.BinSize;
            if isempty(eventTimes)
                coactivity = double.empty(0, 1, nBins);
            elseif (size(spiketimes1, 2) == 0) || (size(spiketimes2, 2) == 0)
                coactivity = double.empty(length(eventTimes), 0, nBins);
            else
                firingRate1 = obj.convolveAndBin(start, duration, obj.BinSize, spiketimes1, gaussianWindowWidth);
                firingRate2 = obj.convolveAndBin(start, duration, obj.BinSize, spiketimes2, gaussianWindowWidth);
                coactivity = min(firingRate1, firingRate2);
            end
        end
                
        
        % Z-score
        function firingRate = zScore(~, firingRate, reference)
            
            % Returns firing rates z-scored with respect to some array of reference data
            % 
            % Inputs
            %   - firingRate:  array of firing rates. Can be of size N x B where N is number of neurons and B is number of time bins,
            %                      or T x N x B where T is number of behavioural timestamps (e.g. trials)
            %   - reference:   array of firing rates from which to calculate the mean and standard deviation for z-scoring.
            %                      Can be of size 1 x B where B is number of time bins, or N x B where N is number of neurons
            % 
            % Outputs
            %   - firingRate:  z-scored firing rate, of same size as input firingRate
            
            % Calculate mean and standard deviation of reference data
            if (ndims(firingRate)~=ndims(reference)) || (mean(size(firingRate)==size(reference)) < 1)
                % Calculate mu and sigma for each cell and repeat for all
                % trials corresponding to the cell
                mu = nanmean(reference, 2);
                sigma = nanstd(reference, 0, 2);
                if ndims(firingRate)==2
                    mu = repmat(mu', size(firingRate, 1), size(firingRate, 2));
                    sigma = repmat(sigma', size(firingRate, 1), size(firingRate, 2));
                elseif ndims(firingRate)==3
                    mu = repmat(mu', size(firingRate, 1), 1, size(firingRate, 3));
                    sigma = repmat(sigma', size(firingRate, 1), 1, size(firingRate, 3));
                end
            elseif min(size(firingRate))==1 && min(size(reference))==1
                mu = nanmean(reference);
                sigma = nanstd(reference);
            elseif (ndims(firingRate)==ndims(reference)) || (mean(size(firingRate)==size(reference)) == 1)
                % Calculate mu and sigma for each time window
                mu = nanmean(reference, 3);
                sigma = nanstd(reference, 0, 3);
                mu = repmat(mu, 1, 1, size(firingRate, 3));
                sigma = repmat(sigma, 1, 1, size(firingRate, 3));
            end
            
            % Z-score
            firingRate = (firingRate - mu) ./ sigma;
            
            % Cap at 100
            firingRate = min(firingRate, 100);
                        
        end
        
    end
    
end
