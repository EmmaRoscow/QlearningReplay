
classdef TestUtilities
    
    methods(Static)

        % Mock spiketimes where only first and third cells are responsive
        function [spiketimes, CPEntryTimes, rewardArrivalTimes, CPExitTimes] = createRewardResponsiveData()
            % Cell 1 fires 1 spike in 15 bins before CP entry, and 5 spikes
            % in 20 bins before reward arrival
            % Cell 2 fires 6 spikes in 15 bins before CP entry, and 0 spikes
            % around reward arrival
            % Cell 3 fires 1 spike in 15 bins before CP entry, and 5 spikes
            % in 20 bins after reward arrival
            spiketimes = {...
                [0.3, 1.9, 4.9, 7.9, 10.9, 13.9],...
                [0.5, 3.6, 6.7, 9.8, 12.8, 12.9],...
                [0.9, 2.7, 5.7, 8.7, 11.7, 14.7]...
                };
            CPEntryTimes = [1.0, 4.0, 7.0, 10.0, 13.0];
            rewardArrivalTimes = CPEntryTimes + 1;
            CPExitTimes = CPEntryTimes + 0.5;
        end

        % Create and save mock data
        function saveMockEphysData(varargin)

            p = inputParser;
            addRequired(p, 'directory', @(x) ischar(x) || isstring(x));
            addParameter(p, 'spiketimes', {}, @(x) iscell(x));
            addParameter(p, 'CPentryTimes', [], @(x) ismatrix(x));
            addParameter(p, 'RarrivalTimes', [], @(x) ismatrix(x));
            addParameter(p, 'CPexitTimes', [], @(x) ismatrix(x));
            addParameter(p, 'nNAcUnits', 0, @(x) isnumeric(x));
            addParameter(p, 'binnedfr', [], @(x) ismatrix(x));
            addParameter(p, 'startEndTimes', [], @(x) ismatrix(x));
            addParameter(p, 'runningSpeed', struct('speed', [], 'timestamps', []), @(x) isstruct(x) && isfield(x, 'speed') && isfield(x, 'timestamps'));
            addParameter(p, 'ripples3std', struct('rippleStart', [], 'rippleStop', [], 'ripplePeak', []), @(x) isstruct(x) && isfield(x, 'rippleStart') && isfield(x, 'rippleStop') && isfield(x, 'ripplePeak'));

            parse(p, varargin{:});
            dir = p.Results.directory;
            dataFiles = rmfield(p.Results, 'directory');
            mkdir(dir);

            variableNames = fieldnames(dataFiles);
            nVariables = numel(variableNames);
            for file = 1:nVariables
                field = variableNames{file};
                filename = [field '.mat'];
                if isstruct(dataFiles.(field))
                    subStruct = dataFiles.(field);
                    save(fullfile(dir, filename), '-struct', 'subStruct', '-v7.3');
                else
                    save(fullfile(dir, filename), '-struct', 'dataFiles', sprintf(field), '-v7.3');
                end
            end

        end

        % Create and save mock behavioural data
        function saveMockBehaviouralData(varargin)

            p = inputParser;
            addRequired(p, 'directory', @(x) ischar(x) || isstring(x));
            addRequired(p, 'ratName', @(x) ischar(x) || isstring(x));
            addParameter(p, 'n_trials', [], @(x) ismatrix(x));
            addParameter(p, 'actions', {}, @(x) iscell(x));
            addParameter(p, 'rewarded', {}, @(x) iscell(x));
            addParameter(p, 'arm_values', {}, @(x) iscell(x));
            addParameter(p, 'reward_probs', {}, @(x) iscell(x));
            addParameter(p, 'proportion_optimal', [], @(x) ismatrix(x));

            parse(p, varargin{:});
            dir = p.Results.directory;
            results = rmfield(p.Results, 'directory');
            filename = ['behaviour_' lower(results.ratName) '.mat'];
            dataFiles = struct(...
                results.ratName, rmfield(results, 'ratName'));
            mkdir(dir);

            save(fullfile(dir, filename), '-struct', 'dataFiles', sprintf(results.ratName), '-v7.3');

        end
        
       % Verify that at least one of multiple plotted horizontal lines represents a particular value (such as a median) in the right colour
       function [valueFound, matchesExpectedColour] = findPlottedMedian(lines, expectedValue, expectedColour)
            
            valueFound = false;
            matchesExpectedColour = false;
            i = 0;
            while ~valueFound && i < numel(lines)
                i = i + 1;
                lineData = lines(i).YData;
                if lineData == [expectedValue, expectedValue];
                    valueFound = true;
                end
            end
            if isnan(expectedColour)
                matchesExpectedColour = true;
            else
                medianLineColor = lines(i).Color;
                if medianLineColor==expectedColour
                    matchesExpectedColour = true;
                end
            end

        end

        
    end

end
