
classdef ExplainedVarianceReactivationTest < matlab.unittest.TestCase
    
    properties
        ExplainedVarianceReactivationObj
        TestDataDir
    end
    
    
    methods(Access = private)
        
        function [binnedfrPRE, binnedfrPOST, binnedfrTASK, EV, REV, preCC, taskCC, postCC] = createBinnedFR(~, varargin)

                function [binnedfr, CC] = createFiringRates(rho, seed)

                    covariance_matrix = [...
                        1, rho(1), rho(2);
                        rho(1), 1, rho(3);
                        rho(2), rho(3), 1];
                    cc(1) = inf; cc(2) = inf; cc(3) = inf;
                    tolerance = 0.01;
                    % In case the provided random seed doesn't work well, search over others
                    i = seed;
                    while cc(1) < (rho(1)-tolerance) || cc(1) > (rho(1)+tolerance) || cc(2) < (rho(2)-tolerance) || cc(2) > (rho(2)+tolerance) || cc(3) < (rho(3)-tolerance) || cc(3) > (rho(3)+tolerance)
                        rng(i)
                        binnedfr = mvnrnd(zeros(1, 3), covariance_matrix, 10)';
                        CC = corrcoef(binnedfr');
                        cc = [CC(1, 2), CC(1, 3), CC(2, 3)];
                        i = i + 1;
                    end

                end
                
            % Check if correlations are between areas
            if (nargin > 1) && ischar(varargin{1}) && strcmp(varargin, 'inter-area')
                interarea = true;
            else
                interarea = false;
            end

            % PRE activity for 3 units, with normally distributed firing rates, and correlations of -0.1, 0 and 0.2 respectively
            rho = [-0.1, 0, 0.2];
            [binnedfrPRE, preCC] = createFiringRates(rho, 77299);

            % TASK activity for 3 units, with normally distributed firing rates, and correlations of -0.2, 0.3 and 0.5 respectively
            rho = [-0.2, 0.3, 0.5];
            [binnedfrTASK, taskCC] = createFiringRates(rho, 9726);

            % POST activity for 3 units, with normally distributed firing rates, and correlations of -0.2, 0.1 and 0.6 respectively
            rho = [-0.2, 0.1, 0.6];
            [binnedfrPOST, postCC] = createFiringRates(rho,  38326);
            
            % If inter-area, double these and take first 2 cells as accumbens cells
            if interarea
                binnedfrPRE = repmat(binnedfrPRE, 2, 1);
                binnedfrPOST = repmat(binnedfrPOST, 2, 1);
                binnedfrTASK = repmat(binnedfrTASK, 2, 1);
                preCC = repmat(preCC, 2, 2);
                taskCC = repmat(taskCC, 2, 2);
                postCC = repmat(postCC, 2, 2);
                preCC = preCC(1:2, 3:6);
                taskCC = taskCC(1:2, 3:6);
                postCC = postCC(1:2, 3:6);
            end

            % Calculate EV and REV
            prePostCC = corr2(preCC, postCC);
            taskPreCC = corr2(taskCC, preCC);
            taskPostCC = corr2(taskCC, postCC);
            EV = ( (taskPostCC - taskPreCC * prePostCC) / sqrt( (1 - taskPreCC.^2) * (1 - prePostCC.^2) ) ).^2;
            REV = ( (taskPreCC - taskPostCC * prePostCC) / sqrt( (1 - taskPostCC.^2) * (1 - prePostCC.^2) ) ).^2;

        end
        
        function [spiketimes, binnedfr, startEndTimes, EV, REV, rippleStart, rippleStop] = createSpiketimes(obj, varargin)
                        
            % Process variables
            p = inputParser;
            addParameter(p, 'nNAcUnits', 3, @(x) isnumeric(x) & floor(x)==x);
            addParameter(p, 'nCA1Units', 0, @(x) isnumeric(x) & floor(x)==x);
            addParameter(p, 'interarea', false, @(x) islogical(x));
            addParameter(p, 'ripples', false, @(x)  islogical(x));
            addParameter(p, 'seed', 0, @(x) isnumeric(x));
            parse(p, varargin{:});
            nNAcUnits = p.Results.nNAcUnits;
            nCA1Units = p.Results.nCA1Units;
            interarea = p.Results.interarea;
            ripples = p.Results.ripples;
            seed = p.Results.seed;
            
            if ~interarea
                assert(nNAcUnits==0 || nCA1Units==0, 'For intraarea analysis, one of nNAcUnits or nCA1Units must be 0. Please revise input arguments.')
            end
            
            function spiketimes = spikes(meanFR, nNeurons, maxDuration)
                
                nSpikes = 10000;
                goodCoverage = false;
                idx = [];
                for i = 1:length(nNeurons)
                    idx = [idx 0:nNeurons(i)-1];
                end
                while ~goodCoverage
                    % Create exponentially distributed random spiketimes,
                    % with interspike intervals defined by mean firing rates
                    spiketimes = arrayfun(@(x) cumsum(random('exponential', 1/meanFR(rem(x, length(meanFR))+1), 1, nSpikes)), idx, 'UniformOutput', false);
                    % Ensure there are enough spikes to cover the full desired duration
                    if all(cellfun(@max, spiketimes) >= maxDuration)
                        goodCoverage = true;
                    end
                    nSpikes = nSpikes * 10;
                end
                
                % Trim the spikes to the desired duration
                spiketimes = cellfun(@(x) x(x<= maxDuration), spiketimes, 'UniformOutput', false);
                
            end
                
            
            % Create exponentially distributed spiketimes for 3 neurons with
            % different, constant mean firing rates lasting 300 seconds
            rng(seed);
            startEndTimes = [0, 100; 100, 200; 200, 300];
            maxDuration = startEndTimes(3, 2);
            meanFR = [5, 10, 20, 50, 100];
            rippleMeanFR = flip(meanFR);
            spiketimes_NAc = spikes(meanFR, nNAcUnits, maxDuration);
            spiketimes_CA1 = spikes(meanFR, nCA1Units, maxDuration);
            spiketimes = [spiketimes_NAc, spiketimes_CA1];
            
            if ripples
                % Create some ripple times lasting 10 seconds
                rippleDuration = 10;
                rippleStartPRE = [20, 40, 60];
                rippleStopPRE = rippleStartPRE + rippleDuration;
                rippleStartPOST = [220, 240, 260];
                rippleStopPOST = rippleStartPOST + rippleDuration;
                rippleStart = [rippleStartPRE rippleStartPOST];
                rippleStop = [rippleStopPRE, rippleStopPOST];
                % Create spiketimes during ripples with different mean firing rates from the rest
                spiketimesRipples = spikes(rippleMeanFR, [nNAcUnits, nCA1Units], maxDuration);
                % Replace spiking during ripples in PRE with the new ripple spiketimes
                for iRip = 1:length(rippleStartPRE)
                    % Remove original spiketimes from ripple period
                    spiketimes = cellfun(@(x) x( x < rippleStartPRE(iRip) |  x > rippleStopPRE(iRip) ), spiketimes, 'UniformOutput', false);
                    % Combine with the ripple spiking
                    spiketimes = cellfun(@(x, y) [x y], spiketimes, spiketimesRipples, 'UniformOutput', false);
                    spiketimes = cellfun(@sort, spiketimes, 'UniformOutput', false);
                end
            else
                rippleStart = [];
                rippleStop = [];
            end
                            
            % Bin the firing rates using binless method
            binnedfrPRE = squeeze(obj.ExplainedVarianceReactivationObj.convolveAndBin(startEndTimes(1, 1), diff(startEndTimes(1, :)), obj.ExplainedVarianceReactivationObj.BinSize, spiketimes));
            binnedfrTASK = squeeze(obj.ExplainedVarianceReactivationObj.convolveAndBin(startEndTimes(2, 1), diff(startEndTimes(2, :)), obj.ExplainedVarianceReactivationObj.BinSize, spiketimes));
            binnedfrPOST = squeeze(obj.ExplainedVarianceReactivationObj.convolveAndBin(startEndTimes(3, 1), diff(startEndTimes(3, :)), obj.ExplainedVarianceReactivationObj.BinSize, spiketimes));
            binnedfr = [binnedfrPRE binnedfrTASK binnedfrPOST];
            
            % Create correlation matrices between them
            if interarea
                preCC = zeros(nNAcUnits, nCA1Units);
                taskCC = zeros(nNAcUnits, nCA1Units);
                postCC = zeros(nNAcUnits, nCA1Units);
                for iNAc = 1:nNAcUnits
                    for iHPC = 1:nCA1Units
                        indexHPC = iHPC + nNAcUnits; 
                        c = corrcoef(binnedfrPRE(iNAc,:), binnedfrPRE(indexHPC,:));
                        preCC(iNAc, iHPC) = c(2);
                        c = corrcoef(binnedfrTASK(iNAc,:), binnedfrTASK(indexHPC,:));
                        taskCC(iNAc, iHPC) = c(2);
                        c = corrcoef(binnedfrPOST(iNAc,:), binnedfrPOST(indexHPC,:));
                        postCC(iNAc, iHPC) = c(2);
                    end
                end

            else
                preCC = corrcoef(binnedfrPRE');
                taskCC = corrcoef(binnedfrTASK');
                postCC = corrcoef(binnedfrPOST');
            end
            taskPreCC = corr2(taskCC,preCC);
            taskPostCC = corr2(taskCC,postCC);
            prePostCC = corr2(preCC,postCC);

            % Calculate EV and REV
            EV = ( (taskPostCC - taskPreCC * prePostCC) / sqrt( (1 - taskPreCC.^2) * (1 - prePostCC.^2) ) ).^2;
            REV = ( (taskPreCC - taskPostCC * prePostCC) / sqrt( (1 - taskPostCC.^2) * (1 - prePostCC.^2) ) ).^2;
            
        end
        
        function [nNAcUnits] = getIntraAreaNAc(~, relevantNeurons, iFunc, nRelevantNeurons, nIrrelevantNeurons)
            if strcmp(relevantNeurons{iFunc}, 'nNAcUnits')
                nNAcUnits = nRelevantNeurons;
            elseif strcmp(relevantNeurons{iFunc}, 'nCA1Units')
                nNAcUnits = nIrrelevantNeurons;
            else
                error('Unexpected value of relevantNeurons: %s', relevantNeurons{iFun})
            end
        end
        
        function [] = findPlottedMedian(testCase, lines, expectedMedian, expectedColour)
            
            [medianFound, matchesExpectedColour] = TestUtilities.findPlottedMedian(lines, expectedMedian, expectedColour);
            testCase.verifyTrue(medianFound, 'Figure does not contain a line indicating the median of the data.');
            testCase.verifyTrue(matchesExpectedColour, 'Figure contains a line indicating the median of the data, but it is not the expected colour.');
            
            
        end
                    
    end
    
    
    methods(TestMethodSetup)
        
        % Initialise ExplainedVarianceReactivation object with temporary data path
        function createExplainedVarianceReactivation(testCase)
            testCase.ExplainedVarianceReactivationObj = ExplainedVarianceReactivation();
            % Create a temporary directory for test data
            testCase.TestDataDir = tempname;
            mkdir(testCase.TestDataDir);
        end
        
        % Suppress warnings during testing to not overload the command window
        function disableWarnings(testCase)
            warning('off','all')
        end
        

    end
    
    
    methods(TestMethodTeardown)
        
        % Remove the temporary test data directory after tests
        function removeTestData(testCase)
            if exist(testCase.TestDataDir, 'dir') == 7
                rmdir(testCase.TestDataDir, 's');
            end
        end
        
        % Close figures
        function closeFigures(testCase)
            close('all', 'force')
        end
        
        % Reenable warnings
        function enableWarnings(testCase)
            warning('off','all')
        end
        
    end

    
    methods (Test)
        
        % ========================= initiating ========================= %
        
        % Test that the required data subdirectories are on the path
        function checkDataDirectories(testCase)
            
            requiredPaths = {...
                fullfile(pwd, '.', 'data'),...
                fullfile(pwd, '.', 'data', 'behavioural_data'),...
                fullfile(pwd, '.', 'data', 'ephys_data')};
            
            for requiredPath = requiredPaths
                testCase.verifyTrue(~isempty(strfind(path, requiredPath{1})),...
                sprintf('Required directory "%s" is not on the MATLAB path. Make sure the current directory is the root', requiredPath{1}));
            end
            
        end
        
        % Test that the required data subdirectories are not empty
        function checkDataSubdirectories(testCase)
            
            % Behavioural data
            directory = fullfile('.', 'data', 'behavioural_data', '**');
            testCase.verifyFalse(isempty(dir(directory)),...
                sprintf('Behavioural data directory contains no data files. Add data or refer to the README'))
            
            % Ephys data
            directory = fullfile('.', 'data', 'ephys_data', '**');
            testCase.verifyFalse(isempty(dir(directory)),...
                sprintf('Ephys data directory contains no data files. Add data or refer to the README'))
            ephysSubDirectories = dir('data/ephys_data/**');
            ephysSubDirectories = {ephysSubDirectories.name};
            for subdir = ephysSubDirectories(3:end)
                directory = fullfile('.', 'data', 'ephys_data', subdir{1});
                testCase.verifyFalse(isempty(dir(directory)),...
                    sprintf('Directory "%s" contains no data files. Add data or refer to the README', directory))
            end
            
        end
        
        
        % Test that the main fields for the inclusion criteria are
        % created, even if no inclusion criteria are passed
        function createsInclusionCriteria(testCase)
            obj = ExplainedVarianceReactivation();
            testCase.verifyTrue(isstruct(obj.inclusionCriteria));
            expected_fields = {'rats', 'sessions', 'cells'};
            for iField = 1:length(expected_fields)
                testCase.verifyTrue(isfield(obj.inclusionCriteria, expected_fields{iField})),
            end
        end
        
        % Test that the subfields (for the second field, sessions) are created,
        % even if no inclusion criteria are passed
        function createsSessionInclusionCriteria(testCase)
            obj = ExplainedVarianceReactivation();
            expectedFields = {'initialLearningOnly', 'skipEarliest', 'significantPerformanceThreshold', 'minimumNCA1Cells', 'minimumNVStrCells'};
            for iField = 1:length(expectedFields)
                testCase.verifyTrue(isfield(obj.inclusionCriteria.sessions, expectedFields{iField}), sprintf('Expected field "%s" is missing from inclusion criteria.', expectedFields{iField}));
            end
        end
        
        % Test that if non-default inclusion criterion is passed, it gets
        % included correctly (for the first field, rats)
        function includes_nondefault_rats(testCase)
            obj = ExplainedVarianceReactivation(struct('rats', {{'Severus'}}));
            testCase.verifyEqual(obj.inclusionCriteria.rats{1}, 'Severus');
        end
        
        % Test that if non-default inclusion criterion is passed, it gets
        % included correctly (for the first subfield, sessions.initialLearningOnly)
        function includesInitialLearning_true(testCase)
            obj = ExplainedVarianceReactivation(struct('sessions', struct('initialLearningOnly', true)));
            expectedFields = {'initialLearningOnly', 'skipEarliest', 'significantPerformanceThreshold', 'minimumNCA1Cells', 'minimumNVStrCells'};
            for iField = 1:length(expectedFields)
                testCase.verifyTrue(isfield(obj.inclusionCriteria.sessions, expectedFields{iField}), sprintf('Expected field "%s" is missing from inclusion criteria.', expectedFields{iField}));
            end
            testCase.verifyTrue(obj.inclusionCriteria.sessions.initialLearningOnly);
        end
        
        % Test that if non-default inclusion criterion is passed, it gets
        % included correctly (for the first subfield, sessions.initialLearningOnly)
        function includesInitialLearning_false(testCase)
            obj = ExplainedVarianceReactivation(struct('sessions', struct('initialLearningOnly', false)));
            expectedFields = {'initialLearningOnly', 'skipEarliest', 'significantPerformanceThreshold', 'minimumNCA1Cells', 'minimumNVStrCells'};
            for iField = 1:length(expectedFields)
                testCase.verifyTrue(isfield(obj.inclusionCriteria.sessions, expectedFields{iField}), sprintf('Expected field "%s" is missing from inclusion criteria.', expectedFields{iField}));
            end
            testCase.verifyFalse(obj.inclusionCriteria.sessions.initialLearningOnly);
        end
        
        
        % ===================== selecting sessions ===================== %
        
        % Test that the initial learning sessions contain the right number
        % of sessions if 'skipEarliest' parameter is 'true'
        function testSelectSessions_skipsInitialSessions(testCase)
            
            % Create mock ephys  data
            for session = 1:20
                dataPath = fullfile(testCase.TestDataDir, 'Rat_Q', sprintf('Session%s', num2str(session)));
                TestUtilities.saveMockEphysData(dataPath, 'nNAcUnits', 5, 'spiketimes', cell(1, 10));
            end

            % Specify inclusion criteria (one rat, skip initial true)
            inclusionCriteria = struct(...
                'rats', {{'Quirinius'}},...
                'sessions', struct(...
                    'initialLearningOnly', true,...
                    'skipEarliest', true,...
                    'significantPerformanceThreshold', 0));
            
            testCase.ExplainedVarianceReactivationObj = ExplainedVarianceReactivation(inclusionCriteria);
            output = testCase.ExplainedVarianceReactivationObj.selectSessions('ephysDataPath', testCase.TestDataDir);
            testCase.verifyEqual(output.sessionsForSpecificAnalysis.Quirinius, [3:12]);
            
        end
        
         % Test that the initial learning sessions contain the right number
        % of sessions if 'skipEarliest' parameter is 'false'
        function testSelectSessions_includesAllInitialSessions(testCase)
            
            % Create mock ephys  data
            for session = 1:20
                ephysDataPath = fullfile(testCase.TestDataDir, 'Rat_Q', sprintf('Session%s', num2str(session)));
                TestUtilities.saveMockEphysData(ephysDataPath, 'nNAcUnits', 5, 'spiketimes', cell(1, 10));
            end

            % Specify inclusion criteria (one rat, skip initial false)
            inclusionCriteria = struct(...
                'rats', {{'Quirinius'}},...
                'sessions', struct(...
                    'initialLearningOnly', true,...
                    'skipEarliest', false,...
                    'significantPerformanceThreshold', 0));
            
            testCase.ExplainedVarianceReactivationObj = ExplainedVarianceReactivation(inclusionCriteria);
            output = testCase.ExplainedVarianceReactivationObj.selectSessions('ephysDataPath', testCase.TestDataDir);
            testCase.verifyEqual(output.sessionsForSpecificAnalysis.Quirinius, [1:12]);
        end
        
        % Test that the initial learning sessions contain the right number
        % of sessions if 'significantPerformanceTthreshold' parameter is
        % higher than 0
        function testSelectSessions_skipsLowPerformanceSessions(testCase)
                        
            % Create mock ephys  data (20 sessions with 5 hippocampal and 5 striatal neurons)
            for session = 1:20
                dataPath = fullfile(testCase.TestDataDir, 'Rat_Q', sprintf('Session%s', num2str(session)));
                TestUtilities.saveMockEphysData(dataPath, 'nNAcUnits', 5, 'spiketimes', cell(1, 10));
            end
            
            % Create mock behavioural data (1000 trials, first 4 and last 4 sessions below 1/3 threshold)
            behaviouralDataPath = fullfile(testCase.TestDataDir, 'test_session'); 
            proportion_optimal = [0:0.1:1.0, 0.8:-0.1:0.0];
            n_trials = ones(1, 20) * 1000;
            TestUtilities.saveMockBehaviouralData(behaviouralDataPath, 'Quirinius', 'proportion_optimal', proportion_optimal, 'n_trials', n_trials);

            % Specify inclusion criteria (one rat, skip initial false, significance threshold 1/3)
            inclusionCriteria = struct(...
                'rats', {{'Quirinius'}},...
                'sessions', struct(...
                    'initialLearningOnly', true,...
                    'skipEarliest', false,...
                    'significantPerformanceThreshold', 1/3));
            
            testCase.ExplainedVarianceReactivationObj = ExplainedVarianceReactivation(inclusionCriteria);
            output = testCase.ExplainedVarianceReactivationObj.selectSessions('behaviouralDataPath', behaviouralDataPath, 'ephysDataPath', testCase.TestDataDir);
            testCase.verifyEqual(output.sessionsForSpecificAnalysis.Quirinius, find(proportion_optimal(1:12) > 1/3));
        end

        
        
        % ======================= selecting cells ======================= %
        
        % % Test that all cells are included when inclusionCriteria.rewardResponsive is false
        function testSelectRewardResponsiveCells_allCellsIncluded(testCase)
                        
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.cells.rewardResponsive = false;
            [spiketimes, CPEntryTimes, rewardArrivalTimes, ~] = TestUtilities.createRewardResponsiveData();         
            expected = boolean(ones(1, 3));
            actual = testCase.ExplainedVarianceReactivationObj.selectRewardResponsiveCells(spiketimes, rewardArrivalTimes, CPEntryTimes);
            testCase.verifyEqual(actual, expected, 'All cells should be included when rewardResponsive is false');

        end
        
        % Test that only reward-responsive cells are included when
        % inclusionCriteria.rewardResponsive is true
        function testSelectRewardResponsiveCells_someCellsIncluded(testCase)
            
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.cells.rewardResponsive = true;
            [spiketimes, CPEntryTimes, rewardArrivalTimes, ~] = TestUtilities.createRewardResponsiveData();
            expected = [true, false, true];
            actual = testCase.ExplainedVarianceReactivationObj.selectRewardResponsiveCells(spiketimes, rewardArrivalTimes, CPEntryTimes);
            testCase.verifyEqual(actual, expected, 'Only first and third cells should be responsive');
            
        end
        
        
        % ======================= analysing firing ======================= %
        
        % Test that data is correctly z-scored when referenced to itself
        function test_zscoring_autoreferenced(testCase)
            
            % Dummy data: 100 trials, 10 cells, 20 timestamps
            firing = random('Normal', 5, 2, 100, 10, 20);
            mu = mean(firing, 3);
            mu = repmat(mu, 1, 1, 20);
            sigma = std(firing, [], 3);
            sigma = repmat(sigma, 1, 1, 20);
            expected_result = (firing - mu) ./ sigma;
            
            obj = ExplainedVarianceReactivation();
            zscored_firing = obj.zScore(firing, firing);
            testCase.verifyEqual(expected_result, zscored_firing);
            
        end
        
        
        % Test that function raises an error if reference data is the wrong
        % size
        function test_zscoring_fails_wrong_dimensions(testCase)
        
            % Dummy firing data: 100 trials, 10 cells, 20 timestamps
            firing = random('Normal', 5, 2, 100, 10, 20);
            
            % Dummy reference data: 1 trial, 5 cells, 50 timestamps
            reference = random('Normal', 5, 2, 1, 5, 50);
            obj = ExplainedVarianceReactivation();
            errorCondition = @() obj.zScore(firing, reference);
            testCase.verifyError(errorCondition, ?MException);
            
        end
        
        % Test that data is correctly z-scored when referenced to other
        % data
        function test_zscoring_referenced(testCase)
            
            % Dummy firing data: 100 trials, 10 cells, 20 timestamps
            firing = random('Normal', 5, 2, 100, 10, 20);
            
            % Dummy reference data: 10 cells, 50 timestamps
            reference = random('Normal', 5, 2, 10, 50);
            mu = mean(reference, 2)';
            mu = repmat(mu, 100, 1, 20);
            sigma = std(reference, [], 2)';
            sigma = repmat(sigma, 100, 1, 20);
            expected_result = (firing - mu) ./ sigma;
            
            obj = ExplainedVarianceReactivation();
            zscored_firing = obj.zScore(firing, reference);
            testCase.verifyEqual(expected_result, zscored_firing);
            
        end
        
        
        % ======================== binning activity ====================== %
        
        % Test basic functionality of binActivity
        function testBinActivity_basicFunctionality(testCase)
            
            % Mock data: PRE, TASK, and POST lasting 10 seconds each;
            % 3 cells with different but constant firing rates
            startEndTimes = [0, 10; 10, 20; 20, 30];
            meanFR = [10, 20, 50];
            spiketimes = {{...
                0 : 1/meanFR(1) : 30,...
                0 : 1/meanFR(2) : 30,...
                0 : 1/meanFR(3) : 30}};
            
            % Run EV-REV analysis
            data = struct(...
                'spiketimes', spiketimes);
            timestamps = struct('startEndTimes', startEndTimes);
            analysisConfig = struct(...
                'postDuration', NaN,...
                'rippleActivityMode', 'none');
            
            % Call binActivity
            [binnedfrPRE, binnedfrPOST, binnedfrTASK] = testCase.ExplainedVarianceReactivationObj.binActivity(data, timestamps, analysisConfig);
            expectedNNeurons = 3;
            expectedNBins = 10 / testCase.ExplainedVarianceReactivationObj.BinSize;
            testCase.verifyEqual(size(binnedfrPRE), [expectedNNeurons, expectedNBins], 'Output binnedfrPRE does not contain the expected number of bins')
            testCase.verifyEqual(size(binnedfrPOST), [expectedNNeurons, expectedNBins], 'Output binnedfrPOST does not contain the expected number of bins')
            testCase.verifyEqual(size(binnedfrTASK), [expectedNNeurons, expectedNBins], 'Output binnedfrTASK does not contain the expected number of bins')
            
        end
        
        
        % Test binActivity's processing of equalised time around ripples
        function testBinActivity_equalisedRipples(testCase)
            
            % Mock data: PRE, TASK, and POST lasting 10 seconds each;
            % 4 ripples of which 2 during PRE and 2 during POST, with duration 0;
            % 3 cells with different but constant firing rates with spiketimes only around ripples
            startEndTimes = [0, 10; 10, 20; 20, 30];
            ripplePeak = [2, 5, 22, 25];
            rippleStart = ripplePeak;
            rippleStop = ripplePeak;
            meanFR = [10, 20, 50];
            spiketimes = cell(1, length(meanFR));
            for iRip = 1:length(ripplePeak)
                for iCell = 1:length(meanFR)
                    spiketimes{iCell} = [spiketimes{iCell} ripplePeak(iRip) : 1/meanFR(iCell) : ripplePeak(iRip) + 0.2];
                end
            end
            
            % Run EV-REV analysis
            data = struct(...
                'spiketimes', {spiketimes});
            timestamps = struct(...
                'startEndTimes', startEndTimes,...
                'ripplePeak', ripplePeak,...
                'rippleStart', rippleStart,...
                'rippleStop', rippleStop);
            analysisConfig = struct(...
                'postDuration', NaN,...
                'rippleActivityMode', 'equalise');
            
            % Call binActivity
            [binnedfrPRE, binnedfrPOST, binnedfrTASK] = testCase.ExplainedVarianceReactivationObj.binActivity(data, timestamps, analysisConfig);

            % Check that the number of bins is as expected
            expectedNNeurons = 3;
            expectedNBinsRest = 0.4 / testCase.ExplainedVarianceReactivationObj.BinSize;
            expectedNBinsTask = 10 / testCase.ExplainedVarianceReactivationObj.BinSize;
            testCase.verifyEqual(size(binnedfrPRE), [expectedNNeurons, expectedNBinsRest], 'Output binnedfrPRE does not contain the expected number of bins')
            testCase.verifyEqual(size(binnedfrPOST), [expectedNNeurons, expectedNBinsRest], 'Output binnedfrPOST does not contain the expected number of bins')
            testCase.verifyEqual(size(binnedfrTASK), [expectedNNeurons, expectedNBinsTask], 'Output binnedfrTASK does not contain the expected number of bins')
            
            % Check that the values are as expected
            testCase.verifyFalse(any(any(binnedfrPRE==0)));
            testCase.verifyFalse(any(any(binnedfrPOST==0)));
            testCase.verifyTrue(all(all(binnedfrTASK==0)));
            testCase.verifyFalse(any(any(diff(binnedfrPRE)==0)));
            testCase.verifyFalse(any(any(diff(binnedfrPOST)==0)));
            
        end

        
        % Test binActivity's processing of variable time around ripples
        function testBinActivity_variableRipples(testCase)
            
            % Mock data: PRE, TASK, and POST lasting 10 seconds each;
            % 4 ripples of which 2 during PRE and 2 during POST with varying duration;
            % 3 cells with spikes only during ripple times
            startEndTimes = [0, 10; 10, 20; 20, 30];
            rippleStart = [2, 5, 22, 25];
            rippleLength = [0.04, 0.05, 0.15, 0.22];
            rippleStop = rippleStart + rippleLength;
            ripplePeak = rippleStart;
            meanFR = [20, 50, 100];
            spiketimes = cell(1, length(meanFR));
            for iRip = 1:length(ripplePeak)
                for iCell = 1:length(meanFR)
                    spiketimes{iCell} = [spiketimes{iCell} rippleStart(iRip) : 1/meanFR(iCell) : rippleStop(iRip)];
                end
            end
            
            % Run EV-REV analysis
            data = struct(...
                'spiketimes', {spiketimes});
            timestamps = struct(...
                'startEndTimes', startEndTimes,...
                'ripplePeak', ripplePeak,...
                'rippleStart', rippleStart,...
                'rippleStop', rippleStop);
            analysisConfig = struct(...
                'postDuration', NaN,...
                'rippleActivityMode', 'variable');
            
            % Call binActivity
            [binnedfrPRE, binnedfrPOST, binnedfrTASK] = testCase.ExplainedVarianceReactivationObj.binActivity(data, timestamps, analysisConfig);

            % Check that the number of bins is as expected
            expectedNNeurons = 3;
            expectedNBinsPre = sum(ceil(rippleLength(1:2) / testCase.ExplainedVarianceReactivationObj.BinSize));
            expectedNBinsPost = sum(ceil(rippleLength(3:4) / testCase.ExplainedVarianceReactivationObj.BinSize));
            expectedNBinsTask = 10 / testCase.ExplainedVarianceReactivationObj.BinSize;
            testCase.verifyEqual(size(binnedfrPRE), [expectedNNeurons, expectedNBinsPre], 'Output binnedfrPRE does not contain the expected number of bins')
            testCase.verifyEqual(size(binnedfrPOST), [expectedNNeurons, expectedNBinsPost], 'Output binnedfrPOST does not contain the expected number of bins')
            testCase.verifyEqual(size(binnedfrTASK), [expectedNNeurons, expectedNBinsTask], 'Output binnedfrTASK does not contain the expected number of bins')
            
            % Check that the values are as expected
            testCase.verifyFalse(any(any(binnedfrPRE==0)));
            testCase.verifyFalse(any(any(binnedfrPOST==0)));
            testCase.verifyTrue(all(all(binnedfrTASK==0)));
            testCase.verifyFalse(any(any(diff(binnedfrPRE)==0)));
            testCase.verifyFalse(any(any(diff(binnedfrPOST)==0)));
            
        end
        
        % Test binActivity's processing of duration of POST-task rest, unrestricted to ripple times
        function testBinActivity_postDuration(testCase)
            
            % Mock data: 12 seconds each of PRE, TASK and POST, with a gap of 3 seconds between each;
            % 3 cells with different but constant firing rates
            startEndTimes = [0, 12; 15, 27; 30, 42];
            meanFR = [10, 20, 50];
            spiketimes = {{...
                0 : 1/meanFR(1) : 42,...
                0 : 1/meanFR(2) : 42,...
                0 : 1/meanFR(3) : 42}};
            
            % Run EV-REV analysis
            data = struct(...
                'spiketimes', spiketimes);
            timestamps = struct('startEndTimes', startEndTimes);
            analysisConfig = struct(...
                'postDuration', 5,...
                'rippleActivityMode', 'none');
            
            % Call binActivity
            [binnedfrPRE, binnedfrPOST, binnedfrTASK] = testCase.ExplainedVarianceReactivationObj.binActivity(data, timestamps, analysisConfig);
            expectedNNeurons = 3;
            expectedNBinsPre = 12 / testCase.ExplainedVarianceReactivationObj.BinSize;
            expectedNBinsPost = 5 / testCase.ExplainedVarianceReactivationObj.BinSize;
            expectedNBinsTask = 12 / testCase.ExplainedVarianceReactivationObj.BinSize;
            testCase.verifyEqual(size(binnedfrPRE), [expectedNNeurons, expectedNBinsPre], 'Output binnedfrPRE does not contain the expected number of bins')
            testCase.verifyEqual(size(binnedfrPOST), [expectedNNeurons, expectedNBinsPost], 'Output binnedfrPOST does not contain the expected number of bins')
            testCase.verifyEqual(size(binnedfrTASK), [expectedNNeurons, expectedNBinsTask], 'Output binnedfrTASK does not contain the expected number of bins')
            
        end
        
        % Test binActivity's processing of duration of POST-task rest, with equalised time around ripples
        function testBinActivity_postDurationEqualisedRipples(testCase)
            
            % Mock data: 5 ripples, of which 2 are in PRE, 1 is in TASK,
            % 2 are in the first 5 seconds of POST, and 1 is in the last 5 seconds of POST;
            % spiketimes only during ripples
            startEndTimes = [0, 10; 10, 20; 20, 30];
            ripplePeak = [2, 5, 11, 22, 26, 28];
            rippleStart = ripplePeak;
            rippleStop = ripplePeak;
            meanFR = [10, 20, 50];
            spiketimes = cell(1, length(meanFR));
            for iRip = 1:length(ripplePeak)
                for iCell = 1:length(meanFR)
                    spiketimes{iCell} = [spiketimes{iCell} ripplePeak(iRip) : 1/meanFR(iCell) : ripplePeak(iRip) + 0.2];
                end
            end
            
            % Run EV-REV analysis
            data = struct(...
                'spiketimes', {spiketimes});
            timestamps = struct(...
                'startEndTimes', startEndTimes,...
                'ripplePeak', ripplePeak,...
                'rippleStart', rippleStart,...
                'rippleStop', rippleStop);
            analysisConfig = struct(...
                'postDuration', 5,...
                'rippleActivityMode', 'equalise');
            
            % Call binActivity
            [binnedfrPRE, binnedfrPOST, binnedfrTASK] = testCase.ExplainedVarianceReactivationObj.binActivity(data, timestamps, analysisConfig);

            % Check that the number of bins is as expected
            expectedNNeurons = 3;
            expectedNBinsPre = 0.4 / testCase.ExplainedVarianceReactivationObj.BinSize;
            expectedNBinsPost = 0.2 / testCase.ExplainedVarianceReactivationObj.BinSize;
            expectedNBinsTask = 10 / testCase.ExplainedVarianceReactivationObj.BinSize;
            testCase.verifyEqual(size(binnedfrPRE), [expectedNNeurons, expectedNBinsPre], 'Output binnedfrPRE does not contain the expected number of bins')
            testCase.verifyEqual(size(binnedfrPOST), [expectedNNeurons, expectedNBinsPost], 'Output binnedfrPOST does not contain the expected number of bins')
            testCase.verifyEqual(size(binnedfrTASK), [expectedNNeurons, expectedNBinsTask], 'Output binnedfrTASK does not contain the expected number of bins')
            
            % Check that the values are as expected
            testCase.verifyFalse(any(any(binnedfrPRE==0)));
            testCase.verifyFalse(any(any(binnedfrPOST==0)));
            testCase.verifyFalse(any(any(diff(binnedfrPRE)==0)));
            testCase.verifyFalse(any(any(diff(binnedfrPOST)==0)));
            
        end
        
        % Test binActivity's processing of duration of POST-task rest, with variable time around ripples
        function testBinActivity_postDurationVariableRipples(testCase)
            
            % Mock data: 6 ripples of variable length, of which 3 are in PRE (first 10 seconds), 
            % 2 are in the first 5 seconds of POST, and 1 is in the last 5 seconds of POST; 
            % spiketimes only during ripples
            startEndTimes = [0, 10; 10, 20; 20, 30];
            rippleStart = [2, 5, 9, 22, 23.8, 25];
            rippleLength = [0.04, 0.05, 0.06, 0.15, 0.2, 0.22];
            rippleStop = rippleStart + rippleLength;
            ripplePeak = rippleStart;
            meanFR = [20, 50, 100];
            spiketimes = cell(1, length(meanFR));
            for iRip = 1:length(ripplePeak)
                for iCell = 1:length(meanFR)
                    spiketimes{iCell} = [spiketimes{iCell} rippleStart(iRip) : 1/meanFR(iCell) : rippleStop(iRip)];
                end
            end
            
            % Run EV-REV analysis
            data = struct(...
                'spiketimes', {spiketimes});
            timestamps = struct(...
                'startEndTimes', startEndTimes,...
                'ripplePeak', ripplePeak,...
                'rippleStart', rippleStart,...
                'rippleStop', rippleStop);
            analysisConfig = struct(...
                'postDuration', 4,...
                'rippleActivityMode', 'variable');
            
            % Call binActivity
            [binnedfrPRE, binnedfrPOST, binnedfrTASK] = testCase.ExplainedVarianceReactivationObj.binActivity(data, timestamps, analysisConfig);

            % Check that the number of bins is as expected
            expectedNNeurons = 3;
            expectedNBinsPre = sum(ceil(rippleLength(1:3) / testCase.ExplainedVarianceReactivationObj.BinSize));
            expectedNBinsPost = sum(ceil(rippleLength(4:5) / testCase.ExplainedVarianceReactivationObj.BinSize));
            expectedNBinsTask = 10 / testCase.ExplainedVarianceReactivationObj.BinSize;
            testCase.verifyEqual(size(binnedfrPRE), [expectedNNeurons, expectedNBinsPre], 'Output binnedfrPRE does not contain the expected number of bins')
            testCase.verifyEqual(size(binnedfrPOST), [expectedNNeurons, expectedNBinsPost], 'Output binnedfrPOST does not contain the expected number of bins')
            testCase.verifyEqual(size(binnedfrTASK), [expectedNNeurons, expectedNBinsTask], 'Output binnedfrTASK does not contain the expected number of bins')
            
            % Check that the values are as expected
            testCase.verifyFalse(any(any(binnedfrPRE==0)));
            testCase.verifyFalse(any(any(binnedfrPOST==0)));
            testCase.verifyTrue(all(all(binnedfrTASK==0)));
            testCase.verifyFalse(any(any(diff(binnedfrPRE)==0)));
            testCase.verifyFalse(any(any(diff(binnedfrPOST)==0)));
            
        end


        
        % ========================== EV-REV ========================== %
        
        % Test that intra-area EV-REV works correctly with specified binned
        % firing rates for each epoch
        function testExplainedVariance_basicFunctionalityBinnedFR(testCase)
            
            % Create firing rates and calculate their EV-REV
            [binnedfrPRE, binnedfrPOST, binnedfrTASK, expectedEV, expectedREV, ~, ~, ~] = testCase.createBinnedFR();
            
            % Run EV-REV analysis
            data = struct(...
                'binnedfrPRE', binnedfrPRE,...
                'binnedfrTASK', binnedfrTASK,...
                'binnedfrPOST', binnedfrPOST,...
                'nNAcUnits', []);
            timestamps = struct();
            analysisConfig = struct();
                
            [EV, REV, ~, ~, ~, ~, ~,  ] = testCase.ExplainedVarianceReactivationObj.explainedVariance(data, timestamps, analysisConfig);
            testCase.verifyEqual(EV-REV, expectedEV-expectedREV, 'AbsTol', 0.001);

        end
        
        % Test that intra-area EV-REV returns binned firing rates for each epoch unchanged
        function testExplainedVariance_correctOutputs(testCase)
            
            % Create firing rates and calculate their EV-REV
            [expectedBinnedfrPRE, expectedBinnedfrPOST, expectedBinnedfrTASK, ~, ~, expectedPreCC, expectedTaskCC, expectedPostCC] = testCase.createBinnedFR();
            
            % Run EV-REV analysis
            data = struct(...
                'binnedfrPRE', expectedBinnedfrPRE,...
                'binnedfrTASK', expectedBinnedfrTASK,...
                'binnedfrPOST', expectedBinnedfrPOST,...
                'nNAcUnits', []);
            timestamps = struct();
            analysisConfig = struct();
                
            [~, ~, binnedfrPRE, binnedfrPOST, binnedfrTASK, preCC, taskCC, postCC] = testCase.ExplainedVarianceReactivationObj.explainedVariance(data, timestamps, analysisConfig);
            testCase.verifyEqual(binnedfrPRE, expectedBinnedfrPRE, 'RelTol', 0.001);
            testCase.verifyEqual(binnedfrTASK, expectedBinnedfrTASK, 'RelTol', 0.001);
            testCase.verifyEqual(binnedfrPOST, expectedBinnedfrPOST, 'RelTol', 0.001);
            testCase.verifyEqual(preCC, expectedPreCC, 'RelTol', 0.001);
            testCase.verifyEqual(taskCC, expectedTaskCC, 'RelTol', 0.001);
            testCase.verifyEqual(postCC, expectedPostCC, 'RelTol', 0.001);

        end
        
        % Test that EV-REV returns NaN if empty data matrices are passed
        function testExplainedVariance_noData(testCase)
                       
            % Run EV-REV analysis
            data = struct(...
                'binnedfrPRE', [],...
                'binnedfrTASK', [],...
                'binnedfrPOST', [],...
                'nNAcUnits', []);
            timestamps = struct();
            analysisConfig = struct();
                
            % Call explainedVariance
            [EV, REV, ~, ~, ~, ~, ~, ~] = testCase.ExplainedVarianceReactivationObj.explainedVariance(data, timestamps, analysisConfig);
            testCase.verifyEqual(EV, NaN);
            testCase.verifyEqual(REV, NaN);

        end
        
        % Test EV-REV between brain areas
        function testExplainedVariance_interRegionEVREV(testCase)
            
            % Create firing rates and calculate their EV-REV
            [binnedfrPRE, binnedfrPOST, binnedfrTASK, expectedEV, expectedREV, ~, ~, ~] = testCase.createBinnedFR('inter-area');
            
            % Create inputs to method
            data = struct(...
                'binnedfrPRE', binnedfrPRE,...
                'binnedfrTASK', binnedfrTASK,...
                'binnedfrPOST', binnedfrPOST,...
                'nNAcUnits', 2);
            timestamps = struct();
            analysisConfig = struct();
                
            % Call explainedVariance
            [EV, REV, ~, ~, ~, ~, ~,  ] = testCase.ExplainedVarianceReactivationObj.explainedVariance(data, timestamps, analysisConfig);
            testCase.verifyEqual(EV-REV, expectedEV-expectedREV, 'AbsTol', 0.001);

        end
        
        % Test EV-REV between brain areas with one cell-pair excluded
        function testExplainedVariance_excludeOnePair(testCase)
            
            % Create firing rates and calculate their EV-REV
            [binnedfrPRE, binnedfrPOST, binnedfrTASK, expectedEV, expectedREV, ~, ~, ~] = testCase.createBinnedFR('inter-area');
            
            % Create inputs to method
            data = struct(...
                'binnedfrPRE', binnedfrPRE,...
                'binnedfrTASK', binnedfrTASK,...
                'binnedfrPOST', binnedfrPOST,...
                'nNAcUnits', 2);
            timestamps = struct();
            analysisConfig = struct('exclude', [1, 2]);
                
            % Call explainedVariance
            [EV, REV, ~, ~, ~, ~, ~,  ] = testCase.ExplainedVarianceReactivationObj.explainedVariance(data, timestamps, analysisConfig);
            testCase.verifyNotEqual(EV-REV, expectedEV-expectedREV);

        end
        
        % Test that explainedVariance returns an error if binned firing
        % rate matrices have incompatible sizes
        function testExplainedVariance_incompatibleBinnedFringRates(testCase)
            
            % Create firing rates and calculate their EV-REV
            [binnedfrPRE, binnedfrPOST, binnedfrTASK, ~, ~, ~, ~, ~] = testCase.createBinnedFR('inter-area');
            
            % Create inputs to method, removing one neuron's data from binnedfrPOST
            data = struct(...
                'binnedfrPRE', binnedfrPRE,...
                'binnedfrTASK', binnedfrTASK,...
                'binnedfrPOST', binnedfrPOST(1:end-1, :),...
                'nNAcUnits', 2);
            timestamps = struct();
            analysisConfig = struct();
            
            expectedErrorId = 'ExplainedVarianceReactivation:invalidData';
            testCase.verifyError(@() testCase.ExplainedVarianceReactivationObj.explainedVariance(data, timestamps, analysisConfig), expectedErrorId);
            
        end
        
        % Test that binned firing rates are created correctly if not provided
        function testExplainedVariance_basicFunctionalitySpiketimes(testCase)
            
            % Create 3 cells of spiketimes with constant firing rates for
            % 30 seconds, and 10 seconds each of PRE, TASK and POST
            meanFR = [10, 20, 50];
            spiketimes = {{...
                0 : 1/meanFR(1) : 30,...
                0 : 1/meanFR(2) : 30,...
                0 : 1/meanFR(3) : 30}};
            startEndTimes = [0, 10; 10, 20; 20, 30];
            
            % Run EV-REV analysis
            data = struct(...
                'spiketimes', spiketimes,...
                'nNAcUnits', []);
            timestamps = struct('startEndTimes', startEndTimes);
            analysisConfig = struct(...
                'postDuration', NaN,...
                'rippleActivityMode', 'none');
                
            [EV, REV, ~, ~, ~, ~, ~,  ] = testCase.ExplainedVarianceReactivationObj.explainedVariance(data, timestamps, analysisConfig);
            testCase.verifyEqual(EV-REV, 0, 'AbsTol', 0.001);

        end
        
        
        % ================ per-cell contributions to EV-REV ================ %
        
        % Test basic functionality of calculatePerCellPairContributions
        function testCalculatePerCellPairContributions_basicFunctionality(testCase)

            % Mock data: 5 cells with 100 bins of random firing rates
            % for each of PRE, TASK and POST
            rng(0); binnedfrPRE = randn(5, 100);
            rng(1); binnedfrTASK = randn(5, 100);
            rng(2); binnedfrPOST = randn(5, 100);
            preCC = corr(binnedfrPRE', binnedfrPRE');
            taskCC = corr(binnedfrTASK', binnedfrTASK');
            postCC = corr(binnedfrPOST', binnedfrPOST');
            spiketimes = {repmat({{}}, 1, 5)};

            data = struct(...
                'spiketimes', spiketimes,...
                'binnedfrPRE', binnedfrPRE,...
                'binnedfrTASK', binnedfrTASK,...
                'binnedfrPOST', binnedfrPOST,...
                'preCC', preCC,...
                'taskCC', taskCC,...
                'postCC', postCC,...
                'nUnits', 5,...
                'nNAcUnits', 5);
            includeCells = 1:5;
            intraarea = true;
            EV_REV0 = 0;

            % Call calculatePerCellPairContributions
            perPairContribution = testCase.ExplainedVarianceReactivationObj.calculatePerCellPairContributions(data, includeCells, intraarea, EV_REV0);
            testCase.verifyTrue(isstruct(perPairContribution))
            expectedFields = {'cells', 'contribution', 'correlations'};
            for field = expectedFields
                testCase.verifyEqual(isfield(perPairContribution, field{1}), true, sprintf('Expected field "%s" missing.', field{1}));
            end

        end
        
        % Test that calculatePerCellPairContributions returns the correct data
        % when calculating contributions to EV-REV within one brain area
        function testCalculatePerCellPairContributions_checkDataIntraArea(testCase)

            % Mock data: 5 cells with 100 bins of random firing rates
            % for each of PRE, TASK and POST
            nNeurons = 5;
            rng(0); binnedfrPRE = randn(nNeurons, 100);
            rng(1); binnedfrTASK = randn(nNeurons, 100);
            rng(2); binnedfrPOST = randn(nNeurons, 100);
            preCC = corr(binnedfrPRE', binnedfrPRE');
            taskCC = corr(binnedfrTASK', binnedfrTASK');
            postCC = corr(binnedfrPOST', binnedfrPOST');
            spiketimes = {repmat({{}}, 1, nNeurons)};

            data = struct(...
                'spiketimes', spiketimes,...
                'binnedfrPRE', binnedfrPRE,...
                'binnedfrTASK', binnedfrTASK,...
                'binnedfrPOST', binnedfrPOST,...
                'preCC', preCC,...
                'taskCC', taskCC,...
                'postCC', postCC,...
                'nUnits', nNeurons,...
                'nNAcUnits', nNeurons);
            includeCells = 1:nNeurons;
            intraarea = true;
            EV_REV0 = 0;

            % Call calculatePerCellPairContributions
            perPairContribution = testCase.ExplainedVarianceReactivationObj.calculatePerCellPairContributions(data, includeCells, intraarea, EV_REV0);
                        
            % Test that the data is the expected size given number of neurons
            actualSize = size(perPairContribution);
            expectedSize = [1, nNeurons^2];
            testCase.verifyEqual(actualSize, expectedSize);
            
            % Test that perPairContribution.cells contains the right data
            expectedSubfieldSize = [1, 2];
            for iField = 1:expectedSize(2)
                actualSubfieldSize = size(perPairContribution(iField).cells);
                testCase.verifyEqual(actualSubfieldSize, expectedSubfieldSize);
            end
            actualFirstCell = cellfun(@(x) x(1), {perPairContribution.cells});
            expectedFirstCell = repelem(1:nNeurons, nNeurons);
            testCase.verifyEqual(actualFirstCell, expectedFirstCell);
            actualSecondCell = cellfun(@(x) x(2), {perPairContribution.cells});
            expectedSecondCell = repmat(1:nNeurons, 1, nNeurons);
            testCase.verifyEqual(actualSecondCell, expectedSecondCell);
            
            % Test that perPairContribution.contribution contains the right data
            for iField = 1:expectedSize(2)
                actualSubfieldSize = numel(perPairContribution(iField).contribution);
                testCase.verifyEqual(actualSubfieldSize, 1);
                testCase.verifyNotEqual(actualSubfieldSize, NaN);
            end
            % They should be different for each cell-pair, and not NaN
            testCase.verifyFalse(any(diff([perPairContribution.contribution])==0));
            testCase.verifyFalse(any(isnan([perPairContribution.contribution])));
            
            % Test that perPairContribution.correlations contains the right data
            expectedSubfieldSize = [3, 1];
            for iField = 1:expectedSize(2)
                actualSubfieldSize = size(perPairContribution(iField).correlations);
                testCase.verifyEqual(actualSubfieldSize, expectedSubfieldSize);
                if perPairContribution(iField).cells(1) == perPairContribution(iField).cells(2)
                    % The correlations between a neuron's activity and its own activity should be 1
                    testCase.verifyEqual(perPairContribution(iField).correlations, ones(3, 1), 'RelTol', 0.001),
                else
                    % The correlations between a neuron's activity and another neuron's activity
                    % should be different in the three epochs
                    testCase.verifyNotEqual(perPairContribution(iField).correlations(1), perPairContribution(iField).correlations(2));
                    testCase.verifyNotEqual(perPairContribution(iField).correlations(2), perPairContribution(iField).correlations(3));
                    testCase.verifyNotEqual(perPairContribution(iField).correlations(1), perPairContribution(iField).correlations(3));
                    % None should be NaN
                    testCase.verifyFalse(any(isnan(perPairContribution(iField).correlations)));
                end
            end
            
        end
        
        
        % Test that calculatePerCellPairContributions returns the correct data
        % when calculating contributions to EV-REV between brain areas
        function testCalculatePerCellPairContributions_checkDataInterArea(testCase)

            % Mock data: 5 cells with 100 bins of random firing rates
            % for each of PRE, TASK and POST
            nNeurons = 5;
            nNAcUnits = 2;
            nCA1Units = nNeurons - nNAcUnits;
            rng(0); binnedfrPRE = randn(nNeurons, 100);
            rng(1); binnedfrTASK = randn(nNeurons, 100);
            rng(2); binnedfrPOST = randn(nNeurons, 100);
            preCC = corr(binnedfrPRE(1:nNAcUnits, :)', binnedfrPRE(nNAcUnits+1:end, :)');
            taskCC = corr(binnedfrTASK(1:nNAcUnits, :)', binnedfrTASK(nNAcUnits+1:end, :)');
            postCC = corr(binnedfrPOST(1:nNAcUnits, :)', binnedfrPOST(nNAcUnits+1:end, :)');
            spiketimes = {repmat({{}}, 1, nNeurons)};

            data = struct(...
                'spiketimes', spiketimes,...
                'binnedfrPRE', binnedfrPRE,...
                'binnedfrTASK', binnedfrTASK,...
                'binnedfrPOST', binnedfrPOST,...
                'preCC', preCC,...
                'taskCC', taskCC,...
                'postCC', postCC,...
                'nUnits', nNeurons,...
                'nNAcUnits', nNAcUnits);
            includeCells = 1:nNAcUnits;
            intraarea = false;
            EV_REV0 = 0;

            % Call calculatePerCellPairContributions
            perPairContribution = testCase.ExplainedVarianceReactivationObj.calculatePerCellPairContributions(data, includeCells, intraarea, EV_REV0);
            
            % Check that all expected fields are returned
            testCase.verifyTrue(isstruct(perPairContribution))
            expectedFields = {'cells', 'contribution', 'correlations'};
            for field = expectedFields
                testCase.verifyEqual(isfield(perPairContribution, field{1}), true, sprintf('Expected field "%s" missing.', field{1}));
            end
            
            % Test that the data is the expected size given number of neurons
            actualSize = size(perPairContribution);
            expectedSize = [1, nCA1Units*nNAcUnits];
            testCase.verifyEqual(actualSize, expectedSize);
            
            % Test that perPairContribution.cells contains the right data
            expectedSubfieldSize = [1, 2];
            for iField = 1:expectedSize(2)
                actualSubfieldSize = size(perPairContribution(iField).cells);
                testCase.verifyEqual(actualSubfieldSize, expectedSubfieldSize);
            end
            actualFirstCell = cellfun(@(x) x(1), {perPairContribution.cells});
            expectedFirstCell = repelem(1:nNAcUnits, nCA1Units); 
            testCase.verifyEqual(actualFirstCell, expectedFirstCell);
            actualSecondCell = cellfun(@(x) x(2), {perPairContribution.cells});
            expectedSecondCell = repmat(1:nCA1Units, 1, nNAcUnits);
            testCase.verifyEqual(actualSecondCell, expectedSecondCell);
            
            % Test that perPairContribution.contribution contains the right data
            for iField = 1:expectedSize(2)
                actualSubfieldSize = numel(perPairContribution(iField).contribution);
                testCase.verifyEqual(actualSubfieldSize, 1);
                testCase.verifyNotEqual(actualSubfieldSize, NaN);
            end
            % They should be different for each cell-pair, and not NaN
            testCase.verifyFalse(any(diff([perPairContribution.contribution])==0));
            testCase.verifyFalse(any(isnan([perPairContribution.contribution])));
            
            % Test that perPairContribution.correlations contains the right data
            expectedSubfieldSize = [3, 1];
            for iField = 1:expectedSize(2)
                actualSubfieldSize = size(perPairContribution(iField).correlations);
                testCase.verifyEqual(actualSubfieldSize, expectedSubfieldSize);
                % The correlations between a neuron's activity and another neuron's activity
                % should be different in the three epochs
                testCase.verifyNotEqual(perPairContribution(iField).correlations(1), perPairContribution(iField).correlations(2));
                testCase.verifyNotEqual(perPairContribution(iField).correlations(2), perPairContribution(iField).correlations(3));
                testCase.verifyNotEqual(perPairContribution(iField).correlations(1), perPairContribution(iField).correlations(3));
                % None should be NaN
                testCase.verifyFalse(any(isnan(perPairContribution(iField).correlations)));
            end

        end
        
        
        % Test that calculatePerCellPairContributions returns data only for included cells
        % when calculating contributions to EV-REV within one brain area
        function testCalculatePerCellPairContributions_intraAreaIncludeCells(testCase)

            % Mock data: 5 cells with 100 bins of random firing rates
            % for each of PRE, TASK and POST
            nNeurons = 5;
            rng(0); binnedfrPRE = randn(nNeurons, 100);
            rng(1); binnedfrTASK = randn(nNeurons, 100);
            rng(2); binnedfrPOST = randn(nNeurons, 100);
            preCC = corr(binnedfrPRE', binnedfrPRE');
            taskCC = corr(binnedfrTASK', binnedfrTASK');
            postCC = corr(binnedfrPOST', binnedfrPOST');
            spiketimes = {repmat({{}}, 1, nNeurons)};

            data = struct(...
                'spiketimes', spiketimes,...
                'binnedfrPRE', binnedfrPRE,...
                'binnedfrTASK', binnedfrTASK,...
                'binnedfrPOST', binnedfrPOST,...
                'preCC', preCC,...
                'taskCC', taskCC,...
                'postCC', postCC,...
                'nUnits', nNeurons,...
                'nNAcUnits', nNeurons);
            includeCells = [2, 3, 5];
            intraarea = true;
            EV_REV0 = 0;

            % Call calculatePerCellPairContributions
            perPairContribution = testCase.ExplainedVarianceReactivationObj.calculatePerCellPairContributions(data, includeCells, intraarea, EV_REV0);
                        
            % Test that the data is the expected size given number of neurons
            actualSize = size(perPairContribution);
            expectedSize = [1, length(includeCells)^2];
            testCase.verifyEqual(actualSize, expectedSize);
            
            % Test that perPairContribution.cells contains the right data
            expectedSubfieldSize = [1, 2];
            for iField = 1:expectedSize(2)
                actualSubfieldSize = size(perPairContribution(iField).cells);
                testCase.verifyEqual(actualSubfieldSize, expectedSubfieldSize);
            end
            actualFirstCell = cellfun(@(x) x(1), {perPairContribution.cells});
            expectedFirstCell = repelem(includeCells, length(includeCells));
            testCase.verifyEqual(actualFirstCell, expectedFirstCell);
            actualSecondCell = cellfun(@(x) x(2), {perPairContribution.cells});
            expectedSecondCell = repmat(includeCells, 1, length(includeCells));
            testCase.verifyEqual(actualSecondCell, expectedSecondCell);
            
            % Test that perPairContribution.contribution contains the right data
            for iField = 1:expectedSize(2)
                actualSubfieldSize = numel(perPairContribution(iField).contribution);
                testCase.verifyEqual(actualSubfieldSize, 1);
                testCase.verifyNotEqual(actualSubfieldSize, NaN);
            end
            % They should be different for each cell-pair, and not NaN
            testCase.verifyFalse(any(diff([perPairContribution.contribution])==0));
            testCase.verifyFalse(any(isnan([perPairContribution.contribution])));
            
            % Test that perPairContribution.correlations contains the right data
            expectedSubfieldSize = [3, 1];
            for iField = 1:expectedSize(2)
                actualSubfieldSize = size(perPairContribution(iField).correlations);
                testCase.verifyEqual(actualSubfieldSize, expectedSubfieldSize);
                if perPairContribution(iField).cells(1) == perPairContribution(iField).cells(2)
                    % The correlations between a neuron's activity and its own activity should be 1
                    testCase.verifyEqual(perPairContribution(iField).correlations, ones(3, 1), 'RelTol', 0.001),
                else
                    % The correlations between a neuron's activity and another neuron's activity
                    % should be different in the three epochs
                    testCase.verifyNotEqual(perPairContribution(iField).correlations(1), perPairContribution(iField).correlations(2));
                    testCase.verifyNotEqual(perPairContribution(iField).correlations(2), perPairContribution(iField).correlations(3));
                    testCase.verifyNotEqual(perPairContribution(iField).correlations(1), perPairContribution(iField).correlations(3));
                    % None should be NaN
                    testCase.verifyFalse(any(isnan(perPairContribution(iField).correlations)));
                end
            end
            
        end

        % Test that calculatePerCellPairContributions returns data only for included cells
        % when calculating contributions to EV-REV between brain areas
        function testCalculatePerCellPairContributions_interAreaIncludeCells(testCase)
            
            % Mock data: 10 cells with 100 bins of random firing rates
            % for each of PRE, TASK and POST
            nNeurons = 10;
            nNAcUnits = 6;
            nCA1Units = nNeurons - nNAcUnits;
            rng(0); binnedfrPRE = randn(nNeurons, 100);
            rng(1); binnedfrTASK = randn(nNeurons, 100);
            rng(2); binnedfrPOST = randn(nNeurons, 100);
            preCC = corr(binnedfrPRE(1:nNAcUnits, :)', binnedfrPRE(nNAcUnits+1:end, :)');
            taskCC = corr(binnedfrTASK(1:nNAcUnits, :)', binnedfrTASK(nNAcUnits+1:end, :)');
            postCC = corr(binnedfrPOST(1:nNAcUnits, :)', binnedfrPOST(nNAcUnits+1:end, :)');
            spiketimes = {repmat({{}}, 1, nNeurons)};

            data = struct(...
                'spiketimes', spiketimes,...
                'binnedfrPRE', binnedfrPRE,...
                'binnedfrTASK', binnedfrTASK,...
                'binnedfrPOST', binnedfrPOST,...
                'preCC', preCC,...
                'taskCC', taskCC,...
                'postCC', postCC,...
                'nUnits', nNeurons,...
                'nNAcUnits', nNAcUnits);
            includeCells = [1, 4, 5];
            intraarea = false;
            EV_REV0 = 0;

            % Call calculatePerCellPairContributions
            perPairContribution = testCase.ExplainedVarianceReactivationObj.calculatePerCellPairContributions(data, includeCells, intraarea, EV_REV0);
                        
            % Test that the data is the expected size given number of neurons
            actualSize = size(perPairContribution);
            expectedSize = [1, nCA1Units * length(includeCells)];
            testCase.verifyEqual(actualSize, expectedSize);
            
            % Test that perPairContribution.cells contains the right data
            expectedSubfieldSize = [1, 2];
            for iField = 1:expectedSize(2)
                actualSubfieldSize = size(perPairContribution(iField).cells);
                testCase.verifyEqual(actualSubfieldSize, expectedSubfieldSize);
            end
            actualFirstCell = cellfun(@(x) x(1), {perPairContribution.cells});
            expectedFirstCell = repelem(includeCells, nCA1Units); 
            testCase.verifyEqual(actualFirstCell, expectedFirstCell);
            actualSecondCell = cellfun(@(x) x(2), {perPairContribution.cells});
            expectedSecondCell = repmat(1:nCA1Units, 1, length(includeCells));
            testCase.verifyEqual(actualSecondCell, expectedSecondCell);
            
            % Test that perPairContribution.contribution contains the right data
            for iField = 1:expectedSize(2)
                actualSubfieldSize = numel(perPairContribution(iField).contribution);
                testCase.verifyEqual(actualSubfieldSize, 1);
                testCase.verifyNotEqual(actualSubfieldSize, NaN);
            end
            % They should be different for each cell-pair, and not NaN
            testCase.verifyFalse(any(diff([perPairContribution.contribution])==0));
            testCase.verifyFalse(any(isnan([perPairContribution.contribution])));
            
            % Test that perPairContribution.correlations contains the right data
            expectedSubfieldSize = [3, 1];
            for iField = 1:expectedSize(2)
                actualSubfieldSize = size(perPairContribution(iField).correlations);
                testCase.verifyEqual(actualSubfieldSize, expectedSubfieldSize);
                % The correlations between a neuron's activity and another neuron's activity
                % should be different in the three epochs
                testCase.verifyNotEqual(perPairContribution(iField).correlations(1), perPairContribution(iField).correlations(2));
                testCase.verifyNotEqual(perPairContribution(iField).correlations(2), perPairContribution(iField).correlations(3));
                testCase.verifyNotEqual(perPairContribution(iField).correlations(1), perPairContribution(iField).correlations(3));
                % None should be NaN
                testCase.verifyFalse(any(isnan(perPairContribution(iField).correlations)));
            end
            
        end
        
        % Test that calculatePerCellPairContributions returns an error if
        % binned firing rates are missing
        function testCalculatePerCellPairContributions_missingBinnedFRTask(testCase)
            
            % Mock data: 10 cells with 100 bins of random firing rates
            % for each of PRE, TASK and POST
            nNeurons = 10;
            nNAcUnits = 6;
            nCA1Units = nNeurons - nNAcUnits;
            rng(0); binnedfrPRE = randn(nNeurons, 100);
            rng(2); binnedfrPOST = randn(nNeurons, 100);
            preCC = corr(binnedfrPRE(1:nNAcUnits, :)', binnedfrPRE(nNAcUnits+1:end, :)');
            taskCC = preCC;
            postCC = corr(binnedfrPOST(1:nNAcUnits, :)', binnedfrPOST(nNAcUnits+1:end, :)');
            spiketimes = {repmat({{}}, 1, nNeurons)};

            data = struct(...
                'spiketimes', spiketimes,...
                'binnedfrPRE', binnedfrPRE,...
                'binnedfrPOST', binnedfrPOST,...
                'preCC', preCC,...
                'taskCC', taskCC,...
                'postCC', postCC,...
                'nUnits', nNeurons,...
                'nNAcUnits', nNAcUnits);
            includeCells = 1:nNAcUnits;
            intraarea = false;
            EV_REV0 = 0;

            % Call calculatePerCellPairContributions
            expectedErrorId = 'ExplainedVarianceReactivation:dataNotAssigned';
            testCase.verifyError(@() testCase.ExplainedVarianceReactivationObj.calculatePerCellPairContributions(data, includeCells, intraarea, EV_REV0), expectedErrorId);

        end
        
        % Test that calculatePerCellPairContributions returns an error if
        % binned firing rates are empty
        function testCalculatePerCellPairContributions_emptyBinnedFRPost(testCase)
            
            % Mock data: 10 cells with 100 bins of random firing rates
            % for each of PRE, TASK and POST
            nNeurons = 10;
            nNAcUnits = 6;
            nCA1Units = nNeurons - nNAcUnits;
            rng(0); binnedfrPRE = randn(nNeurons, 100);
            rng(1); binnedfrTASK = randn(nNeurons, 100);
            preCC = corr(binnedfrPRE(1:nNAcUnits, :)', binnedfrPRE(nNAcUnits+1:end, :)');
            taskCC = corr(binnedfrTASK(1:nNAcUnits, :)', binnedfrTASK(nNAcUnits+1:end, :)');
            postCC = taskCC;
            spiketimes = {repmat({{}}, 1, nNeurons)};

            data = struct(...
                'spiketimes', spiketimes,...
                'binnedfrPRE', binnedfrPRE,...
                'binnedfrTASK', binnedfrTASK,...
                'binnedfrPOST', [],...
                'preCC', preCC,...
                'taskCC', taskCC,...
                'postCC', postCC,...
                'nUnits', nNeurons,...
                'nNAcUnits', nNAcUnits);
            includeCells = 1:nNAcUnits;
            intraarea = false;
            EV_REV0 = 0;

            % Call calculatePerCellPairContributions
            expectedErrorId = 'ExplainedVarianceReactivation:invalidData';
            testCase.verifyError(@() testCase.ExplainedVarianceReactivationObj.calculatePerCellPairContributions(data, includeCells, intraarea, EV_REV0), expectedErrorId);

        end
        
        % Test that calculatePerCellPairContributions returns an error if
        % neural activity is for inter-area but intraarea is specified
        function testCalculatePerCellPairContributions_misspecifiedIntraarea(testCase)
            
            % Mock data: 10 cells with 100 bins of random firing rates
            % for each of PRE, TASK and POST
            nNeurons = 10;
            nNAcUnits = 6;
            nCA1Units = nNeurons - nNAcUnits;
            rng(0); binnedfrPRE = randn(nNeurons, 100);
            rng(1); binnedfrTASK = randn(nNeurons, 100);
            rng(2); binnedfrPOST = randn(nNeurons, 100);
            preCC = corr(binnedfrPRE(1:nNAcUnits, :)', binnedfrPRE(nNAcUnits+1:end, :)');
            taskCC = corr(binnedfrTASK(1:nNAcUnits, :)', binnedfrTASK(nNAcUnits+1:end, :)');
            postCC = corr(binnedfrPOST(1:nNAcUnits, :)', binnedfrPOST(nNAcUnits+1:end, :)');
            spiketimes = {repmat({{}}, 1, nNeurons)};

            data = struct(...
                'spiketimes', spiketimes,...
                'binnedfrPRE', binnedfrPRE,...
                'binnedfrTASK', binnedfrTASK,...
                'binnedfrPOST', binnedfrPOST,...
                'preCC', preCC,...
                'taskCC', taskCC,...
                'postCC', postCC,...
                'nUnits', nNeurons,...
                'nNAcUnits', nNAcUnits);
            includeCells = 1:nNAcUnits;
            intraarea = true;
            EV_REV0 = 0;

            % Call calculatePerCellPairContributions
            expectedErrorId = 'ExplainedVarianceReactivation:invalidData';
            testCase.verifyError(@() testCase.ExplainedVarianceReactivationObj.calculatePerCellPairContributions(data, includeCells, intraarea, EV_REV0), expectedErrorId);

        end
        
        % Test that calculatePerCellPairContributions returns an error if
        % there are too many includeCells
        function testCalculatePerCellPairContributions_tooManyIncludeCells(testCase)
            
            % Mock data: 10 cells with 100 bins of random firing rates
            % for each of PRE, TASK and POST
            nNeurons = 10;
            nNAcUnits = 6;
            nCA1Units = nNeurons - nNAcUnits;
            rng(0); binnedfrPRE = randn(nNeurons, 100);
            rng(1); binnedfrTASK = randn(nNeurons, 100);
            rng(2); binnedfrPOST = randn(nNeurons, 100);
            preCC = corr(binnedfrPRE(1:nNAcUnits, :)', binnedfrPRE(nNAcUnits+1:end, :)');
            taskCC = corr(binnedfrTASK(1:nNAcUnits, :)', binnedfrTASK(nNAcUnits+1:end, :)');
            postCC = corr(binnedfrPOST(1:nNAcUnits, :)', binnedfrPOST(nNAcUnits+1:end, :)');
            spiketimes = {repmat({{}}, 1, nNeurons)};

            data = struct(...
                'spiketimes', spiketimes,...
                'binnedfrPRE', binnedfrPRE,...
                'binnedfrTASK', binnedfrTASK,...
                'binnedfrPOST', binnedfrPOST,...
                'preCC', preCC,...
                'taskCC', taskCC,...
                'postCC', postCC,...
                'nUnits', nNeurons,...
                'nNAcUnits', nNAcUnits);
            includeCells = 1:nNAcUnits*2;
            intraarea = false;
            EV_REV0 = 0;

            % Call calculatePerCellPairContributions
            expectedErrorId = 'ExplainedVarianceReactivation:invalidData';
            testCase.verifyError(@() testCase.ExplainedVarianceReactivationObj.calculatePerCellPairContributions(data, includeCells, intraarea, EV_REV0), expectedErrorId);

        end
        
        % Test that calculatePerCellPairContributions returns an empty
        % structure for interarea if there are no hippocampal units
        function testCalculatePerCellPairContributions_noCA1Units(testCase)
            
            % Mock data: 10 cells with 100 bins of random firing rates
            % for each of PRE, TASK and POST
            nNeurons = 5;
            nNAcUnits = 5;
            nCA1Units = 0;
            rng(0); binnedfrPRE = randn(nNeurons, 100);
            rng(1); binnedfrTASK = randn(nNeurons, 100);
            rng(2); binnedfrPOST = randn(nNeurons, 100);
            preCC = corr(binnedfrPRE(1:nNAcUnits, :)', binnedfrPRE(nNAcUnits+1:end, :)');
            taskCC = corr(binnedfrTASK(1:nNAcUnits, :)', binnedfrTASK(nNAcUnits+1:end, :)');
            postCC = corr(binnedfrPOST(1:nNAcUnits, :)', binnedfrPOST(nNAcUnits+1:end, :)');
            spiketimes = {repmat({{}}, 1, nNeurons)};

            data = struct(...
                'spiketimes', spiketimes,...
                'binnedfrPRE', binnedfrPRE,...
                'binnedfrTASK', binnedfrTASK,...
                'binnedfrPOST', binnedfrPOST,...
                'preCC', preCC,...
                'taskCC', taskCC,...
                'postCC', postCC,...
                'nUnits', nNeurons,...
                'nNAcUnits', nNAcUnits);
            includeCells = 1:nNAcUnits;
            intraarea = false;
            EV_REV0 = 0;

            % Call calculatePerCellPairContributions
            perPairContribution = testCase.ExplainedVarianceReactivationObj.calculatePerCellPairContributions(data, includeCells, intraarea, EV_REV0);
            
            % Check that the output contains the right fields, but no values
            testCase.verifyTrue(isstruct(perPairContribution), 'Expected a structure to be returned');
            testCase.verifyEmpty(perPairContribution, 'Returned structure expected to be empty');
            expectedFields = {'cells', 'contribution', 'correlations'};
            for field = expectedFields
                testCase.verifyEqual(isfield(perPairContribution, field{1}), true, sprintf('Expected field "%s" missing.', field{1}));
            end
            
        end
            
        % Test that calculatePerCellPairContributions returns an empty
        % structure for intraarea if there are no accumbens units
        function testCalculatePerCellPairContributions_noNAcUnits(testCase)
            
            % Mock data: 10 cells with 100 bins of random firing rates
            % for each of PRE, TASK and POST
            nNeurons = 5;
            nNAcUnits = 0;
            nCA1Units = 5;
            rng(0); binnedfrPRE = randn(nNeurons, 100);
            rng(1); binnedfrTASK = randn(nNeurons, 100);
            rng(2); binnedfrPOST = randn(nNeurons, 100);
            preCC = [];
            taskCC = [];
            postCC = [];
            spiketimes = {repmat({{}}, 1, nNeurons)};

            data = struct(...
                'spiketimes', spiketimes,...
                'binnedfrPRE', binnedfrPRE,...
                'binnedfrTASK', binnedfrTASK,...
                'binnedfrPOST', binnedfrPOST,...
                'preCC', preCC,...
                'taskCC', taskCC,...
                'postCC', postCC,...
                'nUnits', nNeurons,...
                'nNAcUnits', nNAcUnits);
            includeCells = [];
            intraarea = true;
            EV_REV0 = 0;

            % Call calculatePerCellPairContributions
            perPairContribution = testCase.ExplainedVarianceReactivationObj.calculatePerCellPairContributions(data, includeCells, intraarea, EV_REV0);
            
            % Check that the output contains the right fields, but no values
            testCase.verifyTrue(isstruct(perPairContribution), 'Expected a structure to be returned');
            testCase.verifyEmpty(perPairContribution, 'Returned structure expected to be empty');
            expectedFields = {'cells', 'contribution', 'correlations'};
            for field = expectedFields
                testCase.verifyEqual(isfield(perPairContribution, field{1}), true, sprintf('Expected field "%s" missing.', field{1}));
            end
                        
        end

        
        % =================== getting immobility periods =================== %
        
        % Test basic functionality of getFixedDurationImmobility
        function testGetFixedDurationImmobility_basicFunctionality(testCase)
            
            % Mock data: 2 rest periods of 1 second each in PRE, and 2 in POST
            testCase.ExplainedVarianceReactivationObj.analysisConfig.postDuration = NaN;
            data = struct(...
                'startEndTimes', [0, 10; 10, 20; 20, 30],...
                'restStart', [1, 5, 21, 25],...
                'restStop', [2, 6, 22, 26]);
            
            
            % Call getFixedDurationImmobility
            [restStartPRE, restStopPRE, restStartPOST, restStopPOST] = testCase.ExplainedVarianceReactivationObj.getFixedDurationImmobility(data);
            testCase.verifyEqual(size(restStartPRE), [1, 2], 'Expected restStartPRE to be a row vector of length 2');
            testCase.verifyEqual(size(restStopPRE), [1, 2], 'Expected restStopPRE to be a row vector of length 2');
            testCase.verifyEqual(size(restStartPOST), [1, 2], 'Expected restStartPOST to be a row vector of length 2');
            testCase.verifyEqual(size(restStopPOST), [1, 2], 'Expected restStopPOST to be a row vector of length 2');
            
        end
        
        % Test that getFixedDurationImmobility includes only immobility during PRE and POST
        function testGetFixedDurationImmobility_excludesImmobilityOutsideEpochs(testCase)
            
            % Mock data: 2 rest periods of 1 second each in PRE, and 2 in POST,
            % plus 1 each before PRE, during TASK and after POST which should be ignored
            testCase.ExplainedVarianceReactivationObj.analysisConfig.postDuration = NaN;
            data = struct(...
                'startEndTimes', [0, 10; 10, 20; 20, 30],...
                'restStart', [-3, 1, 5, 15, 21, 25, 31],...
                'restStop', [-2, 2, 6, 16, 22, 26, 32]);
            
            
            % Call getFixedDurationImmobility
            [restStartPRE, restStopPRE, restStartPOST, restStopPOST] = testCase.ExplainedVarianceReactivationObj.getFixedDurationImmobility(data);
            testCase.verifyEqual(restStartPRE, [1, 5]);
            testCase.verifyEqual(restStopPRE, [2, 6]);
            testCase.verifyEqual(restStartPOST, [21, 25]);
            testCase.verifyEqual(restStopPOST, [22, 26]);
            
        end
        
        % Test that getFixedDurationImmobility uses postDuration parameter correctly
        function testGetFixedDurationImmobility_postDuration(testCase)
            
            % Mock data: 2 rest periods of 1 second each in PRE, and 2 in POST,
            % plus 1 each before PRE, during TASK and after POST which should be ignored
            testCase.ExplainedVarianceReactivationObj.analysisConfig.postDuration = 2;
            data = struct(...
                'startEndTimes', [0, 10; 10, 20; 20, 30],...
                'restStart', [1, 5, 9,  21, 25, 29],...
                'restStop', [2, 6, 10, 22, 26, 30]);
            
            
            % Call getFixedDurationImmobility
            [restStartPRE, restStopPRE, restStartPOST, restStopPOST] = testCase.ExplainedVarianceReactivationObj.getFixedDurationImmobility(data);
            testCase.verifyEqual(restStartPRE, [1, 5, 9]);
            testCase.verifyEqual(restStopPRE, [2, 6, 10]);
            testCase.verifyEqual(restStartPOST, [21, 25]);
            testCase.verifyEqual(restStopPOST, [22, 26]);
            
        end
        
        % Test that getFixedDurationImmobility uses postDuration parameter
        % correctly to split immobiility bouts where necessary
        function testGetFixedDurationImmobility_postDurationSplitBouts(testCase)
            
            % Mock data: 2 rest periods of 1 second each in PRE, and 2 in POST,
            % plus 1 each before PRE, during TASK and after POST which should be ignored
            testCase.ExplainedVarianceReactivationObj.analysisConfig.postDuration = 2.5;
            data = struct(...
                'startEndTimes', [0, 10; 10, 20; 20, 30],...
                'restStart', [1, 5, 8,  21, 25, 29],...
                'restStop', [2, 6, 9, 22, 26, 30]);
            
            
            % Call getFixedDurationImmobility
            [restStartPRE, restStopPRE, restStartPOST, restStopPOST] = testCase.ExplainedVarianceReactivationObj.getFixedDurationImmobility(data);
            testCase.verifyEqual(restStartPRE, [1, 5, 8]);
            testCase.verifyEqual(restStopPRE, [2, 6, 9]);
            testCase.verifyEqual(restStartPOST, [21, 25, 29]);
            testCase.verifyEqual(restStopPOST, [22, 26, 29.5]);
            
        end
        
        % Test that getFixedDurationImmobility uses postDuration parameter
        % correctly when greater than the total immobility duration
        function testGetFixedDurationImmobility_postDurationTooLong(testCase)
            
            % Mock data: 2 rest periods of 1 second each in PRE, and 2 in POST,
            % plus 1 each before PRE, during TASK and after POST which should be ignored
            testCase.ExplainedVarianceReactivationObj.analysisConfig.postDuration = 5;
            data = struct(...
                'startEndTimes', [0, 10; 10, 20; 20, 30],...
                'restStart', [1, 5, 8,  21, 25, 29],...
                'restStop', [2, 6, 9, 22, 26, 30]);
            
            
            % Call getFixedDurationImmobility
            [restStartPRE, restStopPRE, restStartPOST, restStopPOST] = testCase.ExplainedVarianceReactivationObj.getFixedDurationImmobility(data);
            testCase.verifyEqual(restStartPRE, [1, 5, 8]);
            testCase.verifyEqual(restStopPRE, [2, 6, 9]);
            testCase.verifyEqual(restStartPOST, [21, 25, 29]);
            testCase.verifyEqual(restStopPOST, [22, 26, 30]);
            
        end
        
        % Test that getFixedDurationImmobility handles there being no immobility
        function testGetFixedDurationImmobility_noData(testCase)
            
            % Mock data: 2 rest periods of 1 second each in PRE, and 2 in POST,
            % plus 1 each before PRE, during TASK and after POST which should be ignored
            testCase.ExplainedVarianceReactivationObj.analysisConfig.postDuration = 5;
            data = struct(...
                'startEndTimes', [0, 10; 10, 20; 20, 30],...
                'restStart', [],...
                'restStop', []);
            
            
            % Call getFixedDurationImmobility
            [restStartPRE, restStopPRE, restStartPOST, restStopPOST] = testCase.ExplainedVarianceReactivationObj.getFixedDurationImmobility(data);
            testCase.verifyEmpty(restStartPRE);
            testCase.verifyEmpty(restStopPRE);
            testCase.verifyEmpty(restStartPOST);
            testCase.verifyEmpty(restStopPOST);
            
        end
        
        % Test that getFixedDurationImmobility handles there being one rest bout
        function testGetFixedDurationImmobility_singleBout(testCase)
            
            % Mock data: 2 rest periods of 1 second each in PRE, and 2 in POST,
            % plus 1 each before PRE, during TASK and after POST which should be ignored
            testCase.ExplainedVarianceReactivationObj.analysisConfig.postDuration = 2;
            data = struct(...
                'startEndTimes', [0, 10; 10, 20; 20, 30],...
                'restStart', [28],...
                'restStop', [30]);
            
            
            % Call getFixedDurationImmobility
            [restStartPRE, restStopPRE, restStartPOST, restStopPOST] = testCase.ExplainedVarianceReactivationObj.getFixedDurationImmobility(data);
            testCase.verifyEmpty(restStartPRE);
            testCase.verifyEmpty(restStopPRE);
            testCase.verifyEqual(restStartPOST, [28]);
            testCase.verifyEqual(restStopPOST, [30]);
            
        end
        
        
        % ================ reactivation between brain areas ================ %
        
        % Test basic functionality
        function testCalculateInterAreaEVREV_basicFunctionality(testCase)
            
            % Mock inclusion criteria and sessions: 1 rat with 1 session
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rewardResponsive = false;
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.minimumFiringRate = 0;
            testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
            testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;
            
            % Mock analysis configuration
            testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'none';
            
            % Create mock behavioural data
            behavDirectory = fullfile(testCase.TestDataDir, 'test_session');
            ratName = 'A_Rat';
            TestUtilities.saveMockBehaviouralData(behavDirectory, ratName);
            
            % Mock the rest of the data: 3 vStr and 3 CA1 neurons and EV<REV
            ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
            nNAcUnits = 3;
            nCA1Units = 3;
            [spiketimes, binnedfr, startEndTimes, EV, REV] = testCase.createSpiketimes('interarea', true, 'nNAcUnits', nNAcUnits, 'nCA1Units', nCA1Units, 'seed', 1);
            assert(EV<REV, 'Mock data EV-REV should be negative for simpler testing in this function')
            TestUtilities.saveMockEphysData(ephysDirectory,...
                'spiketimes', spiketimes,...
                'startEndTimes', startEndTimes,...
                'nNAcUnits', nNAcUnits,...
                'binnedfr', binnedfr);
            
            % Call calculateInterAreaEVREV
            output = testCase.ExplainedVarianceReactivationObj.calculateInterAreaEVREV('dataPath', testCase.TestDataDir);
            
            % Verify output structure
            testCase.verifyNotEmpty(output);
            expectedFields = {'CA1_vStr', 'perCellContributions'};
            actualFields = fieldnames(output);
            for field = expectedFields
                testCase.verifyEqual(ismember(field{1}, actualFields), true, sprintf('Expected field "%s" missing.', field{1}));
                testCase.verifyEqual(ismember((field{1}), actualFields), true, sprintf('Field "%s" expected to be a struct.', field{1}));
            end
            expectedSubfields = {'EV', 'REV'};
            for subfield = expectedSubfields
                testCase.verifyEqual(isfield(output.CA1_vStr, subfield{1}), true, sprintf('Expected subfield "%s" missing.', subfield{1}));
                testCase.verifyEqual(isfield(output.CA1_vStr.(subfield{1}), 'A_rat'), true, sprintf('Expected subfield "A_rat" missing from field "%s".', subfield{1}));
                testCase.verifyEqual(length(output.CA1_vStr.(subfield{1}).A_rat), 1, sprintf('Expected 1 %s value for A_rat.', subfield{1}));
            end
            
            % Verify returned EV and REV values
            testCase.verifyEqual(output.CA1_vStr.EV.A_rat, EV, 'AbsTol', 1e-3, 'Returned value of EV does not match the expected value.');
            testCase.verifyEqual(output.CA1_vStr.REV.A_rat, REV, 'AbsTol', 1e-3, 'Returned value of REV does not match the expected value.');
            
        end
        
        % Test calculateInterAreaEVREV with one session of positive EV-REV
        % and one session of negative EV-REV
        function testCalculateInterAreaEVREV_twoSessions(testCase)
            
            % Mock inclusion criteria and sessions: 1 rat with 1 session
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rewardResponsive = false;
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.minimumFiringRate = 0;
            testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [1, 2]);
            testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;
            
            % Mock analysis configuration
            testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'none';
            
            % Create mock behavioural data
            behavDirectory = fullfile(testCase.TestDataDir, 'test_session');
            ratName = 'A_Rat';
            TestUtilities.saveMockBehaviouralData(behavDirectory, ratName);
            
            % Mock ephys data for first session: 3 vStr and 3 CA1 neurons and EV<REV
            ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
            nNAcUnits = 3;
            nCA1Units = 3;
            [spiketimes, binnedfr, startEndTimes, EV_session1, REV_session1] = testCase.createSpiketimes('interarea', true, 'nNAcUnits', nNAcUnits, 'nCA1Units', nCA1Units, 'seed', 0);
            assert(EV_session1<REV_session1, 'Mock data EV-REV should be negative for session 1')
            TestUtilities.saveMockEphysData(ephysDirectory,...
                'spiketimes', spiketimes,...
                'startEndTimes', startEndTimes,...
                'nNAcUnits', nNAcUnits,...
                'binnedfr', binnedfr);
            
            % Mock ephys data for second session: 3 vStr and 4 CA1 neurons and EV<REV
            ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session2');
            nNAcUnits = 3;
            nCA1Units = 4;
            [spiketimes, binnedfr, startEndTimes, EV_session2, REV_session2] = testCase.createSpiketimes('interarea', true, 'nNAcUnits', nNAcUnits, 'nCA1Units', nCA1Units, 'seed', 1);
            assert(EV_session2>REV_session2, 'Mock data EV-REV should be positive for session 2')
            TestUtilities.saveMockEphysData(ephysDirectory,...
                'spiketimes', spiketimes,...
                'startEndTimes', startEndTimes,...
                'nNAcUnits', nNAcUnits,...
                'binnedfr', binnedfr);
            
            % Call calculateInterAreaEVREV
            output = testCase.ExplainedVarianceReactivationObj.calculateInterAreaEVREV('dataPath', testCase.TestDataDir);
            
            % Verify output structure
            testCase.verifyNotEmpty(output);
            expectedFields = {'CA1_vStr', 'perCellContributions'};
            actualFields = fieldnames(output);
            for field = expectedFields
                testCase.verifyEqual(ismember(field{1}, actualFields), true, sprintf('Expected field "%s" missing.', field{1}));
                testCase.verifyEqual(ismember((field{1}), actualFields), true, sprintf('Field "%s" expected to be a struct.', field{1}));
            end
            expectedSubfields = {'EV', 'REV'};
            for subfield = expectedSubfields
                testCase.verifyEqual(isfield(output.CA1_vStr, subfield{1}), true, sprintf('Expected subfield "%s" missing.', subfield{1}));
                testCase.verifyEqual(isfield(output.CA1_vStr.(subfield{1}), 'A_rat'), true, sprintf('Expected subfield "A_rat" missing from field "%s".', subfield{1}));
                testCase.verifyEqual(length(output.CA1_vStr.(subfield{1}).A_rat), 2, sprintf('Expected 1 %s value for A_rat.', subfield{1}));
            end
            testCase.verifyEqual(output.perCellContributions.A_rat{1}, NaN, 'Expected the per-cell contributions for session 1 to be NaN.')
            testCase.verifyEqual(isstruct(output.perCellContributions.A_rat{2}), true, 'Expected the per-cell contributions for session 2 to be a struct.')
            testCase.verifyEqual(length(output.perCellContributions.A_rat{2}), nNAcUnits*nCA1Units, 'Expected the per-cell contributions for session 2 to match the number of cell pairs.')
            expectedSubfields = {'cells', 'contribution', 'correlations'};
            for subfield = expectedSubfields
                testCase.verifyEqual(isfield(output.perCellContributions.A_rat{2}, subfield{1}), true, sprintf('Expected subfield "%s" missing from per-cell contributions.', subfield{1}));
            end
            
            % Verify returned EV and REV values
            for session = 1:2
                eval(['EV = EV_session' num2str(session) ';']);
                eval(['REV = REV_session' num2str(session) ';']);
                testCase.verifyEqual(output.CA1_vStr.EV.A_rat(session), EV, 'AbsTol', 1e-5, sprintf('Returned value of EV does not match the expected value for session %i.', session));
                testCase.verifyEqual(output.CA1_vStr.REV.A_rat(session), REV, 'AbsTol', 1e-5, sprintf('Returned value of REV does not match the expected value for session %i.', session));
            end
            
        end
        
        % Test that calculateInterAreaEVREV correctly implements the
        % minimum firing rate threshold
        function testCalculateInterAreaEVREV_testMinimumFiringRate(testCase)
            
            % Mock inclusion criteria and sessions: 1 rat with 1 session
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rewardResponsive = false;
            testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
            testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;
            
            % Mock analysis configuration
            testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'none';
            
            % Create mock behavioural data
            behavDirectory = fullfile(testCase.TestDataDir, 'test_session');
            ratName = 'A_Rat';
            TestUtilities.saveMockBehaviouralData(behavDirectory, ratName);
            
            % Mock the rest of the data: 5 vStr of which 3 are above the minimum firing rate and 3 CA1 neurons
            ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
            nNAcUnits = 5;
            nCA1Units = 4;
            seed = 5;
            [spiketimes, binnedfr, startEndTimes, EV, REV, ~, ~] = testCase.createSpiketimes('interarea', true, 'ripples', false, 'nNAcUnits', nNAcUnits, 'nCA1Units', nCA1Units, 'seed', seed);
            minFR = 15;
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.minimumFiringRate = minFR;
            assert(sum(mean(binnedfr(1:nNAcUnits, :), 2) > minFR) == 3, 'Mock data should include three NAc neurons above the minimum firing rate threshold for testing in this function')
            TestUtilities.saveMockEphysData(ephysDirectory,...
                'spiketimes', spiketimes,...
                'startEndTimes', startEndTimes,...
                'nNAcUnits', nNAcUnits,...
                'binnedfr', binnedfr);
            
            % Call calculateInterAreaEVREV
            output = testCase.ExplainedVarianceReactivationObj.calculateInterAreaEVREV('dataPath', testCase.TestDataDir);
            
            % Verify output structure
            testCase.verifyNotEmpty(output);
            expectedFields = {'CA1_vStr', 'perCellContributions'};
            actualFields = fieldnames(output);
            for field = expectedFields
                testCase.verifyEqual(ismember(field{1}, actualFields), true, sprintf('Expected field "%s" missing.', field{1}));
                testCase.verifyEqual(ismember((field{1}), actualFields), true, sprintf('Field "%s" expected to be a struct.', field{1}));
            end
            expectedSubfields = {'EV', 'REV'};
            for subfield = expectedSubfields
                testCase.verifyEqual(isfield(output.CA1_vStr, subfield{1}), true, sprintf('Expected subfield "%s" missing.', subfield{1}));
                testCase.verifyEqual(isfield(output.CA1_vStr.(subfield{1}), 'A_rat'), true, sprintf('Expected subfield "A_rat" missing from field "%s".', subfield{1}));
                testCase.verifyEqual(length(output.CA1_vStr.(subfield{1}).A_rat), 1, sprintf('Expected 1 %s value for A_rat.', subfield{1}));
            end
            
            % Check values of EV and REV
            testCase.verifyNotEqual(output.CA1_vStr.EV.A_rat(1), EV, 'Returned value of EV was expected to be different from the value calculated using all neurons in the session');
            testCase.verifyNotEqual(output.CA1_vStr.REV.A_rat(1), REV, 'Returned value of REV was expected to be different from the value calculated using all neurons in the session');
            
            % Check structure of perCellContributions
            assert(output.CA1_vStr.EV.A_rat(1) > output.CA1_vStr.REV.A_rat(1), 'EV-REV should be positive for this test case')
            testCase.verifyEqual(isstruct(output.perCellContributions.A_rat{1}), true, 'Expected the per-cell contributions for session 1 to be a struct.')
            testCase.verifyEqual(length(output.perCellContributions.A_rat{1}), nNAcUnits*nCA1Units, 'Expected the per-cell contributions to match the number of cell pairs.')
            expectedSubfields = {'cells', 'contribution', 'correlations'};
            for subfield = expectedSubfields
                testCase.verifyEqual(isfield(output.perCellContributions.A_rat{1}, subfield{1}), true, sprintf('Expected subfield "%s" missing from per-cell contributions.', subfield{1}));
            end
            
        end
                
        % Test that calculateInterAreaEVREV returns different values depending on the ripple mode
        function testCalculateInterAreaEVREV_differentRippleModes(testCase)
            
            % Mock inclusion criteria and sessions: 1 rat with 1 session
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rewardResponsive = false;
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.minimumFiringRate = 0;
            testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
            testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;
            
            % Create mock behavioural data: 1 session with 1 rewarded trial
            behavDirectory = fullfile(testCase.TestDataDir, 'test_session');
            ratName = 'A_Rat';
            TestUtilities.saveMockBehaviouralData(behavDirectory, ratName);
            
            % Mock the rest of the data: 3 vStr and 3 CA1 neurons and EV<REV
            ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
            nNAcUnits = 3;
            nCA1Units = 3;
            seed = 1;
            [spiketimes, binnedfr, startEndTimes, EV, REV, rippleStart, rippleStop] = testCase.createSpiketimes('interarea', true, 'ripples', true, 'nNAcUnits', nNAcUnits, 'nCA1Units', nCA1Units, 'seed', seed);
            assert(EV<REV, 'Mock data EV-REV using ripple times should be negative for testing in this function')
            ripples3std = struct('rippleStart', rippleStart, 'rippleStop', rippleStop, 'ripplePeak', rippleStart);
            TestUtilities.saveMockEphysData(ephysDirectory,...
                'spiketimes', spiketimes,...
                'startEndTimes', startEndTimes,...
                'nNAcUnits', nNAcUnits,...
                'binnedfr', binnedfr,...
                'ripples3std', ripples3std);
            
            % Call calculateInterAreaEVREV with rippleActivityMode "none"
            testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'none';
            output_noRipples = testCase.ExplainedVarianceReactivationObj.calculateInterAreaEVREV('dataPath', testCase.TestDataDir);
            testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'variable';
            output_variableRipples = testCase.ExplainedVarianceReactivationObj.calculateInterAreaEVREV('dataPath', testCase.TestDataDir);
            testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'equalise';
            output_equalisedRipples = testCase.ExplainedVarianceReactivationObj.calculateInterAreaEVREV('dataPath', testCase.TestDataDir);
            
            % Verify output structure
            testCase.verifyNotEmpty(output_noRipples);
            testCase.verifyNotEmpty(output_variableRipples);
            testCase.verifyNotEmpty(output_equalisedRipples);
            expectedFields = {'CA1_vStr', 'perCellContributions'};
            for field = expectedFields
                testCase.verifyEqual(ismember(field{1}, fieldnames(output_noRipples)), true, sprintf('Expected field "%s" missing from output ignoring ripple times.', field{1}));
                testCase.verifyEqual(ismember((field{1}), fieldnames(output_noRipples)), true, sprintf('Field "%s" expected to be a struct (output ignoring ripple times).', field{1}));
                testCase.verifyEqual(ismember(field{1}, fieldnames(output_variableRipples)), true, sprintf('Expected field "%s" missing from output with variable ripple duration.', field{1}));
                testCase.verifyEqual(ismember((field{1}), fieldnames(output_variableRipples)), true, sprintf('Field "%s" expected to be a struct (output with variable ripple duration).', field{1}));
                testCase.verifyEqual(ismember(field{1}, fieldnames(output_equalisedRipples)), true, sprintf('Expected field "%s" missing from output with equalised ripple activity.', field{1}));
                testCase.verifyEqual(ismember((field{1}), fieldnames(output_equalisedRipples)), true, sprintf('Field "%s" expected to be a struct (output with equalised ripple activity).', field{1}));
            end
            expectedSubfields = {'EV', 'REV'};
            for subfield = expectedSubfields
                for output = {output_noRipples, output_variableRipples, output_equalisedRipples}
                    testCase.verifyEqual(isfield(output{1}.CA1_vStr, subfield{1}), true, sprintf('Expected subfield "%s" missing.', subfield{1}));
                    testCase.verifyEqual(isfield(output{1}.CA1_vStr.(subfield{1}), 'A_rat'), true, sprintf('Expected subfield "A_rat" missing from field "%s".', subfield{1}));
                    testCase.verifyEqual(length(output{1}.CA1_vStr.(subfield{1}).A_rat), 1, sprintf('Expected 1 %s value for A_rat.', subfield{1}));
                end
            end
            
            % Verify that returned EV and REV values are different
            testCase.verifyNotEqual(output_noRipples.CA1_vStr.EV.A_rat, output_variableRipples.CA1_vStr.EV.A_rat, 'Function calculated the same explained variance with rippleActivityMode = "none" or "variable".')
            testCase.verifyNotEqual(output_noRipples.CA1_vStr.REV.A_rat, output_variableRipples.CA1_vStr.REV.A_rat, 'Function calculated the same reverse explained variance with rippleActivityMode = "none" or "variable".')
            testCase.verifyNotEqual(output_noRipples.CA1_vStr.EV.A_rat, output_equalisedRipples.CA1_vStr.EV.A_rat, 'Function calculated the same explained variance with rippleActivityMode = "none" or "equalise".')
            testCase.verifyNotEqual(output_noRipples.CA1_vStr.REV.A_rat, output_equalisedRipples.CA1_vStr.REV.A_rat, 'Function calculated the same reverse explained variance with rippleActivityMode = "none" or "equalise".')
            testCase.verifyNotEqual(output_variableRipples.CA1_vStr.EV.A_rat, output_equalisedRipples.CA1_vStr.EV.A_rat, 'Function calculated the same explained variance with rippleActivityMode = "variable" or "equalise".')
            testCase.verifyNotEqual(output_variableRipples.CA1_vStr.REV.A_rat, output_equalisedRipples.CA1_vStr.REV.A_rat, 'Function calculated the same reverse explained variance with rippleActivityMode = "variable" or "equalise".')

        end
       
        % Test that calculateInterAreaEVREV returns different values
        % depending on the number of ripples
        function testCalculateInterAreaEVREV_nRipples(testCase)
            
            % Mock inclusion criteria and sessions: 1 rat with 1 session
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rewardResponsive = false;
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.minimumFiringRate = 0;
            testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
            testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;
            
            % Mock analysis configuration
            testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'variable';
            
            % Create mock behavioural data
            behavDirectory = fullfile(testCase.TestDataDir, 'test_session');
            ratName = 'A_Rat';
            TestUtilities.saveMockBehaviouralData(behavDirectory, ratName);
            
            % Mock the rest of the data: 3 vStr and 3 CA1 neurons and EV<REV
            ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
            nNAcUnits = 3;
            nCA1Units = 3;
            seed = 1;
            [spiketimes, binnedfr, startEndTimes, EV, REV, rippleStart, rippleStop] = testCase.createSpiketimes('interarea', true, 'ripples', true, 'nNAcUnits', nNAcUnits, 'nCA1Units', nCA1Units, 'seed', seed);
            assert(EV<REV, 'Mock data EV-REV using ripple times should be negative for testing in this function')
            ripples3std = struct('rippleStart', rippleStart, 'rippleStop', rippleStop, 'ripplePeak', rippleStart);
            TestUtilities.saveMockEphysData(ephysDirectory,...
                'spiketimes', spiketimes,...
                'startEndTimes', startEndTimes,...
                'nNAcUnits', nNAcUnits,...
                'binnedfr', binnedfr,...
                'ripples3std', ripples3std);
            
            % Call calculateInterAreaEVREV with rippleActivityMode "none", 1, 2, and 3
            testCase.ExplainedVarianceReactivationObj.analysisConfig.nRipples = NaN;
            output_noNRipples = testCase.ExplainedVarianceReactivationObj.calculateInterAreaEVREV('dataPath', testCase.TestDataDir);
            testCase.ExplainedVarianceReactivationObj.analysisConfig.nRipples = 1;
            output_1Ripple = testCase.ExplainedVarianceReactivationObj.calculateInterAreaEVREV('dataPath', testCase.TestDataDir);
            testCase.ExplainedVarianceReactivationObj.analysisConfig.nRipples = 2;
            output_2Ripples = testCase.ExplainedVarianceReactivationObj.calculateInterAreaEVREV('dataPath', testCase.TestDataDir);
            testCase.ExplainedVarianceReactivationObj.analysisConfig.nRipples = 3;
            output_3Ripples = testCase.ExplainedVarianceReactivationObj.calculateInterAreaEVREV('dataPath', testCase.TestDataDir);
            
            % Verify output structure
            testCase.verifyNotEmpty(output_noNRipples);
            testCase.verifyNotEmpty(output_1Ripple);
            testCase.verifyNotEmpty(output_2Ripples);
            testCase.verifyNotEmpty(output_3Ripples);
            expectedFields = {'CA1_vStr', 'perCellContributions'};
            for field = expectedFields
                testCase.verifyEqual(ismember(field{1}, fieldnames(output_noNRipples)), true, sprintf('Expected field "%s" missing from output with nRipples = NaN.', field{1}));
                testCase.verifyEqual(ismember((field{1}), fieldnames(output_noNRipples)), true, sprintf('Field "%s" expected to be a struct with nRipples = NaN.', field{1}));
                testCase.verifyEqual(ismember(field{1}, fieldnames(output_1Ripple)), true, sprintf('Expected field "%s" missing from output with nRipples = 1.', field{1}));
                testCase.verifyEqual(ismember((field{1}), fieldnames(output_1Ripple)), true, sprintf('Field "%s" expected to be a struct with nRipples = 1.', field{1}));
                testCase.verifyEqual(ismember(field{1}, fieldnames(output_2Ripples)), true, sprintf('Expected field "%s" missing from output with nRipples = 2.', field{1}));
                testCase.verifyEqual(ismember((field{1}), fieldnames(output_2Ripples)), true, sprintf('Field "%s" expected to be a struct output with nRipples = 2.', field{1}));
                testCase.verifyEqual(ismember(field{1}, fieldnames(output_3Ripples)), true, sprintf('Expected field "%s" missing from output with nRipples = 3.', field{1}));
                testCase.verifyEqual(ismember((field{1}), fieldnames(output_3Ripples)), true, sprintf('Field "%s" expected to be a struct output with nRipples = 3.', field{1}));
            end
            expectedSubfields = {'EV', 'REV'};
            for subfield = expectedSubfields
                for output = {output_noNRipples, output_1Ripple, output_2Ripples, output_3Ripples}
                    testCase.verifyEqual(isfield(output{1}.CA1_vStr, subfield{1}), true, sprintf('Expected subfield "%s" missing.', subfield{1}));
                    testCase.verifyEqual(isfield(output{1}.CA1_vStr.(subfield{1}), 'A_rat'), true, sprintf('Expected subfield "A_rat" missing from field "%s".', subfield{1}));
                    testCase.verifyEqual(length(output{1}.CA1_vStr.(subfield{1}).A_rat), 1, sprintf('Expected 1 %s value for A_rat.', subfield{1}));
                end
            end
            
            % Verify that returned EV and REV values are different for different numbers of ripples
            testCase.verifyNotEqual(output_noNRipples.CA1_vStr.EV.A_rat, output_1Ripple.CA1_vStr.EV.A_rat, 'Function calculated the same explained variance with nRipples = NaN or 1.')
            testCase.verifyNotEqual(output_noNRipples.CA1_vStr.REV.A_rat, output_1Ripple.CA1_vStr.REV.A_rat, 'Function calculated the same reverse explained variance with nRipples = NaN or 1.')
            testCase.verifyNotEqual(output_noNRipples.CA1_vStr.EV.A_rat, output_2Ripples.CA1_vStr.EV.A_rat, 'Function calculated the same explained variance with nRipples = NaN or 2.')
            testCase.verifyNotEqual(output_noNRipples.CA1_vStr.REV.A_rat, output_2Ripples.CA1_vStr.REV.A_rat, 'Function calculated the same reverse explained variance with nRipples = NaN or 2.')
            testCase.verifyNotEqual(output_1Ripple.CA1_vStr.EV.A_rat, output_2Ripples.CA1_vStr.EV.A_rat, 'Function calculated the same explained variance with nRipples = 1 or 2.')
            testCase.verifyNotEqual(output_1Ripple.CA1_vStr.REV.A_rat, output_2Ripples.CA1_vStr.REV.A_rat, 'Function calculated the same reverse explained variance with nRipples = 1 or 2.')
            
            % Verify that returned EV and REV values are the same for NaN and 3 ripples
            testCase.verifyEqual(output_noNRipples.CA1_vStr.EV.A_rat, output_3Ripples.CA1_vStr.EV.A_rat, 'RelTol', 1e-5, 'Function calculated different explained variance with nRipples = NaN or 3, for a total of 3 ripples.')
            testCase.verifyEqual(output_noNRipples.CA1_vStr.REV.A_rat, output_3Ripples.CA1_vStr.REV.A_rat, 'RelTol', 1e-5, 'Function calculated different reverseexplained variance with nRipples = NaN or 3, for a total of 3 ripples.')

        end
        
        % Test that calculateInterAreaEVREV returns different values
        % depending on the number of ripples (separate number for PRE and POST)
        function testCalculateInterAreaEVREV_nRipplesPrePost(testCase)
            
            % Mock inclusion criteria and sessions: 1 rat with 1 session
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rewardResponsive = false;
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.minimumFiringRate = 0;
            testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
            testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;
            
            % Mock analysis configuration
            testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'variable';
            
            % Create mock behavioural data
            behavDirectory = fullfile(testCase.TestDataDir, 'test_session');
            ratName = 'A_Rat';
            TestUtilities.saveMockBehaviouralData(behavDirectory, ratName);
            
            % Mock the rest of the data: 3 vStr and 3 CA1 neurons and EV<REV
            ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
            nNAcUnits = 3;
            nCA1Units = 3;
            seed = 1;
            [spiketimes, binnedfr, startEndTimes, EV, REV, rippleStart, rippleStop] = testCase.createSpiketimes('interarea', true, 'ripples', true, 'nNAcUnits', nNAcUnits, 'nCA1Units', nCA1Units, 'seed', seed);
            assert(EV<REV, 'Mock data EV-REV using ripple times should be negative for testing in this function')
            ripples3std = struct('rippleStart', rippleStart, 'rippleStop', rippleStop, 'ripplePeak', rippleStart);
            TestUtilities.saveMockEphysData(ephysDirectory,...
                'spiketimes', spiketimes,...
                'startEndTimes', startEndTimes,...
                'nNAcUnits', nNAcUnits,...
                'binnedfr', binnedfr,...
                'ripples3std', ripples3std);
            
            % Call calculateInterAreaEVREV with rippleActivityMode "none", 1, 2, and 3
            testCase.ExplainedVarianceReactivationObj.analysisConfig.nRipples = NaN;
            output_noNRipples = testCase.ExplainedVarianceReactivationObj.calculateInterAreaEVREV('dataPath', testCase.TestDataDir);
            testCase.ExplainedVarianceReactivationObj.analysisConfig.nRipples = [1, 1];
            output_1Ripple = testCase.ExplainedVarianceReactivationObj.calculateInterAreaEVREV('dataPath', testCase.TestDataDir);
            testCase.ExplainedVarianceReactivationObj.analysisConfig.nRipples = [2, 1];
            output_2Ripples = testCase.ExplainedVarianceReactivationObj.calculateInterAreaEVREV('dataPath', testCase.TestDataDir);
            testCase.ExplainedVarianceReactivationObj.analysisConfig.nRipples = [3, 3];
            output_3Ripples = testCase.ExplainedVarianceReactivationObj.calculateInterAreaEVREV('dataPath', testCase.TestDataDir);
            
            % Verify output structure
            testCase.verifyNotEmpty(output_noNRipples);
            testCase.verifyNotEmpty(output_1Ripple);
            testCase.verifyNotEmpty(output_2Ripples);
            testCase.verifyNotEmpty(output_3Ripples);
            expectedFields = {'CA1_vStr', 'perCellContributions'};
            for field = expectedFields
                testCase.verifyEqual(ismember(field{1}, fieldnames(output_noNRipples)), true, sprintf('Expected field "%s" missing from output with nRipples = NaN.', field{1}));
                testCase.verifyEqual(ismember((field{1}), fieldnames(output_noNRipples)), true, sprintf('Field "%s" expected to be a struct with nRipples = NaN.', field{1}));
                testCase.verifyEqual(ismember(field{1}, fieldnames(output_1Ripple)), true, sprintf('Expected field "%s" missing from output with nRipples = 1.', field{1}));
                testCase.verifyEqual(ismember((field{1}), fieldnames(output_1Ripple)), true, sprintf('Field "%s" expected to be a struct with nRipples = 1.', field{1}));
                testCase.verifyEqual(ismember(field{1}, fieldnames(output_2Ripples)), true, sprintf('Expected field "%s" missing from output with nRipples = 2.', field{1}));
                testCase.verifyEqual(ismember((field{1}), fieldnames(output_2Ripples)), true, sprintf('Field "%s" expected to be a struct output with nRipples = 2.', field{1}));
                testCase.verifyEqual(ismember(field{1}, fieldnames(output_3Ripples)), true, sprintf('Expected field "%s" missing from output with nRipples = 3.', field{1}));
                testCase.verifyEqual(ismember((field{1}), fieldnames(output_3Ripples)), true, sprintf('Field "%s" expected to be a struct output with nRipples = 3.', field{1}));
            end
            expectedSubfields = {'EV', 'REV'};
            for subfield = expectedSubfields
                for output = {output_noNRipples, output_1Ripple, output_2Ripples, output_3Ripples}
                    testCase.verifyEqual(isfield(output{1}.CA1_vStr, subfield{1}), true, sprintf('Expected subfield "%s" missing.', subfield{1}));
                    testCase.verifyEqual(isfield(output{1}.CA1_vStr.(subfield{1}), 'A_rat'), true, sprintf('Expected subfield "A_rat" missing from field "%s".', subfield{1}));
                    testCase.verifyEqual(length(output{1}.CA1_vStr.(subfield{1}).A_rat), 1, sprintf('Expected 1 %s value for A_rat.', subfield{1}));
                end
            end
            
            % Verify that returned EV and REV values are different for different numbers of ripples
            testCase.verifyNotEqual(output_noNRipples.CA1_vStr.EV.A_rat, output_1Ripple.CA1_vStr.EV.A_rat, 'Function calculated the same explained variance with nRipples = NaN or 1.')
            testCase.verifyNotEqual(output_noNRipples.CA1_vStr.REV.A_rat, output_1Ripple.CA1_vStr.REV.A_rat, 'Function calculated the same reverse explained variance with nRipples = NaN or 1.')
            testCase.verifyNotEqual(output_noNRipples.CA1_vStr.EV.A_rat, output_2Ripples.CA1_vStr.EV.A_rat, 'Function calculated the same explained variance with nRipples = NaN or 2.')
            testCase.verifyNotEqual(output_noNRipples.CA1_vStr.REV.A_rat, output_2Ripples.CA1_vStr.REV.A_rat, 'Function calculated the same reverse explained variance with nRipples = NaN or 2.')
            testCase.verifyNotEqual(output_1Ripple.CA1_vStr.EV.A_rat, output_2Ripples.CA1_vStr.EV.A_rat, 'Function calculated the same explained variance with nRipples = 1 or 2.')
            testCase.verifyNotEqual(output_1Ripple.CA1_vStr.REV.A_rat, output_2Ripples.CA1_vStr.REV.A_rat, 'Function calculated the same reverse explained variance with nRipples = 1 or 2.')
            
            % Verify that returned EV and REV values are the same for NaN and 3 ripples
            testCase.verifyEqual(output_noNRipples.CA1_vStr.EV.A_rat, output_3Ripples.CA1_vStr.EV.A_rat, 'RelTol', 1e-5, 'Function calculated different explained variance with nRipples = NaN or 3, for a total of 3 ripples.')
            testCase.verifyEqual(output_noNRipples.CA1_vStr.REV.A_rat, output_3Ripples.CA1_vStr.REV.A_rat, 'RelTol', 1e-5, 'Function calculated different reverseexplained variance with nRipples = NaN or 3, for a total of 3 ripples.')

        end
        
        % Test error handling of calculateInterAreaEVREV when there are no
        % ripples
        function testCalculateInterAreaEVREV_noRipples(testCase)
            
            % Mock inclusion criteria and sessions: 1 rat with 1 session
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rewardResponsive = false;
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.minimumFiringRate = 0;
            testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
            testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;
            
            % Mock analysis configuration
            testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'variable';
            testCase.ExplainedVarianceReactivationObj.analysisConfig.nRipples = NaN;
            
            % Create mock behavioural data
            behavDirectory = fullfile(testCase.TestDataDir, 'test_session');
            ratName = 'A_Rat';
            TestUtilities.saveMockBehaviouralData(behavDirectory, ratName);
            
            % Mock the rest of the data: 3 vStr and 3 CA1 neurons and EV<REV
            ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
            nNAcUnits = 3;
            nCA1Units = 3;
            seed = 1;
            [spiketimes, binnedfr, startEndTimes, EV, REV, ~, ~] = testCase.createSpiketimes('interarea', true, 'ripples', false, 'nNAcUnits', nNAcUnits, 'nCA1Units', nCA1Units, 'seed', seed);
            assert(EV<REV, 'Mock data EV-REV using ripple times should be negative for testing in this function')
            ripples3std = struct('rippleStart', [], 'rippleStop', [], 'ripplePeak', []);
            TestUtilities.saveMockEphysData(ephysDirectory,...
                'spiketimes', spiketimes,...
                'startEndTimes', startEndTimes,...
                'nNAcUnits', nNAcUnits,...
                'binnedfr', binnedfr,...
                'ripples3std', ripples3std);
            
            % Call calculateInterAreaEVREV with rippleActivityMode "none", 1, 2, and 3
            expectedErrorId = 'ExplainedVarianceReactivation:invalidData';
            testCase.verifyError(@() testCase.ExplainedVarianceReactivationObj.calculateInterAreaEVREV('dataPath', testCase.TestDataDir), expectedErrorId);
                        
        end
        
        % Test calculateInterAreaEVREV when no rats are specified
        function testCalculateInterAreaEVREV_noRats(testCase)
            
            % Mock inclusion criteria and sessions: 1 rat with 1 session
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {};
            
            % Call calculateInterAreaEVREV
            expectedErrorId = 'ExplainedVarianceReactivation:sessions:noSessions';
            testCase.verifyError(@() testCase.ExplainedVarianceReactivationObj.calculateInterAreaEVREV(), expectedErrorId);
             
        end
        
         % Test error handling for calculateInterAreaEVREV with a missing rat
        function testCalculateInterAreaEVREV_missingRat(testCase)
            
            % Mock inclusion criteria and sessions: 1 rat with 1 session
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'missing_rat'};
            
            % Call calculateInterAreaEVREV
            expectedErrorId = 'ExplainedVarianceReactivation:sessions:noSessions';
            testCase.verifyError(@() testCase.ExplainedVarianceReactivationObj.calculateInterAreaEVREV(), expectedErrorId);
            
        end
        
        % Test that calculateInterAreaEVREV skips sessions that are not to be included
        function testCalculateInterAreaEVREV_sessionSelection(testCase)
            
            % Mock inclusion criteria and sessions: 1 rat with 1 session
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rewardResponsive = false;
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.minimumFiringRate = 0;
            testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [2, 3]);
            testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;
            
            % Mock analysis configuration
            testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'none';
            
            % Create mock behavioural data
            behavDirectory = fullfile(testCase.TestDataDir, 'test_session');
            ratName = 'A_Rat';
            TestUtilities.saveMockBehaviouralData(behavDirectory, ratName);
            
            % Mock ephys data for session 1: 3 vStr and 3 CA1 neurons and EV>REV
            ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
            nNAcUnits = 3;
            nCA1Units = 3;
            [spiketimes, binnedfr, startEndTimes, EV, REV] = testCase.createSpiketimes('interarea', true, 'nNAcUnits', nNAcUnits, 'nCA1Units', nCA1Units, 'seed', 2);
            assert(EV>REV, 'Mock data EV-REV should be positive for simpler testing this function')
            TestUtilities.saveMockEphysData(ephysDirectory,...
                'spiketimes', spiketimes,...
                'startEndTimes', startEndTimes,...
                'nNAcUnits', nNAcUnits,...
                'binnedfr', binnedfr);
            
            % Mock ephys data for session 1: 3 vStr and 4 CA1 neurons and EV>REV
            ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session2');
            nNAcUnits_session2 = 3;
            nCA1Units_session2 = 4;
            [spiketimes, binnedfr, startEndTimes, EV, REV] = testCase.createSpiketimes('interarea', true, 'nNAcUnits', nNAcUnits_session2, 'nCA1Units', nCA1Units_session2, 'seed', 3);
            assert(EV>REV, 'Mock data EV-REV should be positive for testing in this function')
            TestUtilities.saveMockEphysData(ephysDirectory,...
                'spiketimes', spiketimes,...
                'startEndTimes', startEndTimes,...
                'nNAcUnits', nNAcUnits_session2,...
                'binnedfr', binnedfr);
            
            % Mock ephys data for session 1: 4 vStr and 5 CA1 neurons and EV>REV
            ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session3');
            nNAcUnits_session3 = 4;
            nCA1Units_session3 = 5;
            [spiketimes, binnedfr, startEndTimes, EV, REV] = testCase.createSpiketimes('interarea', true, 'nNAcUnits', nNAcUnits_session3, 'nCA1Units', nCA1Units_session3, 'seed', 5);
            assert(EV>REV, 'Mock data EV-REV should be positive for testing in this function')
            TestUtilities.saveMockEphysData(ephysDirectory,...
                'spiketimes', spiketimes,...
                'startEndTimes', startEndTimes,...
                'nNAcUnits', nNAcUnits_session3,...
                'binnedfr', binnedfr);
            
            % Call calculateInterAreaEVREV
            output = testCase.ExplainedVarianceReactivationObj.calculateInterAreaEVREV('dataPath', testCase.TestDataDir);
                        
            % Verify output structure
            testCase.verifyNotEmpty(output);
            expectedFields = {'CA1_vStr', 'perCellContributions'};
            actualFields = fieldnames(output);
            for field = expectedFields
                testCase.verifyEqual(ismember(field{1}, actualFields), true, sprintf('Expected field "%s" missing.', field{1}));
                testCase.verifyEqual(ismember((field{1}), actualFields), true, sprintf('Field "%s" expected to be a struct.', field{1}));
            end
            expectedSubfields = {'EV', 'REV'};
            for subfield = expectedSubfields
                testCase.verifyEqual(isfield(output.CA1_vStr, subfield{1}), true, sprintf('Expected subfield "%s" missing.', subfield{1}));
                testCase.verifyEqual(isfield(output.CA1_vStr.(subfield{1}), 'A_rat'), true, sprintf('Expected subfield "A_rat" missing from field "%s".', subfield{1}));
                testCase.verifyEqual(length(output.CA1_vStr.(subfield{1}).A_rat), 2, sprintf('Expected 1 %s value for A_rat.', subfield{1}));
            end
            testCase.verifyEmpty(output.perCellContributions.A_rat{1}, 'Expected the per-cell contributions for session 1 to be empty.')
            testCase.verifyEqual(isstruct(output.perCellContributions.A_rat{2}), true, 'Expected the per-cell contributions for session 2 to be a struct.')
            testCase.verifyEqual(isstruct(output.perCellContributions.A_rat{3}), true, 'Expected the per-cell contributions for session 3 to be a struct.')
            testCase.verifyEqual(length(output.perCellContributions.A_rat{2}), nNAcUnits_session2*nCA1Units_session2, 'Expected the per-cell contributions for session 2 to match the number of cell pairs.')
            testCase.verifyEqual(length(output.perCellContributions.A_rat{3}), nNAcUnits_session3*nCA1Units_session3, 'Expected the per-cell contributions for session 3 to match the number of cell pairs.')
            expectedSubfields = {'cells', 'contribution', 'correlations'};
            for subfield = expectedSubfields
                for i = 2:3
                    testCase.verifyEqual(isfield(output.perCellContributions.A_rat{2}(i), subfield{1}), true, sprintf('Expected subfield "%s" missing from per-cell contributions for session %i.', subfield{1}, i));
                end
            end
            
        end
        
        
        % ================ reactivation within brain areas ================ %
        
        % Test basic functionality
        function testCalculateIntraAreaEVREV_basicFunctionality(testCase)
            
            functions = {...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateVStrEVREV('dataPath', x),...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateCA1EVREV('dataPath', x)};
            relevantNeurons = {'nNAcUnits', 'nCA1Units'};
            irrelevantNeurons = {'nCA1Units', 'nNAcUnits'};
            outputField = {'vStr', 'CA1'};
            
            for iFunc = 1:numel(functions)
                
                % Mock inclusion criteria and sessions: 1 rat with 1 session
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rewardResponsive = false;
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.minimumFiringRate = 0;
                testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
                testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;

                % Mock analysis configuration
                testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'none';

                % Create mock behavioural data
                behavDirectory = fullfile(testCase.TestDataDir, 'test_session');
                ratName = 'A_Rat';
                TestUtilities.saveMockBehaviouralData(behavDirectory, ratName);

                % Mock the rest of the data: 3 relevant and 0 irrelevant neurons
                ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
                nRelevantNeurons = 3;
                nIrrelevantNeurons = 0;
                [spiketimes, binnedfr, startEndTimes, EV, REV] = testCase.createSpiketimes('interarea', false, relevantNeurons{iFunc}, nRelevantNeurons, irrelevantNeurons{iFunc}, nIrrelevantNeurons, 'seed', 0);
                assert(EV<REV, 'Mock data EV-REV should be negative for simpler testing in this function')
                nNAcUnits = testCase.getIntraAreaNAc(relevantNeurons, iFunc, nRelevantNeurons, nIrrelevantNeurons);
                TestUtilities.saveMockEphysData(ephysDirectory,...
                    'spiketimes', spiketimes,...
                    'startEndTimes', startEndTimes,...
                    'nNAcUnits', nNAcUnits,...
                    'binnedfr', binnedfr);
            
                % Call method
                func = functions{iFunc};
                output = func(testCase.TestDataDir);
            
                % Verify output structure
                testCase.verifyNotEmpty(output);
                expectedFields = outputField(iFunc);
                actualFields = fieldnames(output);
                for field = expectedFields
                    testCase.verifyEqual(ismember(field{1}, actualFields), true, sprintf('Expected field "%s" missing.', field{1}));
                    testCase.verifyEqual(ismember((field{1}), actualFields), true, sprintf('Field "%s" expected to be a struct.', field{1}));
                end
                expectedSubfields = {'EV', 'REV'};
                for subfield = expectedSubfields
                    testCase.verifyEqual(isfield(output.(outputField{iFunc}), subfield{1}), true, sprintf('Expected subfield "%s" missing from field "%s".', subfield{1}, outputField{iFunc}));
                    testCase.verifyEqual(isfield(output.(outputField{iFunc}).(subfield{1}), 'A_rat'), true, sprintf('Expected subfield "A_rat" missing from field "%s" (%s).', subfield{1}, outputField{iFunc}));
                    testCase.verifyEqual(length(output.(outputField{iFunc}).(subfield{1}).A_rat), 1, sprintf('Expected 1 %s value for A_rat.', subfield{1}));
                end

                % Verify returned EV and REV values
                testCase.verifyEqual(output.(outputField{iFunc}).EV.A_rat, EV, 'AbsTol', 1e-3, sprintf('Returned value of EV does not match the expected value (%s).', outputField{iFunc}));
                testCase.verifyEqual(output.(outputField{iFunc}).REV.A_rat, REV, 'AbsTol', 1e-3, sprintf('Returned value of REV does not match the expected value (%s).', outputField{iFunc}));
                
            end
            
        end
        
        % Test ignoring spiking/firing data from irrelevant brain area:
        % including or excluding spiketimes and binned firing rates from the other
        % brain area should make no difference to the output
        function testCalculateIntraAreaEVREV_relevantBrainAreaOnly(testCase)
            
            functions = {...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateVStrEVREV('dataPath', x),...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateCA1EVREV('dataPath', x)};
            relevantNeurons = {'nNAcUnits', 'nCA1Units'};
            outputField = {'vStr', 'CA1'};
            
            for iFunc = 1:numel(functions)
            
                % Mock inclusion criteria and sessions: 1 rat with 1 session
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rewardResponsive = false;
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.minimumFiringRate = 0;
                testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
                testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;
            
                % Mock analysis configuration
                testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'none';
            
                % Create mock behavioural data
                behavDirectory = fullfile(testCase.TestDataDir, 'test_session');
                ratName = 'A_Rat';
                TestUtilities.saveMockBehaviouralData(behavDirectory, ratName);
            
                % Mock the rest of the data: 5 relevant and 5 irrelevant neurons
                nNAcUnits = 5;
                nCA1Units = 5;
                [spiketimes_NAc, binnedfr_NAc, startEndTimes, ~, ~] = testCase.createSpiketimes('interarea', false, 'nNAcUnits', nNAcUnits, 'nCA1Units', 0, 'seed', 0);
                [spiketimes_CA1, binnedfr_CA1, ~, ~, ~] = testCase.createSpiketimes('interarea', false, 'nNAcUnits', 0, 'nCA1Units', nCA1Units, 'seed', 1);
                spiketimes_combined = [spiketimes_NAc, spiketimes_CA1];
                binnedfr_combined = [binnedfr_NAc; binnedfr_CA1];
                ephysDirectory_singleAreaData = fullfile(testCase.TestDataDir, 'singleAreaOnly', 'Rat_A', 'Session1');
                ephysDirectory_combined = fullfile(testCase.TestDataDir, 'combined', 'Rat_A', 'Session1');
                if strcmp(relevantNeurons{iFunc}, 'nNAcUnits')
                    spiketimes_singleArea = spiketimes_NAc;
                    binnedfr_singleArea = binnedfr_NAc;
                    nNAcUnits_singleArea = nNAcUnits;
                elseif strcmp(relevantNeurons{iFunc}, 'nCA1Units')
                    spiketimes_singleArea = spiketimes_CA1;
                    binnedfr_singleArea = binnedfr_CA1;
                    nNAcUnits_singleArea = 0;
                else
                    error('Unexpected value of relevantNeurons: %s', relevantNeurons{iFun})
                end
                TestUtilities.saveMockEphysData(ephysDirectory_singleAreaData,...
                    'spiketimes', spiketimes_singleArea,...
                    'startEndTimes', startEndTimes,...
                    'nNAcUnits', nNAcUnits_singleArea,...
                    'binnedfr', binnedfr_singleArea);
                TestUtilities.saveMockEphysData(ephysDirectory_combined,...
                    'spiketimes', spiketimes_combined,...
                    'startEndTimes', startEndTimes,...
                    'nNAcUnits', nNAcUnits,...
                    'binnedfr', binnedfr_combined);

                % Call method
                func = functions{iFunc};
                output_singleArea = func(fullfile(testCase.TestDataDir, 'singleAreaOnly'));
                output_combined = func(fullfile(testCase.TestDataDir, 'combined'));

                % Verify that the outputs are the same
                expectedSubfields = {'EV', 'REV'};
                for subfield = expectedSubfields
                    testCase.verifyEqual(output_singleArea.(outputField{iFunc}).(subfield{1}), output_combined.(outputField{iFunc}).(subfield{1}), sprintf('Subfield "%s" is not identical between outputs (%s).', subfield{1}, outputField{iFunc}));
                end

            end
            
        end
             
        % Test intra-area EV-REV calculation with one session of positive EV-REV
        % and one session of negative EV-REV
        function testCalculateIntraAreaEVREV_twoSessions(testCase)
            
            functions = {...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateVStrEVREV('dataPath', x),...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateCA1EVREV('dataPath', x)};
            relevantNeurons = {'nNAcUnits', 'nCA1Units'};
            irrelevantNeurons = {'nCA1Units', 'nNAcUnits'};
            outputField = {'vStr', 'CA1'};
            
            for iFunc = 1:numel(functions)
            
                % Mock inclusion criteria and sessions: 1 rat with 1 session
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rewardResponsive = false;
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.minimumFiringRate = 0;
                testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [1, 2]);
                testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;

                % Mock analysis configuration
                testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'none';

                % Create mock behavioural data
                behavDirectory = fullfile(testCase.TestDataDir, 'test_session');
                ratName = 'A_Rat';
                TestUtilities.saveMockBehaviouralData(behavDirectory, ratName);
            
                % Mock ephys data for first session: 3 relevant and 0 irrelevant neurons and EV<REV
                ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
                nRelevantNeurons = 3;
                nIrrelevantNeurons = 0;
                [spiketimes, binnedfr, startEndTimes, EV_session1, REV_session1] = testCase.createSpiketimes('interarea', false, relevantNeurons{iFunc}, nRelevantNeurons, irrelevantNeurons{iFunc}, nIrrelevantNeurons, 'seed', 0);
                assert(EV_session1<REV_session1, 'Mock data EV-REV should be negative for simpler testing in this function')
                nNAcUnits = testCase.getIntraAreaNAc(relevantNeurons, iFunc, nRelevantNeurons, nIrrelevantNeurons);
                TestUtilities.saveMockEphysData(ephysDirectory,...
                    'spiketimes', spiketimes,...
                    'startEndTimes', startEndTimes,...
                    'nNAcUnits', nNAcUnits,...
                    'binnedfr', binnedfr);
            
                % Mock ephys data for second session: 3 vStr and 3 CA1 neurons and EV<REV
                ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session2');
                nNAcUnits = 3;
                nCA1Units = 3;
                [spiketimesNAc, binnedfrNAc, startEndTimesNAc, EVNAc_session2, REVNAc_session2] = testCase.createSpiketimes('interarea', false, 'nNAcUnits', nNAcUnits, 'nCA1Units', 0, 'seed', 1);
                [spiketimesCA1, binnedfrCA1, startEndTimesCA1, EVCA1_session2, REVCA1_session2] = testCase.createSpiketimes('interarea', false, 'nNAcUnits', 0, 'nCA1Units', nCA1Units, 'seed', 1);
                spiketimes = [spiketimesNAc, spiketimesCA1];
                binnedfr = [binnedfrNAc; binnedfrCA1];
                if strcmp(relevantNeurons{iFunc}, 'nNAcUnits')
                    nNAcUnits = nRelevantNeurons;
                    startEndTimes = startEndTimesNAc;
                    EV_session2 = EVNAc_session2;
                    REV_session2 = REVNAc_session2;
                elseif strcmp(relevantNeurons{iFunc}, 'nCA1Units')
                    nNAcUnits = 0;
                    startEndTimes = startEndTimesCA1;
                    EV_session2 = EVCA1_session2;
                    REV_session2 = REVCA1_session2;
                else
                    error('Unexpected value of relevantNeurons: %s', relevantNeurons{iFun})
                end
                assert(EV_session2>REV_session2, 'Mock data EV-REV should be positive for session 2')
                TestUtilities.saveMockEphysData(ephysDirectory,...
                    'spiketimes', spiketimes,...
                    'startEndTimes', startEndTimes,...
                    'nNAcUnits', nNAcUnits,...
                    'binnedfr', binnedfr);
            
                % Call method
                func = functions{iFunc};
                output = func(testCase.TestDataDir);

                 % Verify output structure
                testCase.verifyNotEmpty(output);
                expectedFields = outputField(iFunc);
                actualFields = fieldnames(output);
                for field = expectedFields
                    testCase.verifyEqual(ismember(field{1}, actualFields), true, sprintf('Expected field "%s" missing.', field{1}));
                    testCase.verifyEqual(ismember((field{1}), actualFields), true, sprintf('Field "%s" expected to be a struct.', field{1}));
                end
                expectedSubfields = {'EV', 'REV'};
                for subfield = expectedSubfields
                    testCase.verifyEqual(isfield(output.(outputField{iFunc}), subfield{1}), true, sprintf('Expected subfield "%s" missing from field "%s".', subfield{1}, outputField{iFunc}));
                    testCase.verifyEqual(isfield(output.(outputField{iFunc}).(subfield{1}), 'A_rat'), true, sprintf('Expected subfield "A_rat" missing from field "%s" (%s).', subfield{1}));
                    testCase.verifyEqual(length(output.(outputField{iFunc}).(subfield{1}).A_rat), 2, sprintf('Expected 1 %s value for A_rat (%s).', subfield{1}, outputField{iFunc}));
                end

                % Verify returned EV and REV values
                for session = 1:2
                    eval(['EV = EV_session' num2str(session) ';']);
                    eval(['REV = REV_session' num2str(session) ';']);
                    testCase.verifyEqual(output.(outputField{iFunc}).EV.A_rat(session), EV, 'AbsTol', 1e-5, sprintf('Returned value of EV does not match the expected value for session %i (%s).', session, outputField{iFunc}));
                    testCase.verifyEqual(output.(outputField{iFunc}).REV.A_rat(session), REV, 'AbsTol', 1e-5, sprintf('Returned value of REV does not match the expected value for session %i (%S).', session, outputField{iFunc}));
                end
            end
        end
       
        % Test that intra-area EV-REV calculation correctly implements the
        % minimum firing rate threshold
        function testCcalculateIntraAreaEVREV_testMinimumFiringRate(testCase)
            
            functions = {...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateVStrEVREV('dataPath', x),...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateCA1EVREV('dataPath', x)};
            relevantNeurons = {'nNAcUnits', 'nCA1Units'};
            irrelevantNeurons = {'nCA1Units', 'nNAcUnits'};
            outputField = {'vStr', 'CA1'};
            
            for iFunc = 1:numel(functions)
            
                % Mock inclusion criteria and sessions: 1 rat with 1 session
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rewardResponsive = false;
                testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
                testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;

                % Mock analysis configuration
                testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'none';

                % Create mock behavioural data
                behavDirectory = fullfile(testCase.TestDataDir, 'test_session');
                ratName = 'A_Rat';
                TestUtilities.saveMockBehaviouralData(behavDirectory, ratName);
            
                % Mock the rest of the data: 5 relevant neurons of which 3 are above the minimum firing rate
                ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
                nRelevantNeurons = 5;
                nIrrelevantNeurons = 0;
                seed = 0;
                [spiketimes, binnedfr, startEndTimes, EV, REV, ~, ~] = testCase.createSpiketimes('interarea', false, 'ripples', false, relevantNeurons{iFunc}, nRelevantNeurons, irrelevantNeurons{iFunc}, nIrrelevantNeurons, 'seed', seed);
                minFR = 15;
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.minimumFiringRate = minFR;
                assert(sum(mean(binnedfr(1:nRelevantNeurons, :), 2) > minFR) == 3, 'Mock data should include three NAc neurons above the minimum firing rate threshold for testing in this function')
                nNAcUnits = testCase.getIntraAreaNAc(relevantNeurons, iFunc, nRelevantNeurons, nIrrelevantNeurons);
                TestUtilities.saveMockEphysData(ephysDirectory,...
                    'spiketimes', spiketimes,...
                    'startEndTimes', startEndTimes,...
                    'nNAcUnits', nNAcUnits,...
                    'binnedfr', binnedfr);

                % Call method
                func = functions{iFunc};
                output = func(testCase.TestDataDir);

                % Verify output structure
                testCase.verifyNotEmpty(output);
                expectedFields = outputField(iFunc);
                actualFields = fieldnames(output);
                for field = expectedFields
                    testCase.verifyEqual(ismember(field{1}, actualFields), true, sprintf('Expected field "%s" missing.', field{1}));
                    testCase.verifyEqual(ismember((field{1}), actualFields), true, sprintf('Field "%s" expected to be a struct.', field{1}));
                end
                expectedSubfields = {'EV', 'REV'};
                for subfield = expectedSubfields
                    testCase.verifyEqual(isfield(output.(outputField{iFunc}), subfield{1}), true, sprintf('Expected subfield "%s" missing from field "%s".', subfield{1}, outputField{iFunc}));
                    testCase.verifyEqual(isfield(output.(outputField{iFunc}).(subfield{1}), 'A_rat'), true, sprintf('Expected subfield "A_rat" missing from field "%s" (%s).', subfield{1}, outputField{iFunc}));
                    testCase.verifyEqual(length(output.(outputField{iFunc}).(subfield{1}).A_rat), 1, sprintf('Expected 1 %s value for A_rat (%s).', subfield{1}, outputField{iFunc}));
                end

                % Check values of EV and REV
                testCase.verifyNotEqual(output.(outputField{iFunc}).EV.A_rat(1), EV, sprintf('Returned value of EV was expected to be different from the value calculated using all neurons in the session (%s).', outputField{iFunc}));
                testCase.verifyNotEqual(output.(outputField{iFunc}).REV.A_rat(1), REV, sprintf('Returned value of REV was expected to be different from the value calculated using all neurons in the session (%s).', outputField{iFunc}));
            
            end
        end
                
        % Test that intra-area EV-REV calculation returns different values depending on the ripple mode
        function testCalculateIntraAreaEVREV_differentRippleModes(testCase)
            
            functions = {...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateVStrEVREV('dataPath', x),...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateCA1EVREV('dataPath', x)};
            relevantNeurons = {'nNAcUnits', 'nCA1Units'};
            irrelevantNeurons = {'nCA1Units', 'nNAcUnits'};
            outputField = {'vStr', 'CA1'};
            
            for iFunc = 1:numel(functions)
            
                % Mock inclusion criteria and sessions: 1 rat with 1 session
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rewardResponsive = false;
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.minimumFiringRate = 0;
                testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
                testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;

                % Create mock behavioural data: 1 session with 1 rewarded trial
                behavDirectory = fullfile(testCase.TestDataDir, 'test_session');
                ratName = 'A_Rat';
                TestUtilities.saveMockBehaviouralData(behavDirectory, ratName);
            
                % Mock the rest of the data: 3 relevant and 0 irrelevant neurons
                ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
                nRelevantNeurons = 3;
                seed = 0;
                [spiketimes, binnedfr, startEndTimes, ~, ~, rippleStart, rippleStop] = testCase.createSpiketimes('interarea', true, 'ripples', true, relevantNeurons{iFunc}, nRelevantNeurons, irrelevantNeurons{iFunc}, 0, 'seed', seed);
                ripples3std = struct('rippleStart', rippleStart, 'rippleStop', rippleStop, 'ripplePeak', rippleStart);
                nNAcUnits = testCase.getIntraAreaNAc(relevantNeurons, iFunc, nRelevantNeurons, 0);
                TestUtilities.saveMockEphysData(ephysDirectory,...
                    'spiketimes', spiketimes,...
                    'startEndTimes', startEndTimes,...
                    'nNAcUnits', nNAcUnits,...
                    'binnedfr', binnedfr,...
                    'ripples3std', ripples3std);
            
                % Call method with rippleActivityMode "none", "variable" and "equalise"
                func = functions{iFunc};
                testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'none';
                output_noRipples = func(testCase.TestDataDir);
                testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'variable';
                output_variableRipples = func(testCase.TestDataDir);
                testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'equalise';
                output_equalisedRipples = func(testCase.TestDataDir);
            
                % Verify output structure
                testCase.verifyNotEmpty(output_noRipples);
                testCase.verifyNotEmpty(output_variableRipples);
                testCase.verifyNotEmpty(output_equalisedRipples);
                expectedFields = outputField(iFunc);
                for field = expectedFields
                    testCase.verifyEqual(ismember(field{1}, fieldnames(output_noRipples)), true, sprintf('Expected field "%s" missing from output ignoring ripple times (%s).', field{1}, outputField{iFunc}));
                    testCase.verifyEqual(ismember((field{1}), fieldnames(output_noRipples)), true, sprintf('Field "%s" expected to be a struct (output ignoring ripple times) (%s).', field{1}, outputField{iFunc}));
                    testCase.verifyEqual(ismember(field{1}, fieldnames(output_variableRipples)), true, sprintf('Expected field "%s" missing from output with variable ripple duration (%s).', field{1}, outputField{iFunc}));
                    testCase.verifyEqual(ismember((field{1}), fieldnames(output_variableRipples)), true, sprintf('Field "%s" expected to be a struct (output with variable ripple duration) (%s).', field{1}, outputField{iFunc}));
                    testCase.verifyEqual(ismember(field{1}, fieldnames(output_equalisedRipples)), true, sprintf('Expected field "%s" missing from output with equalised ripple activity (%s).', field{1}, outputField{iFunc}));
                    testCase.verifyEqual(ismember((field{1}), fieldnames(output_equalisedRipples)), true, sprintf('Field "%s" expected to be a struct (output with equalised ripple activity) (%s).', field{1}, outputField{iFunc}));
                end
                expectedSubfields = {'EV', 'REV'};
                for subfield = expectedSubfields
                    for output = {output_noRipples, output_variableRipples, output_equalisedRipples}
                        testCase.verifyEqual(isfield(output{1}.(outputField{iFunc}), subfield{1}), true, sprintf('Expected subfield "%s" missing (%s).', subfield{1}, outputField{iFunc}));
                        testCase.verifyEqual(isfield(output{1}.(outputField{iFunc}).(subfield{1}), 'A_rat'), true, sprintf('Expected subfield "A_rat" missing from field "%s" (%s).', subfield{1}, outputField{iFunc}));
                        testCase.verifyEqual(length(output{1}.(outputField{iFunc}).(subfield{1}).A_rat), 1, sprintf('Expected 1 %s value for A_rat (%s).', subfield{1}, outputField{iFunc}));
                    end
                end

                % Verify that returned EV and REV values are different
                testCase.verifyNotEqual(output_noRipples.(outputField{iFunc}).EV.A_rat, output_variableRipples.(outputField{iFunc}).EV.A_rat, sprintf('Function calculated the same explained variance with rippleActivityMode = "none" or "variable" (%s).', outputField{iFunc}))
                testCase.verifyNotEqual(output_noRipples.(outputField{iFunc}).REV.A_rat, output_variableRipples.(outputField{iFunc}).REV.A_rat, sprintf('Function calculated the same reverse explained variance with rippleActivityMode = "none" or "variable" (%s).', outputField{iFunc}))
                testCase.verifyNotEqual(output_noRipples.(outputField{iFunc}).EV.A_rat, output_equalisedRipples.(outputField{iFunc}).EV.A_rat, sprintf('Function calculated the same explained variance with rippleActivityMode = "none" or "equalise" (%s).', outputField{iFunc}))
                testCase.verifyNotEqual(output_noRipples.(outputField{iFunc}).REV.A_rat, output_equalisedRipples.(outputField{iFunc}).REV.A_rat, sprintf('Function calculated the same reverse explained variance with rippleActivityMode = "none" or "equalise" (%s).', outputField{iFunc}))
                testCase.verifyNotEqual(output_variableRipples.(outputField{iFunc}).EV.A_rat, output_equalisedRipples.(outputField{iFunc}).EV.A_rat, sprintf('Function calculated the same explained variance with rippleActivityMode = "variable" or "equalise" (%s).', outputField{iFunc}))
                testCase.verifyNotEqual(output_variableRipples.(outputField{iFunc}).REV.A_rat, output_equalisedRipples.(outputField{iFunc}).REV.A_rat, sprintf('Function calculated the same reverse explained variance with rippleActivityMode = "variable" or "equalise" (%s).', outputField{iFunc}))

            end
        end
       
       % Test that intra-area EV-REV calculation returns different values
        % depending on the number of ripples
        function testCalculateIntraAreaEVREV_nRipples(testCase)
            
            functions = {...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateVStrEVREV('dataPath', x),...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateCA1EVREV('dataPath', x)};
            relevantNeurons = {'nNAcUnits', 'nCA1Units'};
            irrelevantNeurons = {'nCA1Units', 'nNAcUnits'};
            outputField = {'vStr', 'CA1'};
            
            for iFunc = 1:numel(functions)
            
                % Mock inclusion criteria and sessions: 1 rat with 1 session
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rewardResponsive = false;
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.minimumFiringRate = 0;
                testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
                testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;

                % Mock analysis configuration
                testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'variable';

                % Create mock behavioural data
                behavDirectory = fullfile(testCase.TestDataDir, 'test_session');
                ratName = 'A_Rat';
                TestUtilities.saveMockBehaviouralData(behavDirectory, ratName);
            
                % Mock the rest of the data: 3 relevant neurons and EV<REV
                ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
                nRelevantNeurons = 3;
                nIrrelevantNeurons = 0;
                seed = 1;
                [spiketimes, binnedfr, startEndTimes, ~, ~, rippleStart, rippleStop] = testCase.createSpiketimes('interarea', false, 'ripples', true, relevantNeurons{iFunc}, nRelevantNeurons, irrelevantNeurons{iFunc}, nIrrelevantNeurons, 'seed', seed);
                nNAcUnits = testCase.getIntraAreaNAc(relevantNeurons, iFunc, nRelevantNeurons, nIrrelevantNeurons);
                ripples3std = struct('rippleStart', rippleStart, 'rippleStop', rippleStop, 'ripplePeak', rippleStart);
                TestUtilities.saveMockEphysData(ephysDirectory,...
                    'spiketimes', spiketimes,...
                    'startEndTimes', startEndTimes,...
                    'nNAcUnits', nNAcUnits,...
                    'binnedfr', binnedfr,...
                    'ripples3std', ripples3std);

                % Call EV and REV with rippleActivityMode "none", 1, 2, and 3
                func = functions{iFunc};
                testCase.ExplainedVarianceReactivationObj.analysisConfig.nRipples = NaN;
                output_noNRipples = func(testCase.TestDataDir);
                testCase.ExplainedVarianceReactivationObj.analysisConfig.nRipples = 1;
                output_1Ripple = func(testCase.TestDataDir);
                testCase.ExplainedVarianceReactivationObj.analysisConfig.nRipples = 2;
                output_2Ripples = func(testCase.TestDataDir);
                testCase.ExplainedVarianceReactivationObj.analysisConfig.nRipples = 3;
                output_3Ripples = func(testCase.TestDataDir);
            
                % Verify output structure
                testCase.verifyNotEmpty(output_noNRipples);
                testCase.verifyNotEmpty(output_1Ripple);
                testCase.verifyNotEmpty(output_2Ripples);
                testCase.verifyNotEmpty(output_3Ripples);
                expectedFields = outputField(iFunc);
                for field = expectedFields
                    testCase.verifyEqual(ismember(field{1}, fieldnames(output_noNRipples)), true, sprintf('Expected field "%s" missing from output with nRipples = NaN (%s).', field{1}, outputField{iFunc}));
                    testCase.verifyEqual(ismember((field{1}), fieldnames(output_noNRipples)), true, sprintf('Field "%s" expected to be a struct with nRipples = NaN (%s).', field{1}, outputField{iFunc}));
                    testCase.verifyEqual(ismember(field{1}, fieldnames(output_1Ripple)), true, sprintf('Expected field "%s" missing from output with nRipples = 1 (%s).', field{1}, outputField{iFunc}));
                    testCase.verifyEqual(ismember((field{1}), fieldnames(output_1Ripple)), true, sprintf('Field "%s" expected to be a struct with nRipples = 1 (%s).', field{1}, outputField{iFunc}));
                    testCase.verifyEqual(ismember(field{1}, fieldnames(output_2Ripples)), true, sprintf('Expected field "%s" missing from output with nRipples = 2 (%s).', field{1}, outputField{iFunc}));
                    testCase.verifyEqual(ismember((field{1}), fieldnames(output_2Ripples)), true, sprintf('Field "%s" expected to be a struct output with nRipples = 2 (%s).', field{1}, outputField{iFunc}));
                    testCase.verifyEqual(ismember(field{1}, fieldnames(output_3Ripples)), true, sprintf('Expected field "%s" missing from output with nRipples = 3 (%s).', field{1}, outputField{iFunc}));
                    testCase.verifyEqual(ismember((field{1}), fieldnames(output_3Ripples)), true, sprintf('Field "%s" expected to be a struct output with nRipples = 3 (%s).', field{1}, outputField{iFunc}));
                end
                expectedSubfields = {'EV', 'REV'};
                for subfield = expectedSubfields
                    for output = {output_noNRipples, output_1Ripple, output_2Ripples, output_3Ripples}
                        testCase.verifyEqual(isfield(output{1}.(outputField{iFunc}), subfield{1}), true, sprintf('Expected subfield "%s" missing (%s).', subfield{1}, outputField{iFunc}));
                        testCase.verifyEqual(isfield(output{1}.(outputField{iFunc}).(subfield{1}), 'A_rat'), true, sprintf('Expected subfield "A_rat" missing from field "%s" (%s).', subfield{1}, outputField{iFunc}));
                        testCase.verifyEqual(length(output{1}.(outputField{iFunc}).(subfield{1}).A_rat), 1, sprintf('Expected 1 %s value for A_rat. (%s)', subfield{1}, outputField{iFunc}));
                    end
                end
            
                % Verify that returned EV and REV values are different for different numbers of ripples
                testCase.verifyNotEqual(output_noNRipples.(outputField{iFunc}).EV.A_rat, output_1Ripple.(outputField{iFunc}).EV.A_rat, sprintf('Function calculated the same explained variance with nRipples = NaN or 1 (%s).', outputField{iFunc}));
                testCase.verifyNotEqual(output_noNRipples.(outputField{iFunc}).REV.A_rat, output_1Ripple.(outputField{iFunc}).REV.A_rat, sprintf('Function calculated the same reverse explained variance with nRipples = NaN or 1 (%s).', outputField{iFunc}));
                testCase.verifyNotEqual(output_noNRipples.(outputField{iFunc}).EV.A_rat, output_2Ripples.(outputField{iFunc}).EV.A_rat, sprintf('Function calculated the same explained variance with nRipples = NaN or 2 (%s).', outputField{iFunc}));
                testCase.verifyNotEqual(output_noNRipples.(outputField{iFunc}).REV.A_rat, output_2Ripples.(outputField{iFunc}).REV.A_rat, sprintf('Function calculated the same reverse explained variance with nRipples = NaN or 2 (%s).', outputField{iFunc}));
                testCase.verifyNotEqual(output_1Ripple.(outputField{iFunc}).EV.A_rat, output_2Ripples.(outputField{iFunc}).EV.A_rat, sprintf('Function calculated the same explained variance with nRipples = 1 or 2 (%s).', outputField{iFunc}));
                testCase.verifyNotEqual(output_1Ripple.(outputField{iFunc}).REV.A_rat, output_2Ripples.(outputField{iFunc}).REV.A_rat, sprintf('Function calculated the same reverse explained variance with nRipples = 1 or 2 (%s).', outputField{iFunc}));
            
                % Verify that returned EV and REV values are the same for NaN and 3 ripples
                testCase.verifyEqual(output_noNRipples.(outputField{iFunc}).EV.A_rat, output_3Ripples.(outputField{iFunc}).EV.A_rat, 'RelTol', 1e-5, sprintf('Function calculated different explained variance with nRipples = NaN or 3, for a total of 3 ripples (%s).', outputField{iFunc}));
                testCase.verifyEqual(output_noNRipples.(outputField{iFunc}).REV.A_rat, output_3Ripples.(outputField{iFunc}).REV.A_rat, 'RelTol', 1e-5, sprintf('Function calculated different reverseexplained variance with nRipples = NaN or 3, for a total of 3 ripples (%s).', outputField{iFunc}));
                
            end
        end
        
        % Test that intra-area EV-REV calculation returns different values
        % depending on the number of ripples (separate number for PRE and POST)
        function testCalculateIntraAreaEVREV_nRipplesPrePost(testCase)
            
            functions = {...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateVStrEVREV('dataPath', x),...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateCA1EVREV('dataPath', x)};
            relevantNeurons = {'nNAcUnits', 'nCA1Units'};
            irrelevantNeurons = {'nCA1Units', 'nNAcUnits'};
            outputField = {'vStr', 'CA1'};
            
            for iFunc = 1:numel(functions)
            
                % Mock inclusion criteria and sessions: 1 rat with 1 session
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rewardResponsive = false;
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.minimumFiringRate = 0;
                testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
                testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;

                % Mock analysis configuration
                testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'variable';

                % Create mock behavioural data
                behavDirectory = fullfile(testCase.TestDataDir, 'test_session');
                ratName = 'A_Rat';
                TestUtilities.saveMockBehaviouralData(behavDirectory, ratName);
            
                % Mock the rest of the data: 3 relevant neurons
                ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
                nRelevantNeurons = 3;
                nIrrelevantNeurons = 0;
                seed = 1;
                [spiketimes, binnedfr, startEndTimes, ~, ~, rippleStart, rippleStop] = testCase.createSpiketimes('interarea', false, 'ripples', true, relevantNeurons{iFunc}, nRelevantNeurons, irrelevantNeurons{iFunc}, nIrrelevantNeurons, 'seed', seed);
                nNAcUnits = testCase.getIntraAreaNAc(relevantNeurons, iFunc, nRelevantNeurons, nIrrelevantNeurons);
                ripples3std = struct('rippleStart', rippleStart, 'rippleStop', rippleStop, 'ripplePeak', rippleStart);
                TestUtilities.saveMockEphysData(ephysDirectory,...
                    'spiketimes', spiketimes,...
                    'startEndTimes', startEndTimes,...
                    'nNAcUnits', nNAcUnits,...
                    'binnedfr', binnedfr,...
                    'ripples3std', ripples3std);
            
                % Call calculateInterAreaEVREV with rippleActivityMode "none", 1, 2, and 3
                func = functions{iFunc};
                testCase.ExplainedVarianceReactivationObj.analysisConfig.nRipples = NaN;
                output_noNRipples = func(testCase.TestDataDir);
                testCase.ExplainedVarianceReactivationObj.analysisConfig.nRipples = [1, 1];
                output_1Ripple = func(testCase.TestDataDir);
                testCase.ExplainedVarianceReactivationObj.analysisConfig.nRipples = [2, 1];
                output_2Ripples = func(testCase.TestDataDir);
                testCase.ExplainedVarianceReactivationObj.analysisConfig.nRipples = [3, 3];
                output_3Ripples = func(testCase.TestDataDir);

                % Verify output structure
                testCase.verifyNotEmpty(output_noNRipples);
                testCase.verifyNotEmpty(output_1Ripple);
                testCase.verifyNotEmpty(output_2Ripples);
                expectedFields = outputField(iFunc);
                for field = expectedFields
                    testCase.verifyEqual(ismember(field{1}, fieldnames(output_noNRipples)), true, sprintf('Expected field "%s" missing from output with nRipples = NaN (%s).', field{1}, outputField{iFunc}));
                    testCase.verifyEqual(ismember((field{1}), fieldnames(output_noNRipples)), true, sprintf('Field "%s" expected to be a struct with nRipples = NaN (%s).', field{1}, outputField{iFunc}));
                    testCase.verifyEqual(ismember(field{1}, fieldnames(output_1Ripple)), true, sprintf('Expected field "%s" missing from output with nRipples = 1 (%s).', field{1}, outputField{iFunc}));
                    testCase.verifyEqual(ismember((field{1}), fieldnames(output_1Ripple)), true, sprintf('Field "%s" expected to be a struct with nRipples = 1 (%s).', field{1}, outputField{iFunc}));
                    testCase.verifyEqual(ismember(field{1}, fieldnames(output_2Ripples)), true, sprintf('Expected field "%s" missing from output with nRipples = 2 (%s).', field{1}, outputField{iFunc}));
                    testCase.verifyEqual(ismember((field{1}), fieldnames(output_2Ripples)), true, sprintf('Field "%s" expected to be a struct output with nRipples = 2 (%s).', field{1}, outputField{iFunc}));
                    testCase.verifyEqual(ismember(field{1}, fieldnames(output_3Ripples)), true, sprintf('Expected field "%s" missing from output with nRipples = 3 (%s).', field{1}, outputField{iFunc}));
                    testCase.verifyEqual(ismember((field{1}), fieldnames(output_3Ripples)), true, sprintf('Field "%s" expected to be a struct output with nRipples = 3 (%s).', field{1}, outputField{iFunc}));
                end
                expectedSubfields = {'EV', 'REV'};
                for subfield = expectedSubfields
                    for output = {output_noNRipples, output_1Ripple, output_2Ripples, output_3Ripples}
                        testCase.verifyEqual(isfield(output{1}.(outputField{iFunc}), subfield{1}), true, sprintf('Expected subfield "%s" missing (%s).', subfield{1}, outputField{iFunc}));
                        testCase.verifyEqual(isfield(output{1}.(outputField{iFunc}).(subfield{1}), 'A_rat'), true, sprintf('Expected subfield "A_rat" missing from field "%s" (%s).', subfield{1}, outputField{iFunc}));
                        testCase.verifyEqual(length(output{1}.(outputField{iFunc}).(subfield{1}).A_rat), 1, sprintf('Expected 1 %s value for A_rat (%s).', subfield{1}, outputField{iFunc}));
                    end
                end
            
                % Verify that returned EV and REV values are different for different numbers of ripples
                testCase.verifyNotEqual(output_noNRipples.(outputField{iFunc}).EV.A_rat, output_1Ripple.(outputField{iFunc}).EV.A_rat, sprintf('Function calculated the same explained variance with nRipples = NaN or 1 (%s).', outputField{iFunc}));
                testCase.verifyNotEqual(output_noNRipples.(outputField{iFunc}).REV.A_rat, output_1Ripple.(outputField{iFunc}).REV.A_rat, sprintf('Function calculated the same reverse explained variance with nRipples = NaN or 1 (%s).', outputField{iFunc}));
                testCase.verifyNotEqual(output_noNRipples.(outputField{iFunc}).EV.A_rat, output_2Ripples.(outputField{iFunc}).EV.A_rat, sprintf('Function calculated the same explained variance with nRipples = NaN or 2 (%s).', outputField{iFunc}));
                testCase.verifyNotEqual(output_noNRipples.(outputField{iFunc}).REV.A_rat, output_2Ripples.(outputField{iFunc}).REV.A_rat, sprintf('Function calculated the same reverse explained variance with nRipples = NaN or 2 (%s).', outputField{iFunc}));
                testCase.verifyNotEqual(output_1Ripple.(outputField{iFunc}).EV.A_rat, output_2Ripples.(outputField{iFunc}).EV.A_rat, sprintf('Function calculated the same explained variance with nRipples = 1 or 2 (%s).', outputField{iFunc}));
                testCase.verifyNotEqual(output_1Ripple.(outputField{iFunc}).REV.A_rat, output_2Ripples.(outputField{iFunc}).REV.A_rat, sprintf('Function calculated the same reverse explained variance with nRipples = 1 or 2 (%s).', outputField{iFunc}));

                % Verify that returned EV and REV values are the same for NaN and 3 ripples
                testCase.verifyEqual(output_noNRipples.(outputField{iFunc}).EV.A_rat, output_3Ripples.(outputField{iFunc}).EV.A_rat, 'RelTol', 1e-5, sprintf('Function calculated different explained variance with nRipples = NaN or 3, for a total of 3 ripples (%s).', outputField{iFunc}));
                testCase.verifyEqual(output_noNRipples.(outputField{iFunc}).REV.A_rat, output_3Ripples.(outputField{iFunc}).REV.A_rat, 'RelTol', 1e-5, sprintf('Function calculated different reverseexplained variance with nRipples = NaN or 3, for a total of 3 ripples (%s).', outputField{iFunc}));

            end
        end
        
        % Test error handling of intra-area EV-REV calculation when there are no ripples
        function testCalculateIntraAreaEVREV_noRipples(testCase)
            
            functions = {...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateVStrEVREV('dataPath', x),...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateCA1EVREV('dataPath', x)};
            relevantNeurons = {'nNAcUnits', 'nCA1Units'};
            irrelevantNeurons = {'nCA1Units', 'nNAcUnits'};
            outputField = {'vStr', 'CA1'};
            
            for iFunc = 1:numel(functions)
            
                % Mock inclusion criteria and sessions: 1 rat with 1 session
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rewardResponsive = false;
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.minimumFiringRate = 0;
                testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
                testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;

                % Mock analysis configuration
                testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'variable';
                testCase.ExplainedVarianceReactivationObj.analysisConfig.nRipples = NaN;

                % Create mock behavioural data
                behavDirectory = fullfile(testCase.TestDataDir, 'test_session');
                ratName = 'A_Rat';
                TestUtilities.saveMockBehaviouralData(behavDirectory, ratName);
            
                % Mock the rest of the data: 3 relevant neurons
                ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
                nRelevantNeurons = 3;
                nIrrelevantNeurons = 0;
                seed = 1;
                [spiketimes, binnedfr, startEndTimes, ~, ~] = testCase.createSpiketimes('interarea', false, relevantNeurons{iFunc}, nRelevantNeurons, irrelevantNeurons{iFunc}, nIrrelevantNeurons, 'seed', seed);
                nNAcUnits = testCase.getIntraAreaNAc(relevantNeurons, iFunc, nRelevantNeurons, nIrrelevantNeurons);
                ripples3std = struct('rippleStart', [], 'rippleStop', [], 'ripplePeak', []);
                TestUtilities.saveMockEphysData(ephysDirectory,...
                    'spiketimes', spiketimes,...
                    'startEndTimes', startEndTimes,...
                    'nNAcUnits', nNAcUnits,...
                    'binnedfr', binnedfr,...
                    'ripples3std', ripples3std);

                % Call method with rippleActivityMode "none", 1, 2, and 3
                func = functions{iFunc};
                expectedErrorId = 'ExplainedVarianceReactivation:invalidData';
                testCase.verifyError(@() func(testCase.TestDataDir), expectedErrorId);
                
            end
        end
        
        % Test error handling of intra-area EV-REV when no rats are specified
        function testCalculateIntraAreaEVREV_noRats(testCase)
            
            functions = {...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateVStrEVREV('dataPath', x),...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateCA1EVREV('dataPath', x)};
            
            for iFunc = 1:numel(functions)
            
                % Mock inclusion criteria and sessions: 1 rat with 1 session
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {};
            
                % Call method
                func = functions{iFunc};
                expectedErrorId = 'ExplainedVarianceReactivation:sessions:noSessions';
                testCase.verifyError(@() func(testCase.TestDataDir), expectedErrorId);
                
            end
        end
        
        % Test error handling of intra-area EV-REV calculation with a missing rat
        function testCalculateIntraAreaEVREV_missingRat(testCase)
            
            functions = {...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateVStrEVREV('dataPath', x),...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateCA1EVREV('dataPath', x)};
            
            for iFunc = 1:numel(functions)
                            
                % Mock inclusion criteria and sessions: 1 rat with 1 session
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'missing_rat'};
                
                % Call method
                func = functions{iFunc};
                expectedErrorId = 'ExplainedVarianceReactivation:sessions:noSessions';
                testCase.verifyError(@() func(testCase.TestDataDir), expectedErrorId);
            
            end
        end
        
        % Test that intra-area EV-REV calculation skips sessions that are not to be included
        function testCalculateIntraAreaEVREV_sessionSelection(testCase)
            
            functions = {...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateVStrEVREV('dataPath', x),...
                @(x) testCase.ExplainedVarianceReactivationObj.calculateCA1EVREV('dataPath', x)};
            relevantNeurons = {'nNAcUnits', 'nCA1Units'};
            irrelevantNeurons = {'nCA1Units', 'nNAcUnits'};
            outputField = {'vStr', 'CA1'};
            
            for iFunc = 1:numel(functions)
                            
                % Mock inclusion criteria and sessions: 1 rat with 1 session
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rewardResponsive = false;
                testCase.ExplainedVarianceReactivationObj.inclusionCriteria.minimumFiringRate = 0;
                testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [2, 3]);
                testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;
            
                % Mock analysis configuration
                testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'none';

                % Create mock behavioural data
                behavDirectory = fullfile(testCase.TestDataDir, 'test_session');
                ratName = 'A_Rat';
                TestUtilities.saveMockBehaviouralData(behavDirectory, ratName);
            
                % Mock ephys data for session 1: 3 relevant neurons
                ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
                nRelevantNeurons = 3;
                nIrrelevantNeurons = 0;
                [spiketimes, binnedfr, startEndTimes, ~, ~] = testCase.createSpiketimes('interarea', false, relevantNeurons{iFunc}, nRelevantNeurons, irrelevantNeurons{iFunc}, nIrrelevantNeurons, 'seed', 2);
                nNAcUnits = testCase.getIntraAreaNAc(relevantNeurons, iFunc, nRelevantNeurons, nIrrelevantNeurons);
                TestUtilities.saveMockEphysData(ephysDirectory,...
                    'spiketimes', spiketimes,...
                    'startEndTimes', startEndTimes,...
                    'nNAcUnits', nNAcUnits,...
                    'binnedfr', binnedfr);
            
                % Mock ephys data for session 2: 3 relevant neurons
                ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session2');
                nRelevantNeurons_session2 = 3;
                nIrrelevantNeurons = 0;
                [spiketimes, binnedfr, startEndTimes, ~, ~] = testCase.createSpiketimes('interarea', false, relevantNeurons{iFunc}, nRelevantNeurons_session2, irrelevantNeurons{iFunc}, nIrrelevantNeurons, 'seed', 3);
                nNAcUnits = testCase.getIntraAreaNAc(relevantNeurons, iFunc, nRelevantNeurons, nIrrelevantNeurons);
                TestUtilities.saveMockEphysData(ephysDirectory,...
                    'spiketimes', spiketimes,...
                    'startEndTimes', startEndTimes,...
                    'nNAcUnits', nNAcUnits,...
                    'binnedfr', binnedfr);
            
                % Mock ephys data for session 3: 4 relevant neurons
                ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session3');
                nRelevantNeurons_session3 = 3;
                nIrrelevantNeurons = 0;
                [spiketimes, binnedfr, startEndTimes, EV, REV] = testCase.createSpiketimes('interarea', false, relevantNeurons{iFunc}, nRelevantNeurons_session3, irrelevantNeurons{iFunc}, nIrrelevantNeurons, 'seed', 5);
                nNAcUnits = testCase.getIntraAreaNAc(relevantNeurons, iFunc, nRelevantNeurons, nIrrelevantNeurons);
                TestUtilities.saveMockEphysData(ephysDirectory,...
                    'spiketimes', spiketimes,...
                    'startEndTimes', startEndTimes,...
                    'nNAcUnits', nNAcUnits,...
                    'binnedfr', binnedfr);
                
                % Call method
                func = functions{iFunc};
                output = func(testCase.TestDataDir);
                        
                % Verify output structure
                testCase.verifyNotEmpty(output);
                expectedFields = outputField(iFunc);
                actualFields = fieldnames(output);
                for field = expectedFields
                    testCase.verifyEqual(ismember(field{1}, actualFields), true, sprintf('Expected field "%s" missing (%s).', field{1}, outputField{iFunc}));
                    testCase.verifyEqual(ismember((field{1}), actualFields), true, sprintf('Field "%s" expected to be a struct (%s).', field{1}, outputField{iFunc}));
                end
                expectedSubfields = {'EV', 'REV'};
                for subfield = expectedSubfields
                    testCase.verifyEqual(isfield(output.(outputField{iFunc}), subfield{1}), true, sprintf('Expected subfield "%s" missing (%s).', subfield{1}, outputField{iFunc}));
                    testCase.verifyEqual(isfield(output.(outputField{iFunc}).(subfield{1}), 'A_rat'), true, sprintf('Expected subfield "A_rat" missing from field "%s" (%s).', subfield{1}, outputField{iFunc}));
                    testCase.verifyEqual(length(output.(outputField{iFunc}).(subfield{1}).A_rat), 2, sprintf('Expected 1 %s value for A_rat (%s).', subfield{1}, outputField{iFunc}));
                end
            
            end
        end
        
        % Test that all EV-REV outputs are retained after running all three functions
        function test_EVREVFunctionIntegration(testCase)
            
            % Mock inclusion criteria and sessions: 1 rat with 1 session
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rewardResponsive = false;
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.minimumFiringRate = 0;
            testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
            testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;
            
            % Mock analysis configuration
            testCase.ExplainedVarianceReactivationObj.analysisConfig.rippleActivityMode = 'none';

            % Create mock behavioural data
            behavDirectory = fullfile(testCase.TestDataDir, 'test_session');
            ratName = 'A_Rat';
            TestUtilities.saveMockBehaviouralData(behavDirectory, ratName);
            
            % Mock ephys data: 5 vStr and 5 CA1 neurons and EV>REV between brain areas
            ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
            nNAcUnits = 5;
            nCA1Units = 5;
            [spiketimes, binnedfr, startEndTimes, EV, REV] = testCase.createSpiketimes('interarea', true, 'nNAcUnits', nNAcUnits, 'nCA1Units', nCA1Units, 'seed', 5);
            assert(EV>REV, 'Mock data EV-REV should be positive for testing in this function')
            TestUtilities.saveMockEphysData(ephysDirectory,...
                'spiketimes', spiketimes,...
                'startEndTimes', startEndTimes,...
                'nNAcUnits', nNAcUnits,...
                'binnedfr', binnedfr);
            
            % Call all three methods
            testCase.ExplainedVarianceReactivationObj = testCase.ExplainedVarianceReactivationObj.calculateInterAreaEVREV('dataPath', testCase.TestDataDir);
            testCase.ExplainedVarianceReactivationObj = testCase.ExplainedVarianceReactivationObj.calculateVStrEVREV('dataPath', testCase.TestDataDir);
            output = testCase.ExplainedVarianceReactivationObj.calculateCA1EVREV('dataPath', testCase.TestDataDir);
            
            % Verify output structure for overall EV-REV
            testCase.verifyNotEmpty(output);
            expectedStructFields = {'CA1_vStr', 'vStr', 'CA1'};
            expectedSubfields = {'EV', 'REV'};
            actualFields = fieldnames(output);
            for field = expectedStructFields
                testCase.verifyEqual(ismember(field{1}, actualFields), true, sprintf('Expected field "%s" missing.', field{1}));
                testCase.verifyEqual(ismember((field{1}), actualFields), true, sprintf('Field "%s" expected to be a struct.', field{1}));
                for subfield = expectedSubfields
                    testCase.verifyEqual(isfield(output.(field{1}), subfield{1}), true, sprintf('Expected subfield "%s" missing (%s).', subfield{1}, field{1}));
                    testCase.verifyEqual(isfield(output.(field{1}).(subfield{1}), 'A_rat'), true, sprintf('Expected subfield "A_rat" missing from field "%s" (%s).', subfield{1}, field{1}));
                    testCase.verifyEqual(length(output.(field{1}).(subfield{1}).A_rat), 1, sprintf('Expected 1 %s value for A_rat (%s).', subfield{1}, field{1}));
                    testCase.verifyNotEqual(output.(field{1}).(subfield{1}).A_rat, NaN, sprintf('Expected a non-NaN %s value for A_rat (%s).', subfield{1}, field{1}));
                end

            end
            
            % Verify output structure for per-cell contributions
            testCase.verifyEqual(ismember('perCellContributions', actualFields), true, sprintf('Expected field "%s" missing.', 'perCellContributions'));
            testCase.verifyEqual(isstruct(output.perCellContributions.A_rat{1}), true, 'Expected the per-cell contributions to be a struct.')
            testCase.verifyEqual(length(output.perCellContributions.A_rat{1}), nNAcUnits*nCA1Units, 'Expected the per-cell contributions to match the number of cell pairs.')
            expectedSubfields = {'cells', 'contribution', 'correlations'};
            for subfield = expectedSubfields
                testCase.verifyEqual(isfield(output.perCellContributions.A_rat{1}, subfield{1}), true, sprintf('Expected subfield "%s" missing from per-cell contributions.', subfield{1}));
            end

        end
        
        
        % =============== cell-pair reactivation significance ================ %
        
        % Test basic functionality of getSignificantlyReactivatedCellPairs
        function testGetSignificantlyReactivatedCellPairs_basicFunctionality(testCase)
            
            % Mock inclusion criteria and sessions: 1 rat with 1 session
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
            testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
            testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;
            
            % Mock data on per-cell contributions to EV-REV
            testCase.ExplainedVarianceReactivationObj.perCellContributions = struct('A_rat', {{struct('cells', {[1, 1]}, 'contribution', {[0.1]}, 'correlations', {[-0.1; 0; 0.1]})}});
            
            % Call getSignificantlyReactivatedCellPairs
            output = testCase.ExplainedVarianceReactivationObj.getSignificantlyReactivatedCellPairs();
            
            % Verify that the required outputs are produced
            expectedFields = {'controlCellPairs', 'reactivatedCellPairs'};
            for field = expectedFields
                fieldName = field{1};
                testCase.verifyTrue(isfield(output.perCellContributions, fieldName), sprintf('Expected field "%s" is missing.', fieldName));
                testCase.verifyTrue(isstruct(output.perCellContributions.(fieldName)), sprintf('Output field "%s" should be a struct.', fieldName));
                testCase.verifyTrue(isfield(output.perCellContributions.(fieldName), 'A_rat'), sprintf('Field "%s" is missing a subfield for each rat.', fieldName));
            end
        
        end
        
        % Test implementation of negativeContributions with one session
        function testGetSignificantlyReactivatedCellPairs_negativeContributions(testCase)
            
            % Mock inclusion criteria and sessions: 1 rat with 1 session
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
            testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
            testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;
            
            % Mock analysis config
            testCase.ExplainedVarianceReactivationObj.analysisConfig.controlCriterion = 'negativeContributions';
            testCase.ExplainedVarianceReactivationObj.analysisConfig.reactivationQuantileThreshold = 0.1;
            
            % Mock data on per-cell contributions to EV-REV: 10 cell pairs, 
            % with one cell pair above and one cell pair below the quantile threshold
            contributions = mat2cell([-0.5 linspace(-0.1, 0.1, 8) 0.5], 1, ones(1, 10));
            cells = {{1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {1, 6}, {1, 7}, {1, 8}, {1, 9}, {1, 10}};
            correlations = repmat({[-0.5; 0.0; 0.5]}, 1, 10);   % All cell pairs show increasing correlation from PRE to POST
            testCase.ExplainedVarianceReactivationObj.perCellContributions = struct('A_rat', {{struct('cells', cells, 'contribution', contributions, 'correlations', correlations)}});
            
            % Call getSignificantlyReactivatedCellPairs
            output = testCase.ExplainedVarianceReactivationObj.getSignificantlyReactivatedCellPairs();
            
            % Verify that the outputs are correct
            expectedControls = [1, 1];
            expectedReactivated = [1, 10];
            actualControls = cell2mat(output.perCellContributions.controlCellPairs.A_rat{1});
            actualReactivated = cell2mat(output.perCellContributions.reactivatedCellPairs.A_rat{1});
            testCase.verifyEqual(expectedControls, actualControls, 'Control cell pairs are not as expected.')
            testCase.verifyEqual(expectedReactivated, actualReactivated, 'Control cell pairs are not as expected.')
            
        end
        
         % Test implementation of minimalContributions with one session
        function testGetSignificantlyReactivatedCellPairs_minimalContributions(testCase)
            
            % Mock inclusion criteria and sessions: 1 rat with 1 session
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
            testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
            testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;
            
            % Mock analysis config
            testCase.ExplainedVarianceReactivationObj.analysisConfig.controlCriterion = 'minimalContributions';
            testCase.ExplainedVarianceReactivationObj.analysisConfig.reactivationQuantileThreshold = 0.1;
            
            % Mock data on per-cell contributions to EV-REV: 10 cell pairs, 
            % with one cell pair above and one cell pair below the quantile threshold
            contributions = mat2cell([0 linspace(-0.1, 0.1, 8) 0.5], 1, ones(1, 10));
            cells = {{1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {1, 6}, {1, 7}, {1, 8}, {1, 9}, {1, 10}};
            correlations = repmat({[-0.5; 0.0; 0.5]}, 1, 10);   % All cell pairs show increasing correlation from PRE to POST
            testCase.ExplainedVarianceReactivationObj.perCellContributions = struct('A_rat', {{struct('cells', cells, 'contribution', contributions, 'correlations', correlations)}});
            
            % Call getSignificantlyReactivatedCellPairs
            output = testCase.ExplainedVarianceReactivationObj.getSignificantlyReactivatedCellPairs();
            
            % Verify that the outputs are correct
            expectedControls = [1, 1];
            expectedReactivated = [1, 10];
            actualControls = cell2mat(output.perCellContributions.controlCellPairs.A_rat{1});
            actualReactivated = cell2mat(output.perCellContributions.reactivatedCellPairs.A_rat{1});
            testCase.verifyEqual(expectedControls, actualControls, 'Control cell pairs are not as expected.')
            testCase.verifyEqual(expectedReactivated, actualReactivated, 'Control cell pairs are not as expected.')
            
        end
        
        % Test implementation of negativeContributions with two rats and three sessions
        function testGetSignificantlyReactivatedCellPairs_mutipleRats(testCase)
            
            % Mock inclusion criteria and sessions: 2 rats with 2 and 1 sessions respectively
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat', 'B_rat'};
            testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis = struct('A_rat', [1, 2], 'B_rat', [2]);
            testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = testCase.ExplainedVarianceReactivationObj.sessionsForSpecificAnalysis;
            
            % Mock analysis config
            testCase.ExplainedVarianceReactivationObj.analysisConfig.controlCriterion = 'negativeContributions';
            testCase.ExplainedVarianceReactivationObj.analysisConfig.reactivationQuantileThreshold = 0.1;
            
            % Mock data on per-cell contributions to EV-REV: 10 cell pairs, 
            % with one cell pair above and one cell pair below the quantile threshold
            contributions = mat2cell([-0.5 linspace(-0.1, 0.1, 8) 0.5], 1, ones(1, 10));
            cells = {{1, 1}, {1, 2}, {1, 3}, {1, 4}, {1, 5}, {1, 6}, {1, 7}, {1, 8}, {1, 9}, {1, 10}};
            correlations = repmat({[-0.5; 0.0; 0.5]}, 1, 10);   % All cell pairs show increasing correlation from PRE to POST
            testCase.ExplainedVarianceReactivationObj.perCellContributions = struct(...
                'A_rat', {{...
                    struct('cells', cells, 'contribution', contributions, 'correlations', correlations),...
                    struct('cells', cells, 'contribution', contributions, 'correlations', correlations)}},...
                'B_rat', {{...
                    [],...
                    struct('cells', cells, 'contribution', contributions, 'correlations', correlations)}});
            
            % Call getSignificantlyReactivatedCellPairs
            output = testCase.ExplainedVarianceReactivationObj.getSignificantlyReactivatedCellPairs();
            
            % Verify that the outputs are correct
            expectedControls = [1, 1];
            expectedReactivated = [1, 10];
            testCase.verifyTrue(isfield(output.perCellContributions.controlCellPairs, 'A_rat'), 'Expected field "A_rat" is missing from output.')
            testCase.verifyTrue(isfield(output.perCellContributions.controlCellPairs, 'B_rat'), 'Expected field "B_rat" is missing from output.')
            testCase.verifyEqual(expectedControls, cell2mat(output.perCellContributions.controlCellPairs.A_rat{1}), 'Control cell pairs are not as expected for A_rat session 1.');
            testCase.verifyEqual(expectedControls, cell2mat(output.perCellContributions.controlCellPairs.A_rat{2}), 'Control cell pairs are not as expected for A_rat session 2.');
            testCase.verifyEmpty(cell2mat(output.perCellContributions.controlCellPairs.B_rat{1}), 'Expected output for B_rat session 1 to be empty, because it is not listed in the inclusion criteria.');
            testCase.verifyEqual(expectedControls, cell2mat(output.perCellContributions.controlCellPairs.B_rat{2}), 'Control cell pairs are not as expected for B_rat session 2.');
            testCase.verifyEqual(expectedReactivated, cell2mat(output.perCellContributions.reactivatedCellPairs.A_rat{1}), 'Reactivated cell pairs are not as expected for A_rat session 1.');
            testCase.verifyEqual(expectedReactivated, cell2mat(output.perCellContributions.reactivatedCellPairs.A_rat{2}), 'Reactivated cell pairs are not as expected for A_rat session 2.');
            testCase.verifyEmpty(cell2mat(output.perCellContributions.reactivatedCellPairs.B_rat{1}), 'Expected output for B_rat session 1 to be empty, because it is not listed in the inclusion criteria.');
            testCase.verifyEqual(expectedReactivated, cell2mat(output.perCellContributions.reactivatedCellPairs.B_rat{2}), 'Reactivated cell pairs are not as expected for B_rat session 2.');
            
        end
        
        
        % ================== extract and transform data ================== %
        
        % Test basic functionality of getReactivatedCellPairActivity
        function testGetReactivatedCellPairActivity_basicFunctionality(testCase)
            
            % Mock inclusion criteria and sessions: 1 rat with 1 session
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
            
            % Mock behavioural data
            behaviouralDataPath = fullfile(testCase.TestDataDir, 'test_session'); 
            n_trials = [2];
            actions = {[1, 2]};
            arm_values = {{...
                {'Optimal', 'Legitimate', 'Legitimate'},...
                {'Legitimate', 'Optimal', 'Legitimate'}}};
            reward_probs = {repmat({'high', 'medium', 'low'}, 1, 2)};
            TestUtilities.saveMockBehaviouralData(behaviouralDataPath, 'A_rat', 'n_trials', n_trials, 'actions', actions, 'reward_probs', reward_probs, 'arm_values', arm_values);
            
            % Mock ephys data
            dataPath = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
            spiketimes = {[0.5, 1.5, 2.5], [1.0, 2.5, 2.0]};
            binnedfr = random('Normal', 1, 1, 2, 50);
            CPEntryTimes = [0.5, 2];
            RarrivalTimes = [1.5];
            CPExitTimes = [1, 2.5];
            TestUtilities.saveMockEphysData(dataPath, 'spiketimes', spiketimes, 'binnedfr', binnedfr, 'CPEntryTimes', CPEntryTimes, 'RarrivalTimes', RarrivalTimes, 'CPExitTimes', CPExitTimes);
            
            % Mock per-cell contributions of reactivated and control cell pairs
            testCase.ExplainedVarianceReactivationObj.perCellContributions.reactivatedCellPairs = struct('A_rat', {{[1, 1]}});
            testCase.ExplainedVarianceReactivationObj.perCellContributions.controlCellPairs = struct('A_rat', {{[1, 2]}});
            
            % Call getReactivatedCellPairActivity
            output = testCase.ExplainedVarianceReactivationObj.getReactivatedCellPairActivity('centralPlatform', 'ephysDataPath', testCase.TestDataDir, 'behaviouralDataPath', behaviouralDataPath);
            testCase.verifyTrue(ismember('taskCoactivity', fieldnames(output)), 'Output of getReactivatedCellPairActivity() is missing the field "taskCoactivity".'); 
            testCase.verifyNotEmpty(output.taskCoactivity, 'Output of getReactivatedCellPairActivity() is missing the field "taskCoactivity".');
            
        end
        
        % Test basic functionality of getEventTriggeredActivity
        function testGetEventTriggeredActivity_basicFunctionality(testCase)
            
            % Mock inclusion criteria and sessions: 1 rat with 1 session
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
            
            % Mock event-triggered data
            binSize = testCase.ExplainedVarianceReactivationObj.BinSize;
            rng(0);
            taskCoactivity = struct(...
                'perCellHighRewardExpectation', struct(...
                    'reactivated', struct('A_rat', random('Normal', 1, 1, 10, 200)),...
                    'control', struct('A_rat', random('Normal', 1, 1, 10, 200))),...
                'perCellMediumRewardExpectation', struct(...
                    'reactivated', struct('A_rat', random('Normal', 1, 1, 10, 200)),...
                    'control', struct('A_rat', random('Normal', 1, 1, 10, 200))),...
                'ts', (-5 + binSize/2) : binSize : (5 - binSize/2));
            testCase.ExplainedVarianceReactivationObj.taskCoactivity = taskCoactivity;
            
            % Call getEventTriggeredActivity
            output = testCase.ExplainedVarianceReactivationObj.getEventTriggeredActivity();
            testCase.verifyTrue(ismember('meanEventTriggeredActivity', fieldnames(output)), 'Output of getReactivatedCellPairActivity() is missing the field "meanEventTriggeredActivity".'); 
            testCase.verifyNotEmpty(output.meanEventTriggeredActivity, 'Output of getReactivatedCellPairActivity() is missing the field "meanEventTriggeredActivity".');
            
        end
            
        % Test mixedEffectsNestedAnova with no significant effects
        function testMixedEffectsNestedAnova_noEffects(testCase)
            
            % Mock inclusion criteria and sessions: 1 rat with 1 session
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
            
            % Mock event-triggered data: identical for all conditions, so
            % no significant differences
            rng(0);
            data = {random('Normal', 1, 1, 200, 1)};
            meanEventTriggeredActivity = struct(...
                'high', struct('reactivated', struct('A_rat', data), 'control', struct('A_rat', data)),...
                'medium', struct('reactivated', struct('A_rat', data), 'control', struct('A_rat', data)));
            testCase.ExplainedVarianceReactivationObj.meanEventTriggeredActivity = meanEventTriggeredActivity;
            
            % Call mixedEffectsNestedAnova
            [table, p_interaction, p_medium, p_high] = testCase.ExplainedVarianceReactivationObj.mixedEffectsNestedAnova();
            
            % Verify that all effects are insignificant
            testCase.verifyTrue(isnan(p_interaction));
            testCase.verifyTrue(p_medium > 0.05);
            testCase.verifyTrue(p_high > 0.05);
            
        end
            
        % Test mixedEffectsNestedAnova with an interaction between
        % cell-pair type and arm type
        function testMixedEffectsNestedAnova_interactionEffect(testCase)
            
            % Mock inclusion criteria and sessions: 1 rat with 1 session
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
            
            % Mock event-triggered data: identical for all conditions except reactivted-high,
            % so interaction effect is significant
            rng(0);
            data = {random('Normal', 1, 1, 200, 1)};
            differentData = {data{1}*2};
            meanEventTriggeredActivity = struct(...
                'high', struct('reactivated', struct('A_rat', differentData), 'control', struct('A_rat', data)),...
                'medium', struct('reactivated', struct('A_rat', data), 'control', struct('A_rat', data)));
            testCase.ExplainedVarianceReactivationObj.meanEventTriggeredActivity = meanEventTriggeredActivity;
            
            % Call mixedEffectsNestedAnova
            [table, p_interaction, p_medium, p_high] = testCase.ExplainedVarianceReactivationObj.mixedEffectsNestedAnova();
            
            % Verify a main effect and an interaction effect
            testCase.verifyTrue(p_interaction < 0.05);
            testCase.verifyTrue(p_medium > 0.05);
            testCase.verifyTrue(p_high < 0.05);
            
        end
        
        % Test mixedEffectsNestedAnova with main effect of cell-pair type
        function testMixedEffectsNestedAnova_mainEffectOnly(testCase)
            
            % Mock inclusion criteria and sessions: 1 rat with 1 session
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
            
            % Mock event-triggered data: systematically different between
            % reactivated and control cell pairs, so main effect with no
            % interaction effect
            rng(0);
            data1 = {random('Normal', 1, 1, 200, 1)};
            data2 = {random('Normal', 1, 1, 200, 1)};
            differentData = {data1{1}*2};
            meanEventTriggeredActivity = struct(...
                'high', struct('reactivated', struct('A_rat', differentData), 'control', struct('A_rat', data1)),...
                'medium', struct('reactivated', struct('A_rat', differentData), 'control', struct('A_rat', data2)));
            testCase.ExplainedVarianceReactivationObj.meanEventTriggeredActivity = meanEventTriggeredActivity;
            
            % Call mixedEffectsNestedAnova
            [table, p_interaction, p_medium, p_high] = testCase.ExplainedVarianceReactivationObj.mixedEffectsNestedAnova();
            
            % Verify post-hoc significance and no interaction effect
            testCase.verifyTrue(p_interaction > 0.05);
            testCase.verifyTrue(p_medium < 0.05);
            testCase.verifyTrue(p_high < 0.05);
            
        end

        
        % ========================== plotting ========================== %
        
        % Test basic functionality of plotEVREV
        function testPlotEVREV_basicFunctionality(testCase)
            
            % Mock EV and REV values for within and between brain areas
            % (single rat, single session)
            testCase.ExplainedVarianceReactivationObj.CA1 = struct('EV', struct('A_rat', [0.1]), 'REV', struct('A_rat', [0.1]));
            testCase.ExplainedVarianceReactivationObj.vStr = struct('EV', struct('A_rat', [0.1]), 'REV', struct('A_rat', [0.1]));
            testCase.ExplainedVarianceReactivationObj.CA1_vStr = struct('EV', struct('A_rat', [0.1]), 'REV', struct('A_rat', [0.1]));
            
            % Call plotEVREV
            testCase.ExplainedVarianceReactivationObj.plotEVREV();
            
            % Verify that all objects apear on the figure
            f = gcf;
            patches = findobj(f, 'Type', 'patch');
            lines = findobj(f, 'Type', 'line');
            text = findobj(f, 'Type', 'text');
            expectedNBars = 6;
            expectedNTTests = 3;
            expectedNPatches = expectedNBars;
            expectedNLines = (expectedNBars*2) + (expectedNTTests*3);
            expectedNTexts = expectedNTTests;
            testCase.verifyEqual(expectedNPatches, numel(patches), sprintf('Expected %i patch objects on the figure, but there are %i.', expectedNPatches, numel(patches)));
            testCase.verifyEqual(expectedNLines, numel(lines), sprintf('Expected %i line objects on the figure, but there are %i.', expectedNLines, numel(lines)));
            testCase.verifyEqual(expectedNTexts, numel(text), sprintf('Expected %i text objects on the figure, but there are %i.', expectedNTexts, numel(text)));
            
            % Verify that axes are labelled
            ax = gca;
            testCase.verifyNotEmpty(ax.YLabel.String, 'Figure does not have a y-axis label.');
            testCase.verifyNotEmpty(ax.XTickLabel, 'Figure does not have x-axis labels.');
            
        end
        
        % Test basic functionality of plotEVREV
        function testPlotEVREV_multipleRats(testCase)
            
            colour_scheme;
            
            % Mock EV and REV values for within and between brain areas
            % (two rats, eight sessions)
            EV_values = [0.1, 0.2, 0.2, 0.3, 0.3, 0.4, 0.4, 0.5];
            REV_values = EV_values / 2;
            EV_expectedQuartiles = [0.2, 0.3, 0.4];
            REV_expectedQuartiles = [0.1, 0.15, 0.2];
            testCase.ExplainedVarianceReactivationObj.CA1 = struct(...
                'EV', struct('A_rat', EV_values(1:4), 'B_rat', EV_values(5:end)),...
                'REV', struct('A_rat', REV_values(1:4), 'B_rat', REV_values(5:end)));
            testCase.ExplainedVarianceReactivationObj.vStr = struct(...
                'EV', struct('A_rat', EV_values(1:4) + 0.1, 'B_rat', EV_values(5:end) + 0.1),...
                'REV', struct('A_rat', REV_values(1:4) + 0.1, 'B_rat', REV_values(5:end) + 0.1));
            testCase.ExplainedVarianceReactivationObj.CA1_vStr = struct(...
                'EV', struct('A_rat', EV_values(1:4) + 0.2, 'B_rat', EV_values(5:end) + 0.2),...
                'REV', struct('A_rat', REV_values(1:4) + 0.2, 'B_rat', REV_values(5:end) + 0.2));
            
            % Call plotEVREV
            testCase.ExplainedVarianceReactivationObj.plotEVREV();
            
            % Verify that all objects apear on the figure
            f = gcf;
            patches = findobj(f, 'Type', 'patch');
            lines = findobj(f, 'Type', 'line');
            texts = findobj(f, 'Type', 'text');
            expectedNBars = 6;
            expectedNTTests = 3;
            expectedNPatches = expectedNBars;
            expectedNLines = (expectedNBars*2) + (expectedNTTests*3);
            expectedNTexts = expectedNTTests;
            testCase.verifyEqual(expectedNPatches, numel(patches), sprintf('Expected %i patch objects on the figure, but there are %i.', expectedNPatches, numel(patches)));
            testCase.verifyEqual(expectedNLines, numel(lines), sprintf('Expected %i line objects on the figure, but there are %i.', expectedNLines, numel(lines)));
            testCase.verifyEqual(expectedNTexts, numel(texts), sprintf('Expected %i text objects on the figure, but there are %i.', expectedNTexts, numel(texts)));
            
            % Verify that all medians appear represented by a line with the right colour
            EV_expectedMedians = EV_expectedQuartiles(2) + [0, 0.1, 0.2];
            REV_expectedMedians = REV_expectedQuartiles(2) + [0, 0.1, 0.2];
            REV_expectedColours = {colourscheme.CA1, colourscheme.vStr, [0, 0, 0]};
            for i = 1:3
                testCase.findPlottedMedian(lines, EV_expectedMedians(i), NaN);
                testCase.findPlottedMedian(lines, REV_expectedMedians(i), REV_expectedColours{i});
            end
            
            % Verify that all interquartile range patches for EV appear with the right colour
            tolerance = 1e-5;
            foundCA1Patch = false; foundVStrPatch = false; foundVStr_CA1Patch = false;
            for iPatch = 1:numel(patches)
                patchYData = patches(iPatch).Vertices(:, 2);
                patchColour = patches(iPatch).FaceColor;
                if sum(abs(patchYData - [EV_expectedQuartiles(1); EV_expectedQuartiles(1); EV_expectedQuartiles(3); EV_expectedQuartiles(3)])) < tolerance
                    if patchColour == colourscheme.CA1
                        foundCA1Patch = true;
                    end
                elseif sum(abs(patchYData - ([EV_expectedQuartiles(1); EV_expectedQuartiles(1); EV_expectedQuartiles(3); EV_expectedQuartiles(3)] + 0.1))) < tolerance
                    if patchColour == colourscheme.vStr
                        foundVStrPatch = true;
                    end
                elseif sum(abs(patchYData - ([EV_expectedQuartiles(1); EV_expectedQuartiles(1); EV_expectedQuartiles(3); EV_expectedQuartiles(3)] + 0.2))) < tolerance
                    if patchColour == [0, 0, 0,];
                        foundVStr_CA1Patch = true;
                    end
                end
            end
            testCase.verifyTrue(foundCA1Patch, 'Figure does not contain a patch showing the interquartile range of CA1 EV in the right colour.')
            testCase.verifyTrue(foundVStrPatch, 'Figure does not contain a patch showing the interquartile range of vStr EV in the right colour.')
            testCase.verifyTrue(foundVStr_CA1Patch, 'Figure does not contain a patch showing the interquartile range of CA1-vStr EV in the right colour.')
            
            % Verify that all interquartile range patches for REV appear with the right colour
            tolerance = 1e-5;
            foundCA1Patch = false; foundVStrPatch = false; foundVStr_CA1Patch = false;
            for iPatch = 1:numel(patches)
                patchYData = patches(iPatch).Vertices(:, 2);
                patchColour = patches(iPatch).EdgeColor;
                if sum(abs(patchYData - [REV_expectedQuartiles(1); REV_expectedQuartiles(1); REV_expectedQuartiles(3); REV_expectedQuartiles(3)])) < tolerance
                    if patchColour == colourscheme.CA1
                        foundCA1Patch = true;
                    end
                elseif sum(abs(patchYData - ([REV_expectedQuartiles(1); REV_expectedQuartiles(1); REV_expectedQuartiles(3); REV_expectedQuartiles(3)] + 0.1))) < tolerance
                    if patchColour == colourscheme.vStr
                        foundVStrPatch = true;
                    end
                elseif sum(abs(patchYData - ([REV_expectedQuartiles(1); REV_expectedQuartiles(1); REV_expectedQuartiles(3); REV_expectedQuartiles(3)] + 0.2))) < tolerance
                    if patchColour == [0, 0, 0,];
                        foundVStr_CA1Patch = true;
                    end
                end
            end
            testCase.verifyTrue(foundCA1Patch, 'Figure does not contain a patch showing the interquartile range of CA1 REV in the right colour.')
            testCase.verifyTrue(foundVStrPatch, 'Figure does not contain a patch showing the interquartile range of vStr REV in the right colour.')
            testCase.verifyTrue(foundVStr_CA1Patch, 'Figure does not contain a patch showing the interquartile range of CA1-vStr REV in the right colour.')
            
            % Verify that all t-tests are indicated significant
            for iText = 1:3
                testCase.verifyTrue(strcmp(texts(iText).String, '*'), sprintf('Figure text is expected to be "*" but it is %s.', texts(iText).String));
            end
                        
        end
        
        % Test plotEVREV error handling when there is no data to plot
        function testPlotEVREV_noData(testCase)
            
            % Call plotEVREV without mocking data
            expectedErrorId = 'ExplainedVarianceReactivation:dataNotAssigned';
            testCase.verifyError(@() testCase.ExplainedVarianceReactivationObj.plotEVREV(), expectedErrorId);
            
        end
        
        % Test plotRunningSpeedAndPopulationActivity basic functionality
        function testPlotRunningSpeedAndPopulationActivity_basicFunctionality(testCase)
            
            % Mock inclusion criteria: 1 rat
            testCase.ExplainedVarianceReactivationObj.inclusionCriteria.rats = {'A_rat'};
            testCase.ExplainedVarianceReactivationObj.sessionsForBroadAnalysis = struct('A_rat', [1]);
            
            % Mock running speed data
            dataPath = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
            RarrivalTimes = [5];
            [spiketimes, binnedfr, ~, ~, ~, ~, ~] = testCase.createSpiketimes('nNAcUnits', 3, 'nCA1Units', 3, 'interarea', true);
            runningSpeed = struct('speed', random('Normal', 1, 1, 1, 1000), 'timestamps', linspace(0, 20, 1000));
            TestUtilities.saveMockEphysData(dataPath, 'RarrivalTimes', RarrivalTimes, 'CPEntryTimes', RarrivalTimes, 'CPExitTimes', RarrivalTimes, 'runningSpeed', runningSpeed, 'spiketimes', spiketimes, 'binnedfr', binnedfr, 'nNAcUnits', 3);
            
            % Mock behavioural data
            behaviouralDataPath = fullfile(testCase.TestDataDir, 'test_session'); 
            TestUtilities.saveMockBehaviouralData(behaviouralDataPath, 'A_rat');

            
            % Call plotRunningSpeedAndPopulationActivity
            event = 'rewardArrival';
            scope = 'all';
            testCase.ExplainedVarianceReactivationObj.plotRunningSpeedAndPopulationActivity(event, scope, 'behaviouralDataPath', behaviouralDataPath, 'ephysDataPath', testCase.TestDataDir);
            
            % Verify that a figure is created
            f = gcf;
            axes = findobj(f, 'Type', 'Axes');
            testCase.verifyEqual(numel(axes), 3, 'Expected three axes.');
            
            % Verify that each axis plots a line with shaded region and has
            % a y-label, and the bottom axis has an x-label
            for iAx = 1:3
                lines = findobj(axes(iAx).Children, 'Type', 'line');
                patches = findobj(axes(iAx).Children, 'Type', 'patch');
                testCase.verifyTrue(numel(lines) >= 1, sprintf('No line plotted on axis %i.', iAx));
                testCase.verifyTrue(numel(patches) >= 1, sprintf('No line plotted on axis %i.', iAx));
                testCase.verifyNotEmpty(axes(iAx).YLabel.String, sprintf('Subplot %i has no y-axis label.', iAx));
            end
            testCase.verifyNotEmpty(axes(1).XLabel.String, 'Bottom subplot has no x-axis label.');
            
        end
        
    end
    
    end    
    