
classdef RewardRelatedFiringTest < matlab.unittest.TestCase
    
    properties
        RewardRelatedFiringObj
        TestDataDir
    end
    
    
    methods(Static)
        
        % All expected fields of rewardFiring structure
        function [perCellFields, stackedFields, otherFields] = RewardFiringExpectedFields()
            perCellFields = {'perCellHighRewardExpectation', 'perCellMediumRewardExpectation', 'perCellRewarded', 'perCellUnrewarded'};
            stackedFields = {'highRewardExpectation', 'mediumRewardExpectation', 'rewarded', 'unrewarded'};
            otherFields = {'ts'};
        end
                
        % Create mock data for a session of 4 trials
        function [behavDirectory, n_trials, actions, rewarded] = createMockBehaviouralDataSingleSession(obj, ratName)
    
            behavDirectory = fullfile(obj.TestDataDir, 'test_session');
            n_trials = [4];
            actions = {[1, 2, 1, 2]};
            rewarded = {[1, 1, 0, 0]};
            arm_values = {{...
                {'Optimal', 'Legitimate', 'Legitimate'},...
                {'Legitimate', 'Optimal', 'Legitimate'},...
                {'Optimal', 'Legitimate', 'Legitimate'},...
                {'Legitimate', 'Optimal', 'Legitimate'}}};
            reward_probs = {repmat({'high', 'medium', 'low'}, 1, 4)};
            TestUtilities.saveMockBehaviouralData(...
                behavDirectory,...
                ratName,...
                'n_trials', n_trials,...
                'actions', actions,...
                'rewarded', rewarded,...
                'arm_values', arm_values,...
                'reward_probs', reward_probs);
        end
        
        % Create mock data for 3 sessions of 4 trials
        function [behavDirectory, n_trials, actions, rewarded] = createMockBehaviouralDataThreeSessions(obj, ratName)
    
            behavDirectory = fullfile(obj.TestDataDir, 'test_session');
            n_trials = [2, 4, 6];
            actions = {[1, 3], [1, 2, 1, 2], [1, 3, 1, 2, 3, 2]};
            rewarded = {[1, 0], [1, 1, 0, 0], [1, 0, 0, 1, 0, 0]};
            arm_values = {...
                {{'Optimal', 'Legitimate', 'Legitimate'},...
                {'Legitimate', 'Optimal', 'Legitimate'}},...
                {{'Optimal', 'Legitimate', 'Legitimate'},...
                {'Legitimate', 'Optimal', 'Legitimate'},...
                {'Optimal', 'Legitimate', 'Legitimate'},...
                {'Legitimate', 'Optimal', 'Legitimate'}},...
                {{'Optimal', 'Legitimate', 'Legitimate'},...
                {'Legitimate', 'Optimal', 'Legitimate'},...
                {'Optimal', 'Legitimate', 'Legitimate'},...
                {'Legitimate', 'Optimal', 'Legitimate'},...
                {'Optimal', 'Legitimate', 'Legitimate'},...
                {'Optimal', 'Legitimate', 'Legitimate'}}};
            reward_probs = {...
                repmat({'high', 'medium', 'low'}, 1, 2)...
                repmat({'high', 'medium', 'low'}, 1, 4)...
                repmat({'high', 'medium', 'low'}, 1, 6)};
            TestUtilities.saveMockBehaviouralData(...
                behavDirectory,...
                ratName,...
                'n_trials', n_trials,...
                'actions', actions,...
                'rewarded', rewarded,...
                'arm_values', arm_values,...
                'reward_probs', reward_probs);
        end
        
        % Test the structure of rewardFiring
        function verifyRewardFiringStructure(testCase, output, expectedStructFields, expectedNonStructFields, expectedRatNames)
            testCase.verifyNotEmpty(output);
            for field = expectedStructFields
                % Test that the field exists
                testCase.verifyEqual(isfield(output.rewardFiring, field{1}), true, sprintf('Expected field "%s" missing.', field{1}));
                % Test that the field is itself a struct
                testCase.verifyEqual(isstruct(output.rewardFiring.(field{1})), true, sprintf('Outputted field "%s" is not a structure.', field{1}));
                % Test that the field contains a subfield for all rats being analysed
                for ratName = expectedRatNames
                    testCase.verifyEqual(isfield(output.rewardFiring.(field{1}), ratName{1}), true, sprintf('Outputted structure "%s" does not hold expected field for "%s"', field{1}, ratName{1}));
                end
            end
            for field = expectedNonStructFields
                % Test that the field exists
                testCase.verifyEqual(isfield(output.rewardFiring, field{1}), true, sprintf('Expected field "%s" missing.', field{1}));
            end
        end
        
        function [] = findPlottedMedian(testCase, lines, expectedMedian, expectedColour)
            
            [medianFound, matchesExpectedColour] = TestUtilities.findPlottedMedian(lines, expectedMedian, expectedColour);
            testCase.verifyTrue(medianFound, 'Figure does not contain a line indicating the median of the data.');
            testCase.verifyTrue(matchesExpectedColour, 'Figure contains a line indicating the median of the data, but it is not the expected colour.');
                    
        end


    
    end
    
    
    methods(TestMethodSetup)
        
        % Initialise RewardRelatedFiring object with temporary data path
        function createRewardRelatedFiring(testCase)
            testCase.RewardRelatedFiringObj = RewardRelatedFiring();
            % Create a temporary directory for test data
            testCase.TestDataDir = tempname;
            mkdir(testCase.TestDataDir);
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
            close all
        end
        
    end
    
    methods (Test)
        
        % ========================= iniitalising ========================= %
        
        % Test that the required data subdirectory is on the path
        function testDataSubdirectoryExists(testCase)
            paths = strsplit(path, pathsep);
            testCase.verifyTrue(any(strcmp(fullfile(pwd, './data'), paths)),...
                sprintf('Required "./data" directory is not on the MATLAB path'));
        end
        
        % Test that the main fields for the inclusion criteria are
        % created, even if no inclusion criteria are passed
        function testCreatesInclusionCriteria(testCase)
            obj = RewardRelatedFiring(nan);
            testCase.verifyTrue(isstruct(obj.inclusionCriteria));
            expectedFields = {'rats', 'sessions', 'cells'};
            for iField = 1:length(expectedFields)
                testCase.verifyTrue(isfield(obj.inclusionCriteria, expectedFields{iField})),
            end
        end
        
        % Test that the subfields (for the second field, sessions) are created,
        % even if no inclusion criteria are passed
        function testCreatesSessionInclusionCriteria(testCase)
            obj = RewardRelatedFiring(nan);
            expectedFields = {'skipEarliest', 'significantPerformanceThreshold'};
            for iField = 1:length(expectedFields)
                testCase.verifyTrue(isfield(obj.inclusionCriteria.sessions, expectedFields{iField})),
            end
        end
        
        % Test that if non-default inclusion criterion is passed, it gets
        % included correctly (for the first subfield, sessions.initialLearningOnly)
        function testIncludesInitialLearningTrue(testCase)
            obj = RewardRelatedFiring(struct('sessions', struct('initialLearningOnly', true)));
            testCase.verifyEqual(obj.inclusionCriteria.sessions.initialLearningOnly, true);
        end
        
        % Test that if non-default inclusion criterion is passed, it gets
        % included correctly (for the first subfield, sessions.initialLearningOnly)
        function testIncludesInitialLearningFalse(testCase)
            obj = RewardRelatedFiring(struct('sessions', struct('initialLearningOnly', false)));
            testCase.verifyEqual(obj.inclusionCriteria.sessions.initialLearningOnly, false);
        end
        
        
        % ===================== selecting sessions ===================== %
        
        % Test that the initial learning sessions contain the right number
        % of sessions if this parameter is 'true'
        function testSelectSessions_skipsEarliestSessions(testCase)
            obj = RewardRelatedFiring(struct('sessions', struct('skipEarliest', true)));
            obj = obj.selectSessions();
            testCase.verifyEqual(obj.initialLearningSessions.Quirinius, [3:12]);
        end
        
         % Test that the initial learning sessions contain the right number
        % of sessions if this parameter is 'false'
        function testSelectSessions_includesAllInitialSessions(testCase)
            obj = RewardRelatedFiring(struct('sessions', struct('skipEarliest', false)));
            obj = obj.selectSessions();
            testCase.verifyEqual(obj.initialLearningSessions.Quirinius, [1:12]);
        end
        
        
        % =================== selecting high firing rates =================== %    
        
        % Test that all cells are included when all firing rates are slightly above
        % threshold
        function testSelectFiringRate_allCellsAboveThreshold(testCase)
            
            % Create uniform binned firing rates of 10 Hz
            meanFR = 10;
            nCells = 5;
            binnedfr = repmat(meanFR, nCells, 100);
            
            % Select cells based on minimum firing rate of 9.9 Hz
            testCase.RewardRelatedFiringObj.inclusionCriteria.cells.minimumFiringRate = 9.9; 
            includeCells = testCase.RewardRelatedFiringObj.selectHighFiringRateCells(binnedfr);
            
            expectedOutput = true(1, nCells);
            testCase.verifyEqual(includeCells, expectedOutput)
        end
        
        % Test that no cells are included when all firing rates are
        % slightly below threshold
        function testSelectFiringRate_allCellsBelowThreshold(testCase)
            
            % Create uniform binned firing rates of 10 Hz
            meanFR = 10;
            nCells = 5;
            binnedfr = repmat(meanFR, nCells, 100);
            
            % Select cells based on minimum firing rate of 10.1 Hz
            testCase.RewardRelatedFiringObj.inclusionCriteria.cells.minimumFiringRate = 10.1; 
            includeCells = testCase.RewardRelatedFiringObj.selectHighFiringRateCells(binnedfr);
            
            expectedOutput = false(1, nCells);
            testCase.verifyEqual(includeCells, expectedOutput)
        end
        
        % Test a mix of cells, some above and some below threshold
        function testSelectFiringRate_mixedCells(testCase)
            
            % Create uniform binned firing rates of 10 Hz
            meanFR = [5, 50]';
            nCells = 3;
            binnedfr = repmat(meanFR, nCells, 100);
            
            % Select cells based on minimum firing rate of 10.1 Hz
            testCase.RewardRelatedFiringObj.inclusionCriteria.cells.minimumFiringRate = 10.1; 
            includeCells = testCase.RewardRelatedFiringObj.selectHighFiringRateCells(binnedfr);
            
            expectedOutput = [false, true, false, true, false, true];
            testCase.verifyEqual(includeCells, expectedOutput)
        end
                
        
        % ======================= analysing firing ======================= %    
        
        % Test basic functionality with one rat, one session, one rewarded trial, one cell
        function testGetRewardRelatedFiring_basicFunctionality(testCase)
            
            % Mock inclusion criteria and sessions
            testCase.RewardRelatedFiringObj.inclusionCriteria.rats = {'A_rat'};
            testCase.RewardRelatedFiringObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
            
            % Create mock behavioural data
            behavDirectory = fullfile(testCase.TestDataDir, 'test_session');
            ratName = 'A_Rat';
            n_trials = [1];
            actions = {[1]};
            rewarded = {[1]};
            arm_values = {{{'Legitimate', 'Legitimate', 'Optimal'}}};
            reward_probs = {{{'low', 'medium', 'high'}}};
            TestUtilities.saveMockBehaviouralData(...
                behavDirectory,...
                ratName,...
                'n_trials', n_trials,...
                'actions', actions,...
                'rewarded', rewarded,...
                'arm_values', arm_values,...
                'reward_probs', reward_probs);
            
            % Mock the rest of the data
            ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
            spiketimes = {[0.5, 1.5, 2.5]};
            CPExitTimes = [1];
            RarrivalTimes = [2];
            nNAcUnits = 1;
            binnedfr = ones(1, 100);
            TestUtilities.saveMockEphysData(ephysDirectory,...
                'spiketimes', spiketimes,...
                'CPexitTimes', CPExitTimes,...
                'RarrivalTimes', RarrivalTimes,...
                'nNAcUnits', nNAcUnits,...
                'binnedfr', binnedfr);
            
            % Call getRewardRelatedFiring
            output = testCase.RewardRelatedFiringObj.getRewardRelatedFiring('behaviouralDataPath', behavDirectory, 'ephysDataPath', testCase.TestDataDir);
            
            % Verify rewardFiring structure
            testCase.verifyNotEmpty(output);
            [perCellFields, stackedFields, otherFields] = testCase.RewardFiringExpectedFields();
            expectedStructFields = [perCellFields, stackedFields];
            testCase.verifyRewardFiringStructure(testCase, output, expectedStructFields, otherFields, {'A_rat'});    
        end
        
        % Test with one rat, one session, one trial of each type (high, medium, rewarded, unrewarded), 
        % and two cells of constant firing rates 10 Hz and 20 Hz
        function testGetRewardRelatedFiring_constantFiringRates(testCase)
            
            % Mock inclusion criteria and sessions
            testCase.RewardRelatedFiringObj.inclusionCriteria.cells.reward_responsive = false;
            testCase.RewardRelatedFiringObj.inclusionCriteria.rats = {'A_rat'};
            testCase.RewardRelatedFiringObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
            
            % Create mock behavioural data (1 session, 4 trials)
            [behavDirectory, n_trials, ~, ~] = testCase.createMockBehaviouralDataSingleSession(testCase, 'A_rat');
            
            % Mock the rest of the data (2 striatal neurons)
            ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
            CPExitTimes = [6, 10, 14, 18];
            RarrivalTimes = [7, 11, 15, 19];
            CPEntryTimes = [9, 13, 17, 21];
            periRewardFiringRates = [10, 20];
            periRewardPeriodDuration = 24;
            spiketimes = {...
                [0 : 1/periRewardFiringRates(1) : periRewardPeriodDuration],...
                [0 : 1/periRewardFiringRates(2) : periRewardPeriodDuration]};
            nNAcUnits = 2;
            rng(3);
            variableBinnedFR = randn(nNAcUnits, 100000);                  % Add additional Gaussian firing rates to check z-scored outputs are correct
            meanBinnedFR = mean(variableBinnedFR, 2);
            assert(isequal(meanBinnedFR > 0, ones(nNAcUnits, 1)))  % Ensure mean firing rates are above 0, for both cells to be included
            binnedfr = variableBinnedFR;
            TestUtilities.saveMockEphysData(ephysDirectory,...
                'spiketimes', spiketimes,...
                'CPexitTimes', CPExitTimes,...
                'RarrivalTimes', RarrivalTimes,...
                'CPEntryTimes', CPEntryTimes,...
                'nNAcUnits', nNAcUnits,...
                'binnedfr', binnedfr);
            
            % Call getRewardRelatedFiring
            output = testCase.RewardRelatedFiringObj.getRewardRelatedFiring('behaviouralDataPath', behavDirectory, 'ephysDataPath', testCase.TestDataDir);
            
            % Verify output
            [perCellFields, stackedFields, otherFields] = testCase.RewardFiringExpectedFields();
            expectedNBins = 2 / testCase.RewardRelatedFiringObj.BinSize;
            expectedSizePerCell = [nNAcUnits, expectedNBins];
            expectedSizeStacked = [4, expectedNBins];
            expectedFiringRatePerCell = periRewardFiringRates';
            expectedFiringRateStacked = repelem(expectedFiringRatePerCell, 2);
            expectedStructFields = [perCellFields, stackedFields];
            testCase.verifyRewardFiringStructure(testCase, output, expectedStructFields, otherFields, {'A_rat'});
            for field = perCellFields
                actualData = output.rewardFiring.(field{1}).A_rat;
                actualDataSize = size(actualData);
                % Test that the data is the expected size given number of neurons
                testCase.verifyEqual(actualDataSize, expectedSizePerCell, sprintf('Field "%s" is the wrong shape.', field{1}));
                % Test that the data holds the right values assuming constant, z-scored firing rates
                testCase.verifyEqual(actualData, repmat(expectedFiringRatePerCell, 1, expectedNBins), 'AbsTol', 1e-1, sprintf('Field "%s" does not contain the expected values.', field{1}));
            end
            for field = stackedFields
                actualData = output.rewardFiring.(field{1}).A_rat;
                actualDataSize = size(actualData);
                % Test that the data is the expected size given number of neurons and trials
                testCase.verifyEqual(actualDataSize, expectedSizeStacked, sprintf('Field "%s" is the wrong shape.', field{1}));
                % Test that the data holds the right values assuming constant, z-scored firing rates
                testCase.verifyEqual(actualData, repmat(expectedFiringRateStacked, 1, expectedNBins), 'AbsTol', 1e-1, sprintf('Field "%s" does not contain the expected values.', field{1}));
            end
            
        end
        
        % Test with one trial of each type, 
        % and two cells of constant firing rates, with only one above minimum firing rate
        function testGetRewardRelatedFiring_firingRateBelowThreshold(testCase)
            
            % Mock inclusion criteria and sessions
            testCase.RewardRelatedFiringObj.inclusionCriteria.cells.reward_responsive = false;
            testCase.RewardRelatedFiringObj.inclusionCriteria.rats = {'A_rat'};
            testCase.RewardRelatedFiringObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
            
            % Create mock behavioural data (1 session, 4 trials)
            [behavDirectory, n_trials, ~, ~] = testCase.createMockBehaviouralDataSingleSession(testCase, 'A_rat');
            
            % Mock the rest of the data (2 striatal neurons)
            ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
            CPExitTimes = [6, 10, 14, 18];
            RarrivalTimes = [7, 11, 15, 19];
            CPEntryTimes = [9, 13, 17, 21];
            periRewardFiringRates = [10, 20];
            periRewardPeriodDuration = 24;
            spiketimes = {...
                [0 : 1/periRewardFiringRates(1) : periRewardPeriodDuration],...
                [0 : 1/periRewardFiringRates(2) : periRewardPeriodDuration]};
            nNAcUnits = 2;
            rng(3);
            variableBinnedFR = randn(nNAcUnits, 100000);                  % Add additional Gaussian firing rates to check z-scored outputs are correct
            meanBinnedFR = mean(variableBinnedFR, 2);
            assert(isequal(meanBinnedFR > 0, ones(nNAcUnits, 1)))     % Ensure mean firing rates are above 0, for both cells to be included
            assert(meanBinnedFR(1) > meanBinnedFR(2))                   % Ensure that only the first cell's mean firing rate is above the minimum threshold
            testCase.RewardRelatedFiringObj.inclusionCriteria.cells.minimumFiringRate = mean(meanBinnedFR);
            binnedfr = variableBinnedFR;
            TestUtilities.saveMockEphysData(ephysDirectory,...
                'spiketimes', spiketimes,...
                'CPexitTimes', CPExitTimes,...
                'RarrivalTimes', RarrivalTimes,...
                'CPEntryTimes', CPEntryTimes,...
                'nNAcUnits', nNAcUnits,...
                'binnedfr', binnedfr);
            
            % Call getRewardRelatedFiring
            output = testCase.RewardRelatedFiringObj.getRewardRelatedFiring('behaviouralDataPath', behavDirectory, 'ephysDataPath', testCase.TestDataDir);
            
            % Verify output
            [perCellFields, stackedFields, otherFields] = testCase.RewardFiringExpectedFields();
            expectedNBins = 2 / testCase.RewardRelatedFiringObj.BinSize;
            expectedSizePerCell = [1, expectedNBins];
            expectedSizeStacked = [2, expectedNBins];
            expectedFiringRatePerCell = periRewardFiringRates(1);
            expectedFiringRateStacked = repelem(expectedFiringRatePerCell, 2, 1);
            expectedStructFields = [perCellFields, stackedFields];
            testCase.verifyRewardFiringStructure(testCase, output, expectedStructFields, otherFields, {'A_rat'});
            for field = perCellFields
                actualData = output.rewardFiring.(field{1}).A_rat;
                actualDataSize = size(actualData);
                % Test that the data is the expected size given number of neurons
                testCase.verifyEqual(actualDataSize, expectedSizePerCell, sprintf('Field "%s" is the wrong shape.', field{1}));
                % Test that the data holds the right values assuming constant, z-scored firing rates
                testCase.verifyEqual(actualData, repmat(expectedFiringRatePerCell, 1, expectedNBins), 'AbsTol', 1e-1, sprintf('Field "%s" does not contain the expected values.', field{1}));
            end
            for field = stackedFields
                actualData = output.rewardFiring.(field{1}).A_rat;
                actualDataSize = size(actualData);
                % Test that the data is the expected size given number of neurons and trials
                testCase.verifyEqual(actualDataSize, expectedSizeStacked, sprintf('Field "%s" is the wrong shape.', field{1}));
                % Test that the data holds the right values assuming constant, z-scored firing rates
                testCase.verifyEqual(actualData, repmat(expectedFiringRateStacked, 1, expectedNBins), 'AbsTol', 1e-1, sprintf('Field "%s" does not contain the expected values.', field{1}));
            end
            
        end
            
        % Test with one rat, one session, one trial of each type (high, medium, rewarded, unrewarded), 
        % and two cells of variable firing rates
        function testGetRewardRelatedFiring_variableFiringRates(testCase)
            
            % Mock inclusion criteria and sessions
            testCase.RewardRelatedFiringObj.inclusionCriteria.cells.reward_responsive = false;
            testCase.RewardRelatedFiringObj.inclusionCriteria.rats = {'A_rat'};
            testCase.RewardRelatedFiringObj.sessionsForSpecificAnalysis = struct('A_rat', [1]);
            
            % Create mock behavioural data (1 session, 4 trials)
            [behavDirectory, n_trials, actions, rewarded] = testCase.createMockBehaviouralDataSingleSession(testCase, 'A_rat');
            
            % Mock the rest of the data (2 striatal neurons)
            ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
            CPExitTimes = [6, 16, 26, 36];
            RarrivalTimes = [7, 17, 27, 37];
            CPEntryTimes = [9, 19, 29, 39];
            periRewardFiringRates = [actions{1}*10; rewarded{1}*10 + 5];       % Neuron 1: higher firing when action=2 than action=1. Neuron 2: higher firing when rewarded than not rewarded
            nNAcUnits = 2;
            spiketimes = cell(1, nNAcUnits);
            for iTrial = 1:n_trials(1)
                spikingStart = RarrivalTimes(iTrial) - 7;
                spikingEnd = RarrivalTimes(iTrial) + 5;
                for iCell = 1:nNAcUnits
                    fr = periRewardFiringRates(iCell, iTrial);
                    spiketimes{iCell} = [spiketimes{iCell} spikingStart : 1/fr : spikingEnd];
                end
            end
            rng(3);
            variableBinnedFR = randn(nNAcUnits, 100000);               % Add additional Gaussian firing rates to check z-scored outputs are correct
            meanBinnedFR = mean(variableBinnedFR, 2);
            assert(isequal(meanBinnedFR > 0, ones(nNAcUnits, 1)))  % Ensure mean firing rates are above 0, for both cells to be included
            binnedfr = variableBinnedFR;
            TestUtilities.saveMockEphysData(ephysDirectory,...
                'spiketimes', spiketimes,...
                'CPexitTimes', CPExitTimes,...
                'RarrivalTimes', RarrivalTimes,...
                'CPEntryTimes', CPEntryTimes,...
                'nNAcUnits', nNAcUnits,...
                'binnedfr', binnedfr);
            
            % Call getRewardRelatedFiring
            output = testCase.RewardRelatedFiringObj.getRewardRelatedFiring('behaviouralDataPath', behavDirectory, 'ephysDataPath', testCase.TestDataDir);
            
            % Verify output
            [perCellFields, stackedFields, otherFields] = testCase.RewardFiringExpectedFields();
            expectedNBins = 2 / testCase.RewardRelatedFiringObj.BinSize;
            expectedOutput = struct(...
                'perCellHighRewardExpectation', mean(periRewardFiringRates(:, find(actions{1}==1)), 2),...
                'perCellMediumRewardExpectation', mean(periRewardFiringRates(:, find(actions{1}==2)), 2),...
                'perCellRewarded', mean(periRewardFiringRates(:, find(rewarded{1}==1)), 2),...
                'perCellUnrewarded', mean(periRewardFiringRates(:, find(rewarded{1}==0)), 2),...
                'highRewardExpectation', reshape(periRewardFiringRates(:, find(actions{1}==1))', [], 1),...
                'mediumRewardExpectation', reshape(periRewardFiringRates(:, find(actions{1}==2))', [], 1),...
                'rewarded', reshape(periRewardFiringRates(:, find(rewarded{1}==1))', [], 1),...
                'unrewarded', reshape(periRewardFiringRates(:, find(rewarded{1}==0))', [], 1),...
                'ts', (-2 - testCase.RewardRelatedFiringObj.BinSize/2) : testCase.RewardRelatedFiringObj.BinSize : 0);
            expectedStructFields = [perCellFields, stackedFields];
            testCase.verifyRewardFiringStructure(testCase, output, expectedStructFields, otherFields, {'A_rat'});
            for field = fieldnames(expectedOutput)
                testCase.verifyEqual(...
                    output.rewardFiring.(field{1}).A_rat, ...
                    repmat(expectedOutput.(field{1}), 1, expectedNBins),...
                    'AbsTol', 0.1, sprintf('Field "%s" does not contain the expected values.', field{1}));
            end            
            
        end
        
        % Test with two rats, three sessions, one trial of each type (high, medium, rewarded, unrewarded), 
        % and two cells of variable firing rates
        function testGetRewardRelatedFiring_multipleSessions(testCase)
            
            % Mock inclusion criteria and sessions
            testCase.RewardRelatedFiringObj.inclusionCriteria.cells.reward_responsive = false;
            testCase.RewardRelatedFiringObj.inclusionCriteria.rats = {'A_rat', 'B_rat'};
            testCase.RewardRelatedFiringObj.sessionsForSpecificAnalysis = struct('A_rat', [1], 'B_rat', [2, 3]);
            
            % Create mock behavioural data for rat A (1 session, 4 trials)
            [behavDirectory, n_trials_A, actions_A, rewarded_A] = testCase.createMockBehaviouralDataSingleSession(testCase, 'A_rat');
            
            % Mock the rest of the data (2 striatal neurons)
            ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
            CPExitTimes = [6, 16, 26, 36];
            RarrivalTimes = [7, 17, 27, 37];
            CPEntryTimes = [9, 19, 29, 39];
            periRewardFiringRates = [actions_A{1}*10; rewarded_A{1}*10 + 5];       % Neuron 1: higher firing when action=2 than action=1. Neuron 2: higher firing when rewarded than not rewarded
            nNAcUnits = 2;
            spiketimes = cell(1, nNAcUnits);
            for iTrial = 1:n_trials_A(1)
                spikingStart = RarrivalTimes(iTrial) - 7;
                spikingEnd = RarrivalTimes(iTrial) + 5;
                for iCell = 1:nNAcUnits
                    fr = periRewardFiringRates(iCell, iTrial);
                    spiketimes{iCell} = [spiketimes{iCell} spikingStart : 1/fr : spikingEnd];
                end
            end
            rng(3);
            variableBinnedFR = randn(nNAcUnits, 100000);               % Add additional Gaussian firing rates to check z-scored outputs are correct
            meanBinnedFR = mean(variableBinnedFR, 2);
            assert(isequal(meanBinnedFR > 0, ones(nNAcUnits, 1)))  % Ensure mean firing rates are above 0, for both cells to be included
            binnedfr = variableBinnedFR;
            TestUtilities.saveMockEphysData(ephysDirectory,...
                'spiketimes', spiketimes,...
                'CPexitTimes', CPExitTimes,...
                'RarrivalTimes', RarrivalTimes,...
                'CPEntryTimes', CPEntryTimes,...
                'nNAcUnits', nNAcUnits,...
                'binnedfr', binnedfr);
            
            % Create mock behavioural data for rat B (3 sessions)
            [behavDirectory, n_trials_B, actions_B, rewarded_B] = testCase.createMockBehaviouralDataThreeSessions(testCase, 'B_rat');
            
            % Mock the rest of the data (2 striatal neurons in each session)
            allPeriRewardFiringRates = periRewardFiringRates;
            for iSession = 1:3
                ephysDirectory = fullfile(testCase.TestDataDir, 'Rat_B', ['Session' num2str(iSession)]);
                CPExitTimes = [6, 16, 26, 36, 46, 56]; CPExitTimes = CPExitTimes(1:n_trials_B(iSession));
                RarrivalTimes = [7, 17, 26, 36, 46, 56]; RarrivalTimes = RarrivalTimes(1:n_trials_B(iSession));
                CPEntryTimes = [9, 19, 29, 39, 49, 59]; CPEntryTimes = CPEntryTimes(1:n_trials_B(iSession));
                periRewardFiringRates = [actions_B{iSession}*10; rewarded_B{iSession}*10 + 5];
                if iSession ~= 1
                    allPeriRewardFiringRates = [allPeriRewardFiringRates periRewardFiringRates];
                end
                nNAcUnits = 2;
                spiketimes = cell(1, nNAcUnits);
                for iTrial = 1:n_trials_B(iSession)
                    spikingStart = RarrivalTimes(iTrial) - 7;
                    spikingEnd = RarrivalTimes(iTrial) + 5;
                    for iCell = 1:nNAcUnits
                        fr = periRewardFiringRates(iCell, iTrial);
                        spiketimes{iCell} = [spiketimes{iCell} spikingStart : 1/fr : spikingEnd];
                    end
                end
                rng(3);
                variableBinnedFR = randn(nNAcUnits, 100000);               % Add additional Gaussian firing rates to check z-scored outputs are correct
                meanBinnedFR = mean(variableBinnedFR, 2);
                assert(isequal(meanBinnedFR > 0, ones(nNAcUnits, 1)))  % Ensure mean firing rates are above 0, for both cells to be included
                binnedfr = variableBinnedFR;
                TestUtilities.saveMockEphysData(ephysDirectory,...
                    'spiketimes', spiketimes,...
                    'CPexitTimes', CPExitTimes,...
                    'RarrivalTimes', RarrivalTimes,...
                    'CPEntryTimes', CPEntryTimes,...
                    'nNAcUnits', nNAcUnits,...
                    'binnedfr', binnedfr);
            end
            
            % Call getRewardRelatedFiring
            output = testCase.RewardRelatedFiringObj.getRewardRelatedFiring('behaviouralDataPath', behavDirectory, 'ephysDataPath', testCase.TestDataDir);
            
            % Verify output
            [perCellFields, stackedFields, otherFields] = testCase.RewardFiringExpectedFields();
            expectedNBins = 2 / testCase.RewardRelatedFiringObj.BinSize;
            expectedOutput = struct(...
                'perCellHighRewardExpectation', struct(...
                    'A_rat', mean(allPeriRewardFiringRates(:, find(cell2mat(actions_A)==1)), 2),...
                    'B_rat', mean(allPeriRewardFiringRates(:, find(cell2mat(actions_B(2:3))==1)), 2)),...
                'perCellMediumRewardExpectation', struct(...
                    'A_rat', mean(allPeriRewardFiringRates(:, find(cell2mat(actions_A)==2)), 2),...
                    'B_rat', mean(allPeriRewardFiringRates(:, find(cell2mat(actions_B(2:3))==2)), 2)),...
                'perCellRewarded', struct(...
                    'A_rat', mean(allPeriRewardFiringRates(:, find(cell2mat(rewarded_A)==1)), 2),...
                    'B_rat', mean(allPeriRewardFiringRates(:, find(cell2mat(rewarded_B(2:3))==1)), 2)),...
                'perCellUnrewarded', struct(...
                    'A_rat', mean(allPeriRewardFiringRates(:, find(cell2mat(rewarded_A)==0)), 2),...
                    'B_rat', mean(allPeriRewardFiringRates(:, find(cell2mat(rewarded_B(2:3))==0)), 2)),...
                'highRewardExpectation', struct(...
                    'A_rat', reshape(allPeriRewardFiringRates(:, find(cell2mat(actions_A)==1))', [], 1),...
                    'B_rat', reshape(allPeriRewardFiringRates(:, cell2mat(actions_B(2:3))==1)', [], 1)),...
                'mediumRewardExpectation', struct(...
                    'A_rat', reshape(allPeriRewardFiringRates(:, find(cell2mat(actions_A)==2))', [], 1),...
                    'B_rat', reshape(allPeriRewardFiringRates(:, find(cell2mat(actions_B(2:3))==2))', [], 1)),...
                'rewarded', struct(...
                    'A_rat', reshape(allPeriRewardFiringRates(:, find(cell2mat(rewarded_A)==1))', [], 1),...
                    'B_rat', reshape(allPeriRewardFiringRates(:, find(cell2mat(rewarded_B(2:3))==1))', [], 1)),...
                'unrewarded', struct(...
                    'A_rat', reshape(allPeriRewardFiringRates(:, find(cell2mat(rewarded_A)==0))', [], 1),...
                    'B_rat', reshape(allPeriRewardFiringRates(:, find(cell2mat(rewarded_B(2:3))==0))', [], 1)),...
                'ts', (-2 - testCase.RewardRelatedFiringObj.BinSize/2) : testCase.RewardRelatedFiringObj.BinSize : 0);
            testCase.verifyRewardFiringStructure(testCase, output, [perCellFields, stackedFields], otherFields, {'A_rat', 'B_rat'});
            for field = fieldnames(expectedOutput)
                testCase.verifyEqual(...
                    output.rewardFiring.(field{1}).A_rat, ...
                    repmat(expectedOutput.(field{1}).A_rat, 1, expectedNBins),...
                    'AbsTol', 0.1, sprintf('Field "%s" does not contain the expected values for rat A.', field{1}));
                testCase.verifyEqual(...
                    output.rewardFiring.(field{1}).B_rat, ...
                    repmat(expectedOutput.(field{1}).B_rat, 2, expectedNBins),...
                    'AbsTol', 0.1, sprintf('Field "%s" does not contain the expected values for rat B.', field{1}));
            end    
            
            
        end

        
        
        % ========================== plotting ========================== %
        
        % Test that plotting of reward-related firing produces a figure
        function testRewardPlotting_basicFunctionality(testCase)
            
            % Mock inclusion criteria
            testCase.RewardRelatedFiringObj.inclusionCriteria.rats = {'A_rat'};

            % Mock reward firing data
            rng(0);
            testCase.RewardRelatedFiringObj.rewardFiring = struct('perCellMediumRewardExpectation', struct('A_rat', []), 'perCellHighRewardExpectation', struct('A_rat', []));
            testCase.RewardRelatedFiringObj.rewardFiring.perCellMediumRewardExpectation.A_rat = random('Normal', 1, 1, 100, 200);
            testCase.RewardRelatedFiringObj.rewardFiring.perCellHighRewardExpectation.A_rat = random('Normal', 1, 1, 100, 200);

            % Call plotRewardRelatedFiring
            output = testCase.RewardRelatedFiringObj.plotRewardRelatedFiring();

            % Verify that a figure is created
            testCase.verifyNotEmpty(output.figures, 'No figure handle stored in stored in obj.figures.');
            testCase.verifyClass(output.figures{end}, 'matlab.ui.Figure', 'The stored object is not a Matlab figure.');
        
        end
        
        % Test that plotting of reward-related firing produces a correct figure
        % with data from one rat
        function testRewardPlotting_dataSingleRat(testCase)
            
            % Mock inclusion criteria
            testCase.RewardRelatedFiringObj.inclusionCriteria.rats = {'A_rat'};

            % Mock reward firing data
            rng(0);
            testCase.RewardRelatedFiringObj.rewardFiring = struct(...
                'perCellMediumRewardExpectation', struct('A_rat', random('Normal', 1, 1, 100, 200)),...
                'perCellHighRewardExpectation', struct('A_rat', random('Normal', 1, 1, 100, 200)));
            meanMediumFR = mean(testCase.RewardRelatedFiringObj.rewardFiring.perCellMediumRewardExpectation.A_rat);
            meanHighFR = mean(testCase.RewardRelatedFiringObj.rewardFiring.perCellHighRewardExpectation.A_rat);

            % Call plotRewardRelatedFiring
            output = testCase.RewardRelatedFiringObj.plotRewardRelatedFiring();
            
            % Verify that there is one figure
            nFigures = length(output.figures);
            testCase.verifyEqual(nFigures, 1, sprintf('%i figures were created when there should be 1', nFigures))

            % Verify that the figure has two axes
            axes = findobj(output.figures{1}, 'Type', 'axes');
            nAxes = numel(axes);
            testCase.verifyEqual(nAxes, 2, sprintf('%i axes were created when there should be 2', nAxes));
            
            % Verify that each axis contains two Lines and a Patch
            for iAxis = 1:numel(axes)
                axisCh = get(axes(iAxis), 'Children');
                testCase.verifyEqual(numel(axisCh), 3);
                testCase.verifyEqual(class(axisCh(1)), 'matlab.graphics.primitive.Line', sprintf('First child of axis %i expected to be a line (vertical dashed line at zero) but it wasn''t', iAxis));
                testCase.verifyEqual(class(axisCh(2)), 'matlab.graphics.chart.primitive.Line', sprintf('Second child of axis %i expected to be a line (mean of firing rates) but it wasn''t', iAxis));
                testCase.verifyEqual(class(axisCh(3)), 'matlab.graphics.primitive.Patch', sprintf('Third child of axis %i expected to be a patch (mean +- SEM of firing rates) but it wasn''t', iAxis));
            end
            
            % Verify that both axes have an x-label
            for iAxis = 1:numel(axes)
                label = get(axes(iAxis), 'xlabel');
                testCase.verifyNotEmpty(label.String, sprintf('Axis %i should have an x-label', iAxis));
            end
            
            % Verify that the left axis has a y-label
            label = get(axes(end), 'ylabel');
            testCase.verifyNotEmpty(label.String, 'Axis 1 should have a y-label');
            
            % Check that the plotted mean for medium firing is correct
            children = get(output.figures{1}, 'Children');
            axisCh = get(children(2), 'Children');
            testCase.verifyEqual(axisCh(2).YData, meanMediumFR, 'Mean firing rate plotted on the axis does not correspond to the expected values for medium-reward-expectation')

            % Check that the plotted mean for high firing is correct
            axisCh = get(children(1), 'Children');
            testCase.verifyEqual(axisCh(2).YData, meanHighFR, 'Mean firing rate plotted on the axis does not correspond to the expected values for high-reward-expectation')
        
        end
        
        
        % Test that plotting of reward-related firing produces a correct figure
        % with data from three rats
        function testRewardPlotting_dataMultipleRats(testCase)
            
            % Mock inclusion criteria
            testCase.RewardRelatedFiringObj.inclusionCriteria.rats = {'A_rat', 'B_rat', 'C_rat'};

            % Mock reward firing data
            rng(0);
            testCase.RewardRelatedFiringObj.rewardFiring = struct(...
                'perCellMediumRewardExpectation', struct('A_rat', random('Normal', 1, 1, 100, 200), 'B_rat', random('Normal', 5, 1, 100, 200), 'C_rat', random('Normal', 10, 1, 100, 200)),...
                'perCellHighRewardExpectation', struct('A_rat', random('Normal', 2, 1, 100, 200), 'B_rat', random('Normal', 10, 1, 100, 200), 'C_rat', random('Normal', 20, 1, 100, 200)));
            meanMediumFR = struct(...
                'A_rat', mean(testCase.RewardRelatedFiringObj.rewardFiring.perCellMediumRewardExpectation.A_rat),...
                'B_rat', mean(testCase.RewardRelatedFiringObj.rewardFiring.perCellMediumRewardExpectation.B_rat),...
                'C_rat', mean(testCase.RewardRelatedFiringObj.rewardFiring.perCellMediumRewardExpectation.C_rat));
            meanHighFR = struct(...
                'A_rat', mean(testCase.RewardRelatedFiringObj.rewardFiring.perCellHighRewardExpectation.A_rat),...
                'B_rat', mean(testCase.RewardRelatedFiringObj.rewardFiring.perCellHighRewardExpectation.B_rat),...
                'C_rat', mean(testCase.RewardRelatedFiringObj.rewardFiring.perCellHighRewardExpectation.C_rat));
            
            % Call plotRewardRelatedFiring
            output = testCase.RewardRelatedFiringObj.plotRewardRelatedFiring();
            
            % Verify that there are three figures
            nFigures = length(output.figures);
            testCase.verifyEqual(nFigures, 3, sprintf('%i figure(s) were created when there should be 3', nFigures))
            
            for iFig = 1:nFigures
                
                % Verify that the figure has two axes
                axes = findobj(output.figures{iFig}, 'Type', 'axes');
                nAxes = numel(axes);
                testCase.verifyEqual(nAxes, 2, sprintf('%i axes were created when there should be 2 (figure %i)', nAxes, iFig));

                % Verify that each axis contains two Lines and a Patch
                for iAxis = 1:numel(axes)
                    axisCh = get(axes(iAxis), 'Children');
                    testCase.verifyEqual(numel(axisCh), 3);
                    testCase.verifyEqual(class(axisCh(1)), 'matlab.graphics.primitive.Line', sprintf('First child of axis %i expected to be a line (vertical dashed line at zero) but it wasn''t (figure %i)', iAxis, iFig));
                    testCase.verifyEqual(class(axisCh(2)), 'matlab.graphics.chart.primitive.Line', sprintf('Second child of axis %i expected to be a line (mean of firing rates) but it wasn''t (figure %i)', iAxis, iFig));
                    testCase.verifyEqual(class(axisCh(3)), 'matlab.graphics.primitive.Patch', sprintf('Third child of axis %i expected to be a patch (mean +- SEM of firing rates) but it wasn''t (figure %i)', iAxis, iFig));
                end
                
                % Verify that both axes have an x-label
                for iAxis = 1:numel(axes)
                    label = get(axes(iAxis), 'xlabel');
                    testCase.verifyNotEmpty(label.String, sprintf('Axis %i should have an x-label (figure %i)', iAxis, iFig));
                end

                % Verify that the left axis has a y-label
                label = get(axes(end), 'ylabel');
                testCase.verifyNotEmpty(label.String, sprintf('Axis 1 should have a y-label (figure %i)', iFig));

                % Check that the plotted mean for medium firing is correct
                children = get(output.figures{iFig}, 'Children');
                ratName = testCase.RewardRelatedFiringObj.inclusionCriteria.rats{iFig};
                axisCh = get(children(2), 'Children');
                testCase.verifyEqual(axisCh(2).YData, meanMediumFR.(ratName), sprintf('Mean firing rate plotted on the axis for %s mock data does not correspond to the expected values for medium-reward-expectation', ratName));

                % Check that the plotted mean for high firing is correct
                axisCh = get(children(1), 'Children');
                testCase.verifyEqual(axisCh(2).YData, meanHighFR.(ratName), sprintf('Mean firing rate plotted on the axis for %s mock data does not correspond to the expected values for high-reward-expectation', ratName));
                
            end
            
        end
        
        
        % Test that plotting of reward-related firing runs without producing a figure
        % with data from no rats
        function testRewardPlotting_dataNoRats(testCase)
            
            % Mock inclusion criteria
            testCase.RewardRelatedFiringObj.inclusionCriteria.rats = {};
            
            % Mock reward firing data
            rng(0);
            testCase.RewardRelatedFiringObj.rewardFiring = struct(...
                'perCellMediumRewardExpectation', struct('missing_rat', random('Normal', 1, 1, 100, 200)),...
                'perCellHighRewardExpectation', struct('missing_rat', random('Normal', 1, 1, 100, 200)));
            
            % Call plotRewardRelatedFiring
            output = testCase.RewardRelatedFiringObj.plotRewardRelatedFiring();
            
            % Verify that there are three figures
            nFigures = length(output.figures);
            testCase.verifyEqual(nFigures, 0, sprintf('%i figure(s) were created when there should be none', nFigures))
            
        end
        
        
        % Test that plotting of reward-related firing produces an error
        % with no data from one rat
        function testRewardPlotting_noData(testCase)
            
            % Mock inclusion criteria
            testCase.RewardRelatedFiringObj.inclusionCriteria.rats = {'A_rat'};
            
            % Call plotRewardRelatedFiring
            expectedErrorId = 'RewardRelatedFiring:dataNotAssigned';
            testCase.verifyError(@() testCase.RewardRelatedFiringObj.plotRewardRelatedFiring(), expectedErrorId);
            
        end

        
        % Test basic functionality of plotBar
        function testPlotBar_basicFunctionality(testCase)
            
            % Mock data with known median and upper and lower quartiles           
            data = [1; 3; 5; 7];
            expectedQuartiles = [2, 4, 6];
            expectedMinMax = [1, 7];
            
            % Select a colour from the colour scheme
            colour_scheme;
            colournames = fieldnames(colourscheme);
            primarycolours = colournames(find(cellfun(@(x) ~isempty(x), strfind(colournames, '_primary'))));
            secondarycolours = colournames(find(cellfun(@(x) ~isempty(x), strfind(colournames, '_secondary'))));
            colour1 = NaN; colour2 = NaN;
            iColour = 0;
            while ~strcmp(colour1, colour2)
                iColour = iColour + 1;
                colour1 = strrep(primarycolours{iColour}, '_primary', '');
                colour2 = strrep(secondarycolours{iColour}, '_secondary', '');
            end
            
            % Call method
            f = figure;
            output = testCase.RewardRelatedFiringObj.plotBar(data, colour1, 1, f);
            
            % Assert that there are two lines (one for median, one for range/whiskers) and a patch (for interquartile range) on the figure
            lines = findobj(output, 'Type', 'line');
            patches = findobj(output, 'Type', 'patch');
            testCase.verifyEqual(numel(lines), 2, sprintf('Expected 2 line objects on the figure but there are %i.', numel(lines)));
            testCase.verifyEqual(numel(patches), 1, sprintf('Expected 1 patch object on the figure but there are %i.', numel(patches)));
            
            % Check that there is a line that represents the median in the secondary colour
            testCase.findPlottedMedian(testCase, lines, expectedQuartiles(2), colourscheme.(secondarycolours{iColour}))
            
            % Check that there is a line that represents the whisker 
            % from minimum to maximum of the data in black
            rangeFound = false; i = 0;
            while ~rangeFound && i < numel(lines)
                i = i + 1;
                lineData = lines(i).YData;
                if sort(lineData) == [expectedMinMax(1), expectedMinMax(2)];
                    rangeFound = true;
                end
            end
            testCase.verifyTrue(rangeFound, 'Figure does not contain a line indicating the data range.');
            rangeLineColor = lines(i).Color;
            testCase.verifyEqual(rangeLineColor, [0, 0, 0], 'Figure contains a line indicating the data range, but it is not the expected colour.');
                        
            % Check that there is a patch that represents the interquartile range
            patchYData = patches(1).Vertices(:, 2);
            expectedYData = [expectedQuartiles(1); expectedQuartiles(1); expectedQuartiles(3); expectedQuartiles(3)];
            testCase.verifyEqual(patchYData, expectedYData, 'Figure contains a patch, but it does not represent the upper and lower quartiles.');
            patchXData = patches(1).Vertices(:, 1);
            testCase.verifyEqual(numel(unique(patchXData)), 2, 'Figure contains a patch indicating the interquartile range, but it is not a rectangle.');
            
            % Check that the ylimits are appropriate for the data
            axes = findobj(output, 'Type', 'Axes');
            ylim = axes(1).YLim;
            testCase.verifyTrue(ylim(2) >= expectedMinMax(2), 'Axis ylimit is not extensive enough to show all the data.')
            testCase.verifyTrue(ylim(1) <= expectedMinMax(1), 'Axis ylimit is not extensive enough to show all the data.')
            
        end
            
        % Test that plotBar correctly adds new data to an axis that already
        % contains data
        function testPlotBar_addPlot(testCase)
            
            % Mock data with known median and upper and lower quartiles           
            data = [1; 3; 5; 7];
            expectedQuartiles = [2, 4, 6];
            expectedMinMax = [1, 7];
            
            % Select a colour from the colour scheme
            colour_scheme;
            colournames = fieldnames(colourscheme);
            primarycolours = colournames(find(cellfun(@(x) ~isempty(x), strfind(colournames, '_primary'))));
            secondarycolours = colournames(find(cellfun(@(x) ~isempty(x), strfind(colournames, '_secondary'))));
            colour1 = NaN; colour2 = NaN;
            iColour = 0;
            while ~strcmp(colour1, colour2)
                iColour = iColour + 1;
                colour1 = strrep(primarycolours{iColour}, '_primary', '');
                colour2 = strrep(secondarycolours{iColour}, '_secondary', '');
            end
            
            % Create a figure with pre-existing data much higher on the y-axis than the new data
            f = figure;
            patch([0, 1, 1, 0], [100, 100, 101, 101], 'k')
            line([0.5, 0.5], [0, 1])
            
            % Call method
            output = testCase.RewardRelatedFiringObj.plotBar(data, colour1, 2, f);
            
            % Check that all objects appear on the figure
            lines = findobj(output, 'Type', 'line');
            patches = findobj(output, 'Type', 'patch');
            testCase.verifyEqual(numel(lines), 3, sprintf('Expected 3 line objects on the figure but there are %i.', numel(lines)));
            testCase.verifyEqual(numel(patches), 2, sprintf('Expected 2 patch objects on the figure but there are %i.', numel(patches)));
            
            % Check that the ylimits are appropriate for the data
            axes = findobj(output, 'Type', 'Axes');
            ylim = axes(1).YLim;
            testCase.verifyTrue(ylim(2) >= 101, 'Axis ylimit is not extensive enough to show all the data.')
            testCase.verifyTrue(ylim(1) <= 0, 'Axis ylimit is not extensive enough to show all the data.')
            
        end
            
        % Test basic functionality of plotTTestSignificance with significant differences
        function testPlotTTestSigificance_significant(testCase)
            
            % Mock statistically significant random Gaussiannoise
            rng(0);
            data1 = random('Normal', 1, 1, 1, 100);
            data2 = random('Normal', 100, 1, 1, 100);
            f = figure();
            
            % Call plotTTestSigificance
            [h, p] = testCase.RewardRelatedFiringObj.plotTTestSignificance(data1, data2, [1, 2], f);
            
            % Verify outpout
            testCase.verifyNotEmpty(h, 'Function returned an empty figure handle.');
            lines = findobj(h, 'Type', 'line');
            text = findobj(h, 'Type', 'text');
            testCase.verifyEqual(numel(lines), 3, sprintf('Expected 3 line objects on the figure but there are %i.', numel(lines)));
            testCase.verifyEqual(numel(text), 1, sprintf('Expected 1 text object on the figure but there are %i.', numel(text)));
            testCase.verifyEqual(text(1).String, '*', sprintf('Expected the figure to show the text "*" but it shows %s.', text(1).String));
            
        end
        
        % Test basic functionality of plotTTestSignificance with no significant differences
        function testPlotTTestSigificance_notSignificant(testCase)
            
            % Mock statistically significant random Gaussiannoise
            rng(0);
            data1 = random('Normal', 1, 1, 1, 100);
            data2 = random('Normal', 1, 1, 1, 100);
            f = figure();
            
            % Call plotTTestSigificance
            [h, p] = testCase.RewardRelatedFiringObj.plotTTestSignificance(data1, data2, [1, 2], f);
            
            % Verify outpout
            testCase.verifyNotEmpty(h, 'Function returned an empty figure handle.');
            lines = findobj(h, 'Type', 'line');
            text = findobj(h, 'Type', 'text');
            testCase.verifyEqual(numel(lines), 3, sprintf('Expected 3 line objects on the figure but there are %i.', numel(lines)));
            testCase.verifyEqual(numel(text), 1, sprintf('Expected 1 text object on the figure but there are %i.', numel(text)));
            testCase.verifyEqual(text(1).String, 'n.s.', sprintf('Expected the figure to show the text "*" but it shows %s.', text(1).String));
            
        end
        
        % Test that plotting of peak fring rates produces a figure
        function testPlotPeakPreRewardFiring_basicFunctionality(testCase)
            
            % Mock inclusion criteria
            testCase.RewardRelatedFiringObj.inclusionCriteria.rats = {'A_rat'};

            % Mock reward firing data
            rng(0);
            testCase.RewardRelatedFiringObj.rewardFiring = struct(...
                'perCellMediumRewardExpectation', struct('A_rat', random('Normal', 1, 1, 10, 200)),...
                'perCellHighRewardExpectation', struct('A_rat', random('Normal', 1, 1, 10, 200)),...
                'ts', linspace(-10, 5, 200));

            % Call plotPeakPreRewardFiring
            output = testCase.RewardRelatedFiringObj.plotPeakPreRewardFiring();

            % Verify that a figure is created
            testCase.verifyNotEmpty(output.figures, 'No figure handle stored in stored in obj.figures.');
            testCase.verifyClass(output.figures{end}, 'matlab.ui.Figure', 'The stored object is not a Matlab figure.');
            
        end
        
        
        % Test that plotting of peak fring rates produces a correct figure
        % with data from one rat with no significant differences
        function testPlotPeakPreRewardFiring_dataSingleRatNoSignificance(testCase)
            
            % Mock inclusion criteria
            testCase.RewardRelatedFiringObj.inclusionCriteria.rats = {'A_rat'};
            
            % Mock reward firing data (zeros except for 2 seconds before
            % reward arrival time), 4 cells
            preRewardArrivalData = repmat([1; 3; 5; 7], 1, 40);
            firingRate = zeros(4, 200);
            firingRate(:, 60:99) = preRewardArrivalData;
            testCase.RewardRelatedFiringObj.rewardFiring = struct(...
                'perCellMediumRewardExpectation', struct('A_rat', firingRate),...
                'perCellHighRewardExpectation', struct('A_rat', firingRate*2),...
                'ts', linspace(-5, 5, 200));
            
             % Call plotPeakPreRewardFiring
            output = testCase.RewardRelatedFiringObj.plotPeakPreRewardFiring();
            
            % Verify that one figure is created
            nFigures = length(output.figures);
            testCase.verifyEqual(nFigures, 1, 'One figure should be created for peak pre-reward firing.');
            
            % Verify that the figure contains seven Lines and two Patches
            fig = output.figures{1};
            nAxes = numel(findobj(fig, 'Type', 'axes'));
            nLines = numel(findobj(fig, 'Type', 'line'));
            nPatches = numel(findobj(fig, 'Type', 'patch'));
            testCase.verifyEqual(nAxes, 1);
            testCase.verifyEqual(nLines, 7);
            testCase.verifyEqual(nPatches, 2);
            
            % Verify that the figure has an x-label and y-label
            ax = findobj(fig, 'Type', 'axes');
            xLabel = get(ax, 'xlabel');
            yLabel = get(ax, 'ylabel');
            testCase.verifyNotEmpty(xLabel.String, 'Figure should have an x-label');
            testCase.verifyNotEmpty(yLabel.String, 'Figure should have a y-label');

            % Verify that the patches' height reflects the quartiles of the
            % pre-reward-arrival data
            expectedQuartiles = [4, 12; 2, 6];
            patches = findobj(fig, 'Type', 'patch');
            for iPatch = 1:2
                yData = patches(iPatch).Vertices(:, 2);
                lowerQuartile = min(yData);
                upperQuartile = max(yData);
                testCase.verifyEqual([lowerQuartile, upperQuartile], expectedQuartiles(iPatch, :));
            end
                
            % Verify that the lines' height reflects the min, median and
            % max of the pre-reward-arrival data
            expectedMedians = [8, 4];
            expectedMins = [2, 1];
            expectedMaxes = [14, 7];
            lines = findobj(fig, 'Type', 'line');
            highLineData = [lines(4).YData lines(5).YData];
            mediumLineData = [lines(6).YData lines(7).YData];
            testCase.verifyEqual(highLineData, [expectedMedians(1), expectedMedians(1), expectedMins(1), expectedMaxes(1)]);
            testCase.verifyEqual(mediumLineData, [expectedMedians(2), expectedMedians(2), expectedMins(2), expectedMaxes(2)]);
            
            % Verify that there is an asterisk marking no significance
            texts = findobj(output.figures{1}, 'Type', 'Text');
            significanceText = findobj(texts, 'String', 'n.s.');
            testCase.verifyNotEmpty(significanceText, 'Non-significance marker "n.s." should be present when p >= 0.05.');
            
        end
        
        % Test that plotting of peak fring rates produces a correct figure
        % with data from three rats
        function testPlotPeakPreRewardFiring_multipleRats(testCase)
            
            % Mock inclusion criteria
            testCase.RewardRelatedFiringObj.inclusionCriteria.rats = {'A_rat', 'B_rat', 'C_rat'};
            
            % Mock reward firing data
            testCase.RewardRelatedFiringObj.rewardFiring = struct(...
                'perCellMediumRewardExpectation', struct(...
                'A_rat', ones(10, 200)*1,...
                'B_rat', ones(10, 200)*2,...
                'C_rat', ones(10, 200)*3),...
                'perCellHighRewardExpectation', struct(...
                'A_rat', ones(10, 200)*1,...
                'B_rat', ones(10, 200)*2,...
                'C_rat', ones(10, 200)*3),...
                'ts', linspace(-10, 0, 200));

            % Call plotPeakPreRewardFiring
            output = testCase.RewardRelatedFiringObj.plotPeakPreRewardFiring();

            % Verify that one figure is created
            nFigures = length(output.figures);
            testCase.verifyEqual(nFigures, 1, 'One figure should be created for peak pre-reward firing.');
            
            % Verify that the figure contains seven Lines and two Patches
            fig = output.figures{1};
            nAxes = numel(findobj(fig, 'Type', 'axes'));
            nLines = numel(findobj(fig, 'Type', 'line'));
            nPatches = numel(findobj(fig, 'Type', 'patch'));
            testCase.verifyEqual(nAxes, 1);
            testCase.verifyEqual(nLines, 7);
            testCase.verifyEqual(nPatches, 2);
            
            % Verify that the figure has an x-label and y-label
            ax = findobj(fig, 'Type', 'axes');
            xLabel = get(ax, 'xlabel');
            yLabel = get(ax, 'ylabel');
            testCase.verifyNotEmpty(xLabel.String, 'Figure should have an x-label');
            testCase.verifyNotEmpty(yLabel.String, 'Figure should have a y-label');

            % Verify that the patches' height reflects the quartiles of the
            % pre-reward-arrival data
            expectedQuartiles = [1, 3; 1, 3];
            patches = findobj(fig, 'Type', 'patch');
            for iPatch = 1:2
                yData = patches(iPatch).Vertices(:, 2);
                lowerQuartile = min(yData);
                upperQuartile = max(yData);
                testCase.verifyEqual([lowerQuartile, upperQuartile], expectedQuartiles(iPatch, :));
            end
                
            % Verify that the lines' height reflects the min, median and
            % max of the pre-reward-arrival data
            expectedMedians = [2, 2];
            expectedMins = [1, 1];
            expectedMaxes = [3, 3];
            lines = findobj(fig, 'Type', 'line');
            highLineData = [lines(4).YData lines(5).YData];
            mediumLineData = [lines(6).YData lines(7).YData];
            testCase.verifyEqual(highLineData, [expectedMedians(1), expectedMedians(1), expectedMins(1), expectedMaxes(1)]);
            testCase.verifyEqual(mediumLineData, [expectedMedians(2), expectedMedians(2), expectedMins(2), expectedMaxes(2)]);
            
        end
        
         % Test that plotting of peak fring rates produces correct significance markers
        function testPlotPeakPreRewardFiring_significance(testCase)
            
            % Mock inclusion criteria
            testCase.RewardRelatedFiringObj.inclusionCriteria.rats = {'A_rat'};

            % Mock reward firing data with firing rates for medium higher than for high expectation
            rng(0);
            testCase.RewardRelatedFiringObj.rewardFiring = struct(...
                'perCellMediumRewardExpectation', struct('A_rat', random('Normal', 1, 1, 100, 200)),...
                'perCellHighRewardExpectation', struct('A_rat', random('Normal', 100, 1, 100, 200)),...
                'ts', linspace(-15, 0, 200));

            % Call plotRewardRelatedFiring
            output = testCase.RewardRelatedFiringObj.plotPeakPreRewardFiring();
            
            % Verify that the figure contains seven Lines and two Patches
            fig = output.figures{1};
            nAxes = numel(findobj(fig, 'Type', 'axes'));
            nLines = numel(findobj(fig, 'Type', 'line'));
            nPatches = numel(findobj(fig, 'Type', 'patch'));
            testCase.verifyEqual(nAxes, 1);
            testCase.verifyEqual(nLines, 7);
            testCase.verifyEqual(nPatches, 2);
            
            % Verify that there is an asterisk marking significance
            texts = findobj(output.figures{1}, 'Type', 'Text');
            significanceText = findobj(texts, 'String', '*');
            testCase.verifyNotEmpty(significanceText, 'Significance marker "*" should be present when p < 0.05.');
        
        end
        
        
        % Test that plotting of peak fring rates produces an error
        % with data from no rats
        function testPlotPeakPreRewardFiring_noRats(testCase)
            
            % Mock inclusion criteria
            testCase.RewardRelatedFiringObj.inclusionCriteria.rats = {};
            
            % Mock reward firing data
            rng(0);
            testCase.RewardRelatedFiringObj.rewardFiring = struct(...
                'perCellMediumRewardExpectation', struct('missing_rat', random('Normal', 1, 1, 100, 200)),...
                'perCellHighRewardExpectation', struct('missing_rat', random('Normal', 1, 1, 100, 200)));
            
            % Call plotRewardRelatedFiring
            expectedErrorId = 'RewardRelatedFiring:dataNotAssigned';
            testCase.verifyError(@() testCase.RewardRelatedFiringObj.plotPeakPreRewardFiring(), expectedErrorId);
                        
        end



    end
            
end

