
classdef EphysAnalysisTest < matlab.unittest.TestCase

    properties
        EphysObj
        TestDataDir
    end
    

    methods(TestMethodSetup)
        
        % Initialise EphysAnalysis object with mock inclusionCriteria and
        % temporary data path
        function createEphysAnalysis(testCase)
            testCase.EphysObj = EphysAnalysis();
            testCase.EphysObj.inclusionCriteria.cells.rewardResponsive = true;
            % Create a temporary directory for test data
            testCase.TestDataDir = tempname;
            mkdir(testCase.TestDataDir);
        end
        

    end
    
    methods(TestMethodTeardown)
        
        % Remove the temporary test data directory after tests
        function removeTestData(testCase)
            rmdir(testCase.TestDataDir, 's');
        end
        
    end


    methods(Test)
        
        % Test that loadData works
        function testLoadData_success(testCase)
            
            % Define mock data
            [spiketimes, CPEntryTimes, rewardArrivalTimes, CPExitTimes] = TestUtilities.createRewardResponsiveData();
            nNAcUnits = 2;
            rng(0);
            binnedfr = rand(10, length(spiketimes));
            startEndTimes = [0, 10; 20, 30];
            rippleStart = [1.5, 2.5];
            rippleStop = [1.7, 2.7];
            ripplePeak = [1.6, 2.6];
            ripples3std = struct('rippleStart', rippleStart, 'rippleStop', rippleStop, 'ripplePeak', ripplePeak);
            directory = fullfile(testCase.TestDataDir, 'test_session');
            
            TestUtilities.saveMockEphysData(directory,...
                'spiketimes', spiketimes,...
                'CPentryTimes', CPEntryTimes,...
                'RarrivalTimes', rewardArrivalTimes,...
                'CPexitTimes', CPExitTimes,...
                'nNAcUnits', nNAcUnits,...
                'binnedfr', binnedfr,...
                'startEndTimes', startEndTimes,...
                'ripples3std', ripples3std);
            
            % Call loadData
            data = testCase.EphysObj.loadData(directory);
            
            % Verify data struct fields
            testCase.verifyEqual(mean(isfield(data, {'spiketimes', 'nNAcUnits', 'binnedfr', 'startEndTimes', ...
                'rippleStart', 'rippleStop', 'ripplePeak', 'nUnits', 'nCA1Units'})),...
                1, 'Data struct does not contain all expected fields.');
            
            % Verify data content
            testCase.verifyEqual(data.spiketimes, spiketimes);
            testCase.verifyEqual(data.nNAcUnits, nNAcUnits);
            testCase.verifyEqual(data.startEndTimes, startEndTimes);
            testCase.verifyEqual(data.rippleStart, rippleStart);
            testCase.verifyEqual(data.rippleStop, rippleStop);
            testCase.verifyEqual(data.ripplePeak, ripplePeak);
            testCase.verifyEqual(data.binnedfr, binnedfr);
            
            % Verify computed fields
            testCase.verifyEqual(data.nUnits, length(spiketimes), 'nUnits does not match.');
            testCase.verifyEqual(data.nCA1Units, data.nUnits - nNAcUnits, 'nCA1Units does not match.');
            
        end
        
        % Test error handling when filepath does not exist
        function testLoadData_nonExistentFilepath(testCase)
            
            % Call loadData
            expectedErrorId = 'EphysAnalysis:load:couldNotFindPath';
            testCase.verifyError(@() testCase.EphysObj.loadData('non-existent/path'), expectedErrorId);
            
        end
            
        
        % Test error handling when one file is missing
        function testLoadData_missingDataFile(testCase)
            
            % Create mock data, but delete spiketimes file
            directory = fullfile(testCase.TestDataDir, 'test_session');
            TestUtilities.saveMockEphysData(directory);
            delete(fullfile(directory, 'spiketimes.mat'))
            
            % Call loadData
            expectedErrorId = 'EphysAnalysis:load:couldNotReadFile';
            testCase.verifyError(@() testCase.EphysObj.loadData(directory), expectedErrorId);
            
        end
        
        % Test error handling when one variable is missing
        function testLoadData_missingDataVariable(testCase)
            
            % Create mock data
            directory = fullfile(testCase.TestDataDir, 'test_session');
            TestUtilities.saveMockEphysData(directory);
            
            % Replace one data file with a similar one missing a variable
            delete(fullfile(directory, 'ripples3std.mat'))
            ripples3std = struct('rippleStart', [], 'rippleStop', []);
            save(fullfile(directory, 'ripples3std.mat'), '-struct', 'ripples3std', '-v7.3');
            
            % Call loadData
            expectedErrorId = 'EphysAnalysis:load:missingData';
            testCase.verifyError(@() testCase.EphysObj.loadData(directory), expectedErrorId);
            
        end
        
        
        % Test successful loading of ventral striatal spikes
        function testLoadVentralStriatalSpikes_success(testCase)
            
            % Create mock data
            [spiketimes, ~, ~, ~] = TestUtilities.createRewardResponsiveData();
            nNAcUnits = 2;
            directory = fullfile(testCase.TestDataDir, 'test_session');
            TestUtilities.saveMockEphysData(directory,...
                'spiketimes', spiketimes,...
                'nNAcUnits', nNAcUnits)
            
            % Call loadVentralStriatalSpikes
            loadedSpiketimes = testCase.EphysObj.loadVentralStriatalSpikes(directory);
            
            testCase.verifyEqual(loadedSpiketimes, spiketimes(1:nNAcUnits), 'Ventral striatal spiketimes data does not match.')
            
        end
        
        % Test error handling when filepath does not exist
        function testLoadVentralStriatalSpikes_nonExistentFilepath(testCase)
            
            % Call loadVentralStriatalSpikes
            expectedErrorId = 'EphysAnalysis:load:couldNotFindPath';
            testCase.verifyError(@() testCase.EphysObj.loadVentralStriatalSpikes('non-existent/path'), expectedErrorId);
            
        end
            
        
        % Test error handling when spiketimes file is missing
        function testLoadVentralStriatalSpikes_missingSpiketimes(testCase)
            
            % Create mock data, but delete spiketimes file
            directory = fullfile(testCase.TestDataDir, 'test_session');
            TestUtilities.saveMockEphysData(directory);
            delete(fullfile(directory, 'spiketimes.mat'))
            
            % Call loadVentralStriatalSpikes
            expectedErrorId = 'EphysAnalysis:load:couldNotReadFile';
            testCase.verifyError(@() testCase.EphysObj.loadVentralStriatalSpikes(directory), expectedErrorId);
            
        end
        
        % Test error handling when nNAcUnits is invalid 
        function testLoadVentralStriatalSpikesn_NAcUnitsArray(testCase)
            
            % Create mock data, but make nNAcUnits an array
            directory = fullfile(testCase.TestDataDir, 'test_session');
            TestUtilities.saveMockEphysData(directory,...
                'nNAcUnits', [1, 2, 3])
            
            % Call loadData
            expectedErrorId = 'EphysAnalysis:load:invalidData';
            testCase.verifyError(@() testCase.EphysObj.loadVentralStriatalSpikes(directory), expectedErrorId);
            
        end
        
        % Test error handling when nNAcUnits is invalid 
        function testLoadVentralStriatalSpikes_nNAcUnitsTooBig(testCase)
            
            % Create mock data, with nNAcUnits too big
            [spiketimes, ~, ~, ~] = TestUtilities.createRewardResponsiveData();
            nNAcUnits = numel(spiketimes)+1;
            directory = fullfile(testCase.TestDataDir, 'test_session');
            TestUtilities.saveMockEphysData(directory,...
                'spiketimes', spiketimes,...
                'nNAcUnits', nNAcUnits)
            
            % Call loadData
            expectedErrorId = 'EphysAnalysis:load:invalidData';
            testCase.verifyError(@() testCase.EphysObj.loadVentralStriatalSpikes(directory), expectedErrorId);
            
        end
        
        % Test error handling when nNAcUnits is the same length as spiketimes 
        function testLoadVentralStriatalSpikes_onlyNAc(testCase)
            
            % Create mock data, with nNAcUnits too big
            [spiketimes, ~, ~, ~] = TestUtilities.createRewardResponsiveData();
            nNAcUnits = length(spiketimes);
            directory = fullfile(testCase.TestDataDir, 'test_session');
            TestUtilities.saveMockEphysData(directory,...
                'spiketimes', spiketimes,...
                'nNAcUnits', nNAcUnits)
            
            % Call loadVentralStriatalSpikes
            loadedSpiketimes = testCase.EphysObj.loadVentralStriatalSpikes(directory);
            
            testCase.verifyEqual(loadedSpiketimes, spiketimes(1:nNAcUnits), 'Ventral striatal spiketimes data does not match.')
            
        end
        
        % Test error handling when nNAcUnits is 0 
        function testLoadVentralStriatalSpikes_NAc0(testCase)
            
            % Create mock data, with nNAcUnits too big
            [spiketimes, ~, ~, ~] = TestUtilities.createRewardResponsiveData();
            directory = fullfile(testCase.TestDataDir, 'test_session');
            TestUtilities.saveMockEphysData(directory,...
                'spiketimes', spiketimes,...
                'nNAcUnits', 0)
            
            % Call loadVentralStriatalSpikes
            loadedSpiketimes = testCase.EphysObj.loadVentralStriatalSpikes(directory);
            
            testCase.verifyEmpty(loadedSpiketimes, 'Spiketimes should be empty.');
            
        end

                
        % Test successful loading of data files with behavioural timestamps
        function testLoadRewardCPTimes_success(testCase)
            
            % Create mock behavioural timestamp data
            RarrivalTimes = [4.5, 9.6, 15.7];
            CPentryTimes = [1.3, 5.4, 12.0];
            CPexitTimes = [3.3, 8.7, 13.9];
            
            % Save mock data to .mat files
            save(fullfile(testCase.TestDataDir, 'RarrivalTimes.mat'), 'RarrivalTimes');
            save(fullfile(testCase.TestDataDir, 'CPentryTimes.mat'), 'CPentryTimes');
            save(fullfile(testCase.TestDataDir, 'CPexitTimes.mat'), 'CPexitTimes');
            
            % Call loadRewardCPTimes
            [loadedRarrival, loadedCPentry, loadedCPexit] = testCase.EphysObj.loadRewardCPTimes(testCase.TestDataDir);
            
            testCase.verifyEqual(loadedRarrival, RarrivalTimes, 'RarrivalTimes do not match.');
            testCase.verifyEqual(loadedCPentry, CPentryTimes, 'CPentryTimes do not match.');
            testCase.verifyEqual(loadedCPexit, CPexitTimes, 'CPexitTimes do not match.');
            
        end
        
        
        % Test unsuccessful loading of data files where one is missing
        function testLoadRewardCPTimes_missingFile(testCase)
            
            % Create mock behavioural timestamp data
            RarrivalTimes = [4.5, 9.6, 15.7];
            CPentryTimes = [1.3, 5.4, 12.0];
            
            % Save mock data to .mat files
            save(fullfile(testCase.TestDataDir, 'RarrivalTimes.mat'), 'RarrivalTimes');
            save(fullfile(testCase.TestDataDir, 'CPentryTimes.mat'), 'CPentryTimes');
            
             % Call loadRewardCPTimes
            expectedErrorId = 'MATLAB:importdata:FileNotFound';
            testCase.verifyError(@() testCase.EphysObj.loadRewardCPTimes(testCase.TestDataDir), expectedErrorId);
            
        end
        
        
        % Test unsuccessful loading of data files where timestamps are
        % inconsistent across trials (too many central platform exits)
        function testLoadRewardCPTimes_inconsistentDataLength1(testCase)
            
            % Define file paths
            RarrivalTimes = [4.5, 9.6, 15.7];
            CPentryTimes = [1.3, 5.4, 12.0];
            CPexitTimes = [3.3, 8.7, 13.9, 16.0, 17.2];
            
            % Save mock data to .mat files
            save(fullfile(testCase.TestDataDir, 'RarrivalTimes.mat'), 'RarrivalTimes');
            save(fullfile(testCase.TestDataDir, 'CPentryTimes.mat'), 'CPentryTimes');
            save(fullfile(testCase.TestDataDir, 'CPexitTimes.mat'), 'CPexitTimes');
            
            % Call loadRewardCPTimes
            expectedError = 'EphysAnalysis:load:invalidData';
            testCase.verifyError(@() testCase.EphysObj.loadRewardCPTimes(testCase.TestDataDir), expectedError);
            
        end
        
        
        % Test unsuccessful loading of data files where timestamps are
        % inconsistent across trials (not enough central platform exits)
        function testLoadRewardCPTimes_inconsistentDataLength2(testCase)
            
            % Define file paths
            RarrivalTimes = [4.5, 9.6, 15.7];
            CPentryTimes = [1.3, 5.4, 12.0];
            CPexitTimes = [3.3, 8.7];
            
            % Save mock data to .mat files
            save(fullfile(testCase.TestDataDir, 'RarrivalTimes.mat'), 'RarrivalTimes');
            save(fullfile(testCase.TestDataDir, 'CPentryTimes.mat'), 'CPentryTimes');
            save(fullfile(testCase.TestDataDir, 'CPexitTimes.mat'), 'CPexitTimes');
            
            % Call loadRewardCPTimes
            expectedError = 'EphysAnalysis:load:invalidData';
            testCase.verifyError(@() testCase.EphysObj.loadRewardCPTimes(testCase.TestDataDir), expectedError);
            
        end
        
        
        % Test loading when timestamp files are empty
        function testLoadRewardCPTimes_emptyFiles(testCase)
            
            % Define empty arrays
            RarrivalTimes = [];
            CPentryTimes = [];
            CPexitTimes = [];
            
            % Save mock data to .mat files
            save(fullfile(testCase.TestDataDir, 'RarrivalTimes.mat'), 'RarrivalTimes');
            save(fullfile(testCase.TestDataDir, 'CPentryTimes.mat'), 'CPentryTimes');
            save(fullfile(testCase.TestDataDir, 'CPexitTimes.mat'), 'CPexitTimes');
            
            % Call loadRewardCPTimes
            [loadedRarrival, loadedCPentry, loadedCPexit] = testCase.EphysObj.loadRewardCPTimes(testCase.TestDataDir);
            
            % Verify the outputs are empty
            testCase.verifyEmpty(loadedRarrival, 'RarrivalTimes should be empty.');
            testCase.verifyEmpty(loadedCPentry, 'CPentryTimes should be empty.');
            testCase.verifyEmpty(loadedCPexit, 'CPexitTimes should be empty.');
            
        end
        
        
        % Test that all cells are included when inclusionCriteria.rewardResponsive is false
        function testSelectRewardResponsiveCells_allIncluded(testCase)
            testCase.EphysObj.inclusionCriteria.cells.rewardResponsive = false;
            [spiketimes, CPEntryTimes, rewardArrivalTimes, ~] = TestUtilities.createRewardResponsiveData();         
            expected = boolean(ones(1, 3));
            actual = testCase.EphysObj.selectRewardResponsiveCells(spiketimes, rewardArrivalTimes, CPEntryTimes);
            testCase.verifyEqual(actual, expected, 'All cells should be included when rewardResponsive is false');
        end
        
        
        % Test that only reward-responsive cells are included when
        % inclusionCriteria.rewardResponsive is true
        function testSelectRewardResponsiveCells_withCriteria(testCase)
            testCase.EphysObj.inclusionCriteria.cells.rewardResponsive = true;
            [spiketimes, CPEntryTimes, rewardArrivalTimes, ~] = TestUtilities.createRewardResponsiveData();
            expected = [true, false, true];
            actual = testCase.EphysObj.selectRewardResponsiveCells(spiketimes, rewardArrivalTimes, CPEntryTimes);
            testCase.verifyEqual(actual, expected, 'Only first and third cells should be responsive');
        end
        
        
        % Test that empty spiketimes returns no reward-responsive cells
        function testSelectRewardResponsiveCells_emptySpiketimes(testCase)
            testCase.EphysObj.inclusionCriteria.cells.rewardResponsive = true;
            [~, CPEntryTimes, rewardArrivalTimes, ~] = TestUtilities.createRewardResponsiveData();
            spiketimes = {};
            expected = false(1, 0);
            actual = testCase.EphysObj.selectRewardResponsiveCells(spiketimes, rewardArrivalTimes, CPEntryTimes);
            
            testCase.verifyEqual(actual, expected, 'Empty spiketimes should return an empty array');
        end

        
        % Test countRewardResponsiveCells with mock data for multiple rats and sessions
        function testCountRewardResponsiveCells_basicFunctionality(testCase)
                        
            % Mock inclusion criteria and sessions
            testCase.EphysObj.inclusionCriteria.rats = {'A_rat', 'B_rat'};
            testCase.EphysObj.sessionsForBroadAnalysis = struct('A_rat', [1, 2], 'B_rat', [2]);
            
            % Define RatA Session1 with 3 cells, of which 2 are NAc cells,
            % of which 1 is reward-responsive
            [spiketimes_ratA_sess1, CPEntryTimes_ratA_sess1, rewardArrivalTimes_ratA_sess1, CPExitTimes_ratA_sess1] = TestUtilities.createRewardResponsiveData();
            nNAcUnits_ratA_sess1 = 2;
            ratA_sess1_dir = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
            TestUtilities.saveMockEphysData(ratA_sess1_dir,...
                'spiketimes', spiketimes_ratA_sess1,...
                'CPentryTimes', CPEntryTimes_ratA_sess1,...
                'RarrivalTimes', rewardArrivalTimes_ratA_sess1,...
                'CPexitTimes', CPExitTimes_ratA_sess1,...
                'nNAcUnits', nNAcUnits_ratA_sess1)
            
            % Define RatA Session2 with 6 cells, of which 4 are NAc cells,
            % of which 3 are reward-responsive
            spiketimes_ratA_sess2 = [spiketimes_ratA_sess1(:); spiketimes_ratA_sess1(:)];
            CPEntryTimes_ratA_sess2 = [CPEntryTimes_ratA_sess1(:); CPEntryTimes_ratA_sess1(:)];
            rewardArrivalTimes_ratA_sess2 = [rewardArrivalTimes_ratA_sess1(:); rewardArrivalTimes_ratA_sess1(:)];
            CPExitTimes_ratA_sess2 = [CPExitTimes_ratA_sess1(:); CPExitTimes_ratA_sess1(:)];
            nNAcUnits_ratA_sess2 = 4;
            ratA_sess2_dir = fullfile(testCase.TestDataDir, 'Rat_A', 'Session2');
            TestUtilities.saveMockEphysData(ratA_sess2_dir,...
                'spiketimes', spiketimes_ratA_sess2,...
                'CPentryTimes', CPEntryTimes_ratA_sess2,...
                'RarrivalTimes', rewardArrivalTimes_ratA_sess2,...
                'CPexitTimes', CPExitTimes_ratA_sess2,...
                'nNAcUnits', nNAcUnits_ratA_sess2)
            
            % Skip RatB Session 1
            
            % Define RatB Session 2 with 2 cells, of which 1 is a NAc cell,
            % of which 1 is reward-responsive
            spiketimes_ratB_sess2 = cellfun(@(x) x + 1, spiketimes_ratA_sess2(1:2), 'UniformOutput', false);
            CPEntryTimes_ratB_sess2 = CPEntryTimes_ratA_sess2 + 1;
            rewardArrivalTimes_ratB_sess2 = rewardArrivalTimes_ratA_sess2 + 1;
            CPExitTimes_ratB_sess2 = CPExitTimes_ratA_sess2 + 1;
            nNAcUnits_ratB_sess2 = 1;
            ratB_sess2_dir = fullfile(testCase.TestDataDir, 'Rat_B', 'Session2');
            TestUtilities.saveMockEphysData(ratB_sess2_dir,...
                'spiketimes', spiketimes_ratB_sess2,...
                'CPentryTimes', CPEntryTimes_ratB_sess2,...
                'RarrivalTimes', rewardArrivalTimes_ratB_sess2,...
                'CPexitTimes', CPExitTimes_ratB_sess2,...
                'nNAcUnits', nNAcUnits_ratB_sess2)
            
            % Call the method
            perRatRate = testCase.EphysObj.countRewardResponsiveCells('DataPath', testCase.TestDataDir);
            
            % Define expected perRatRewardResponsiveRate
            expectedPerRatRate = struct('A_rat', 2/3, 'B_rat', 1.0);
            
            % Verify perRatRewardResponsiveRate
            testCase.verifyEqual(perRatRate, expectedPerRatRate, 'Per-rat responsive rates do not match expected values.');
            
        end
        
        % Test countRewardResponsiveCells with no rats in inclusion criteria
        function testCountRewardResponsiveCells_noRats(testCase)
            
            % Mock inclusion criteria with empty rats
            testCase.EphysObj.inclusionCriteria.rats = {};
            testCase.EphysObj.sessionsForBroadAnalysis = struct();
            
            % Call the method
            perRatRate = testCase.EphysObj.countRewardResponsiveCells('DataPath', testCase.TestDataDir);
            
            % Verify the output is an empty struct
            testCase.verifyEmpty(fieldnames(perRatRate), 'Per-rat responsive rate struct should be empty when no rats are specified.');
            
        end
        
        
        % Test countRewardResponsiveCells when a rat has no sessions
        function testCountRewardResponsiveCells_ratWithNoSessions(testCase)
            
            % Mock inclusion criteria and sessions
            testCase.EphysObj.inclusionCriteria.rats = {'A_rat', 'B_rat'};
            testCase.EphysObj.sessionsForBroadAnalysis = struct('A_rat', [], 'B_rat', [1]);
            
            % Skip RatA
            
            % Define RatB Session1 with 3 cells, of which 2 are NAc cells,
            % of which 1 is reward-responsive
            [spiketimes, CPEntryTimes, rewardArrivalTimes, CPExitTimes] = TestUtilities.createRewardResponsiveData();
            nNAcUnits = 2;
            ratB_sess1_dir = fullfile(testCase.TestDataDir, 'Rat_B', 'Session1');
            TestUtilities.saveMockEphysData(ratB_sess1_dir,...
                'spiketimes', spiketimes,...
                'CPentryTimes', CPEntryTimes,...
                'RarrivalTimes', rewardArrivalTimes,...
                'CPexitTimes', CPExitTimes,...
                'nNAcUnits', nNAcUnits)
                                    
            % Call the method
            perRatRate = testCase.EphysObj.countRewardResponsiveCells('DataPath', testCase.TestDataDir);
            
            % Define expected perRatRewardResponsiveRate
            expectedPerRatRate = struct('A_rat', NaN, 'B_rat', 0.5);
            
            % Verify perRatRewardResponsiveRate
            testCase.verifyEqual(perRatRate, expectedPerRatRate, 'Per-rat responsive rates do not match expected values.');
            
        end
        
        % Test countRewardResponsiveCells when session data files are missing
        function testCountRewardResponsiveCells_missingSessionData(testCase)
            
            % Mock inclusion criteria and sessions
            testCase.EphysObj.inclusionCriteria.rats = {'A_rat'};
            testCase.EphysObj.sessionsForBroadAnalysis = struct('A_rat', [1]);
                        
            % Call the method
            expectedErrorID = 'EphysAnalysis:load:couldNotFindPath';
            testCase.verifyError(@() testCase.EphysObj.countRewardResponsiveCells('DataPath', testCase.TestDataDir), expectedErrorID);
            
        end
        
        % Test countRewardResponsiveCells with no reward-responsive cells
        function testCountRewardResponsiveCells_noCellsResponsive(testCase)
            
            % Mock inclusion criteria and sessions
            testCase.EphysObj.inclusionCriteria.rats = {'A_rat'};
            testCase.EphysObj.sessionsForBroadAnalysis = struct('A_rat', [1]);
                        
            % Define RatA Session1 with 1 cell, of which 1 is a NAc cell,
            % of which 0 are reward-responsive
            [spiketimes, CPEntryTimes, rewardArrivalTimes, CPExitTimes] = TestUtilities.createRewardResponsiveData();
            spiketimes = {spiketimes{2}};
            nNAcUnits = 1;
            ratA_sess1_dir = fullfile(testCase.TestDataDir, 'Rat_A', 'Session1');
            TestUtilities.saveMockEphysData(ratA_sess1_dir,...
                'spiketimes', spiketimes,...
                'CPentryTimes', CPEntryTimes,...
                'RarrivalTimes', rewardArrivalTimes,...
                'CPexitTimes', CPExitTimes,...
                'nNAcUnits', nNAcUnits)
                                    
            % Call the method
            perRatRate = testCase.EphysObj.countRewardResponsiveCells('DataPath', testCase.TestDataDir);
            
            % Define expected perRatRewardResponsiveRate
            expectedPerRatRate = struct('A_rat', 0);
            
            % Verify perRatRewardResponsiveRate
            testCase.verifyEqual(perRatRate, expectedPerRatRate, 'Per-rat responsive rates do not match expected values.');
            
        end
        
        % Test convolveAndBin with constant firing rate of 100 Hz (1 bin)
        function testConvolveAndBin_constantFiringRateSingleBin(testCase)
            
            % Create a single spike train with a constant firing rate of
            % 100 Hz, firing for 11 seconds
            spiketimes = {0:0.01:11};
            
            % Convolve to get a single 1-second bin
            popV = testCase.EphysObj.convolveAndBin(5, 1, 1, spiketimes);
            
            testCase.verifyEqual(popV, 100, 'RelTol', 1e-10);
            
        end
        
        % Test convolveAndBin with a constant firing rate of 50 Hz
        % (multiple bins)
        function testConvolveAndBin_constantFiringRateMultipleBins(testCase)
            
            % Create a single spike train with a constant firing rate of
            % 100 Hz, firing for 11 seconds
            spiketimes = {0:0.02:11};
            
            % Convolve to get 10 bins of 100ms duration
            popV = testCase.EphysObj.convolveAndBin(5, 1, 0.1, spiketimes);
            
            expectedOutput = ones(1, 1, 10) * 50;
            testCase.verifyEqual(popV, expectedOutput, 'RelTol', 1e-10);
            
        end
        
        % Test convolveAndBin with a constant firing rate of 50 Hz
        % (irregular bins)
        function testConvolveAndBin_constantFiringRateIrregularBins(testCase)
            
            % Create a single spike train with a constant firing rate of
            % 100 Hz, firing for 11 seconds
            spiketimes = {0:0.02:11};
            
            % Convolve to get 1 second in 150ms bins (6 full bins and 1
            % smaller bin)
            popV = testCase.EphysObj.convolveAndBin(5, 1, 0.15, spiketimes);
            
            expectedOutput = ones(1, 1, 7) * 50;
            testCase.verifyEqual(popV, expectedOutput, 'RelTol', 1e-10);
            
        end
        
        % Test convolveAndBin with changing, symmetric firing rate
        function testConvolveAndBin_symmetricFiringRate(testCase)
            
            % Create a single spike train that is symmetrical around 0
            rng(0);
            interSpikeInterval = abs(randn(1, 10));
            spiketimes = {[cumsum(interSpikeInterval) - sum(interSpikeInterval) cumsum(flip(interSpikeInterval))]};
            
            % Convolve to get 10 seconds in 100 100ms bins
            popV = testCase.EphysObj.convolveAndBin(-5, 10, 0.1, spiketimes);
            
            % Verify that first 50 bins mirror last 50 bins
            first50 = popV(1, 1, 1:50);
            last50 = popV(1, 1, 51:100);
            testCase.verifyEqual(first50, flip(last50), 'RelTol', 1e-10);
            
        end
        
        % Test convolveAndBin with multiple trials
        function testConvolveAndBin_multipleTrials(testCase)
            
            % Create a single spike train with a different constant firing
            % rate on each 15-second trial (10 Hz, 50 Hz, 100 Hz)
            trial1 = 0:0.1:15;
            trial2 = 15.02:0.02:30;
            trial3 = 30.01:0.01:45;
            spiketimes = {[trial1, trial2, trial3]};
            
            % Convolve to get 1 second of activity in 4 250 ms bins
            startTimes = [5; 20; 35];
            popV = testCase.EphysObj.convolveAndBin(startTimes, 1, 0.25, spiketimes);
            
            expectedOutput = repmat([10; 50; 100], [1, 1, 4]);
            testCase.verifyEqual(popV, expectedOutput, 'RelTol', 1e-10);

        end
        
        % Test convolveAndBin with multiple spike trains
        function testConvolveAndBin_multipleCells(testCase)
            
            % Create a two spike trains, with a different constant firing
            % rate on each of 2 15-second trials
            cell1_trial1 = 0:0.1:15;           % 10 Hz
            cell1_trial2 = 15.02:0.02:30;   % 50 Hz
            cell2_trial1 = 0:0.05:15;         % 20 Hz
            cell2_trial2 = 15.01:0.01:30;   % 100 Hz
            spiketimes = {[cell1_trial1, cell1_trial2], [cell2_trial1, cell2_trial2]};
            
            % Convolve to get 1 second of activity in 4 250 ms bins
            startTimes = [5; 20];
            popV = testCase.EphysObj.convolveAndBin(startTimes, 1, 0.25, spiketimes);
            
            expectedOutput = repmat([10, 20; 50, 100], [1, 1, 4]);
            testCase.verifyEqual(popV, expectedOutput, 'RelTol', 1e-10);

        end


        % Test getEventTriggeredFiring basic functionality
        function testGetEventTriggeredFiring_basicFunctionality(testCase)
            
            eventTimes = [7,8];
            spiketimes = {[0:0.01:15]};
            timeBefore = 1;
            timeAfter = 1;
            
            % Call the method
            testCase.EphysObj.analysisConfig.gaussianWindowWidth = 0.1;
            firingRate = testCase.EphysObj.getEventTriggeredFiring(eventTimes, spiketimes, timeBefore, timeAfter);
            
            expectedOutput = ones(2, 1, 40) * 100;
            testCase.verifyEqual(firingRate, expectedOutput, 'RelTol', 1e-10);
            
        end
        
        % Test getEventTriggeredRunningSpeed basic functionality
        function testGetEventTriggeredRunningSpeed_basicFunctionality(testCase)
            
            eventTimes = [1, 3, 5];
            timestamps = 0:0.02:10;
            speedData = ones(1, numel(timestamps)) * 10;
            timeBefore = 1;
            timeAfter = 1;
            
            % Call the method
            [speed, ts] = testCase.EphysObj.getEventTriggeredRunningSpeed(eventTimes, speedData, timestamps, timeBefore, timeAfter);
            
            expectedOutput_speed = ones(3, 100) * 10;
            expectedOutput_ts = -1:0.02:(1-0.02);
            testCase.verifyEqual(speed, expectedOutput_speed, 'RelTol', 1e-10);
            testCase.verifyEqual(ts, expectedOutput_ts, 'RelTol', 1e-10);
            
        end
        
        % Test getEventTriggeredCoactivity basic functionality
        function testGetEventTriggeredCoactivity_basicFunctionality(testCase)
            
            spiketimes1 = {0:0.01:15};
            spiketimes2 = spiketimes1;
            eventTimes = [5; 10];
            timeBefore = 1;
            timeAfter = 1;
            gaussianWindowWidth = 0.1;
            
            % Call the method
            coactivity = testCase.EphysObj.getEventTriggeredCoactivity(eventTimes, spiketimes1, spiketimes2, timeBefore, timeAfter, gaussianWindowWidth);
            
            expectedOutput = ones(2, 1, 40) * 100;
            testCase.verifyEqual(coactivity, expectedOutput, 'RelTol', 1e-10);
            
        end
        
        % Test getEventTriggeredCoactivity with one silent spike train
        function testGetEventTriggeredCoactivity_silentSpikeTrain(testCase)
            
            spiketimes1 = {[0:0.01:15]};
            spiketimes2 = {[]};
            eventTimes = [5; 10];
            timeBefore = 1;
            timeAfter = 1;
            gaussianWindowWidth = 0.1;
            
            % Call the method
            coactivity = testCase.EphysObj.getEventTriggeredCoactivity(eventTimes, spiketimes1, spiketimes2, timeBefore, timeAfter, gaussianWindowWidth);
            
            expectedOutput = ones(2, 1, 40) * 0;
            testCase.verifyEqual(coactivity, expectedOutput, 'RelTol', 1e-10);
            
        end
        
         % Test getEventTriggeredCoactivity with spike trains missing
        function testGetEventTriggeredCoactivity_missingSpikeTrains(testCase)
            
            spiketimes1 = {};
            spiketimes2 = {};
            eventTimes = [5; 10];
            timeBefore = 1;
            timeAfter = 1;
            gaussianWindowWidth = 0.1;
            
            % Call the method
            coactivity = testCase.EphysObj.getEventTriggeredCoactivity(eventTimes, spiketimes1, spiketimes2, timeBefore, timeAfter, gaussianWindowWidth);
            
            expectedOutput = double.empty(2, 0, 40);
            testCase.verifyEqual(coactivity, expectedOutput, 'RelTol', 1e-10);
            
        end

        % Test zScore for one cell
        function testZScore_singleCell(testCase)
            
            firingRate = [3, 3, 5, 7, 7];
            reference = firingRate;
            mu = 5;
            sigma = 2;

            % Call the method
            zScoredFiringRate = testCase.EphysObj.zScore(firingRate, reference);
            
            expectedOutput = (firingRate - mu) ./ sigma;
            testCase.verifyEqual(zScoredFiringRate, expectedOutput, 'RelTol', 1e-10);
            
        end
        
        % Test zScore for two cells and a separate reference
        function testZScore_twoCellsReference(testCase)
            
            firingRate = [1, 1, 1, 1, 1; 2, 2, 2, 2, 2];
            reference = [3, 3, 5, 7, 7];
            mu = 5;
            sigma = 2;

            % Call the method
            zScoredFiringRate = testCase.EphysObj.zScore(firingRate, reference);
            
            expectedOutput = (firingRate - mu) ./ sigma;
            testCase.verifyEqual(zScoredFiringRate, expectedOutput, 'RelTol', 1e-10);
            
        end
        
        % Test zScore for three dimensions
        function testZScore_threeDimensions(testCase)
            
            rng(0);
            firingRate = randn(3, 2, 5);
            reference = [3, 3, 5, 7, 7; 7, 7, 10, 13, 13];
            mu = [5, 10];
            sigma = [2, 3];

            % Call the method
            zScoredFiringRate = testCase.EphysObj.zScore(firingRate, reference);
            
            expectedOutput = [...
                (firingRate(:, 1, :) - mu(1)) ./ sigma(1),...
                (firingRate(:, 2, :) - mu(2)) ./ sigma(2)];
            testCase.verifyEqual(zScoredFiringRate, expectedOutput, 'RelTol', 1e-10);
            
        end


    end
    
end