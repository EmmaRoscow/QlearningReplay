
rat_names = {'Harry', 'Ignotus', 'Justin', 'Kingsley', 'Lucius', 'Marvolo'};
policies = {'random', 'expectedReward', 'RPE', 'RPEProb'};

for iPol = 1:4
    for iRat = 1:length(rat_names)
    
        rat = rat_names{iRat};
        policy = policies{iPol};

        load(['C:\Users\emmar\OneDrive\Documents\MATLAB\outputs\shuff_Results\' rat '_Results_' policy '_5of12.mat'])
        ksamples5 = ksamples;
        result5 = result;
        test_score5 = test_score;
        training_score5 = training_score;

        load(['C:\Users\emmar\OneDrive\Documents\MATLAB\outputs\shuff_Results\' rat '_Results_' policy '_12of12.mat'])
        result(1:5) = result5;
        test_score(1:5) = test_score5;
        training_score(1:5) = training_score5;

        assert(sum(cellfun(@isempty, result))==0);
        assert(sum(test_score==0)==0);
        assert(sum(training_score==0)==0);
        
        save(['C:\Users\emmar\OneDrive\Documents\MATLAB\rpe_replay\data\optimised_parameters\shuffled\' rat '_' policy '.mat'], 'ksamples', 'result', 'test_score', 'training_score')

        clearvars -except rat_names policies iRat iPol
        
    end
end