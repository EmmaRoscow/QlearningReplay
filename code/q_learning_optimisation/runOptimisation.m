
function [] = runOptimisation(rat_name)

    % Set parameters
    ksamples = [0 1 3 5 10 15 20 30 40 50 75 100];
    label = rat_name;
    solverName = 'bads';
    lb = [0 0 0];
    ub = [1 1 10];

    % Load behavioural data
    behav_data = importdata(['./data/behavioural_data/behaviour_' lower(rat_name) '.mat']);

    % Specify replay policies
    policies = {'RPEProb', 'RPE', 'expectedReward', 'random'};

    for iPolicy = 1:length(policies)

        params.policy = policies{iPolicy};

        % Set file to write to
        fname = sprintf('%s_LiveOutput_%s.txt', label, params.policy);
        fileID = fopen(fname,'a');
        fprintf(fileID,'%s - Beginning optimisation\n', datetime('now'));
        fclose(fileID);

        for k = 6:length(ksamples)

            rng(1, 'Twister');

            % Set changing parameter value
            params.ksamples = ksamples(k);

            % Ensure the right number of parameters
            if ( strcmp(params.policy, 'random') || strcmp(params.policy, 'expectedReward')  || strcmp(params.policy, 'RPE') || strcmp(params.policy, 'RPEProb') )...
                    && params.ksamples > 0 && length(lb)==3
                lb(4) = 0;
                ub(4) = 10;
            end
            if ( strcmp(params.policy, 'RPE') || strcmp(params.policy, 'RPEProb') ) && params.ksamples > 0 && length(lb)<5
                lb(5) = 0;
                ub(5) = 10;
            end

            switch solverName

                case 'bads'

                    % Specify optimisation problem
                    fun = @(x)run25times(x, behav_data, params);

                    % Run optimisation
                    result_interim = cell(1, 16);
                    training_score_interim = NaN(1, 16);
                    x0 = [];

                    parfor ix = 1:16
                        rng(ix);
                        x0 = (rand(1, length(lb)) + lb) .*(ub-lb);
                        [result_interim{ix}, training_score_interim(ix)] = bads(fun, x0, lb, ub);
                    end

                    % Record best
                    [~, idx_best] = min(training_score_interim);
                    result{k} = result_interim{idx_best};
                    training_score(k) = training_score_interim(idx_best);

            end

            % Test with 1000 random seeds
            test_score(k) = test1000times(result{k}, behav_data, params);

             % Write to file
            fileID = fopen(fname,'a');
            fprintf(fileID,'%s optimisation for ksamples = %i\n',datetime('now'), ksamples(k));
            fprintf(fileID,'resulting parameter values: %f %f %f %f %f\n', result{k});
            fprintf(fileID,'training score: %f\n', training_score(k));
            fprintf(fileID,'test score: %f\n\n', test_score(k));
            fclose(fileID);

            % Save results
            resultsfname = sprintf('%s_Results_%s_%iof%i.mat', label, params.policy, k, length(ksamples));
            save(resultsfname, 'ksamples', 'result', 'training_score', 'test_score');

        end

    end
    
end