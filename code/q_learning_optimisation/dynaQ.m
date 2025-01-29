

function [predictive_accuracy, output] = dynaQ(behav_data, params, display)

%% =========== DESCRIPTION ============ %%

% Inputs
%     behav_data:                 struct with the following fields
%       n_trials                      vector of number of trials taken in each learning session
%       actions                      cell array of actions taken on each trial, one cell per session
%    
%     params:                       struct with the following fields
%       params.seed               random seed for reproducibility (optional)
%       params.alpha              learning rate
%       params.discount          discount factor
%       params.epsilson          factor for exploration vs exploitation in action selection
%       params.policy             replay policy: 'none', 'random', 'RPE', 'RPEProb', 'expectedReward'
%       params.ksamples        number of samples to replay between each session (integer)
%       params.recency          factor for weighting more recent trials more highly when selecting samples to replay
%       params.rpe_recency    factor for weighting more recent trials more highly when averaging RPE for state
%    
%     display:                        boolean for whether to create plots or not

% Outputs
%     predictive_accuracy:     a measure of accuracy in predicting the next action, for each trial
%
%     output:                         struct with the following fields
%       Q_history:                  the learned Q-values for each state-action pair, updated after every trial
%       s_history:                  the state for every trial
%       a_history:                  the (true) action taken for every trial
%       rpe_history:               the reward-prediction error resulting from every trial
%       replay_history_abs:    the trial from which every replay event was taken, e.g. replay_history_abs = 12 means the 12th trial was replayed, correponding to s_history(12) and a_history(12)
%       replay_history_ago:    the number of trials elapsed since the trial from which every replay event was taken, e.g. replay_history_ago = 12 means the replayed trial occurred 12 trials ago
%       action_pred:               predicted action for every trial according to action_prob
%       action_prob:               predicted action probabilities for every trial
%       obs_bar:                    observed action probabilities over all learning

%% Initialise variables

% Set random seed
if isfield(params, 'seed')
    rng(params.seed);
else
    rng(0);
end

Q_init = 'biased';
action_selection = 'epsilon';

% Arrays of replay info: which samples replayed and RPE
n_sessions = length(behav_data.n_trials);
n_trials = sum(behav_data.n_trials-1);
trialsReplayed = [];
replay_history_ago = cell(n_sessions, 1);
replay_history_abs = cell(n_sessions, 1);

% Arrays of action predictions
action_prob = NaN(n_trials, 3);
action_pred = NaN(n_trials, 1);

% Q-values: 'biased' promotes bias towards switching
switch Q_init
    case 'zeros'
        Q = zeros(3,3);
    case 'biased'
        Q = [0 0.7 0.7; 0.7 0 0.7; 0.7 0.7 0];
end
Q_history = zeros(n_trials,3,3);

% Set up for RPS reliability score (see Murphy, 1973, "A new vector partition of the probability score", J. Appl. Meteor.)
if ~isfield(params, 'learning_stage')
    [obs_bar, N_state] = get_obs(behav_data, 'all');
else
    [obs_bar, N_state] = get_obs(behav_data, params.learning_stage);
end

% State-action pairs
sa_pairs = reshape([1:9],3,3);

% Memory buffer of past trials from which to replay (separate arrays for each state-action pair)
samples = struct();
for i = 1:9
    samples(i).samples = zeros(10000,6);
    samples(i).dn = 1;
end

% List of all past trials (single array for all state-action pairs)
sa_list = NaN(10000,1);

% List of all rewarded past trials (single array)
r_list = NaN(10000,1);

% Number of trials experienced
dn_tot = 1;

% Reward-prediction error on every trial
rpe_history = NaN(n_trials,1);
    
 %% Run prediction
 
 trial_counter = 0;
 
 for iSession = 1:n_sessions
     
     for iTrial = 2:length(behav_data.actions{iSession})
         
         trial_counter = trial_counter + 1;
         
        % =================================== %
        % GET STATE, ACTION, STATE-ACTION      %
        % INDEX, AND REWARD FROM                 %
        % BEHAVIOURAL DATA                             %
        
        % state, based on action taken on previous trial
        s = behav_data.actions{iSession}(iTrial-1);

        % action
        a = behav_data.actions{iSession}(iTrial);
        
        % state-action pair
        sa = sa_pairs(s,a);
        
        % resulting state (state on next trial)
        s_prime = a;
               
        % reward
        r = behav_data.rewarded{iSession}(iTrial);
        
        
        % =================================== %
        % PREDICT ACTION
        
        % Remove zeros from Q-values
        Qvals = Q(sa_pairs(s,:)) + exp(-10);
        
        % Action selection
        switch action_selection

            case 'softmax'
                action_prob(trial_counter, :) = exp(Qvals) ./ sum(exp(Qvals));
                [~, action_pred(trial_counter)] = max(action_prob(trial_counter, :));

            case 'epsilon'
                action_prob(trial_counter, :) = exp(params.epsilon*Qvals) / sum(exp(params.epsilon*Qvals));
                [~, action_pred(trial_counter)] = find(rand < cumsum(action_prob(trial_counter, :)), 1);

        end
        
        
        % =================================== %
        % EVALUATE PREDICTION
                
        % RPS reliability score
        RPS_rel(trial_counter) = N_state(s) * sum((action_prob(trial_counter, :) - obs_bar(s,:)).^2);
           
        
        % =================================== %
        % UPDATE Q-VALUES                               %
        
        [Q, rpe] = updateQValues(Q, sa, r, s_prime, sa_pairs, params);
        
        % Store history of Q-values, RPEs, states and actions
        rpe_history(trial_counter) = rpe;
        Q_history(trial_counter,:,:) = Q;
        a_history(trial_counter) = a;
        s_history(trial_counter) = s;
        
        
        % =================================== %
        % STORE SAMPLES TO MEMORY BUFFER %
        
        % Add current trial to bank of samples
        samples(sa).samples(samples(sa).dn,:) = [s, a, r, s_prime, rpe, trial_counter];

        % Record how many times this state-action pair has been experienced
        samples(sa).dn = samples(sa).dn + 1;

        % Add current state-action pair to list
        sa_list(dn_tot) = sa;

        % Add current reward to list
        r_list(dn_tot) = r;

        % List of all rewarded state-action pairs experienced
        r_sa = unique(sa_list(find(r_list)));
        
     end
        
        
         % ===================================== %
         % REPLAY                                                  %
         if ~strcmp(params.policy, 'none')
             
             % If ksamples is not an integer, use it as a probability
             ksamples = floor(params.ksamples);
             ksamples_frac = params.ksamples - ksamples;
             if rand < ksamples_frac
                 ksamples = ksamples + 1;
             end
                 
             % Replay the right number of times accordingly
             for k = 1:ksamples

                % Select sample to replay according to replay policy
                [sample, rep, sampleInd] = chooseSampleToReplay(samples, Q, params);

                % Replay (if applicable)
                if rep

                    % Store
                    replay_history_abs{iSession} = [replay_history_abs{iSession} sample(6)];
                    replay_history_ago{iSession} = [replay_history_ago{iSession} trial_counter - sample(6)];

                    % Update Q
                    [Q, rpe] = updateQValues(Q, sa_pairs(sample(1),sample(2)), sample(3), sample(4), sa_pairs, params);

                    % If replay policy based on RPE, replace the sample's RPE with the new one
                    if strcmp(params.policy,'RPE') || strcmp(params.policy,'RPEProb')
                        samples(sampleInd(1)).samples(sampleInd(2), 5) = rpe;
                    end

                    % Store replay info
                    trialsReplayed = [trialsReplayed; iSession samples(sampleInd(1)).samples(sampleInd(2),6) sampleInd(1) rpe];

                end
                
             end

         end
         
 end
     
    % Output results as a structure
    output.Q_history = Q_history;
    output.s_history = s_history;
    output.a_history = a_history;
    output.rpe_history = rpe_history;
    output.replay_history_abs = replay_history_abs;
    output.replay_history_ago = replay_history_ago;
    output.action_pred = action_pred;
    output.action_prob = action_prob;
    output.obs_bar = obs_bar;
    predictive_accuracy = RPS_rel;
     
     if display
        plotQhistory(Q_history, behav_data.rewarded);
        plotReplayHistory(replay_history_abs, behav_data.n_trials)
        plotPredictiveAccuracy(predictive_accuracy, s_history)
     end
      
end

function [Q, rpe] = updateQValues(Q, sa, r, s_prime, sa_pairs, params)
        
        % Q-values for possible next state-action pairs
        available_sa_prime = sa_pairs(s_prime,:);
        Qvals = Q(available_sa_prime);
        
        % on-policy most likely next action
        [maxval, ~] = max(Qvals);
        [bestA] = find(Qvals==maxval);
        if length(bestA) > 1
            a_prime_star = bestA(1);
        else
            a_prime_star = bestA;
        end
        sa_prime_star = sa_pairs(s_prime,a_prime_star);
        
        % reward prediction error for reward on this trial
        rpe = r + params.discount*Q(sa_prime_star) - Q(sa);
        
        % update Q-value for state-action pair
        Q(sa) = (1-params.alpha)*Q(sa) + params.alpha*(r + params.discount*Q(sa_prime_star));

        % cap Q-values at 1
        if max(Q(:)) > 1
            Q = Q/max(Q(:)); Q(isnan(Q)) = 0;
        end
        
end


function [sample, rep, sampleInd] = chooseSampleToReplay(samples, Q, params)

    % Return 0 if there is nothing to replay
    rep = 1;

    % List state-action pairs available in the buffer
    u_sa = find([samples.dn] > 1);

    % Select sample to replay, according to replay policy
    switch params.policy

        case 'none'
            sample = [];
            rep = 0;
            sampleInd = [];
            w = [];
            return                

        case 'random'
            % Randomly choose a state-action transition
            sa = u_sa(randi(length(u_sa)));
            sample_list = samples(sa).samples(find(samples(sa).samples(:,1)),:);
%             sample_list = flip(sample_list,1);
            
        case 'RPE'
                        
            % Get a measure of RPE for each state-action pair
            rpe_weighting = NaN(1, 9);
            for i = 1:length(u_sa)
                % Samples from this SA pair
                sample_list = samples(u_sa(i)).samples(find(samples(u_sa(i)).samples(:,1)),:);
                % Take recency-weighted mean of RPE
                n_samples = size(sample_list, 1);
                rpe_weighting(u_sa(i)) = mean(abs(sample_list(:,5)).*params.rpe_recency.^(1./(1:n_samples)'));
            end
            
            % Choose transition with greatest recent RPE to replay
            [~, sa] = max(abs(rpe_weighting));
            sample_list = samples(sa).samples(find(samples(sa).samples(:,1)),:);
            
        case 'RPEProb'
            
            % Get a measure of RPE for each state-action pair
            rpe_weighting = NaN(1,9);
            for i = 1:length(u_sa)
                % Samples from this SA pair
                sample_list = samples(u_sa(i)).samples(find(samples(u_sa(i)).samples(:,1)),:);
                % Take recency-weighted mean of RPE
                n_samples = size(sample_list, 1);
                rpe_weighting(u_sa(i)) = mean(abs(sample_list(:,5)).*params.rpe_recency.^(1./(1:n_samples)'));
            end
            
            % Choose transition to replay in proportion to RPE
            % Remove any RPE values which are NaN ot inf
            rpe_weighting(isnan(rpe_weighting)) = 0;
            rpe_weighting(find(rpe_weighting==Inf)) = max(rpe_weighting(rpe_weighting~=Inf))*100;
            sa = find(rand < cumsum(rpe_weighting/sum(rpe_weighting)), 1);
            if isempty(sa)
                sa = u_sa(randi(length(u_sa)));
            end
            sample_list = samples(sa).samples(find(samples(sa).samples(:,1)),:);
            
        case 'expectedReward'
            
            % If Q-values are all 0, choose randomly
            if sum(Q(:)) == 0
                sa = randi(length(Q));
                
            % If only one state-action pair in the buffer, select it
            elseif length(u_sa)==1
                sa = u_sa;
                
            % Otherwise probability of selection state-action pair in
            % proportion to Q-value
            else
                prob = Q(u_sa)./sum(Q(u_sa));
                sa = u_sa(find(rand < cumsum(prob),1));
            end
            sample_list = samples(sa).samples(find(samples(sa).samples(:,1)),:);
            
        otherwise
            error('Error: replay policy does not exist')
            
    end
    
    % Choose sample of experience from this state-action transition, biased by recency
    l = size(sample_list,1);
    pdf = sample_list(:, 6) .^ params.recency / sum(sample_list(:, 6) .^ params.recency);
    
    % Probabilitistically choose a low value for n, to create recency replay bias
    i = find(rand < cumsum(pdf),1);
    sample = sample_list(i,:);
    sampleInd = [sa, i];
    
    if length(sample) < 6
        warning('length sample too small')
    end

end

function [] = plotQhistory(Q_history, r)

    % Plot Q-values on first trial, first rewarded trial, last trial, and
    % intermediate trial(s) at roughly equal intervals
    n_subplots = 5;
    idx(1) = 1;
    idx(2:n_subplots) = round(linspace(find(cell2mat(r) > 0, 1), size(Q_history, 1), n_subplots-1));
    idx(3:n_subplots-1) = round(idx(3:n_subplots-1), -1);

    % Plot Q-values on sample trials
    figure;
    for i = 1:n_subplots
        subplot(1, n_subplots, i);
        imagesc(squeeze(Q_history(idx(i), :, :)));
        set(gca, 'CLim', [0, ceil(max(Q_history(:)))])
        title(['Trial no. ' num2str(idx(i))]);
        if i == 1 set(gca, 'YTick', [1:3]); ylabel('State')
        else set(gca, 'YTick', []); end
    end
    subplot(1, n_subplots, 1)
    title({'First trial', ['Trial no. ' num2str(idx(1))]});
    subplot(1, n_subplots, 2)
    title({'First rewarded trial', ['Trial no. ' num2str(idx(2))]});
    subplot(1, n_subplots, n_subplots)
    title({'Last trial', ['Trial no. ' num2str(idx(n_subplots))]});
    cb = colorbar; cb.Orientation = 'Horizontal'; cb.Location = 'southoutside'; cb.Label.String = 'Q-value'; cb.Ticks = [0 0.5 1];
    ch = get(gcf, 'Children');
    x1 = ch(end).Position(1);
    position = get(cb, 'Position');
    set(cb, 'Position', [x1 position(2) 0.775 position(4)])
    subplot(1, n_subplots, ceil(n_subplots/2)); xlabel('Action'); drawnow
    for i = 2:n_subplots+1
        ch(i).Position = [ch(i).Position(1) 0.4 ch(i).Position(3) ch(i).Position(4)-0.29];
    end
    drawnow
    
    % Plot evolution of Q-values for all state-action pairs
    figure;
    for i = 1:3
        subplot(3, 1, i); hold on
        for j = 1:3
            plot(Q_history(:, i, j), 'LineWidth', 2)
        end
        ylim([-0.05 ceil(max(Q_history(:))) + 0.05])
        xlim([0 size(Q_history, 1)+1])
        title(['State = ' num2str(i)])
    end
    l = legend('1', '2', '3');
    l.Box = 'off';
    l.Title.String = 'Action';
    xlabel('Trial')
    subplot(3, 1, 2); ylabel('Q-value')
    
end

function [] = plotReplayHistory(replay_history_abs, n_trials)

    figure;
    
    % Plot the trial from which sample is taken
    subplot(1, 5, 1:4); hold on
    for i = 1:length(replay_history_abs)
        xjitter = linspace(-0.1, 0.1, length(replay_history_abs{i}))+i; xjitter = xjitter(randperm(length(xjitter)));
        scatter(xjitter, replay_history_abs{i}, 3, 'k')
    end
    stairs(0.5:1:length(n_trials)+0.5, [cumsum(n_trials-1) sum(n_trials-1)], 'Color', [0.5 0.5 0.5])
    xlabel({'Replay following session number', ' '})
    ylabel('Trial from which replay sample is taken')
    xlim([0.5 length(n_trials)+0.5])
    ylim([0 sum(n_trials-1)+1])
    
    % Plot number of times each trial is replayed
    subplot(1, 5, 5);
    histogram(cell2mat(replay_history_abs'), 1:sum(n_trials-1), 'Orientation', 'horizontal', 'FaceColor', 'k', 'EdgeColor', 'none')
    set(gca, 'YTickLabel', [])
    xlabel({'Total number', 'time replayed'})
    axis tight; ylim([0 sum(n_trials-1)+1])
    box off
    
end        

function [] = plotPredictiveAccuracy(predictive_accuracy, s_history)

    figure; hold on
    x = find(s_history==1); plot(x, predictive_accuracy(s_history==1), 'LineWidth', 1.5)
    x = find(s_history==2); plot(x, predictive_accuracy(s_history==2), 'LineWidth', 1.5)
    x = find(s_history==3); plot(x, predictive_accuracy(s_history==3), 'LineWidth', 1.5)
    l = legend('s=1', 's=2', 's=3'); l.Box = 'off'; l.Orientation = 'Horizontal'; l.Location = 'southwest';
    xlabel('Trial'); ylabel('Brier score')
    
end

function [obs_bar, N_state] = get_obs(behav_data, section)

        all_actions = cell2mat(cellfun(@(x) x(2:end), behav_data.actions, 'UniformOutput', false));
        all_states = cell2mat(cellfun(@(x) x(1:end-1), behav_data.actions, 'UniformOutput', false));

    switch section
    
        case 'all'
            % Nothing
            
        case 'initial'
            all_actions = cell2mat(cellfun(@(x) x(2:end), behav_data.actions(1:15), 'UniformOutput', false));
            all_states = cell2mat(cellfun(@(x) x(1:end-1), behav_data.actions(1:15), 'UniformOutput', false));
            
        case 'revaluation'
            all_actions = cell2mat(cellfun(@(x) x(2:end), behav_data.actions(16:20), 'UniformOutput', false));
            all_states = cell2mat(cellfun(@(x) x(1:end-1), behav_data.actions(16:20), 'UniformOutput', false));
            
        case 'reversal'
            all_actions = cell2mat(cellfun(@(x) x(2:end), behav_data.actions(21:22), 'UniformOutput', false));
            all_states = cell2mat(cellfun(@(x) x(1:end-1), behav_data.actions(21:22), 'UniformOutput', false));
            
        otherwise
            error('params.learning_stage unknown')
            
    end
    
    % Sense-check
    assert(length(all_actions) == length(all_states), 'Wrong number of actions and/or states when computing observed outcomes')
    
    for iState = 1:3
        state = find(all_states==iState);
        N_state(iState) = numel(state);
        for iAction = 1:3
            obs_bar(iState,iAction) = numel(find(all_actions(state)==iAction))/N_state(iState);
        end
    end
    
    obs_bar(isnan(obs_bar)) = 0;
    
end


    