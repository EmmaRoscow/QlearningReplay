
function [result] = run25times(vector, behav_data, params)

    params.alpha = vector(1);
    params.discount = vector(2);
    params.epsilon = vector(3);
    if length(vector) >=4
        params.recency = vector(4);
    end
    if length(vector)==5
        params.rpe_recency = vector(5);
    end
   
    for i = 1:25
        params.seed = i;
        result(i) = mean(dynaQ(behav_data, params, false));
    end
    
    result = mean(result);
        
end

