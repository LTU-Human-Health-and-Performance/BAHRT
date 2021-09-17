function [position_array] = get_position_array(qtm_data,labelarray)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
    sum_position_array = 0;
    for label = labelarray
        index_tracker = find(ismember(qtm_data.Trajectories.Labeled.Labels, label));
        if isempty(index_tracker)
            continue;
        end
        
        position_array = qtm_data.Trajectories.Labeled.Data(index_tracker, :, :);
        position_array = squeeze(position_array);
        position_array = position_array(1:3,:);
        sum_position_array = sum_position_array + position_array;
    end
    position_array = sum_position_array/length(labelarray);
    
end

