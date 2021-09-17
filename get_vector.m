function [vector] = get_vector(qtm_datastruct, tracker1, tracker2)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

index_tracker1 = find(ismember(qtm_datastruct.Trajectories.Labeled.Labels, tracker1));
index_tracker2 = find(ismember(qtm_datastruct.Trajectories.Labeled.Labels, tracker2));

tracker1_positions = qtm_datastruct.Trajectories.Labeled.Data(index_tracker1, :, :);
tracker1_positions = squeeze(tracker1_positions);
tracker1_positions = tracker1_positions(1:3,:);

tracker2_positions = qtm_datastruct.Trajectories.Labeled.Data(index_tracker2, :, :);
tracker2_positions = squeeze(tracker2_positions);
tracker2_positions = tracker2_positions(1:3,:);

vector = tracker2_positions - tracker1_positions;

end

