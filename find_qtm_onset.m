function [onset_index_filtered] = find_qtm_onset(emg_data, mvic_window_inseconds, framerate)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%     framerate = 3000;
%     mvic_window_inseconds = 0.1;

    mvic_window = mvic_window_inseconds * framerate;
    peak_removal_window = 2;
    
    onset_filtered = filter((1/mvic_window)*ones(1,mvic_window) , 1, emg_data');
    onset_filtered = circshift(onset_filtered, -mvic_window) - onset_filtered;
    onset_index_filtered = [0 0 0];
    onset_max_filtered = [0 0 0];
    for i = 1:3
        [onset_max_filtered(i), onset_index_filtered(i)] = max(onset_filtered);
        onset_filtered([onset_index_filtered(i)-framerate*peak_removal_window:onset_index_filtered(i)+framerate*peak_removal_window]) = 0;
    end
    onset_index_filtered = onset_index_filtered/framerate;
end

