function [mvic] = get_mvic(emg_data, framerate, mvic_window_inseconds)
%UNTITLED8 Summary of this function goes here
%   Detailed explanation goes here

%framerate = 3000;
%mvic_window_inseconds = 1;
frame_cutoff = framerate * 200;

mvic_window = mvic_window_inseconds * framerate;
tib_filtered = filter((1/mvic_window)*ones(1,mvic_window) , 1, emg_data(1:min([frame_cutoff length(emg_data)]))');
mvic = max(tib_filtered);

end

