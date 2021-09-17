function [qtm_data] = load_qtm_data(path_datafile)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    qtm_data = load(path_datafile);
    qtm_data = struct2cell(qtm_data);
    qtm_data = [qtm_data{:}];
end

