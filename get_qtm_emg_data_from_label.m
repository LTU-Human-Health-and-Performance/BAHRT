function [output_data] = get_qtm_emg_data_from_label(qtm_data, labels)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
    
    labelmatrix1 = repmat(qtm_data.Analog(2).Labels, length(labels), 1);
    labelmatrix2 = repmat(labels', 1, length(qtm_data.Analog(2).Labels));
    label_index = strcmp(labelmatrix1, labelmatrix2);
    
    emg_data = qtm_data.Analog(2).Data;
    output_data = emg_data(logical(sum(label_index, 1)), :);
    
end

