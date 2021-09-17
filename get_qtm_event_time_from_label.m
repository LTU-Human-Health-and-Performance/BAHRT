function [output] = get_qtm_event_time_from_label(qtm_data, labels)
%UNTITLED4 Summary of this function goes here
%   Detailed explanation goes here

    event_struct = struct2cell(qtm_data.Events);
    event_labels = string(event_struct(1,:));

    labelmatrix1 = repmat(event_labels, length(labels), 1);
    labelmatrix2 = repmat(labels', 1, length(event_labels));
    label_index = strcmp(labelmatrix1, labelmatrix2);
    
    if not (any(label_index)) 
        output = 0;
        return
    end
    
    output = event_struct{3,logical(sum(label_index, 1))};
	
end

