function [output] = process_and_filter_emg_data(data_vector)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
nan_indexes = or(isnan(data_vector), isinf(data_vector));
data_vector(nan_indexes) = 0;
data_vector = bandpass(data_vector, [20 500], 3000);
data_vector = envelope(data_vector, 150, 'rms');
data_vector(nan_indexes) = NaN;
output = data_vector;
return;
end

