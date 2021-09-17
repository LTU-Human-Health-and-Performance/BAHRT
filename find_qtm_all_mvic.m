function [max_tib_left, max_gas_left, max_tib_right, max_gas_right] = find_qtm_all_mvic(participant, mvic_window_inseconds)


    
    labels_left_tib = ["LTibAnt", "LTA", "LtibAnt", "LTibant"];
    labels_left_gas = ["LGasMed", "LGM", "LgasMed", "LGasmed"];
    labels_right_tib = ["RTibAnt", "RTA", "RtibAnt", "RTibant"];
    labels_right_gas = ["RGasMed", "RGM", "RgasMed", "RGasmed"];
% Find tib_max and gas_max for left side
    
    % Load left data from file and collect some overarching variables
        qtm_data = load_qtm_data(append('Rorelselabb/Sequence6DOF/', int2str(participant), '_footL.mat'));
        framerate_emg = qtm_data.Analog(2).Frequency;        

    % Get maximal left tib
        % Get the emg data for static left tib
            static_tib_left = get_qtm_emg_data_from_label(qtm_data, labels_left_tib);
            if isempty(static_tib_left)
                fprintf(append('No emg data found from current left mvic tib labels, participant ', int2str(participant), '\n'));
                max_tib_left = nan;
            else
                static_tib_left = process_and_filter_emg_data(static_tib_left);      
                max_tib_left = get_mvic(static_tib_left, framerate_emg, mvic_window_inseconds);
            end
    % Get maximal left gas
        % Get the emg data for static left gas
            static_gas_left = get_qtm_emg_data_from_label(qtm_data, labels_left_gas);
            if isempty(static_gas_left)
                fprintf(append('No emg data found from current left mvic gas labels, participant ', int2str(participant), '\n'));
                max_gas_left = nan;
            else
                static_gas_left = process_and_filter_emg_data(static_gas_left);
                max_gas_left = get_mvic(static_gas_left, framerate_emg, mvic_window_inseconds);
            end
% Find tib_max and gas_max for right side    
    % Load right data from file and collect some overarching variables
        qtm_data = load_qtm_data(append('Rorelselabb/Sequence6DOF/', int2str(participant), '_footR.mat'));
        framerate_emg = qtm_data.Analog(2).Frequency;        

    % Get maximal right tib
        % Get the emg data for static left tib
            static_tib_right = get_qtm_emg_data_from_label(qtm_data, labels_right_tib);
            if isempty(static_tib_right)
                fprintf(append('No emg data found from current right mvic tib labels, participant ', int2str(participant), '\n'));
                max_tib_right = nan;
            else
                static_tib_right = process_and_filter_emg_data(static_tib_right);
                max_tib_right = get_mvic(static_tib_right, framerate_emg, mvic_window_inseconds);
            end
    % Get maximal right gas
        % Get the emg data for static right gas
            static_gas_right = get_qtm_emg_data_from_label(qtm_data, labels_right_gas);
            if isempty(static_gas_right)
                fprintf(append('No emg data found from current right mvic gas labels, participant ', int2str(participant), '\n'));
                max_gas_right = nan;
            else
                static_gas_right = process_and_filter_emg_data(static_gas_right);
                if participant == 1790
                    static_gas_right(178*framerate_emg:181*framerate_emg) = 0;
                end
                max_gas_right = get_mvic(static_gas_right, framerate_emg, mvic_window_inseconds);
            end
end

