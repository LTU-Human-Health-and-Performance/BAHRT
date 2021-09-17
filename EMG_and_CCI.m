complete_listof_participants = [1002 1011 1012 1035 1039 1041 1043 1056 1060 1092 1097 1128 1137 1139 1153 1217 1263 1606 1661 1663 1670 1677 1686 1700 1710 1711 1723 1724 1731 1742 1744 1758 1763 1776 1784 1788 1790 1791];
% 1044: No file

% ta med? = 1763
% testa onset-metod på 1744
% participant_with_several_sequences = [1700 1790];
  
listof_participants = 1791
%[1002 1011 1012 1035 1039 1041 1043 1056 1060 1092 1097 1128 1137 1139 1153 1217 1263 1606 1670 1677 1686 1700 1723 1731 1744 1758 1776 1784 1788 1790 1791];
%listof_participants = [1056 1791];
% Tested participants:
% 1002 1011 1012 1035 1039 1041 1043 1056 1060 1092
% 1097 1128 1137 1139 1153 1217 1263 1606 1677 1686 1700 1723 1731 1790 1791

% Problem participants and notes:
% 1011: No sequence data for left gas?
% 1044: Corrupted file?
% 1128: Possibly flipped sequence labels, left gas with right tib (manually flipped it in the file)
% 1606: Left mvic is rubbish, no marked onset
% 1661: Rubbish sequence data
% 1663: Rubbish sequence data, only left tib is good
% 1710: Rubbish sequence data
% 1711: Onset couldn't be found
% 1724: Rubbish sequence data, a couple of very high peaks
% 1742: Rubbish sequence and left mvic data
% 1744: Disturbed data, no marked onset, left gas seems to be scaled up
% 1763: Rubbish sequence data

% participant_with_several_sequences = [1700 1790];

figure_cci_rudolp = figure(1);
figure_cci_winter = figure(2);
figure_emg = figure(3);

time_plotwindow = 2;

if exist('output')
    clear output
end
output_data.participant = [];
output_data.t_start_to_onset_seq1 = [0];
output_data.t_start_to_onset_seq2 = [0];
output_data.t_start_to_onset_seq3 = [0];

for participant_index = 1:length(listof_participants)
    participant = listof_participants(participant_index);
    
% Find tib_max and gas_max for left side
    mvic_window_inseconds = 1;
    labels_left_tib = ["LTibAnt", "LTA", "LtibAnt", "LTibant"];
    labels_left_gas = ["LGasMed", "LGM", "LgasMed", "LGasmed"];
    labels_right_tib = ["RTibAnt", "RTA", "RtibAnt", "RTibant"];
    labels_right_gas = ["RGasMed", "RGM", "RgasMed", "RGasmed"];
    % Load left data from file and collect some overarching variables
        qtm_data = load_qtm_data(append('Rorelselabb/Sequence6DOF/', int2str(participant), '_footL.mat'));
        framerate = qtm_data.Analog(2).Frequency;        

    % Get maximal left tib
        % Get the emg data for static left tib
            static_tib_left = get_qtm_emg_data_from_label(qtm_data, labels_left_tib);
            if isempty(static_tib_left)
                fprintf(append('No emg data found from current left mvic tib labels, participant ', int2str(participant), '\n'));
                max_tib_left = nan;
            else
                static_tib_left = process_and_filter_emg_data(static_tib_left);      
                max_tib_left = get_mvic(static_tib_left, framerate, mvic_window_inseconds);
            end
    % Get maximal left gas
        % Get the emg data for static left gas
            static_gas_left = get_qtm_emg_data_from_label(qtm_data, labels_left_gas);
            if isempty(static_gas_left)
                fprintf(append('No emg data found from current left mvic gas labels, participant ', int2str(participant), '\n'));
                max_gas_left = nan;
            else
                static_gas_left = process_and_filter_emg_data(static_gas_left);
                max_gas_left = get_mvic(static_gas_left, framerate, mvic_window_inseconds);
                %max_gas_left = max_gas_left *15;
            end
% Find tib_max and gas_max for right side    
    % Load right data from file and collect some overarching variables
        qtm_data = load_qtm_data(append('Rorelselabb/Sequence6DOF/', int2str(participant), '_footR.mat'));
        framerate = qtm_data.Analog(2).Frequency;        

    % Get maximal right tib
        % Get the emg data for static left tib
            static_tib_right = get_qtm_emg_data_from_label(qtm_data, labels_right_tib);
            if isempty(static_tib_right)
                fprintf(append('No emg data found from current right mvic tib labels, participant ', int2str(participant), '\n'));
                max_tib_right = nan;
            else
                static_tib_right = process_and_filter_emg_data(static_tib_right);
                max_tib_right = get_mvic(static_tib_right, framerate, mvic_window_inseconds);
                %max_tib_right = max_tib_right *15
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
                    static_gas_right(178*framerate:181*framerate) = 0;
                end
                max_gas_right = get_mvic(static_gas_right, framerate, mvic_window_inseconds);
                %max_gas_right = max_gas_right *15
            end
            
% Get sequence data and find CCI
    % Load static data from file and collect some overarching variables
        qtm_data = load_qtm_data(append('Rorelselabb/Sequence6DOF/', int2str(participant), '_sequence.mat'));
        framerate = qtm_data.Analog(2).Frequency;
        
        % Fetch emg data from struct
            sequence_tib_left = get_qtm_emg_data_from_label(qtm_data, labels_left_tib);
            sequence_gas_left = get_qtm_emg_data_from_label(qtm_data, labels_left_gas);
            sequence_tib_right = get_qtm_emg_data_from_label(qtm_data, labels_right_tib);
            sequence_gas_right = get_qtm_emg_data_from_label(qtm_data, labels_right_gas);
            
        % Add if sequence file is split
            if participant == 1790
                qtm_data = load_qtm_data(append('Rorelselabb/Sequence6DOF/', int2str(participant), '_sequence2.mat'));
                sequence_tib_left = [sequence_tib_left get_qtm_emg_data_from_label(qtm_data, labels_left_tib)];
                sequence_gas_left = [sequence_gas_left get_qtm_emg_data_from_label(qtm_data, labels_left_gas)];
                sequence_tib_right = [sequence_tib_right get_qtm_emg_data_from_label(qtm_data, labels_right_tib)];
                sequence_gas_right = [sequence_gas_right get_qtm_emg_data_from_label(qtm_data, labels_right_gas)];
                
                qtm_data = load_qtm_data(append('Rorelselabb/Sequence6DOF/', int2str(participant), '_sequence.mat'));
            end
        
        % Left side
            sequence_tib_left = process_and_filter_emg_data(sequence_tib_left);
            sequence_tib_left = sequence_tib_left/max_tib_left;

            sequence_gas_left = process_and_filter_emg_data(sequence_gas_left);
            sequence_gas_left = sequence_gas_left/max_gas_left;

            
        % Right side
            sequence_tib_right = process_and_filter_emg_data(sequence_tib_right);
            sequence_tib_right = sequence_tib_right/max_tib_right;

            sequence_gas_right = process_and_filter_emg_data(sequence_gas_right);
            sequence_gas_right = sequence_gas_right/max_gas_right;
    
        % Left CCI
            cci_rudolp_left = (sequence_gas_left./sequence_tib_left).^((sequence_tib_left > sequence_gas_left)*2 - 1) .* (sequence_gas_left + sequence_tib_left);
            cci_winter_left = min([sequence_tib_left; sequence_gas_left]) .* 2 ./ (sequence_tib_left + sequence_gas_left);
        % Right CCI
            cci_rudolp_right = (sequence_gas_right./sequence_tib_right).^((sequence_tib_right > sequence_gas_right)*2 - 1) .* (sequence_gas_right + sequence_tib_right);
            cci_winter_right = min([sequence_tib_right; sequence_gas_right]) .* 2 ./ (sequence_tib_right + sequence_gas_right);
        % Combined CCI
            cci_rudolp_combined = (cci_rudolp_left + cci_rudolp_right) / 2;
            cci_winter_combined = (cci_winter_left + cci_winter_right) / 2;            
        % Collect CCI in one matrix
            cci_tot_rudolp = [cci_rudolp_left; cci_rudolp_right; cci_rudolp_combined]';
            cci_tot_winter = [cci_winter_left; cci_winter_right; cci_winter_combined]';
            
            
        labels_onset = [
            "onset1", "Onset1"; 
            "onset2", "Onset2"; 
            "onset3", "Onset3"];
        labels_start = [
            "Start1", "start1", "Start";
            "Start2", "start2", "Placeholder";
            "Start3", "start3", "Placeholder"];
        
    % Fetch start and onset
        timestamp_start = [0 0 0];
        timestamp_onset = [0 0 0];
        for i = 1:3
            timestamp_start(i) = get_qtm_event_time_from_label(qtm_data, labels_start(i,:));
            timestamp_onset(i) = get_qtm_event_time_from_label(qtm_data, labels_onset(i,:));
            if participant == 1790
                if i == 1
                    timestamps_split_file_offset = length(qtm_data.Analog(2).Data(1,:))/framerate;
                    qtm_data = load_qtm_data(append('Rorelselabb/Sequence6DOF/', int2str(participant), '_sequence2.mat'));
                else
                    timestamp_start(i) = timestamp_start(i) + timestamps_split_file_offset;
                    timestamp_onset(i) = timestamp_onset(i) + timestamps_split_file_offset;
                end
            end
        end
        if any(timestamp_onset == 0)
            fprintf(append('No timestamp for onset found using current labels for sequence ', int2str(i), ', participant ', int2str(participant), '\n'));
            timestamp_onset = sort(find_qtm_onset(sequence_tib_left, 0.1, framerate));
        end
    % Plot CCI for all 3 sequences
        set(0, 'CurrentFigure', figure_cci_rudolp)
        clf
        hold on
        figure_cci_rudolp.WindowState = 'maximized';
        legend_labels_rudolp = ["Left CCI_{Rudolp}", "Right CCI_{Rudolp}", "Mean CCI_{Rudolp}", "Start", "Onset"]';
        
        set(0, 'CurrentFigure', figure_cci_winter)
        clf
        hold on
        figure_cci_winter.WindowState = 'maximized';
        legend_labels_winter = ["Left CCI_{Winter}", "Right CCI_{Winter}", "Mean CCI_{Winter}", "Start", "Onset"]';
        
        timestamp_height = 1.5;
        
        % Plot data around 3 timestamps
        for i = 1:3
            if any(timestamp_start == 0)
                fprintf(append('No timestamp for start found using current labels for sequence ', int2str(i), ', participant ', int2str(participant), '\n'));
                continue
            else
                if any(timestamp_onset == 0)
                    fprintf(append('No timestamp for onset found using current labels for sequence ', int2str(i), ', participant ', int2str(participant), '\n'));
                    timestamp_onset = find_qtm_onset(sequence_tib_left, 0.1, framerate);
                end
                set(0, 'CurrentFigure', figure_cci_rudolp)
                subplot(3,1,i); 
                hold on
                plot([1:length(cci_rudolp_left)]/framerate - timestamp_start(i), cci_tot_rudolp);
                plot([1 1] * timestamp_start(i) - timestamp_start(i), [0 timestamp_height]);
                plot([1 1] * timestamp_onset(i) - timestamp_start(i), [0 timestamp_height]);
                ax = gca;
                ax.XLim = [-1 1] + timestamp_start(i) - timestamp_start(i);
                legend(legend_labels_rudolp);
                title(append("CCI_{Rudolp} - Sequence ", int2str(i),' - particiant ', int2str(participant)));
                grid on

                set(0, 'CurrentFigure', figure_cci_winter)
                subplot(3,1,i); 
                hold on
                plot([1:length(cci_rudolp_left)]/framerate - timestamp_start(i), cci_tot_winter);
                plot([1 1] * timestamp_start(i) - timestamp_start(i), [0 timestamp_height]);
                plot([1 1] * timestamp_onset(i) - timestamp_start(i), [0 timestamp_height]);
                xlim([-time_plotwindow time_plotwindow]);
                legend(legend_labels_winter);
                title(append("CCI_{Winter} - Sequence ", int2str(i),' - particiant ', int2str(participant)));
                grid on
            end
            
            output_data.('participant')(participant_index, 1) = participant;
            output_data.(append('t_start_to_onset_seq', int2str(i)))(participant_index, 1) = (timestamp_onset(i) - timestamp_start(i)) * 1000;
            % Left side
            output_data.(append('cci_rudolp_left_feed_forward_100ms_seq', int2str(i)))(participant_index, 1) = mean(cci_rudolp_left(floor((timestamp_start(i)-0.1)*3000):floor((timestamp_start(i))*3000)));
            output_data.(append('cci_rudolp_left_feed_backward_100ms_seq', int2str(i)))(participant_index, 1) = mean(cci_rudolp_left(floor((timestamp_onset(i))*3000):floor((timestamp_onset(i)+0.1)*3000)));
            output_data.(append('cci_rudolp_left_feed_forward_200ms_seq', int2str(i)))(participant_index, 1) = mean(cci_rudolp_left(floor((timestamp_start(i)-0.2)*3000):floor((timestamp_start(i))*3000)));
            output_data.(append('cci_rudolp_left_feed_backward_200ms_seq', int2str(i)))(participant_index, 1) = mean(cci_rudolp_left(floor((timestamp_onset(i))*3000):floor((timestamp_onset(i)+0.2)*3000)));
            
            output_data.(append('cci_winter_left_feed_forward_100ms_seq', int2str(i)))(participant_index, 1) = mean(cci_winter_left(floor((timestamp_start(i)-0.1)*3000):floor((timestamp_start(i))*3000)));
            output_data.(append('cci_winter_left_feed_backward_100ms_seq', int2str(i)))(participant_index, 1) = mean(cci_winter_left(floor((timestamp_onset(i))*3000):floor((timestamp_onset(i)+0.1)*3000)));
            output_data.(append('cci_winter_left_feed_forward_200ms_seq', int2str(i)))(participant_index, 1) = mean(cci_winter_left(floor((timestamp_start(i)-0.2)*3000):floor((timestamp_start(i))*3000)));
            output_data.(append('cci_winter_left_feed_backward_200ms_seq', int2str(i)))(participant_index, 1) = mean(cci_winter_left(floor((timestamp_onset(i))*3000):floor((timestamp_onset(i)+0.2)*3000)));
            % Right side
            output_data.(append('cci_rudolp_right_feed_forward_100ms_seq', int2str(i)))(participant_index, 1) = mean(cci_rudolp_right(floor((timestamp_start(i)-0.1)*3000):floor((timestamp_start(i))*3000)));
            output_data.(append('cci_rudolp_right_feed_backward_100ms_seq', int2str(i)))(participant_index, 1) = mean(cci_rudolp_right(floor((timestamp_onset(i))*3000):floor((timestamp_onset(i)+0.1)*3000)));
            output_data.(append('cci_rudolp_right_feed_forward_200ms_seq', int2str(i)))(participant_index, 1) = mean(cci_rudolp_right(floor((timestamp_start(i)-0.2)*3000):floor((timestamp_start(i))*3000)));
            output_data.(append('cci_rudolp_right_feed_backward_200ms_seq', int2str(i)))(participant_index, 1) = mean(cci_rudolp_right(floor((timestamp_onset(i))*3000):floor((timestamp_onset(i)+0.2)*3000)));
            
            output_data.(append('cci_winter_right_feed_forward_100ms_seq', int2str(i)))(participant_index, 1) = mean(cci_winter_right(floor((timestamp_start(i)-0.1)*3000):floor((timestamp_start(i))*3000)));
            output_data.(append('cci_winter_right_feed_backward_100ms_seq', int2str(i)))(participant_index, 1) = mean(cci_winter_right(floor((timestamp_onset(i))*3000):floor((timestamp_onset(i)+0.1)*3000)));
            output_data.(append('cci_winter_right_feed_forward_200ms_seq', int2str(i)))(participant_index, 1) = mean(cci_winter_right(floor((timestamp_start(i)-0.2)*3000):floor((timestamp_start(i))*3000)));
            output_data.(append('cci_winter_right_feed_backward_200ms_seq', int2str(i)))(participant_index, 1) = mean(cci_winter_right(floor((timestamp_onset(i))*3000):floor((timestamp_onset(i)+0.2)*3000)));
            % Combined
            output_data.(append('cci_rudolp_combined_feed_forward_100ms_seq', int2str(i)))(participant_index, 1) = mean(cci_rudolp_combined(floor((timestamp_start(i)-0.1)*3000):floor((timestamp_start(i))*3000)));
            output_data.(append('cci_rudolp_combined_feed_backward_100ms_seq', int2str(i)))(participant_index, 1) = mean(cci_rudolp_combined(floor((timestamp_onset(i))*3000):floor((timestamp_onset(i)+0.1)*3000)));
            output_data.(append('cci_rudolp_combined_feed_forward_200ms_seq', int2str(i)))(participant_index, 1) = mean(cci_rudolp_combined(floor((timestamp_start(i)-0.2)*3000):floor((timestamp_start(i))*3000)));
            output_data.(append('cci_rudolp_combined_feed_backward_200ms_seq', int2str(i)))(participant_index, 1) = mean(cci_rudolp_combined(floor((timestamp_onset(i))*3000):floor((timestamp_onset(i)+0.2)*3000)));
            
            output_data.(append('cci_winter_combined_feed_forward_100ms_seq', int2str(i)))(participant_index, 1) = mean(cci_winter_combined(floor((timestamp_start(i)-0.1)*3000):floor((timestamp_start(i))*3000)));
            output_data.(append('cci_winter_combined_feed_backward_100ms_seq', int2str(i)))(participant_index, 1) = mean(cci_winter_combined(floor((timestamp_onset(i))*3000):floor((timestamp_onset(i)+0.1)*3000)));
            output_data.(append('cci_winter_combined_feed_forward_200ms_seq', int2str(i)))(participant_index, 1) = mean(cci_winter_combined(floor((timestamp_start(i)-0.2)*3000):floor((timestamp_start(i))*3000)));
            output_data.(append('cci_winter_combined_feed_backward_200ms_seq', int2str(i)))(participant_index, 1) = mean(cci_winter_combined(floor((timestamp_onset(i))*3000):floor((timestamp_onset(i)+0.2)*3000)));
            
        end
        saveas(figure_cci_rudolp, append('Rorelselabb/Sequence6DOF_Output/',int2str(participant), '_CCI_rudolp_tibgas', '.png'));
        saveas(figure_cci_winter, append('Rorelselabb/Sequence6DOF_Output/',int2str(participant), '_CCI_winter_tibgas', '.png'));
        
        
    % Plot emg for all 3 sequences
        set(0, 'CurrentFigure', figure_emg)
        clf
        hold on
        figure_emg.WindowState = 'maximized';
        legend_labels = ["Left TibAnt", "Left GasMed", "Right TibAnt"; "Right GasMed", "Start", "Onset"]';
        timestamp_height = 2;
        
        emg_tot = [sequence_tib_left; sequence_gas_left; sequence_tib_right; sequence_gas_right]';

         
        % Plot data around 3 timestamps
        for i = 1:3
            subplot(3,1,i); 
            hold on
            if any(timestamp_start == 0)
                fprintf(append('No timestamp for start found using current labels for sequence ', int2str(i), ', participant ', int2str(participant), '\n'));
                continue
            else
                plot([1:length(cci_rudolp_left)]/framerate - timestamp_start(i), emg_tot);
                plot([1 1] * timestamp_start(i) - timestamp_start(i), [0 timestamp_height]);
            end
            
            if any(timestamp_onset == 0)
                fprintf(append('No timestamp for onset found using current labels for sequence ', int2str(i), ', participant ', int2str(participant), '\n'));
                timestamp_onset = find_qtm_onset(sequence_tib_left, 0.1, framerate);
            end
            plot([1 1] * timestamp_onset(i) - timestamp_start(i), [0 timestamp_height]);
            
            xlim([-time_plotwindow time_plotwindow] + timestamp_start(i) - timestamp_start(i));
            legend(legend_labels);
            title(append("Sequence ", int2str(i),', particiant ', int2str(participant)));
            grid on
        end
        saveas(figure_emg, append('Rorelselabb/Sequence6DOF_Output/',int2str(participant), '_EMG_tibgas', '.png'));
 
        
        
        % Only plotting the plot interval data
        %plot([1:2*framerate+1]/framerate - 1, cci_tot(floor((timestamp-1)*framerate):floor((timestamp+1)*framerate),:));
        %plot([1 1] * timestamp - timestamp, [0 1.5]);
        
end
            
% Export to excell
    if exist('cci_table')
        clear cci_table
    end
    cci_table = struct2table(output_data); %table(participant, t_start_to_onset_seq1, cci_left_feed_forward_100ms_seq1, cci_right_feed_forward_100ms_seq1, cci_combined_feed_forward_100ms_seq1, cci_left_feed_forward_200ms_seq1, cci_right_feed_forward_200ms_seq1, cci_combined_feed_forward_200ms_seq1, t_start_to_onset_seq2, cci_left_feed_forward_100ms_seq2, cci_right_feed_forward_100ms_seq2, cci_combined_feed_forward_100ms_seq2, cci_left_feed_forward_200ms_seq2, cci_right_feed_forward_200ms_seq2, cci_combined_feed_forward_200ms_seq2, t_start_to_onset_seq3, cci_left_feed_forward_100ms_seq3, cci_right_feed_forward_100ms_seq3, cci_combined_feed_forward_100ms_seq3, cci_left_feed_forward_200ms_seq3, cci_right_feed_forward_200ms_seq3, cci_combined_feed_forward_200ms_seq3);
    if isfile('Rorelselabb/cci_table.xlsx')
        delete('Rorelselabb/cci_table.xlsx');
    end
    writetable(cci_table,'Rorelselabb/cci_table.xlsx');

    
    %clear
    
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Removed code, that might be useful later
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%         % Get timestamp for max for static left tib in temporary variable
%             timestamp = get_qtm_event_time_from_label(qtm_data, ["Maxtib", "LTib"]);
%             if isempty(timestamp)
%                 fprintf(append('No timestamp found from static tib labels, participant ', int2str(participant), '\n'));
%             end
%         % Find maximal tib from the second around the timestamp
%             max_tib = mean(static_tib_left((timestamp-0.5)*framerate:(timestamp+0.5)*framerate));
%         clear timestamp


%         % Get timestamp for max for static left gas in temporary variable
%             timestamp = get_qtm_event_time_from_label(qtm_data, ["MaxGas", "LGas"]);
%             if isempty(timestamp)
%                 fprintf(append('No timestamp found from static gas labels, participant ', int2str(participant), '\n'));
%             end
%         % Find maximal gas from the second around the timestamp
%             max_gas = mean(static_gas_left((timestamp-0.5)*framerate:(timestamp+0.5)*framerate));
%         clear timestamp
