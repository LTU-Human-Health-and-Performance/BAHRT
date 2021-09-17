clear
delete('Rorelselabb/Sequence6DOF_Output/*.png')
%% Participant list

%complete_listof_participants = [1002 1011 1012 1035 1039 1041 1043 1056 1060 1092 1097 1128 1137 1139 
%    1153 1217 1263 1606 1661 1663 1670 1677 1686 1700 1710 1711 1723 1724 1731 1742 1744 1758 1763 1776 
%   1784 1788 1790 1791];
% 1044: No file

% ta med? = 1763
% testa onset-metod på 1744
% participant_with_several_sequences = [1700 1790];
  
%listof_participants = [1002 1011 1012 1035 1039 1041 1043 1056 1060 1092 1097 1128 1137 1139 1153 1217
%    1263 1606 1670 1677 1686 1700 1723 1731 1744 1758 1776 1784 1788 1790 1791];
% 1661 1663 1710 1711 1724 1742 1763

% Tested participants:
% 1002 1011 1012 1035 1039 1041 1043 1056 1060 1092
% 1097 1128 1137 1139 1153 1217 1263 1606 1677 1686 1700 1723 1731 1790 1791

participants_in_process = [1002 1011 1012 1035 1039 1041 1043 1056 1060 1092 1097 1128 1137 1139 1153 1217 1263 1606 1661 1670 1677 1686 1700 1723 1731 1776 1784 1788 1791 1710 1711 1742 1744 1758 1790 1791];
% Works: 
% 1002 1011 1012 1035 1039 1041 1043 1056 1060 1092 1097 1128 1137 1139 1153
% 1217 1263 1606 1661 1670 1677 1686 1700 1723 1731 1776 1784 1788 1791
% 1710 1711 1742 1744 1758 1790 1791
% Not working:

% 1663 1724 1763 corrupt data

% Other notes
% 1039 - no foot angles
% 1041 1139 - no hip and knee angles
% 1056 - no foot angle seq 1
% 1670 - no COP
% 1700 - no foot angle seq 3, no back angle, no COP
% 1661 1710 1711 1742 - No CCI or COP
% 1744, spiky cop
% 
participant = 1791;
participant_index = 1;

output_data.participant = [];
output_data.t_sequence_seq1 = [0];
output_data.t_sequence_seq2 = [0];
output_data.t_sequence_seq3 = [0];
%
figure_angles = figure(1);
figure_cci_cop = figure(2);


for participant_index = 1:length(participants_in_process)
    participant = participants_in_process(participant_index);
    
    %
    
    %% MVIC from foot files for CCI

    mvic_window_inseconds = 1;
    [max_tib_left, max_gas_left, max_tib_right, max_gas_right] = find_qtm_all_mvic(participant, mvic_window_inseconds);

    %

    %% Get static angles as a reference

    qtm_data = load_qtm_data(append('Rorelselabb/Sequence6DOF/', int2str(participant), '_static.mat'));
    framerate_trackers = qtm_data.FrameRate;

    vector1 = get_position_array(qtm_data, ["LFM1" "LFM2" "LFM5"]) - get_position_array(qtm_data, "LFCC");
    vector2 = get_position_array(qtm_data, "LTTC") - get_position_array(qtm_data, ["LTAM" "LFAL"]);
    angle.static.foot_left = mean(find_angle(vector1, vector2));

    vector1 = get_position_array(qtm_data, ["RFM1" "RFM2" "RFM5"]) - get_position_array(qtm_data, "RFCC");
    vector2 = get_position_array(qtm_data, "RTTC") - get_position_array(qtm_data, ["RTAM" "RFAL"]);
    angle.static.foot_right = mean(find_angle(vector1, vector2));

    vector1 = get_position_array(qtm_data, ["CV7"]) - get_position_array(qtm_data, ["RIPS" "LIPS"]);
    vector2 = [zeros(1, length(vector1)); -ones(1, length(vector1)); zeros(1, length(vector1))];
    angle.static.back_z = mean(find_angle(vector1, vector2));

    vector2 = get_position_array(qtm_data, ["RIPS" "LIPS"]) - get_position_array(qtm_data, ["LIAS" "RIAS"]);
    angle.static.back = mean(find_angle(vector1, vector2));

    vector1 = get_position_array(qtm_data, ["LFEM-1" "LFEM-4"]) - get_position_array(qtm_data, ["LFEM-2" "LFEM-3"]);
    vector2 = get_position_array(qtm_data, ["LTAM" "LFAL"]) - get_position_array(qtm_data, "LTTC");
    angle.static.knee_left = mean(find_angle(vector1, vector2));

    vector1 = get_position_array(qtm_data, ["RFEM-1" "RFEM-4"]) - get_position_array(qtm_data, ["RFEM-2" "RFEM-3"]);
    vector2 = get_position_array(qtm_data, ["RTAM" "RFAL"]) - get_position_array(qtm_data, "RTTC");
    angle.static.knee_right = mean(find_angle(vector1, vector2));

    vector1 = get_position_array(qtm_data, ["LFEM-1" "LFEM-2"]) - get_position_array(qtm_data, ["LFEM-3" "LFEM-4"]);
    vector2 = get_position_array(qtm_data, ["RIAS" "LIAS"]) - get_position_array(qtm_data, ["LIPS" "RIPS"]);
    angle.static.hip_left = mean(find_angle(vector1, vector2));

    vector1 = get_position_array(qtm_data, ["RFEM-1" "RFEM-2"]) - get_position_array(qtm_data, ["RFEM-3" "RFEM-4"]);
    vector2 = get_position_array(qtm_data, ["RIAS" "LIAS"]) - get_position_array(qtm_data, ["LIPS" "RIPS"]);
    angle.static.hip_right = mean(find_angle(vector1, vector2));

    cop.static = mean(qtm_data.Force(2).COP, 2);

    %

    %% Sequence data extraction

        %% Load sequence data
        qtm_data = load_qtm_data(append('Rorelselabb/Sequence6DOF/', int2str(participant), '_sequence.mat'));
        framerate_force = qtm_data.Force.Frequency;    
        %

        %% Get start times
        labels_start = [
                "Start1", "start1", "Start";
                "Start2", "start2", "Placeholder";
                "Start3", "start3", "Placeholder"];
        timestamp_start = [0 0 0];
        for i = [1 2 3]
            timestamp_start(i) = get_qtm_event_time_from_label(qtm_data, labels_start(i,:));
            if participant == 1790
                if i == 1
                    timestamps_split_file_offset = length(qtm_data.Analog(2).Data(1,:))/framerate_force;
                    qtm_data_seq1 = qtm_data;
                    qtm_data = load_qtm_data(append('Rorelselabb/Sequence6DOF/', int2str(participant), '_sequence2.mat'));
                    qtm_data.Trajectories.Labeled.Data = cat(3, qtm_data_seq1.Trajectories.Labeled.Data, qtm_data.Trajectories.Labeled.Data);
                    %qtm_data.Analog(2).Data = cat(2, qtm_data_seq1.Analog(2).Data, qtm_data.Analog(2).Data);
                    qtm_data.Force(2).COP = cat(2, qtm_data_seq1.Force(2).COP, qtm_data.Force(2).COP);
                else
                    timestamp_start(i) = timestamp_start(i) + timestamps_split_file_offset;
                end
            end
        end
        %
        timestamp_start = sort(timestamp_start);

        %% Find sequence end timestamp, in seconds and chronological order
        temp = get_position_array(qtm_data, ["platform - 2", "Platform - 2"]);
        temp = temp(1,:);
        timestamp_end = [0 0 0];
        for i = [1 2 3]
            [temp2, timestamp_end(i)] = max(temp);
            temp((timestamp_end(i)-framerate_trackers):(timestamp_end(i)+framerate_trackers)) = min(temp);
        end
        timestamp_end = sort(timestamp_end/framerate_trackers);
         %


        %% Foot  
        % Left foot
            vector1 = get_position_array(qtm_data, ["LFM1" "LFM2" "LFM5"]) - get_position_array(qtm_data, "LFCC");
            vector2 = get_position_array(qtm_data, "LTTC") - get_position_array(qtm_data, ["LTAM" "LFAL"]);
            angle.foot.left =  find_angle(vector1, vector2) - angle.static.foot_left;

        % Right foot
            vector1 = get_position_array(qtm_data, ["RFM1" "RFM2" "RFM5"]) - get_position_array(qtm_data, "RFCC");
            vector2 = get_position_array(qtm_data, "RTTC") - get_position_array(qtm_data, ["RTAM" "RFAL"]);
            angle.foot.right = find_angle(vector1, vector2) - angle.static.foot_right;

        % Combined for both feet
            angle.foot.combined = (angle.foot.left + angle.foot.right)/2;
        %

        %% Back
        % Back compared to z
            vector1 = get_position_array(qtm_data, ["CV7"]) - get_position_array(qtm_data, ["RIPS" "LIPS"]);
            vector2 = [zeros(1, length(vector1)); -ones(1, length(vector1)); zeros(1, length(vector1))];
            angle.back.z = find_angle(vector1, vector2)-angle.static.back_z;

        % Back compared to pelvis
            vector2 = get_position_array(qtm_data, ["RIPS" "LIPS"]) - get_position_array(qtm_data, ["LIAS" "RIAS"]);
            angle.back.back_to_hip = find_angle(vector1, vector2)-angle.static.back;    
        %   

        %% Knee
        % Left knee
            vector1 = get_position_array(qtm_data, ["LFEM-1" "LFEM-4"]) - get_position_array(qtm_data, ["LFEM-2" "LFEM-3"]);
            vector2 = get_position_array(qtm_data, ["LTAM" "LFAL"]) - get_position_array(qtm_data, "LTTC");
            angle.knee.left = find_angle(vector1, vector2) - angle.static.knee_left;

        % Right knee
            vector1 = get_position_array(qtm_data, ["RFEM-1" "RFEM-4"]) - get_position_array(qtm_data, ["RFEM-2" "RFEM-3"]);
            vector2 = get_position_array(qtm_data, ["RTAM" "RFAL"]) - get_position_array(qtm_data, "RTTC");
            angle.knee.right = find_angle(vector1, vector2) - angle.static.knee_right;

        % Knee combined
            angle.knee.combined = (angle.knee.left + angle.knee.right)/2;    
        %

        %% Hip
        % Left hip
            vector1 = get_position_array(qtm_data, ["LFEM-1" "LFEM-2"]) - get_position_array(qtm_data, ["LFEM-3" "LFEM-4"]);
            vector2 = get_position_array(qtm_data, ["RIAS" "LIAS"]) - get_position_array(qtm_data, ["LIPS" "RIPS"]);
            angle.hip.left = find_angle(vector1, vector2) - angle.static.hip_left;

        % Right hip
            vector1 = get_position_array(qtm_data, ["RFEM-1" "RFEM-2"]) - get_position_array(qtm_data, ["RFEM-3" "RFEM-4"]);
            vector2 = get_position_array(qtm_data, ["RIAS" "LIAS"]) - get_position_array(qtm_data, ["LIPS" "RIPS"]);
            angle.hip.right = find_angle(vector1, vector2) - angle.static.hip_left;

        % Hip combined
            angle.hip.combined = (angle.hip.left + angle.hip.right)/2;
        %

        %% CCI
        labels.tib.left = ["LTibAnt", "LTA", "LtibAnt", "LTibant"];
        labels.gas.left = ["LGasMed", "LGM", "LgasMed", "LGasmed"];
        labels.tib.right = ["RTibAnt", "RTA", "RtibAnt", "RTibant"];
        labels.gas.right = ["RGasMed", "RGM", "RgasMed", "RGasmed"];

        % Fetch emg data from struct
            sequence.tib.left = get_qtm_emg_data_from_label(qtm_data, labels.tib.left);
            sequence.gas.left = get_qtm_emg_data_from_label(qtm_data, labels.gas.left);
            sequence.tib.right = get_qtm_emg_data_from_label(qtm_data, labels.tib.right);
            sequence.gas.right = get_qtm_emg_data_from_label(qtm_data, labels.gas.right);

        % Add if sequence file is split
            if participant == 1790
                %qtm_data = load_qtm_data(append('Rorelselabb/Sequence6DOF/', int2str(participant), '_sequence2.mat'));
                sequence.tib.left = [get_qtm_emg_data_from_label(qtm_data_seq1, labels.tib.left) sequence.tib.left];
                sequence.gas.left = [get_qtm_emg_data_from_label(qtm_data_seq1, labels.gas.left) sequence.gas.left];
                sequence.tib.right = [get_qtm_emg_data_from_label(qtm_data_seq1, labels.tib.right) sequence.tib.right];
                sequence.gas.right = [get_qtm_emg_data_from_label(qtm_data_seq1, labels.gas.right) sequence.gas.right];

                %qtm_data = load_qtm_data(append('Rorelselabb/Sequence6DOF/', int2str(participant), '_sequence.mat'));
            end

        % Left side
            sequence.tib.left = process_and_filter_emg_data(sequence.tib.left);
            sequence.tib.left = sequence.tib.left/max_tib_left;

            sequence.gas.left = process_and_filter_emg_data(sequence.gas.left);
            sequence.gas.left = sequence.gas.left/max_gas_left;


        % Right side
            sequence.tib.right = process_and_filter_emg_data(sequence.tib.right);
            sequence.tib.right = sequence.tib.right/max_tib_right;

            sequence.gas.right = process_and_filter_emg_data(sequence.gas.right);
            sequence.gas.right = sequence.gas.right/max_gas_right;

        % Calculate CCI
            cci_rudolp_left = (sequence.gas.left./sequence.tib.left).^((sequence.tib.left > sequence.gas.left)*2 - 1) .* (sequence.gas.left + sequence.tib.left);
            cci_rudolp_right = (sequence.gas.right./sequence.tib.right).^((sequence.tib.right > sequence.gas.right)*2 - 1) .* (sequence.gas.right + sequence.tib.right);
            cci_rudolp_combined = (cci_rudolp_left + cci_rudolp_right) / 2;

        %% COP


        subplot(5,1,5)
        %plotgrid1 = [2,1];
        %subplot(plotgrid1(1),plotgrid1(2),1);
        hold on
        grid on

        % COP: 1 = x = höger-vänster, 2 = y = fram-bak
        % Positiv x är vänster för deltagaren, possitiv y är framåt

        platform = get_position_array(qtm_data, ["platform - 2", "Platform - 2"]);
        cop.sequence = qtm_data.Force(2).COP;



    %% Plotting angles
        set(0, 'CurrentFigure', figure_angles)
        clf

        plotwindow_before = 1;
        plotwindow_after = 3;
        calcwindow_before = 0.1;

        output_data.participant(participant_index, 1) = participant;
        output_data.t_sequence_seq1(participant_index, 1) = (timestamp_end(1)-timestamp_start(1))*1000;
        output_data.t_sequence_seq2(participant_index, 1) = (timestamp_end(2)-timestamp_start(2))*1000;
        output_data.t_sequence_seq3(participant_index, 1) = (timestamp_end(3)-timestamp_start(3))*1000;
        %
        
        % Pad ends with nan to avoid error from stopping measurements early
        % (too short arrays)
        padding = NaN(1, 700);
        angle.foot.left = [angle.foot.left, padding];
        angle.foot.right = [angle.foot.right, padding];
        angle.foot.combined = [angle.foot.combined, padding];
        angle.knee.left = [angle.knee.left, padding];
        angle.knee.right = [angle.knee.right, padding];
        angle.knee.combined = [angle.knee.combined, padding];
        angle.hip.left = [angle.hip.left, padding];
        angle.hip.right = [angle.hip.right, padding];
        angle.hip.combined = [angle.hip.combined, padding];
        angle.back.back_to_hip = [angle.back.back_to_hip, padding];
        angle.back.z = [angle.back.z, padding];
        cci_rudolp_combined = [cci_rudolp_combined, NaN(1, 4000)];
        cop.sequence = [cop.sequence, NaN(3, 4000)];
        

        %% Plotting foot
        for seq_num = [1 2 3]
            %seq_num = 1;
            subplot(5,3,0+seq_num)
            hold on
            grid on

            % All angles
            plotwindow_data = angle.foot.left(round((timestamp_start(seq_num)-plotwindow_before)*framerate_trackers):round((timestamp_start(seq_num)+plotwindow_after)*framerate_trackers));
            plot((1:length(plotwindow_data))/framerate_trackers - 1, plotwindow_data)
            plotwindow_data = angle.foot.right(round((timestamp_start(seq_num)-plotwindow_before)*framerate_trackers):round((timestamp_start(seq_num)+plotwindow_after)*framerate_trackers));
            plot((1:length(plotwindow_data))/framerate_trackers - 1, plotwindow_data)
            plotwindow_data = angle.foot.combined(round((timestamp_start(seq_num)-plotwindow_before)*framerate_trackers):round((timestamp_start(seq_num)+plotwindow_after)*framerate_trackers));
            plot((1:length(plotwindow_data))/framerate_trackers - 1, plotwindow_data)

            % Timestamps
            plot([0 0 nan [1 1]*(timestamp_end(seq_num)-timestamp_start(seq_num))], [-100 100 nan -100 100])

            % Calculate mean from interval and standard deviation
            angle.foot.mean.(append('seq', int2str(seq_num))) = mean(angle.foot.combined(round(timestamp_start(seq_num)*framerate_trackers):round(timestamp_end(seq_num)*framerate_trackers)));
            angle.foot.std.(append('seq', int2str(seq_num))) = std(angle.foot.combined(round(timestamp_start(seq_num)*framerate_trackers):round(timestamp_end(seq_num)*framerate_trackers)));

            angle.foot.mean_before.(append('seq', int2str(seq_num))) = mean(angle.foot.combined(round((timestamp_start(seq_num)-calcwindow_before)*framerate_trackers):round(timestamp_start(seq_num)*framerate_trackers)));
            angle.foot.std_before.(append('seq', int2str(seq_num))) = std(angle.foot.combined(round((timestamp_start(seq_num)-calcwindow_before)*framerate_trackers):round(timestamp_start(seq_num)*framerate_trackers)));

            % Save mean and std to output
            output_data.(append('foot_mean_seq',int2str(seq_num)))(participant_index, 1) = angle.foot.mean.(append('seq', int2str(seq_num)));
            output_data.(append('foot_std_seq',int2str(seq_num)))(participant_index, 1) = angle.foot.std.(append('seq', int2str(seq_num)));
            output_data.(append('foot_mean_before_seq',int2str(seq_num)))(participant_index, 1) = angle.foot.mean_before.(append('seq', int2str(seq_num)));
            output_data.(append('foot_std_before_seq',int2str(seq_num)))(participant_index, 1) = angle.foot.std_before.(append('seq', int2str(seq_num)));
            [output_data.(append('foot_max_angle_seq',int2str(seq_num)))(participant_index, 1), output_data.(append('foot_max_time_seq',int2str(seq_num)))(participant_index, 1)] = max(plotwindow_data(round(plotwindow_before*framerate_trackers):round((timestamp_end(seq_num)-timestamp_start(seq_num)+plotwindow_before)*framerate_trackers)));
            [output_data.(append('foot_min_angle_seq',int2str(seq_num)))(participant_index, 1), output_data.(append('foot_min_time_seq',int2str(seq_num)))(participant_index, 1)] = min(plotwindow_data(round(plotwindow_before*framerate_trackers):round((timestamp_end(seq_num)-timestamp_start(seq_num)+plotwindow_before)*framerate_trackers)));
            output_data.(append('foot_max_time_seq',int2str(seq_num)))(participant_index, 1) = (output_data.(append('foot_max_time_seq',int2str(seq_num)))(participant_index, 1))/framerate_trackers*1000;
            output_data.(append('foot_min_time_seq',int2str(seq_num)))(participant_index, 1) = (output_data.(append('foot_min_time_seq',int2str(seq_num)))(participant_index, 1))/framerate_trackers*1000;

            % Plot mean and standard deviation
            plot([0 timestamp_end(seq_num)-timestamp_start(seq_num)], [1 1]*angle.foot.mean.(append('seq', int2str(seq_num))))
            plot([0 timestamp_end(seq_num)-timestamp_start(seq_num) nan 0 timestamp_end(seq_num)-timestamp_start(seq_num)], [1 1 nan -1 -1]*angle.foot.std.(append('seq', int2str(seq_num))) + angle.foot.mean.(append('seq', int2str(seq_num))))
            plot([-calcwindow_before 0], [1 1]*angle.foot.mean_before.(append('seq', int2str(seq_num))))

            % Plot formating
            xlim([-plotwindow_before plotwindow_after])
            if sum(isnan(angle.foot.combined)-1)
                ylim([floor(min(angle.foot.combined)/5)*5 ceil(max(angle.foot.combined)/5)*5])
            end
            legend(["Ankle left" "Ankle right" "Ankle comb." "Sequence" "Mean" "Std." "Mean before"])
            title(append("Angle of ankle, sequence ", int2str(seq_num)))
        end

        %

        %% Plotting knee
        for seq_num = [1 2 3]
            %seq_num = 1;
            subplot(5,3,3+seq_num)
            hold on
            grid on

            % All angles
            plotwindow_data = angle.knee.left(round((timestamp_start(seq_num)-plotwindow_before)*framerate_trackers):round((timestamp_start(seq_num)+plotwindow_after)*framerate_trackers));
            plot((1:length(plotwindow_data))/framerate_trackers - 1, plotwindow_data)
            plotwindow_data = angle.knee.right(round((timestamp_start(seq_num)-plotwindow_before)*framerate_trackers):round((timestamp_start(seq_num)+plotwindow_after)*framerate_trackers));
            plot((1:length(plotwindow_data))/framerate_trackers - 1, plotwindow_data)
            plotwindow_data = angle.knee.combined(round((timestamp_start(seq_num)-plotwindow_before)*framerate_trackers):round((timestamp_start(seq_num)+plotwindow_after)*framerate_trackers));
            plot((1:length(plotwindow_data))/framerate_trackers - 1, plotwindow_data)

            % Timestamps
            plot([0 0 nan [1 1]*(timestamp_end(seq_num)-timestamp_start(seq_num))], [-100 100 nan -100 100])

            % Calculate mean from interval and standard deviation
            angle.knee.mean.(append('seq', int2str(seq_num))) = mean(angle.knee.combined(round(timestamp_start(seq_num)*framerate_trackers):round(timestamp_end(seq_num)*framerate_trackers)));
            angle.knee.std.(append('seq', int2str(seq_num))) = std(angle.knee.combined(round(timestamp_start(seq_num)*framerate_trackers):round(timestamp_end(seq_num)*framerate_trackers)));

            angle.knee.mean_before.(append('seq', int2str(seq_num))) = mean(angle.knee.combined(round((timestamp_start(seq_num)-calcwindow_before)*framerate_trackers):round(timestamp_start(seq_num)*framerate_trackers)));
            angle.knee.std_before.(append('seq', int2str(seq_num))) = std(angle.knee.combined(round((timestamp_start(seq_num)-calcwindow_before)*framerate_trackers):round(timestamp_start(seq_num)*framerate_trackers)));

            % Save mean and std to output
            output_data.(append('knee_mean_seq',int2str(seq_num)))(participant_index, 1) = angle.knee.mean.(append('seq', int2str(seq_num)));
            output_data.(append('knee_std_seq',int2str(seq_num)))(participant_index, 1) = angle.knee.std.(append('seq', int2str(seq_num)));
            output_data.(append('knee_mean_before_seq',int2str(seq_num)))(participant_index, 1) = angle.knee.mean_before.(append('seq', int2str(seq_num)));
            output_data.(append('knee_std_before_seq',int2str(seq_num)))(participant_index, 1) = angle.knee.std_before.(append('seq', int2str(seq_num)));
            [output_data.(append('knee_max_angle_seq',int2str(seq_num)))(participant_index, 1), output_data.(append('knee_max_time_seq',int2str(seq_num)))(participant_index, 1)] = max(plotwindow_data(round(plotwindow_before*framerate_trackers):round((timestamp_end(seq_num)-timestamp_start(seq_num)+plotwindow_before)*framerate_trackers)));
            [output_data.(append('knee_min_angle_seq',int2str(seq_num)))(participant_index, 1), output_data.(append('knee_min_time_seq',int2str(seq_num)))(participant_index, 1)] = min(plotwindow_data(round(plotwindow_before*framerate_trackers):round((timestamp_end(seq_num)-timestamp_start(seq_num)+plotwindow_before)*framerate_trackers)));
            output_data.(append('knee_max_time_seq',int2str(seq_num)))(participant_index, 1) = (output_data.(append('knee_max_time_seq',int2str(seq_num)))(participant_index, 1))/framerate_trackers*1000;
            output_data.(append('knee_min_time_seq',int2str(seq_num)))(participant_index, 1) = (output_data.(append('knee_min_time_seq',int2str(seq_num)))(participant_index, 1))/framerate_trackers*1000;

            % Plot mean and standard deviation
            plot([0 timestamp_end(seq_num)-timestamp_start(seq_num)], [1 1]*angle.knee.mean.(append('seq', int2str(seq_num))))
            plot([0 timestamp_end(seq_num)-timestamp_start(seq_num) nan 0 timestamp_end(seq_num)-timestamp_start(seq_num)], [1 1 nan -1 -1]*angle.knee.std.(append('seq', int2str(seq_num))) + angle.knee.mean.(append('seq', int2str(seq_num))))
            plot([-calcwindow_before 0], [1 1]*angle.knee.mean_before.(append('seq', int2str(seq_num))))

            % Plot formating
            xlim([-plotwindow_before plotwindow_after])
            if sum(isnan(angle.knee.combined)-1)
                ylim([floor(min(angle.knee.combined)/5)*5 ceil(max(angle.knee.combined)/5)*5])
            end
            legend(["Knee left" "Knee right" "Knee comb." "Sequence" "Mean" "Std." "Mean before"])
            title(append("Angle of knee, sequence ", int2str(seq_num)))
        end


        %

        %% Plotting hip
        for seq_num = [1 2 3]
            %seq_num = 1;
            subplot(5,3,6+seq_num)
            hold on
            grid on

            % All angles
            plotwindow_data = angle.hip.left(round((timestamp_start(seq_num)-plotwindow_before)*framerate_trackers):round((timestamp_start(seq_num)+plotwindow_after)*framerate_trackers));
            plot((1:length(plotwindow_data))/framerate_trackers - 1, plotwindow_data)
            plotwindow_data = angle.hip.right(round((timestamp_start(seq_num)-plotwindow_before)*framerate_trackers):round((timestamp_start(seq_num)+plotwindow_after)*framerate_trackers));
            plot((1:length(plotwindow_data))/framerate_trackers - 1, plotwindow_data)
            plotwindow_data = angle.hip.combined(round((timestamp_start(seq_num)-plotwindow_before)*framerate_trackers):round((timestamp_start(seq_num)+plotwindow_after)*framerate_trackers));
            plot((1:length(plotwindow_data))/framerate_trackers - 1, plotwindow_data)

            % Timestamps
            plot([0 0 nan [1 1]*(timestamp_end(seq_num)-timestamp_start(seq_num))], [-100 100 nan -100 100])

            % Calculate mean from interval and standard deviation
            angle.hip.mean.(append('seq', int2str(seq_num))) = mean(angle.hip.combined(round(timestamp_start(seq_num)*framerate_trackers):round(timestamp_end(seq_num)*framerate_trackers)));
            angle.hip.std.(append('seq', int2str(seq_num))) = std(angle.hip.combined(round(timestamp_start(seq_num)*framerate_trackers):round(timestamp_end(seq_num)*framerate_trackers)));

            angle.hip.mean_before.(append('seq', int2str(seq_num))) = mean(angle.hip.combined(round((timestamp_start(seq_num)-calcwindow_before)*framerate_trackers):round(timestamp_start(seq_num)*framerate_trackers)));
            angle.hip.std_before.(append('seq', int2str(seq_num))) = std(angle.hip.combined(round((timestamp_start(seq_num)-calcwindow_before)*framerate_trackers):round(timestamp_start(seq_num)*framerate_trackers)));

            % Save mean and std to output
            output_data.(append('hip_mean_seq',int2str(seq_num)))(participant_index, 1) = angle.hip.mean.(append('seq', int2str(seq_num)));
            output_data.(append('hip_std_seq',int2str(seq_num)))(participant_index, 1) = angle.hip.std.(append('seq', int2str(seq_num)));
            output_data.(append('hip_mean_before_seq',int2str(seq_num)))(participant_index, 1) = angle.hip.mean_before.(append('seq', int2str(seq_num)));
            output_data.(append('hip_std_before_seq',int2str(seq_num)))(participant_index, 1) = angle.hip.std_before.(append('seq', int2str(seq_num)));
            [output_data.(append('hip_max_angle_seq',int2str(seq_num)))(participant_index, 1), output_data.(append('hip_max_time_seq',int2str(seq_num)))(participant_index, 1)] = max(plotwindow_data(round(plotwindow_before*framerate_trackers):round((timestamp_end(seq_num)-timestamp_start(seq_num)+plotwindow_before)*framerate_trackers)));
            [output_data.(append('hip_min_angle_seq',int2str(seq_num)))(participant_index, 1), output_data.(append('hip_min_time_seq',int2str(seq_num)))(participant_index, 1)] = min(plotwindow_data(round(plotwindow_before*framerate_trackers):round((timestamp_end(seq_num)-timestamp_start(seq_num)+plotwindow_before)*framerate_trackers)));
            output_data.(append('hip_max_time_seq',int2str(seq_num)))(participant_index, 1) = (output_data.(append('hip_max_time_seq',int2str(seq_num)))(participant_index, 1))/framerate_trackers*1000;
            output_data.(append('hip_min_time_seq',int2str(seq_num)))(participant_index, 1) = (output_data.(append('hip_min_time_seq',int2str(seq_num)))(participant_index, 1))/framerate_trackers*1000;

            % Plot mean and standard deviation
            plot([0 timestamp_end(seq_num)-timestamp_start(seq_num)], [1 1]*angle.hip.mean.(append('seq', int2str(seq_num))))
            plot([0 timestamp_end(seq_num)-timestamp_start(seq_num) nan 0 timestamp_end(seq_num)-timestamp_start(seq_num)], [1 1 nan -1 -1]*angle.hip.std.(append('seq', int2str(seq_num))) + angle.hip.mean.(append('seq', int2str(seq_num))))
            plot([-calcwindow_before 0], [1 1]*angle.hip.mean_before.(append('seq', int2str(seq_num))))

            % Plot formating
            xlim([-plotwindow_before plotwindow_after])
            if sum(isnan(angle.hip.combined)-1)
                ylim([floor(min(angle.hip.combined)/5)*5 ceil(max(angle.hip.combined)/5)*5])
            end
            legend(["Hip left" "Hip right" "Hip comb." "Sequence" "Mean" "Std." "Mean before"])
            title(append("Angle of hip, sequence ", int2str(seq_num)))
        end   
        %

        %% Plotting back to Z
        for seq_num = [1 2 3]
            %seq_num = 1;
            subplot(5,3,9+seq_num)
            hold on
            grid on

            % All angles

            plotwindow_data = angle.back.z(round((timestamp_start(seq_num)-plotwindow_before)*framerate_trackers):round((timestamp_start(seq_num)+plotwindow_after)*framerate_trackers));
            plot((1:length(plotwindow_data))/framerate_trackers - 1, plotwindow_data)

            % Timestamps
            plot([0 0 nan [1 1]*(timestamp_end(seq_num)-timestamp_start(seq_num))], [-100 100 nan -100 100])

            % Calculate mean from interval and standard deviation
            angle.back.mean_z.(append('seq', int2str(seq_num))) = mean(angle.back.z(round(timestamp_start(seq_num)*framerate_trackers):round(timestamp_end(seq_num)*framerate_trackers)));
            angle.back.std_z.(append('seq', int2str(seq_num))) = std(angle.back.z(round(timestamp_start(seq_num)*framerate_trackers):round(timestamp_end(seq_num)*framerate_trackers)));

            angle.back.mean_z_before.(append('seq', int2str(seq_num))) = mean(angle.back.z(round((timestamp_start(seq_num)-calcwindow_before)*framerate_trackers):round(timestamp_start(seq_num)*framerate_trackers)));
            angle.back.std_z_before.(append('seq', int2str(seq_num))) = std(angle.back.z(round((timestamp_start(seq_num)-calcwindow_before)*framerate_trackers):round(timestamp_start(seq_num)*framerate_trackers)));

            % Save mean and std to output
            output_data.(append('back_mean_z_seq',int2str(seq_num)))(participant_index, 1) = angle.back.mean_z.(append('seq', int2str(seq_num)));
            output_data.(append('back_std_z_seq',int2str(seq_num)))(participant_index, 1) = angle.back.std_z.(append('seq', int2str(seq_num)));
            output_data.(append('back_mean_z_before_seq',int2str(seq_num)))(participant_index, 1) = angle.back.mean_z_before.(append('seq', int2str(seq_num)));
            output_data.(append('back_std_z_before_seq',int2str(seq_num)))(participant_index, 1) = angle.back.std_z_before.(append('seq', int2str(seq_num)));
            [output_data.(append('back_max_angle_z_seq',int2str(seq_num)))(participant_index, 1), output_data.(append('back_max_time_z_seq',int2str(seq_num)))(participant_index, 1)] = max(plotwindow_data(round(plotwindow_before*framerate_trackers):round((timestamp_end(seq_num)-timestamp_start(seq_num)+plotwindow_before)*framerate_trackers)));
            [output_data.(append('back_min_angle_z_seq',int2str(seq_num)))(participant_index, 1), output_data.(append('back_min_time_z_seq',int2str(seq_num)))(participant_index, 1)] = min(plotwindow_data(round(plotwindow_before*framerate_trackers):round((timestamp_end(seq_num)-timestamp_start(seq_num)+plotwindow_before)*framerate_trackers)));
            output_data.(append('back_max_time_z_seq',int2str(seq_num)))(participant_index, 1) = (output_data.(append('back_max_time_z_seq',int2str(seq_num)))(participant_index, 1))/framerate_trackers*1000;
            output_data.(append('back_min_time_z_seq',int2str(seq_num)))(participant_index, 1) = (output_data.(append('back_min_time_z_seq',int2str(seq_num)))(participant_index, 1))/framerate_trackers*1000;

            % Plot mean and standard deviation
            plot([0 timestamp_end(seq_num)-timestamp_start(seq_num)], [1 1]*angle.back.mean_z.(append('seq', int2str(seq_num))))
            plot([0 timestamp_end(seq_num)-timestamp_start(seq_num) nan 0 timestamp_end(seq_num)-timestamp_start(seq_num)], [1 1 nan -1 -1]*angle.back.std_z.(append('seq', int2str(seq_num))) + angle.back.mean_z.(append('seq', int2str(seq_num))))
            plot([-calcwindow_before 0], [1 1]*angle.back.mean_z_before.(append('seq', int2str(seq_num))))

            % Plot formating
            xlim([-plotwindow_before plotwindow_after])
            if sum(isnan(angle.back.z)-1)
                ylim([floor(min(angle.back.z)/5)*5 ceil(max(angle.back.z)/5)*5])
            end
            legend(["Back to Z" "Sequence" "Mean" "Std." "Mean before"])
            title(append("Angle of back to Z, sequence ", int2str(seq_num)))
        end
        %

        %% Plotting back to pelvis
        for seq_num = [1 2 3]
            %seq_num = 1;
            subplot(5,3,12+seq_num)
            hold on
            grid on

            % All angles

            plotwindow_data = angle.back.back_to_hip(round((timestamp_start(seq_num)-plotwindow_before)*framerate_trackers):round((timestamp_start(seq_num)+plotwindow_after)*framerate_trackers));
            plot((1:length(plotwindow_data))/framerate_trackers - 1, plotwindow_data)

            % Timestamps
            plot([0 0 nan [1 1]*(timestamp_end(seq_num)-timestamp_start(seq_num))], [-100 100 nan -100 100])

            % Calculate mean from interval and standard deviation
            angle.back.mean_back.(append('seq', int2str(seq_num))) = mean(angle.back.back_to_hip(round(timestamp_start(seq_num)*framerate_trackers):round(timestamp_end(seq_num)*framerate_trackers)));
            angle.back.std_back.(append('seq', int2str(seq_num))) = std(angle.back.back_to_hip(round(timestamp_start(seq_num)*framerate_trackers):round(timestamp_end(seq_num)*framerate_trackers)));

            angle.back.mean_back_before.(append('seq', int2str(seq_num))) = mean(angle.back.back_to_hip(round((timestamp_start(seq_num)-calcwindow_before)*framerate_trackers):round(timestamp_start(seq_num)*framerate_trackers)));
            angle.back.std_back_before.(append('seq', int2str(seq_num))) = std(angle.back.back_to_hip(round((timestamp_start(seq_num)-calcwindow_before)*framerate_trackers):round(timestamp_start(seq_num)*framerate_trackers)));

            % Save mean and std to output
            output_data.(append('back_mean_back_seq',int2str(seq_num)))(participant_index, 1) = angle.back.mean_back.(append('seq', int2str(seq_num)));
            output_data.(append('back_std_back_seq',int2str(seq_num)))(participant_index, 1) = angle.back.std_back.(append('seq', int2str(seq_num)));
            output_data.(append('back_mean_back_before_seq',int2str(seq_num)))(participant_index, 1) = angle.back.mean_back_before.(append('seq', int2str(seq_num)));
            output_data.(append('back_std_back_before_seq',int2str(seq_num)))(participant_index, 1) = angle.back.std_back_before.(append('seq', int2str(seq_num)));
            [output_data.(append('back_max_angle_seq',int2str(seq_num)))(participant_index, 1), output_data.(append('back_max_time_seq',int2str(seq_num)))(participant_index, 1)] = max(plotwindow_data(round(plotwindow_before*framerate_trackers):round((timestamp_end(seq_num)-timestamp_start(seq_num)+plotwindow_before)*framerate_trackers)));
            [output_data.(append('back_min_angle_seq',int2str(seq_num)))(participant_index, 1), output_data.(append('back_min_time_seq',int2str(seq_num)))(participant_index, 1)] = min(plotwindow_data(round(plotwindow_before*framerate_trackers):round((timestamp_end(seq_num)-timestamp_start(seq_num)+plotwindow_before)*framerate_trackers)));
            output_data.(append('back_max_time_seq',int2str(seq_num)))(participant_index, 1) = (output_data.(append('back_max_time_seq',int2str(seq_num)))(participant_index, 1))/framerate_trackers*1000;
            output_data.(append('back_min_time_seq',int2str(seq_num)))(participant_index, 1) = (output_data.(append('back_min_time_seq',int2str(seq_num)))(participant_index, 1))/framerate_trackers*1000;

            % Plot mean and standard deviation
            plot([0 timestamp_end(seq_num)-timestamp_start(seq_num)], [1 1]*angle.back.mean_back.(append('seq', int2str(seq_num))))
            plot([0 timestamp_end(seq_num)-timestamp_start(seq_num) nan 0 timestamp_end(seq_num)-timestamp_start(seq_num)], [1 1 nan -1 -1]*angle.back.std_back.(append('seq', int2str(seq_num))) + angle.back.mean_back.(append('seq', int2str(seq_num))))
            plot([-calcwindow_before 0], [1 1]*angle.back.mean_back_before.(append('seq', int2str(seq_num))))

            % Plot formating
            xlim([-plotwindow_before plotwindow_after])
            if sum(isnan(angle.back.back_to_hip)-1)
                ylim([floor(min(angle.back.back_to_hip)/5)*5 ceil(max(angle.back.back_to_hip)/5)*5])
            end
            legend(["Back to pelvis" "Sequence" "Mean" "Std." "Mean before"])
            title(append("Angle of back to pelvis, sequence ", int2str(seq_num)))
        end

        
    %% Plotting CCI and COP
        set(0, 'CurrentFigure', figure_cci_cop)
        clf
        %

        %% Plotting CCI
        for seq_num = [1 2 3]
            subplot(3,3,0+seq_num)
            hold on
            grid on

            % CCI data
            plotwindow_data = cci_rudolp_combined(round((timestamp_start(seq_num)-plotwindow_before)*framerate_force):round((timestamp_start(seq_num)+plotwindow_after)*framerate_force));
            plot((1:length(plotwindow_data))/framerate_force - plotwindow_before, plotwindow_data)

            % Timestamps
            plot([0 0 nan [1 1]*(timestamp_end(seq_num)-timestamp_start(seq_num))], [-100 100 nan -100 100])

            % Calculate mean from interval and standard deviation
            cci.mean.(append('seq', int2str(seq_num))) = mean(cci_rudolp_combined(round(timestamp_start(seq_num)*framerate_force):round(timestamp_end(seq_num)*framerate_force)));
            cci.std.(append('seq', int2str(seq_num))) = std(cci_rudolp_combined(round(timestamp_start(seq_num)*framerate_force):round(timestamp_end(seq_num)*framerate_force)));

            cci.mean_before.(append('seq', int2str(seq_num))) = mean(cci_rudolp_combined(round((timestamp_start(seq_num)-calcwindow_before)*framerate_force):round(timestamp_start(seq_num)*framerate_force)));
            cci.std_before.(append('seq', int2str(seq_num))) = std(cci_rudolp_combined(round((timestamp_start(seq_num)-calcwindow_before)*framerate_force):round(timestamp_start(seq_num)*framerate_force)));

            % Save mean and std to output
            output_data.(append('cci_mean_seq',int2str(seq_num)))(participant_index, 1) = cci.mean.(append('seq', int2str(seq_num)));
            output_data.(append('cci_std_seq',int2str(seq_num)))(participant_index, 1) = cci.std.(append('seq', int2str(seq_num)));
            output_data.(append('cci_mean_before_seq',int2str(seq_num)))(participant_index, 1) = cci.mean_before.(append('seq', int2str(seq_num)));
            output_data.(append('cci_std_before_seq',int2str(seq_num)))(participant_index, 1) = cci.std_before.(append('seq', int2str(seq_num)));
            [output_data.(append('cci_max_seq',int2str(seq_num)))(participant_index, 1), output_data.(append('cci_max_time_seq',int2str(seq_num)))(participant_index, 1)] = max(plotwindow_data(round(plotwindow_before*framerate_force):round((timestamp_end(seq_num)-timestamp_start(seq_num)+plotwindow_before)*framerate_force)));
            [output_data.(append('cci_min_seq',int2str(seq_num)))(participant_index, 1), output_data.(append('cci_min_time_seq',int2str(seq_num)))(participant_index, 1)] = min(plotwindow_data(round(plotwindow_before*framerate_force):round((timestamp_end(seq_num)-timestamp_start(seq_num)+plotwindow_before)*framerate_force)));
            output_data.(append('cci_max_time_seq',int2str(seq_num)))(participant_index, 1) = (output_data.(append('cci_max_time_seq',int2str(seq_num)))(participant_index, 1))/framerate_force*1000;
            output_data.(append('cci_min_time_seq',int2str(seq_num)))(participant_index, 1) = (output_data.(append('cci_min_time_seq',int2str(seq_num)))(participant_index, 1))/framerate_force*1000;

            % Plot mean and standard deviation
            plot([0 timestamp_end(seq_num)-timestamp_start(seq_num)], [1 1]*cci.mean.(append('seq', int2str(seq_num))))
            plot([0 timestamp_end(seq_num)-timestamp_start(seq_num) nan 0 timestamp_end(seq_num)-timestamp_start(seq_num)], [1 1 nan -1 -1]*cci.std.(append('seq', int2str(seq_num))) + cci.mean.(append('seq', int2str(seq_num))))
            plot([-calcwindow_before 0], [1 1]*cci.mean_before.(append('seq', int2str(seq_num))))

            % Plot formating
            xlim([-plotwindow_before plotwindow_after])
            ylim([0 2])
            legend(["CCI_{Rudolp}" "Sequence" "Mean" "Std." "Mean before"])
            title(append("CCI, sequence ", int2str(seq_num)))
        end 
        %

        %% Plotting COP X
    %     plot((1:length(cop.sequence(1,:)))/3000, (cop.sequence(1,:)-cop.static(1)))
    %     plot((1:length(cop.sequence(1,:)))/3000, (cop.sequence(2,:)-cop.static(2)))

        for seq_num = [1 2 3]
            subplot(3,3,3+seq_num)
            hold on
            grid on

            % COP data
            plotwindow_data = cop.sequence(1, round((timestamp_start(seq_num)-plotwindow_before)*framerate_force):round((timestamp_start(seq_num)+plotwindow_after)*framerate_force));
            plot((1:length(plotwindow_data))/framerate_force - 1, plotwindow_data)

            % Timestamps
            plot([0 0 nan [1 1]*(timestamp_end(seq_num)-timestamp_start(seq_num))], [-1000 1000 nan -1000 1000])

            % Calculate mean from interval and standard deviation
            cop.mean_x.(append('seq', int2str(seq_num))) = mean(cop.sequence(1, round(timestamp_start(seq_num)*framerate_force):round(timestamp_end(seq_num)*framerate_force)));
            cop.std_x.(append('seq', int2str(seq_num))) = std(cop.sequence(1, round(timestamp_start(seq_num)*framerate_force):round(timestamp_end(seq_num)*framerate_force)));

            cop.mean_x_before.(append('seq', int2str(seq_num))) = mean(cop.sequence(1, round((timestamp_start(seq_num)-calcwindow_before)*framerate_force):round(timestamp_start(seq_num)*framerate_force)));
            cop.std_x_before.(append('seq', int2str(seq_num))) = std(cop.sequence(1, round((timestamp_start(seq_num)-calcwindow_before)*framerate_force):round(timestamp_start(seq_num)*framerate_force)));

            % Save mean and std to output
            output_data.(append('cop_mean_side_seq',int2str(seq_num)))(participant_index, 1) = cop.mean_x.(append('seq', int2str(seq_num)));
            output_data.(append('cop_std_side_seq',int2str(seq_num)))(participant_index, 1) = cop.std_x.(append('seq', int2str(seq_num)));
            output_data.(append('cop_mean_side_before_seq',int2str(seq_num)))(participant_index, 1) = cop.mean_x_before.(append('seq', int2str(seq_num)));
            output_data.(append('cop_std_side_before_seq',int2str(seq_num)))(participant_index, 1) = cop.std_x_before.(append('seq', int2str(seq_num)));
            [output_data.(append('cop_max_side_seq',int2str(seq_num)))(participant_index, 1), output_data.(append('cop_max_time_x_seq',int2str(seq_num)))(participant_index, 1)] = max(plotwindow_data(round(plotwindow_before*framerate_force):round((timestamp_end(seq_num)-timestamp_start(seq_num)+plotwindow_before)*framerate_force)));
            [output_data.(append('cop_min_side_seq',int2str(seq_num)))(participant_index, 1), output_data.(append('cop_min_time_x_seq',int2str(seq_num)))(participant_index, 1)] = min(plotwindow_data(round(plotwindow_before*framerate_force):round((timestamp_end(seq_num)-timestamp_start(seq_num)+plotwindow_before)*framerate_force)));
            output_data.(append('cop_max_time_side_seq',int2str(seq_num)))(participant_index, 1) = (output_data.(append('cop_max_time_x_seq',int2str(seq_num)))(participant_index, 1))/framerate_force*1000;
            output_data.(append('cop_min_time_side_seq',int2str(seq_num)))(participant_index, 1) = (output_data.(append('cop_min_time_x_seq',int2str(seq_num)))(participant_index, 1))/framerate_force*1000;

            % Plot mean and standard deviation
            plot([0 timestamp_end(seq_num)-timestamp_start(seq_num)], [1 1]*cop.mean_x.(append('seq', int2str(seq_num))))
            plot([0 timestamp_end(seq_num)-timestamp_start(seq_num) nan 0 timestamp_end(seq_num)-timestamp_start(seq_num)], [1 1 nan -1 -1]*cop.std_x.(append('seq', int2str(seq_num))) + cop.mean_x.(append('seq', int2str(seq_num))))
            plot([-calcwindow_before 0], [1 1]*cop.mean_x_before.(append('seq', int2str(seq_num))))

            % Plot formating
            xlim([-plotwindow_before plotwindow_after])
            ylim([floor(min(cop.sequence(1,:))/5)*5 ceil(max(cop.sequence(1,:))/5)*5])
            legend(["COP" "Sequence" "Mean" "Std."  "Mean before"])
            title(append("COP side, sequence ", int2str(seq_num)))
        end 
        %

        %% Plotting COP Y
        for seq_num = [1 2 3]
            subplot(3,3,6+seq_num)
            hold on
            grid on

            % COP data
            plotwindow_data = cop.sequence(2, round((timestamp_start(seq_num)-plotwindow_before)*framerate_force):round((timestamp_start(seq_num)+plotwindow_after)*framerate_force));
            plot((1:length(plotwindow_data))/framerate_force - 1, plotwindow_data)

            % Timestamps
            plot([0 0 nan [1 1]*(timestamp_end(seq_num)-timestamp_start(seq_num))], [-1000 1000 nan -1000 1000])

            % Calculate mean from interval and standard deviation
            cop.mean_y.(append('seq', int2str(seq_num))) = mean(cop.sequence(2, round(timestamp_start(seq_num)*framerate_force):round(timestamp_end(seq_num)*framerate_force)));
            cop.std_y.(append('seq', int2str(seq_num))) = std(cop.sequence(2, round(timestamp_start(seq_num)*framerate_force):round(timestamp_end(seq_num)*framerate_force)));

            cop.mean_y_before.(append('seq', int2str(seq_num))) = mean(cop.sequence(2, round((timestamp_start(seq_num)-calcwindow_before)*framerate_force):round(timestamp_start(seq_num)*framerate_force)));
            cop.std_y_before.(append('seq', int2str(seq_num))) = std(cop.sequence(2, round((timestamp_start(seq_num)-calcwindow_before)*framerate_force):round(timestamp_start(seq_num)*framerate_force)));

            % Save mean and std to output
            output_data.(append('cop_mean_posant_seq',int2str(seq_num)))(participant_index, 1) = cop.mean_y.(append('seq', int2str(seq_num)));
            output_data.(append('cop_std_posant_seq',int2str(seq_num)))(participant_index, 1) = cop.std_y.(append('seq', int2str(seq_num)));
            output_data.(append('cop_mean_posant_before_seq',int2str(seq_num)))(participant_index, 1) = cop.mean_y_before.(append('seq', int2str(seq_num)));
            output_data.(append('cop_std_posant_before_seq',int2str(seq_num)))(participant_index, 1) = cop.std_y_before.(append('seq', int2str(seq_num)));
            [output_data.(append('cop_max_posant_seq',int2str(seq_num)))(participant_index, 1), output_data.(append('cop_max_time_y_seq',int2str(seq_num)))(participant_index, 1)] = max(plotwindow_data(round(plotwindow_before*framerate_force):round((timestamp_end(seq_num)-timestamp_start(seq_num)+plotwindow_before)*framerate_force)));
            [output_data.(append('cop_min_posant_seq',int2str(seq_num)))(participant_index, 1), output_data.(append('cop_min_time_y_seq',int2str(seq_num)))(participant_index, 1)] = min(plotwindow_data(round(plotwindow_before*framerate_force):round((timestamp_end(seq_num)-timestamp_start(seq_num)+plotwindow_before)*framerate_force)));
            output_data.(append('cop_max_time_posant_seq',int2str(seq_num)))(participant_index, 1) = (output_data.(append('cop_max_time_y_seq',int2str(seq_num)))(participant_index, 1))/framerate_force*1000;
            output_data.(append('cop_min_time_posant_seq',int2str(seq_num)))(participant_index, 1) = (output_data.(append('cop_min_time_y_seq',int2str(seq_num)))(participant_index, 1))/framerate_force*1000;

            % Plot mean and standard deviation
            plot([0 timestamp_end(seq_num)-timestamp_start(seq_num)], [1 1]*cop.mean_y.(append('seq', int2str(seq_num))))
            plot([0 timestamp_end(seq_num)-timestamp_start(seq_num) nan 0 timestamp_end(seq_num)-timestamp_start(seq_num)], [1 1 nan -1 -1]*cop.std_y.(append('seq', int2str(seq_num))) + cop.mean_y.(append('seq', int2str(seq_num))))
            plot([-calcwindow_before 0], [1 1]*cop.mean_y_before.(append('seq', int2str(seq_num))))

            % Plot formating
            xlim([-plotwindow_before plotwindow_after])
            ylim([floor(min(cop.sequence(2,:))/5)*5 ceil(max(cop.sequence(2,:))/5)*5])
            legend(["COP pos-ant" "Sequence" "Mean" "Std."  "Mean before"])
            title(append("COP pos-ant, sequence ", int2str(seq_num)))
        end

        
        %% Save plots as png
        saveas(figure_angles, append('Rorelselabb/Sequence6DOF_Output/',int2str(participant), '_Angles', '.png'));
        saveas(figure_cci_cop, append('Rorelselabb/Sequence6DOF_Output/',int2str(participant), '_CCI_COP', '.png'));
    
        %
        %% Other plot

        %plot((1:length(vector1))/200, platform(1,:)'+40)
        %plot((1:length(vector1))/200, platform(2,:)'/50-13)

fprintf(append(int2str(participant_index), ' of ', int2str(length(participants_in_process)), ' done\n') )
end    
    
    
    
%% Export output to excel

% Export to excell
    if exist('angles_table')
        clear angles_table
    end
    angles_table = struct2table(output_data); %table(participant, t_start_to_onset_seq1, cci_left_feed_forward_100ms_seq1, cci_right_feed_forward_100ms_seq1, cci_combined_feed_forward_100ms_seq1, cci_left_feed_forward_200ms_seq1, cci_right_feed_forward_200ms_seq1, cci_combined_feed_forward_200ms_seq1, t_start_to_onset_seq2, cci_left_feed_forward_100ms_seq2, cci_right_feed_forward_100ms_seq2, cci_combined_feed_forward_100ms_seq2, cci_left_feed_forward_200ms_seq2, cci_right_feed_forward_200ms_seq2, cci_combined_feed_forward_200ms_seq2, t_start_to_onset_seq3, cci_left_feed_forward_100ms_seq3, cci_right_feed_forward_100ms_seq3, cci_combined_feed_forward_100ms_seq3, cci_left_feed_forward_200ms_seq3, cci_right_feed_forward_200ms_seq3, cci_combined_feed_forward_200ms_seq3);
    %angles_table = rows2vars(angles_table);
    if isfile('Rorelselabb/angles_table.xlsx')
        delete('Rorelselabb/angles_table.xlsx');
    end
    writetable(angles_table,'Rorelselabb/angles_table.xlsx');
    %
    
    
%% Old notes
    

%     labels = qtm_data.Trajectories.Labeled.Labels;
% 
%     vector1 = get_vector(qtm_data, 'LTAM', 'LTTC');
%     vector2 = get_vector(qtm_data, 'LTAM', 'LFM1');
% 
%     angle = find_angle(vector1, vector2);
% 
%     plot(angle)
%     ylabel('Vinkel - fot mot ben')
%     hold on
%     
%     vector1 = get_vector(qtm_data, 'RTAM', 'RTTC');
%     vector2 = get_vector(qtm_data, 'RTAM', 'RFM1');
%     angle = find_angle(vector1, vector2);
%     plot(angle)
% 
%     index_platform = find(ismember(labels, 'platform - 2'));
%     platform = qtm_data.Trajectories.Labeled.Data(index_platform, :, :);
%     platform = squeeze(platform);
%     platform = platform(1:3,:);
% 
%     subplot(plotgrid1(1),plotgrid1(2),2);
%     plot(platform(2,:)/10 + 20);
%     ylabel('Plattformsrörelse i y-led')





%     sqrt(force(1,:).^2+force(2,:).^2+force(3,:).^2);
%     vectorlength = size(force2_mag);
%     subplot(plotgrid1(1),plotgrid1(2),3);
%     plot((1:vectorlength(2))/15, force2_mag/10)
%     ylabel('Force magnitude');