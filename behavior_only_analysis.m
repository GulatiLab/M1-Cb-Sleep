clear;clc;close all;

expected_blocks = 2;
sub_folder = ''; %'\Reach_Vids';
days_path = 'Z:\M1_Stroke\I147\PreStroke';

day_names = dir([days_path, '\Day*']);
for day_i = 1:length(day_names)
    day_path = [days_path, '\', day_names(day_i).name, sub_folder];
    disp(day_names(day_i).name);
    if ~exist([day_path, '\Results'],'dir')
        mkdir([day_path, '\Results']);
    end
    
    data_names = dir([day_path, '\*man*.mat']);
    all_time_strs = '';
    all_times = nan(0,1);
    for data_i = 1:length(data_names)
        time_idx = regexp(data_names(data_i).name,'\d\d?h\d?\dm');
        time_str = data_names(data_i).name(time_idx:time_idx+4);
        time = time_str;
        if time(5) == ')'
            time = ['0', time(1), '0', time(3)];
        elseif time(5) == 'm'
            if time(2) == 'h'
                time = ['0', time(1), time(3:4)];
            else
                time = [time(1:2), '0', time(4)];
            end
        else
            time = [time(1:2), time(4:5)];
        end
        time = str2double(time);
        new_time = true;
        for time_i = 1:size(all_time_strs,1)
            if time == all_times(time_i)
                new_time = false;
                continue
            end
        end
        if new_time
            all_time_strs = cat(1, all_time_strs, time_str);
            all_times = cat(1, all_times, time);
        end
    end
    
    if size(all_time_strs,1) ~= expected_blocks
        [all_times,sort_idx] = sort(all_times);
        all_time_strs = all_time_strs(sort_idx, :);
        all_times_new = zeros(size(all_times));
        for time_i = 1:size(all_time_strs,1)
            all_times_new(time_i) = input(['What block should ', all_time_strs(time_i,:), ' be part of? (Enter 0 to discard)', newline]);
            if all_times_new(time_i) > expected_blocks || all_times_new(time_i) < 0
                error('Group outside expected bounds.')
            end
            all_times_new(time_i) = all_times_new(time_i) * 10000;
            all_times_new(time_i) = all_times_new(time_i) + all_times(time_i);
        end
        [~,sort_idx] = sort(all_times(all_times_new < 10000));
        [~,sort_idx_new] = sort(all_times_new(all_times_new < 10000));
        if sum(sort_idx ~= sort_idx_new) > 0
            warning('Times are not grouped in order.')
            response = input('Proceed anyway? y/n');
            if strcmp(conf, 'y') 
                
            elseif strcmp(conf, 'n')
                disp('Aborting.')
                return
            else
                disp('Unrecognised responce. Aborting.')
                return
            end
        end
        group = all_times_new/10000;
        group = floor(group);
        
    else
        [all_times,sort_idx] = sort(all_times);
        all_time_strs = all_time_strs(sort_idx, :);
        group = 1:expected_blocks;
    end
    
    day_succ_count = 0;
    day_reach_count = 0;
    for cur_group = 1:expected_blocks
        if cur_group == 1
            out_filename = ['D', day_names(day_i).name(4:end), '_PreSleep_GUI.mat'];
        elseif cur_group == 2
            out_filename = ['D', day_names(day_i).name(4:end), '_PostSleep_GUI.mat'];
        else
            out_filename = ['D', day_names(day_i).name(4:end), '_Block', num2str(cur_group), '_GUI.mat'];
        end
        comb_data = cell(0,4);
        comb_idx = 0;
        for time_i = 1:size(all_time_strs,1)
            if cur_group == group(time_i)
                comb_idx = size(comb_data,1);
                time_filename = dir([day_path, '\*', all_time_strs(time_i,:), '*.mat']);
                for data_i = 1:size(time_filename,1)
                    load([day_path, '\', time_filename(data_i).name])
                    for row_i = 1:size(data,1)
                        comb_data(str2double(data{row_i})+comb_idx,:) = data(row_i,:);
                    end
                end
            end
        end
        data = comb_data;
        save([day_path, '\Results\', out_filename],'data');
        
        outcomes = cellfun(@str2double,data(:,3));
        data = data(outcomes < 2,:);
        outcomes = outcomes(outcomes < 2,:);
        
        day_succ_count = day_succ_count + sum(outcomes);
        day_reach_count = day_reach_count + length(outcomes);
    end
    
    fid = fopen([day_path, '\Results\success_rate.txt'],'w');
    fprintf(fid,'%s',[num2str(day_succ_count), newline]);
    fprintf(fid,'%s',num2str(day_reach_count), newline);
    fprintf(fid,'%s',num2str(day_succ_count/day_reach_count));
    fclose(fid);
end


% animals = ['I127'];      %['I111';'I112'];  %['I027';'I031';'I033';'I039';'I073';'I075'];
% days = 10;
% blocks = 1;%2;
% camera_framerate = 87;
% 
% for a = 1:size(animals,1)
%     animal = animals(a,:);
%     shared_data.GUI_data = cell(days, blocks);
%     day_succ_rate = zeros(2,days);
%     for day = 1:days
%         load(['D:\Pierson_Data_Analysis\behavior_only_animals\', animal, '\Day', num2str(day), '\Results\D', num2str(day), '_PreSleep_GUI.mat'])
%         shared_data.GUI_data{day,1} = data;
%         if blocks > 1
%             load(['D:\Pierson_Data_Analysis\behavior_only_animals\', animal, '\Day', num2str(day), '\PostSleepGUI.mat'])
%             shared_data.GUI_data{day,2} = data;
%         end
% 
%         outcomes = cellfun(@str2double,data(:,3));
%         data = data(outcomes < 2,:);
%         outcomes = outcomes(outcomes < 2,:);
%         
%         day_succ_rate(1,day) = sum(outcomes);
%         day_succ_rate(2,day) = length(outcomes);
%         
%     end
%     shared_data.day_success_rate = day_succ_rate;
%     
%     %TODO
%     %shared_data.x_trajectory_consistency = 
%     %shared_data.y_trajectory_consistency = 
%     
%     save(['D:\Pierson_Data_Analysis\behavior_only_animals\',animal,'\Shared_Data.mat'], 'shared_data')
%     param.blocks = blocks;
%     param.Camera_framerate = camera_framerate;
%     save(['D:\Pierson_Data_Analysis\behavior_only_animals\',animal,'\Parameters.mat'], 'param')
% end 

beep
disp 'Analysis Complete.'