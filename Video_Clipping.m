clear;clc;close all;

data_filename = 'Z:\M1_Stroke\I127\Day9\Results\D9_PreSleep_GUI.mat';
video_foldername = 'Z:\M1_Stroke\I127\Day9\';
camera_num = 2;

%Each Row is a different time. Each column is a different camera.
video_filename_prefixes = [{'I127S9R1-(2022_12_9)-(9h22m)-Cam1-'}, {'I127S9R1-(2022_12_9)-(9h22m)-Cam2-'}]; 

load(data_filename);

if ~exist(data_filename(1:end-4), 'dir')
    mkdir(data_filename(1:end-4));
end
new_vid_foldername = [data_filename(1:end-4), '\'];

for camera = 1:camera_num
    data_idx = 0;
    for prefix = 1:size(video_filename_prefixes,1)
        for vid_num = 1:(size(data,1) - data_idx)
            vid_filename = [video_foldername, video_filename_prefixes{prefix, camera}, num2str(vid_num), '.avi'];
            if ~exist(vid_filename, 'file')
                break 
            end
            data_idx = data_idx+1;
            if isempty(data{data_idx,2})
                continue
            end
            
            outcome = str2double(data{data_idx,3});
            framestr  = data{data_idx,2};
            framenum  = str2num(char(framestr));
            if outcome == 3
                continue
            end
            if outcome == 0 || outcome == 4
                if framenum(end) > 1
                    continue
                end
            end
            
            if outcome == 0 || outcome == 4
                reach = framenum(end - 4);
                retract = framenum(end - 1);
            else
                reach = framenum(end - 2);
                retract = framenum(end);
            end
                
            video = VideoReader(vid_filename);
            frames = cell(video.FrameRate * video.Duration,1);
            
            new_vid_filename = [new_vid_foldername, video_filename_prefixes{prefix, camera}, num2str(vid_num), '.avi'];
            vid_out = VideoWriter(new_vid_filename);
            open(vid_out)
            
            idx = 1;
            while idx <= retract
                frame = readFrame(video);
                if idx >= reach
                    writeVideo(vid_out,frame);
                end
                idx = idx + 1;
            end
            
            close(vid_out);
        end
    end
end

beep
disp 'Done.'
    