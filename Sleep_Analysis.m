clear;clc;close all;

%I089 has no sleep blocks on days 4 or 5. I110 have no sleep at all
animal = 'I122';

% HPC/Unix compatability
if isunix
    cd('/common/fleischerp/');
    origin_rootpath = '/common/fleischerp/raw_data/';
    rootpath = '/common/fleischerp/data/';
else
    restoredefaultpath()
    if strcmp(animal, 'I061') || strcmp(animal, 'I064') || strcmp(animal, 'I076') || strcmp(animal, 'I096')
        origin_rootpath = 'X:/LabDataBackup/M1_Cb_Reach/';
    else
        origin_rootpath = 'Z:/M1_Cb_Reach/';
    end
    rootpath = 'D:/Pierson_Data_Analysis/';
end

% Create parameter file or load existing file
if exist([rootpath,animal,'/Parameters.mat'],'file')
    load([rootpath,animal,'/Parameters.mat'])
else
    param = struct();
end
    
% fill in initial parameters
if strcmp(animal, 'I060')
    %I060 Settings
    param.training_block_names = {'TDT_data/I050-200316-135106' 'TDT_data/I050-200316-162539';
                                      'TDT_data/I050-200317-113531' 'TDT_data/I050-200317-144149';
                                      'TDT_data/I050-200318-115517' 'TDT_data/I050-200318-144712';
                                      'TDT_data/I050-200319-112351' 'TDT_data/I050-200319-141221';
                                      'TDT_data/I050-200320-112028' 'TDT_data/I050-200320-140715'};

    param.sleep_block_names = {'I050-200316-115829' 'I050-200316-145422';
                                      'I050-200317-092828' 'I050-200317-122703';
                                      'I050-200318-095036' 'I050-200318-123956';
                                      'I050-200319-092129' 'I050-200319-121123';
                                      'I050-200320-092309' 'I050-200320-120325'};
    
    param.sleep_block_vid_names = {'I060T1S1-(2020_3_16)-(11h59m)-SPONT' 'I060T1S2-(2020_3_16)-(14h54m)-SPONT';
        'I060D2S1-(2020_3_17)-(9h29m)-SPONT'  'I060D2S2-(2020_3_17)-(12h27m)-SPONT';
        'I060D3S1-(2020_3_18)-(9h51m)-SPONT'  'I060D3S2-(2020_3_18)-(12h40m)-SPONT';
        'I060D4S1-(2020_3_19)-(9h21m)-SPONT'  'I060D4S2-(2020_3_19)-(12h11m)-SPONT';
        'I060D5S1-(2020_3_20)-(9h23m)-SPONT'  'I060D5S2-(2020_3_20)-(12h3m)-SPONT'};
    
    param.durFiles = ['Durations_200316';'Durations_200317';'Durations_200318';'Durations_200319';'Durations_200320'];
    param.WAV_video_offset = 0.25;
    param.days = 5;
    
elseif strcmp(animal, 'I061')
    %I061 Settings
    param.training_block_names = {['TDT_data/I061-200511-135832';'TDT_data/I061-200511-142659'] ['TDT_data/I061-200511-162529'];
                                      ['TDT_data/I061-200512-122518'] ['TDT_data/I061-200512-151526'];
                                      ['TDT_data/I061-200513-111020'] ['TDT_data/I061-200513-135834'];
                                      ['TDT_data/I061-200514-120137'] ['TDT_data/I061-200514-145444'];
                                      ['TDT_data/I061-200515-113618'] ['TDT_data/I061-200515-142605']};
                    
    param.sleep_block_names = {['TDT_data/I061-200511-105824';'TDT_data/I061-200511-123320'] ['TDT_data/I061-200511-150921'];
                                      ['TDT_data/I061-200512-102347'] ['TDT_data/I061-200512-131257'];
                                      ['TDT_data/I061-200513-090705'] ['TDT_data/I061-200513-115747'];
                                      ['TDT_data/I061-200514-094708'] ['TDT_data/I061-200514-125245'];
                                      ['TDT_data/I061-200515-093406'] ['TDT_data/I061-200515-122415']};
    
    param.durFiles = ['Durations_200511';'Durations_200512';'Durations_200513';'Durations_200514';'Durations_200515'];
    param.WAV_video_offset = 0.25;
    param.days = 5;
    param.LFP_path = '.';
    
elseif strcmp(animal, 'I064')
    param.training_block_names = {['I064-200713-122227'] ['I064-200713-150141'];
                                  ['I064-200714-112650'] ['I064-200714-142951'];
                                  ['I064-200715-124653'] ['I064-200715-151426'];
                                  ['I064-200716-111934'] ['I064-200716-143134'];
                                  ['I064-200717-105528'] ['I064-200717-135405']};
    
    param.sleep_block_names = {['I064-200713-104904'] ['I064-200713-132826'];
                               ['I064-200714-092202'] ['I064-200714-122631'];
                               ['I064-200715-111430'] ['I064-200715-134154'];
                               ['I064-200716-091742'] ['I064-200716-122853'];
                               ['I064-200717-085341'] ['I064-200717-115259']};

    param.durFiles = ['Durations_200713';'Durations_200714';'Durations_200715';'Durations_200716';'Durations_200717'];
    param.WAV_video_offset = 0.25;
    param.LFP_path = 'RS4_Data';
    param.Wave_path = 'Wave_Channels';
    param.days = 5;
    
elseif strcmp(animal, 'I076')
    param.training_block_names = {['I076-201207-122019'] ['I076-201207-155657'];
                                  ['I076-201208-113010'] ['I076-201208-152911'];
                                  ['I076-201209-115653'] ['I076-201209-150030'];
                                  ['I076-201210-114806'] ['I076-201210-143519'];
                                  ['I076-201211-113827'] ['I076-201211-142237'];
                                  ['I076-201212-111646'] []};
                    
    param.sleep_block_names = {['I076-201207-100346'] ['I076-201207-135540'];
                               ['I076-201208-092637'] ['I076-201208-125538'];
                               ['I076-201209-095419'] ['I076-201209-125751'];
                               ['I076-201210-094626'] ['I076-201210-123308'];
                               ['I076-201211-093156'] ['I076-201211-121942'];
                               [] []};
        
    param.durFiles = ['Durations_201207';'Durations_201208';'Durations_201209';'Durations_201210';'Durations_201211'];
    param.WAV_video_offset = 0.25;
    param.LFP_path = 'Wave_Channels';
    param.Wave_path = 'Wave_Channels';
    param.Spike_path = 'RS4_Data';
    param.days = 5;
    
elseif strcmp(animal, 'I086')
    %I086 Settings
    param.training_block_names = {['I086-210517-125705'] ['I086-210517-155945'];
                                  ['I086-210518-115542'] ['I086-210518-145742'];
                                  ['I086-210519-105747'] ['I086-210519-141210';'I086-210519-144001'];
                                  ['I086-210520-110050'] ['I086-210520-135125'];
                                  ['I086-210521-114421'] ['I086-210521-143545']};
                    
    param.sleep_block_names = {['I086-210517-095239';'I086-210517-115522'] ['I086-210517-143132'];
                               ['I086-210518-094149'] ['I086-210518-125139'];
                               ['I086-210519-085503'] ['I086-210519-120941'];
                               ['I086-210520-094043'] ['I086-210520-115044'];
                               ['I086-210521-094113'] ['I086-210521-123217']};
    param.durFiles = ['Durations_210517';'Durations_210518';'Durations_210519';'Durations_210520';'Durations_210521'];
    param.WAV_video_offset = 0.25;
    param.days = 5;
    param.LFP_path = 'Wave_Channels';
    
elseif strcmp(animal, 'I096')
    %I096 Settings
    param.training_block_names = {['I096-211025-122116'] ['I096-211025-153916'];
                                  ['I096-211026-120626'] ['I096-211026-153106'];
                                  ['I096-211027-120654'] ['I096-211027-151750'];
                                  ['I096-211028-115354'] ['I096-211028-153632'];
                                  ['I096-211029-112213'] ['I096-211029-142122']};
    
    param.sleep_block_names = {['I096-211025-102012'] ['I096-211025-133733'];
                               ['I096-211026-100034'] ['I096-211026-132815'];
                               ['I096-211027-100439'] ['I096-211027-131521'];
                               ['I096-211028-095201'] ['I096-211028-132627'];
                               ['I096-211029-095822'] ['I096-211029-121913']};
    
    param.durFiles = ['Durations_211025';'Durations_211026';'Durations_211027';'Durations_211028';'Durations_211029'];
    param.WAV_video_offset = 0.25;
    param.LFP_path = 'Wave_Channels';
    
elseif strcmp(animal, 'I107')
    %I086 Settings
    param.training_block_names = {['I107-211213-130247'] ['I107-211213-153419'];
                                  ['I107-211214-112244'] ['I107-211214-134922'];
                                  ['I107-211215-111931'] ['I107-211215-133736'];
                                  ['I107-211216-112253'] ['I107-211216-132418'];
                                  ['I107-211217-113515'] ['I107-211217-132608']};
                    
       param.sleep_block_names = {['I107-211213-102326'; 'I107-211213-112822'] ['I107-211213-140250'];
                                  ['I107-211214-095126'] ['I107-211214-121752'];
                                  ['I107-211215-094711'] ['I107-211215-120600'];
                                  ['I107-211216-101622'] ['I107-211216-122235'];
                                  ['I107-211217-100739'; 'I107-211217-110845'] ['I107-211217-122433']};
    param.durFiles = ['Durations_211213';'Durations_211214';'Durations_211215';'Durations_211216';'Durations_211217'];
    param.WAV_video_offset = 0.25;
    param.days = 5;
    param.LFP_path = 'Wave_Channels';
    
elseif strcmp(animal, 'I122')
    %I086 Settings
    param.training_block_names = {['I122-220801-121250'] ['I122-220801-155534'];
                                  ['I122-220802-111151'] ['I122-220802-143241'];
                                  ['I122-220803-112007'] ['I122-220803-143210'];
                                  ['I122-220804-110810'] ['I122-220804-142609'];
                                  ['I122-220805-111043'] ['I122-220805-142332']};
                    
    param.sleep_block_names = {['I122-220801-094046'] ['I122-220801-135342'];
                               ['I122-220802-091025'] ['I122-220802-123043'; 'I122-220802-133300'];
                               ['I122-220803-091853'] ['I122-220803-123008'];
                               ['I122-220804-090709'] ['I122-220804-122036'];
                               ['I122-220805-090917'] ['I122-220805-122214']};
                                  
    param.durFiles = ['Durations_220801';'Durations_220802';'Durations_220803';'Durations_220804';'Durations_220805'];
    param.WAV_video_offset = 0;
    param.LFP_path = 'Wave_Channels';
    param.Wave_path = 'Wave_Channels';
    param.Spike_path = 'RS4_Data';
    param.Camera_framerate = 303; %In hz
    param.dom_hand = 'right';
    param.M1_shant_site_num = 1;
    param.Cb_shant_site_num = 16;
    param.M1_neurons = 19;
    param.Cb_neurons = 16;
    param.M1_spike_wave_Fs = 24414;
    param.Cb_spike_wave_Fs = 24414;
    
else
    error('Unrecognized animal.')
end
    
param.s_block_names = {'Sleep1' 'Sleep2'};
param.sleep_blocks = length(param.s_block_names);
save([rootpath,animal,'/Parameters.mat'], 'param');


param.epoch_length = 6; %In seconds
          %TD: TODO, MF: memory failure, L: long, RM: requires one (or more) of the main analyses 
          % A: Analysis method 1 (do not use), B: Analysis method 2 (3rd-party code)(use this), 2: Analyses that need something specified (a trial, a channel, etc.).
          %     A            B        TD A                         MF?       2  2 TD MF  L    TD  L  L           L     2     2              
          %1 2 [3 4 5 6 7 8][9 10 11] 12 [13 14 15 16 17 18 19] 20 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48
enabled = [0 0  0 0 0 0 0 0  0  0  0   0   0  0  0  0  0  0  0   0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  0  1  0  0  0];

%% Truncate Sleep (1)

if enabled(1)
    disp('Block 1...')
    for day = 1:param.days
        load([origin_rootpath,animal,'/',param.durFiles(day,:),'.mat']);
        dM_idx = 0;
        for block = 1:param.sleep_blocks
            if ~exist([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block}],'dir')
                mkdir([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block}]);
            end
            M1_LFP_t = [];
            Cb_LFP_t = [];
            WAVE1_t = [];
            WAVE2_t = [];
            sbn = param.sleep_block_names{day, block};
            for sub_block = 1:size(sbn,1)
                dM_idx = dM_idx + 1;
                
                load([origin_rootpath,animal,'/',param.LFP_path,'/',sbn(sub_block,:),'/LFP_M1.mat']);
                param.M1_Fs = Fs_lfp;
                load([origin_rootpath,animal,'/',param.LFP_path,'/',sbn(sub_block,:),'/LFP_Cb.mat']);
                param.Cb_Fs = Fs_lfp;
                load([origin_rootpath,animal,'/',param.Wave_path,'/',sbn(sub_block,:),'/WAV.mat']);
                param.Wave_Fs = Fs_wave;
                
                M1_truncation_samp = floor(durMAT(5,dM_idx)*param.M1_Fs);
                Cb_truncation_samp = floor(durMAT(5,dM_idx)*param.Cb_Fs);
                Wave_truncation_samp = floor(durMAT(5,dM_idx)*param.Wave_Fs);
                
                M1_LFP_t = cat(2, M1_LFP_t, LFPs1(:,1:M1_truncation_samp));
                Cb_LFP_t = cat(2, Cb_LFP_t, LFPs2(:,1:Cb_truncation_samp));
                WAVE1_t = cat(2, WAVE1_t, Wave1(1:Wave_truncation_samp));
                WAVE2_t = cat(2, WAVE2_t, Wave2(1:Wave_truncation_samp));
            end
            
            tbn = param.training_block_names{day, block};
            for sub_block = 1:size(tbn,1)
                dM_idx = dM_idx + 1;
            end
            
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/LFP_M1_Truncated.mat'],'M1_LFP_t','-v7.3');
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/LFP_Cb_Truncated.mat'],'Cb_LFP_t','-v7.3');
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/WAV_Truncated.mat'],'WAVE1_t','WAVE2_t');
            save([rootpath,animal,'/Parameters.mat'],'param');
        end
    end
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Select Bad Channels (2)

if enabled(2)
    disp('Block 2...')
    all_M1_sleep_bad_chans = [];
    all_Cb_sleep_bad_chans = [];
    use_existing_bad_chans = true;
    
    if use_existing_bad_chans
        param.M1_sleep_bad_chans = param.M1_bad_chans;
        param.Cb_sleep_bad_chans = param.Cb_bad_chans;
        
    else
        for day = 1:param.days
            for block = 1:param.sleep_blocks
                load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/LFP_M1_Truncated.mat']);
                load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/LFP_Cb_Truncated.mat']);
                
                hold on
                for c = 1:size(M1_LFP_t,1)
                    plot(M1_LFP_t(c,1:min(2000000,size(M1_LFP_t,2)))+(c*0.001))
                end
                
                all_M1_sleep_bad_chans = [all_M1_sleep_bad_chans, input(['M1: Surround multiple bad channels with square brackets. Press return when done.' newline])];
                close all
                hold off
                
                hold on
                for c = 1:size(Cb_LFP_t,1)
                    plot(Cb_LFP_t(c,1:min(2000000,size(Cb_LFP_t,2)))+(c*0.0002))
                end
                
                all_Cb_sleep_bad_chans = [all_Cb_sleep_bad_chans, input(['Cb: Surround multiple bad channels with square brackets. Press return when done.' newline])];
                close all
                hold off
            end
            
        end
        [ct, ch] = hist(all_M1_sleep_bad_chans,unique(all_M1_sleep_bad_chans));
        disp('M1: Surround the channels you wish to exclude with square brackets. Press return when done.')
        disp(ch)
        disp(ct)
        param.M1_sleep_bad_chans = sort(input(''));
        
        [ct, ch] = hist(all_Cb_sleep_bad_chans,unique(all_Cb_sleep_bad_chans));
        disp('Cb: Surround the channels you wish to exclude with square brackets. Press return when done.')
        disp(ch)
        disp(ct)
        param.Cb_sleep_bad_chans = sort(input(''));
    end
    param.M1_sleep_good_chans = 1:param.M1_chans;
    param.M1_sleep_good_chans(param.M1_sleep_bad_chans) = [];
    
    param.Cb_sleep_good_chans = 1:param.Cb_chans;
    param.Cb_sleep_good_chans(param.Cb_sleep_bad_chans) = [];
    
    
    save([rootpath,animal,'/Parameters.mat'],'param');
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Split LFPs into epochs (3)

if enabled(3)
    disp('Block 3...')
    reSort_epochs = false;
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            % Read collated LFP matrix
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/LFP_M1_Truncated.mat']);
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/LFP_Cb_Truncated.mat']);
            
            % Split LFPs into epochs
            M1_epochs = zeros(size(M1_LFP_t,1),floor(param.epoch_length*param.M1_Fs),floor(size(M1_LFP_t,2)/(param.epoch_length*param.M1_Fs)));
            Cb_epochs = zeros(size(Cb_LFP_t,1),floor(param.epoch_length*param.Cb_Fs),floor(size(Cb_LFP_t,2)/(param.epoch_length*param.Cb_Fs)));
            for j=1:floor(size(M1_LFP_t,2)/(param.epoch_length*param.M1_Fs))
                M1_epochs(:,:,j) = M1_LFP_t(:,floor(param.epoch_length*param.M1_Fs)*(j-1)+1:floor(param.epoch_length*param.M1_Fs)*j);
                Cb_epochs(:,:,j) = Cb_LFP_t(:,floor(param.epoch_length*param.Cb_Fs)*(j-1)+1:floor(param.epoch_length*param.Cb_Fs)*j);
            end
            
            % Remove bad Cb epochs
            if ~exist([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/bad_Epochs.mat'],'file') || reSort_epochs
                bad_epochs_M1 = visualizeTrialData(M1_epochs,param.M1_sleep_bad_chans,param.M1_Fs);
                bad_epochs_Cb = visualizeTrialData(Cb_epochs,param.Cb_sleep_bad_chans,param.Cb_Fs);
            else
                load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/bad_Epochs.mat'])
            end
            bad_epochs = union(bad_epochs_M1, bad_epochs_Cb);
            
            bad_epoch_idx = zeros(1,size(M1_epochs,3));
            bad_epoch_idx(bad_epochs) = 1;
            bad_timestamps_M1 = repelem(bad_epoch_idx,size(M1_epochs,2));
            bad_timestamps_Cb = repelem(bad_epoch_idx,size(Cb_epochs,2));
            bad_timestamps_M1(size(M1_LFP_t,2)) = 0;
            bad_timestamps_Cb(size(Cb_LFP_t,2)) = 0;
            
            
            M1_epochs(param.M1_sleep_bad_chans,:,:) = [];
            Cb_epochs(param.Cb_sleep_bad_chans,:,:) = [];
            good_M1_epochs = M1_epochs;
            good_Cb_epochs = Cb_epochs;
            good_M1_epochs(:,:,bad_epochs) = [];
            good_Cb_epochs(:,:,bad_epochs) = [];
            
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_Epochs.mat'], 'M1_epochs', 'Cb_epochs', 'good_M1_epochs', 'good_Cb_epochs');
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/bad_Epochs.mat'], 'bad_epochs_M1', 'bad_epochs_Cb', 'bad_timestamps_M1', 'bad_timestamps_Cb')
        end
    end
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Review Bad Epochs(4)

if enabled(4)
    disp('Block 4...')
    days = [5];
    blocks = [2];
    for day = days
        for block = blocks
            % Read collated LFP matrix
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_Epochs.mat']);
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/bad_Epochs.mat']);

            % Remove bad Cb epochs
            %not_so_bad_epochs_M1 = visualizeTrialData(M1_epochs(:,:,bad_epochs_M1),[],param.M1_Fs);
            not_so_bad_epochs_Cb = visualizeTrialData(Cb_epochs(:,:,bad_epochs_Cb),[],param.Cb_Fs);
            
            %bad_epochs_M1(not_so_bad_epochs_M1) = [];
            bad_epochs_Cb(not_so_bad_epochs_Cb) = [];
            
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/bad_Epochs.mat'], 'bad_epochs_M1', 'bad_epochs_Cb')
        end
    end
end

%% Calculate power spectral density using eeglab (5)

if enabled(5)
    disp('Block 5...')
    addpath(genpath('Z:\Matlab for analysis\eeglab\functions'))
    for day = 1:param.days
        for block = 1:param.sleep_blocks
  
            % Read LFP matrix
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_Epochs.mat']);
  

            % Z-scoring
            M1_epoch_data_n = zeros(size(good_M1_epochs));
            for j=1:length(param.M1_sleep_good_chans)
                tmp = good_M1_epochs(j,:,:);
                tmp = tmp(:);
                M1_epoch_data_n(j,:,:) = (good_M1_epochs(j,:,:)-mean(tmp))./std(tmp);
            end
            
            Cb_epoch_data_n = zeros(size(good_Cb_epochs));
            for j=1:length(param.Cb_sleep_good_chans)
                tmp = good_Cb_epochs(j,:,:);
                tmp = tmp(:);
                Cb_epoch_data_n(j,:,:) = (good_Cb_epochs(j,:,:)-mean(tmp))./std(tmp);
            end
            
            % Median subtraction
            epoch_data_M1 = bsxfun(@minus,M1_epoch_data_n,median(M1_epoch_data_n,1));
            epoch_data_Cb = bsxfun(@minus,Cb_epoch_data_n,median(Cb_epoch_data_n,1));
  
            
            % Get Event-Related Spectral Perturbation (ERSP) using eeglab for M1
            for j=1:size(epoch_data_M1,1)
                data = squeeze(epoch_data_M1(j,:,:));
                [~,~,~,M1_times,M1_freqs,~,~,ersp_data_j] = newtimef(data,size(epoch_data_M1,2),[-size(epoch_data_M1,2)/2 size(epoch_data_M1,2)/2],param.M1_Fs,[0.01 0.1],...
                    'baseline',NaN,'freqs', [0 60],'verbose','off','trialbase','on');
                clf;
                if ~exist('ersp_M1_data', 'var')
                    ersp_M1_data = zeros([length(param.M1_sleep_good_chans) size(ersp_data_j)]);
                end
                ersp_M1_data(j,:,:,:) = ersp_data_j;
            end
            
            for j=1:size(epoch_data_Cb,1)
                data = squeeze(epoch_data_Cb(j,:,:));
                [~,~,~,Cb_times,Cb_freqs,~,~,ersp_data_j] = newtimef(data,size(epoch_data_Cb,2),[-size(epoch_data_Cb,2)/2 size(epoch_data_Cb,2)/2],param.Cb_Fs,[0.01 0.1],...
                    'baseline',NaN,'freqs', [0 60],'verbose','off','trialbase','on');
                clf;
                if ~exist('ersp_Cb_data', 'var')
                    ersp_Cb_data = zeros([length(param.Cb_sleep_good_chans) size(ersp_data_j)]);
                end
                ersp_Cb_data(j,:,:,:) = ersp_data_j;
            end
  
            % Get normalized ESRP data
            data     = abs(ersp_M1_data);
            bl       = mean(data(:,:,M1_times>-1500&M1_times<-500,:),3);
            bl_std   = std(data(:,:,M1_times>-1500&M1_times<-500,:),[],3);
            datanorm = bsxfun(@minus,data,bl);
            ersp_datanorm_M1 = bsxfun(@rdivide,datanorm,bl_std);
            ersp_M1_data = data;
            
            data     = abs(ersp_Cb_data);
            bl       = mean(data(:,:,M1_times>-1500&Cb_times<-500,:),3);
            bl_std   = std(data(:,:,M1_times>-1500&Cb_times<-500,:),[],3);
            datanorm = bsxfun(@minus,data,bl);
            ersp_datanorm_Cb = bsxfun(@rdivide,datanorm,bl_std);
            ersp_Cb_data = data;
  
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/epoch_ERSP.mat'],'ersp_M1_data','ersp_datanorm_M1','ersp_Cb_data','ersp_datanorm_Cb','M1_times','M1_freqs','Cb_times','Cb_freqs','-v7.3');
            clearvars ersp_M1_data ersp_Cb_data
        end
    end
    rmpath(genpath('Z:\Matlab for analysis\eeglab\functions'))
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Classify sleep as NREM and REM (6)

if enabled(6)
    disp('Block 6...')
    check_clusters = false;
    delta_SO = [0.1 4]; %Hz
    gamma = [30 60]; %Hz
    figure('Color','white','Position',[680 494 598 484]);
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_Epochs.mat'])
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/epoch_ERSP.mat'])
  
  
            % Extract SO-D frequency across epochs
            SO_d = ersp_M1_data(:,M1_freqs>=delta_SO(1)&M1_freqs<=delta_SO(2),:,:);
            G = ersp_M1_data(:,M1_freqs>=gamma(1)&M1_freqs<=gamma(2),:,:);

            % Average across channel-time-frequency for every epoch
            SO_d = squeeze(mean(mean(mean(SO_d,1),2),3));
            G    = squeeze(mean(mean(mean(G,1),2),3));

  
            % K-means to classify sleep epochs
            if ~isfield(param,'sleep_C') %Train a classifier if one doesn't exist
                [idx,param.sleep_C] = kmeans([SO_d G],2,'Distance','sqeuclidean'); %,'Distance','cosine' 'hamming'
                idx = idx';
                gscatter(SO_d,G,idx,'bgm');
                hold on;
                plot(param.sleep_C(:,1),param.sleep_C(:,2),'kx')
                legend('Cluster 1','Cluster 2','Cluster Centroid')
                save([rootpath,animal,'/Parameters.mat'],'param');
            else %Use a classifier if it exists
                [~,idx] = pdist2(param.sleep_C,[SO_d G],'euclidean','Smallest',1);
                gscatter(SO_d,G,idx,'bgm');
                hold on;
                plot(param.sleep_C(:,1),param.sleep_C(:,2),'kx')
                legend('Cluster 1','Cluster 2','Cluster Centroid')
            end
            if check_clusters
                input('Did it go alright? (press enter)') %#ok<UNRCH>
            end
            hold off;
            close;
            
            if param.sleep_C(1,1) < param.sleep_C(2,1)
                temp = idx;
                idx(temp==2) = 1;
                idx(temp==1) = 2;
                clear temp
            end
  
            % Assign a different label to less than 5 consecutive sleep epochs
            G_SO_transition = [0 find(diff(idx)~=0) length(idx)];
            transition_diff = diff(G_SO_transition);
            for n=1:length(transition_diff)
                if transition_diff(n) < 5 && idx(G_SO_transition(n+1)) == 1
                    idx(G_SO_transition(n)+1:G_SO_transition(n+1)) = 3;
                end
            end
  
            % Sperate SO-D and Gamma epochs from M1 and Cb LFPs
            SO_d_epochs_M1 = good_M1_epochs(:,:,idx==1);
            SO_d_epochs_Cb = good_Cb_epochs(:,:,idx==1);
            
            G_epochs_M1 = good_M1_epochs(:,:,idx==2);
            G_epochs_Cb = good_Cb_epochs(:,:,idx==2);
            
            % Save
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Epoch_Classification.mat'],...
                'SO_d_epochs_M1','SO_d_epochs_Cb','G_epochs_M1','G_epochs_Cb','idx');

        end
    end
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Detect Peaks, Throughs and Zero-crossing from a filtered average M1-Cb LFP(7)

if enabled(7)
    disp('Block 7...')
    zci = @(v) find(v(:).*circshift(v(:), [-1 0]) <= 0); %zero-crossing detection function
    std_factor = 0.21;
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Epoch_Classification.mat']);
            
            % Average across channels
            SO_d_M1 = squeeze(mean(SO_d_epochs_M1));
            SO_d_Cb = squeeze(mean(SO_d_epochs_Cb));
            
            % Filter M1 in 0.1-4Hz band
            Fn = param.M1_Fs/2;
            CutOff_freqs = 0.1;
            Wn = CutOff_freqs./Fn;
            filterOrder = 2;
            [b,a] = butter(filterOrder,Wn,'high');
            SO_d_M1_filtered = filtfilt(b,a,double(SO_d_M1));
            CutOff_freqs = 4;
            Wn = CutOff_freqs./Fn;
            filterOrder = 5;
            [b,a] = butter(filterOrder,Wn,'low');
            SO_d_M1_filtered = filtfilt(b,a,SO_d_M1_filtered);
            
            % Filter Cb in 0.1-4Hz band
            Fn = param.M1_Fs/2;
            CutOff_freqs = 0.1;
            Wn = CutOff_freqs./Fn;
            filterOrder = 2;
            [b,a] = butter(filterOrder,Wn,'high');
            SO_d_Cb_filtered = filtfilt(b,a,double(SO_d_Cb));
            CutOff_freqs = 4;
            Wn = CutOff_freqs./Fn;
            filterOrder = 5;
            [b,a] = butter(filterOrder,Wn,'low');
            SO_d_Cb_filtered = filtfilt(b,a,SO_d_Cb_filtered);
            
            % Detect Peaks and troughs zero-crossings in M1
            zero_crossings_M1 = cell(1,size(SO_d_M1_filtered,2));
            peaks_M1 = cell(1,size(SO_d_M1_filtered,2));
            troughs_M1 = cell(1,size(SO_d_M1_filtered,2));
            for e=1:size(SO_d_M1_filtered,2)
                epoch = squeeze(SO_d_M1_filtered(:,e));
                std_epoch = std(epoch);
                zc = [1;zci(epoch)];
                p = [];
                t = [];
                for z=1:size(zc,1)-1
                    if max(epoch(zc(z):zc(z+1))) > 0+std_factor*std_epoch
                        p = [p find(epoch==max(epoch(zc(z):zc(z+1))))];
                    elseif min(epoch(zc(z):zc(z+1))) < 0-std_factor*std_epoch
                        t = [t find(epoch==min(epoch(zc(z):zc(z+1))))];
                    end
                end
                zero_crossings_M1{e} = zc(2:end);
                peaks_M1{e} = p;
                troughs_M1{e} = t;
            end
            
            % Detect Peaks and troughs zero-crossings in Cb
            zero_crossings_Cb = cell(1,size(SO_d_Cb_filtered,2));
            peaks_Cb = cell(1,size(SO_d_Cb_filtered,2));
            troughs_Cb = cell(1,size(SO_d_Cb_filtered,2));
            for e=1:size(SO_d_Cb_filtered,2)
                epoch = squeeze(SO_d_Cb_filtered(:,e));
                std_epoch = std(epoch);
                zc = [1;zci(epoch)];
                p = [];
                t = [];
                for z=1:size(zc,1)-1
                    if max(epoch(zc(z):zc(z+1))) > 0+std_factor*std_epoch
                        p = [p find(epoch==max(epoch(zc(z):zc(z+1))))];
                    elseif min(epoch(zc(z):zc(z+1))) < 0-std_factor*std_epoch
                        t = [t find(epoch==min(epoch(zc(z):zc(z+1))))];
                    end
                end
                zero_crossings_Cb{e} = zc(2:end);
                peaks_Cb{e} = p;
                troughs_Cb{e} = t;
            end
            
            % Save
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Events_Sleep.mat'],...
                'SO_d_M1_filtered','SO_d_Cb_filtered','zero_crossings_M1','zero_crossings_Cb',...
                'peaks_M1','peaks_Cb','troughs_M1','troughs_Cb');
        end
    end
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Detect spindles from average M1-Cb LFPs(8)

if enabled(8)
    disp('Block 8...')
    spindle_dur = 0.5; %seconds
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            % Read data
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/LFP_M1_Truncated.mat']);
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/LFP_Cb_Truncated.mat']);
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/bad_Epochs.mat']);
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Epoch_Classification.mat'],'idx');
            
            % Z-score LFPs
            M1_LFP_t = zscore(M1_LFP_t,[],2);
            Cb_LFP_t = zscore(Cb_LFP_t,[],2);
            
            % Average across channels
            M1_LFP_t = squeeze(mean(M1_LFP_t));
            Cb_LFP_t = squeeze(mean(Cb_LFP_t));
            
            % Filter M1 in 10-14Hz band
            Fn = param.M1_Fs/2;
            CutOff_freqs = 10;
            Wn = CutOff_freqs./Fn;
            filterOrder = 6;
            [b,a] = butter(filterOrder,Wn,'high');
            M1_filtered = filtfilt(b,a,double(M1_LFP_t));
            CutOff_freqs = 14;
            Wn = CutOff_freqs./Fn;
            filterOrder = 8;
            [b,a] = butter(filterOrder,Wn,'low');
            M1_filtered = filtfilt(b,a,M1_filtered);
            
            % Filter Cb in 10-14Hz band
            Fn = param.Cb_Fs/2;
            CutOff_freqs = 10;
            Wn = CutOff_freqs./Fn;
            filterOrder = 6;
            [b,a] = butter(filterOrder,Wn,'high');
            Cb_filtered = filtfilt(b,a,double(Cb_LFP_t));
            CutOff_freqs = 14;
            Wn = CutOff_freqs./Fn;
            filterOrder = 8;
            [b,a] = butter(filterOrder,Wn,'low');
            Cb_filtered = filtfilt(b,a,Cb_filtered);
            
            % Extract Envelop using Hilbert transform and Gaussian covolution
            [env_upper_M1,env_lower_M1] = envelope(M1_filtered);
            [env_upper_Cb,env_lower_Cb] = envelope(Cb_filtered);
            
            % Split LFPs into epochs
            epochs_M1 = zeros(floor(param.epoch_length*param.M1_Fs),floor(size(M1_LFP_t,2)/(param.epoch_length*param.M1_Fs)));
            epochs_env_upper_M1 = zeros(floor(param.epoch_length*param.M1_Fs),floor(size(M1_LFP_t,2)/(param.epoch_length*param.M1_Fs)));
            epochs_env_lower_M1 = zeros(floor(param.epoch_length*param.M1_Fs),floor(size(M1_LFP_t,2)/(param.epoch_length*param.M1_Fs)));
            for j=1:floor(size(M1_LFP_t,2)/(param.epoch_length*param.M1_Fs))
                epochs_M1(:,j) = M1_filtered(floor(param.epoch_length*param.M1_Fs)*(j-1)+1:floor(param.epoch_length*param.M1_Fs)*j);
                
                epochs_env_upper_M1(:,j) = env_upper_M1(floor(param.epoch_length*param.M1_Fs)*(j-1)+1:floor(param.epoch_length*param.M1_Fs)*j);
                epochs_env_lower_M1(:,j) = env_lower_M1(floor(param.epoch_length*param.M1_Fs)*(j-1)+1:floor(param.epoch_length*param.M1_Fs)*j);
            end
            
            epochs_Cb = zeros(floor(param.epoch_length*param.Cb_Fs),floor(size(Cb_LFP_t,2)/(param.epoch_length*param.Cb_Fs)));
            epochs_env_upper_Cb = zeros(floor(param.epoch_length*param.Cb_Fs),floor(size(Cb_LFP_t,2)/(param.epoch_length*param.Cb_Fs)));
            epochs_env_lower_Cb = zeros(floor(param.epoch_length*param.Cb_Fs),floor(size(Cb_LFP_t,2)/(param.epoch_length*param.Cb_Fs)));
            for j=1:floor(size(Cb_LFP_t,2)/(param.epoch_length*param.Cb_Fs))
                epochs_Cb(:,j) = Cb_filtered(floor(param.epoch_length*param.Cb_Fs)*(j-1)+1:floor(param.epoch_length*param.Cb_Fs)*j);
                
                epochs_env_upper_Cb(:,j) = env_upper_Cb(floor(param.epoch_length*param.Cb_Fs)*(j-1)+1:floor(param.epoch_length*param.Cb_Fs)*j);
                epochs_env_lower_Cb(:,j) = env_lower_Cb(floor(param.epoch_length*param.Cb_Fs)*(j-1)+1:floor(param.epoch_length*param.Cb_Fs)*j);
            end
            
            % Remove bad epochs
            bad_epochs = union(bad_epochs_M1, bad_epochs_Cb);
            epochs_M1(:,bad_epochs) = [];
            epochs_Cb(:,bad_epochs) = [];
            
            epochs_env_upper_M1(:,bad_epochs) = [];
            epochs_env_lower_M1(:,bad_epochs) = [];
            
            epochs_env_upper_Cb(:,bad_epochs) = [];
            epochs_env_lower_Cb(:,bad_epochs) = [];
            
            % Extract only NREM sleep epochs
            epochs_M1 = epochs_M1(:,idx==1);
            epochs_Cb = epochs_Cb(:,idx==1);
            
            epochs_env_upper_M1 = epochs_env_upper_M1(:,idx==1);
            epochs_env_lower_M1 = epochs_env_lower_M1(:,idx==1);
            
            epochs_env_upper_Cb = epochs_env_upper_Cb(:,idx==1);
            epochs_env_lower_Cb = epochs_env_lower_Cb(:,idx==1);
            
            % Detect spindles
            nEpochs = size(epochs_M1,2);
            %disp(nEpochs)
            
            % Initialize cell arrays for storage
            spindle_ts_M1  = cell(nEpochs,1);
            spindle_amp_M1 = cell(nEpochs,1);
            spindle_ts_Cb  = cell(nEpochs,1);
            spindle_amp_Cb = cell(nEpochs,1);
            
            % Take mean around envelope
            m_env_M1 = mean(epochs_M1,1);
            m_env_Cb = mean(epochs_Cb,1);
            
            % Take std around envelope
            sd_env_M1 = std(epochs_M1,0,1);
            sd_env_Cb = std(epochs_Cb,0,1);
            
            for j=1:nEpochs
                
                % Get epoch
                epoch_M1 = epochs_M1(:,j);
                epoch_Cb = epochs_Cb(:,j);
                len_epoch = 1:length(epoch_M1);
                
%                 % Take mean around envelope
%                 m_env_M1 = mean(epochs_env_upper_M1(:,j),1);
%                 m_env_Cb = mean(epochs_env_upper_Cb(:,j),1);
%                 
%                 % Take std around envelope
%                 sd_env_M1 = std(epochs_env_upper_M1(:,j));
%                 sd_env_Cb = std(epochs_env_upper_Cb(:,j));
                
                % upper threshold
                u_thresh_M1 = m_env_M1(j)+2.5*sd_env_M1;
                u_thresh_Cb = m_env_Cb(j)+2.5*sd_env_Cb;
                
                % lower threshold
                l_thresh_M1 = m_env_M1(j)+1.5*sd_env_M1(j);
                l_thresh_Cb = m_env_Cb(j)+1.5*sd_env_Cb(j);
                
                % Detect indices that fulfill threshold criteria
                u_thresh_bool_M1  = epochs_env_upper_M1(:,j) > u_thresh_M1;
                l_thresh_bool_M1  = epochs_env_upper_M1(:,j) > l_thresh_M1;
                
                u_thresh_bool_Cb  = epochs_env_upper_Cb(:,j) > u_thresh_Cb;
                l_thresh_bool_Cb  = epochs_env_upper_Cb(:,j) > l_thresh_Cb;
                
                % Detect threshold crossings
                l_thresh_xing_M1 = diff(l_thresh_bool_M1);
                l_thresh_xing_Cb = diff(l_thresh_bool_Cb);
                
                % Find lower thresh crossing timestamps
                l_thresh_xing_M1_idx = find(l_thresh_xing_M1);
                l_thresh_xing_Cb_idx = find(l_thresh_xing_Cb);
                
                % Mark lower thresh crossing timestamps that go from below to above
                l_thresh_rising_M1_idx = find(l_thresh_xing_M1 > 0);
                l_thresh_rising_Cb_idx = find(l_thresh_xing_Cb > 0);
                
                if ~isempty(l_thresh_xing_M1_idx)
                    %If the first crossing is not rising, i.e. if the epoch started above threshold, prepend 0 to the rising indices
                    if ~isempty(l_thresh_rising_M1_idx) && l_thresh_rising_M1_idx(1) ~= l_thresh_xing_M1_idx(1)
                        l_thresh_rising_M1_idx = [0; l_thresh_rising_M1_idx];
                    end
                    
                    %If the last crossing is rising, i.e. if the epoch ends above threshold, append the last index of the epoch to the crossings
                    if ~isempty(l_thresh_rising_M1_idx) && l_thresh_rising_M1_idx(end) == l_thresh_xing_M1_idx(end)
                        l_thresh_xing_M1_idx = [l_thresh_xing_M1_idx; length(epochs_env_upper_M1)];
                    end
                    
                    %Detect and save spindles
                    epoch_spindles_ts_M1 = cell(0);
                    epoch_spindles_amp_M1 = cell(0);
                    for r_idx = l_thresh_rising_M1_idx'+1 % the +1 is so the indices are at the first point above the threshold
                        above_l_dur = l_thresh_xing_M1_idx(find(l_thresh_xing_M1_idx >= r_idx,1))-r_idx;
                        if above_l_dur >= spindle_dur*param.M1_Fs && ismember(1,u_thresh_bool_M1(r_idx:(r_idx+above_l_dur-1)))
                            epoch_spindles_ts_M1 = [epoch_spindles_ts_M1; {len_epoch(r_idx:(r_idx+above_l_dur-1))}];
                            epoch_spindles_amp_M1 = [epoch_spindles_amp_M1; {epoch_M1(r_idx:(r_idx+above_l_dur-1))}];
                        end
                    end
                    spindle_ts_M1{j} = epoch_spindles_ts_M1;
                    spindle_amp_M1{j} = epoch_spindles_amp_M1;
                end
                
                if ~isempty(l_thresh_xing_Cb_idx)
                    %If the first crossing is not rising, i.e. if the epoch started above threshold, prepend 0 to the rising indices
                    if ~isempty(l_thresh_rising_Cb_idx) && l_thresh_rising_Cb_idx(1) ~= l_thresh_xing_Cb_idx(1)
                        l_thresh_rising_Cb_idx = [0; l_thresh_rising_Cb_idx];
                    end
                    
                    %If the last crossing is rising, i.e. if the epoch ends above threshold, append the last index of the epoch to the crossings
                    if ~isempty(l_thresh_rising_Cb_idx) && l_thresh_rising_Cb_idx(end) == l_thresh_xing_Cb_idx(end)
                        l_thresh_xing_Cb_idx = [l_thresh_xing_Cb_idx; length(epochs_env_upper_Cb)];
                    end
                    
                    %Detect and save spindles
                    epoch_spindles_ts_Cb = cell(0);
                    epoch_spindles_amp_Cb = cell(0);
                    for r_idx = l_thresh_rising_Cb_idx'+1 % the +1 is so the indices are at the first point above the threshold
                        above_l_dur = l_thresh_xing_Cb_idx(find(l_thresh_xing_Cb_idx >= r_idx,1))-r_idx;
                        if above_l_dur >= spindle_dur*param.Cb_Fs && ismember(1,u_thresh_bool_Cb(r_idx:(r_idx+above_l_dur-1)))
                            epoch_spindles_ts_Cb = [epoch_spindles_ts_Cb; {len_epoch(r_idx:(r_idx+above_l_dur-1))}];
                            epoch_spindles_amp_Cb = [epoch_spindles_amp_Cb; {epoch_Cb(r_idx:(r_idx+above_l_dur-1))}];
                        end
                    end
                    spindle_ts_Cb{j} = epoch_spindles_ts_Cb;
                    spindle_amp_Cb{j} = epoch_spindles_amp_Cb;
                end
                
            end
            % Save
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/spindles_Sleep.mat'],...
                'spindle_amp_M1','spindle_ts_M1','spindle_amp_Cb','spindle_ts_Cb');
            
            %TODO: add spindle counts and total sleep times to shared data
        end
    end
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Use Jaekyung Code to classify time-steps as NREM or REM/awake (use with 10 and 11 as an alternative to 3-8) (9) 

if enabled(9)
    disp('Block 9...')
    addpath(genpath('C:\Users\FleischerP\Documents\MATLAB\power-based sleep detection code'))
    load([rootpath,animal,'/Shared_Data.mat'])
    shared_data.sleep_dur = nan(param.days, param.sleep_blocks);
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/LFP_M1_Truncated.mat']);
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/LFP_Cb_Truncated.mat']);
            if exist([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_data.mat'],'file')
                load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_data.mat'])
            else
                M1_data = struct();
                Cb_data = struct();
            end
            
            if exist([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Channel_Sleep_data.mat'],'file')
                success = false;
                fails = 0;
                while fails < 5 && ~success
                    try
                        load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Channel_Sleep_data.mat'])
                        success = true;
                    catch
                        success = false;
                        fails = fails + 1;
                    end
                end
                if fails == 5
                    error('Too much failure')
                end
            else
                M1_channel_data = struct();
                Cb_channel_data = struct();
            end
            
            M1_data.LFP = double(mean(M1_LFP_t(param.M1_sleep_good_chans,:),1)');
            M1_data.Fs_LFP = param.M1_Fs;
            Cb_data.LFP = double(mean(Cb_LFP_t(param.Cb_sleep_good_chans,:),1)');
            Cb_data.Fs_LFP = param.Cb_Fs;
            
            M1_channel_data.LFP = double(M1_LFP_t(param.M1_sleep_good_chans,:)');
            M1_channel_data.Fs_LFP = param.M1_Fs;
            Cb_channel_data.LFP = double(Cb_LFP_t(param.Cb_sleep_good_chans,:)');
            Cb_channel_data.Fs_LFP = param.Cb_Fs;
            
            [sleep_idx_M1, artifact_idx_M1, M1_data.pwr_out_M1] = sleep_classification2(M1_data.LFP,1,M1_data.Fs_LFP);
            saveas(gcf, [rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_Classification_M1.fig']);
            close all
            [sleep_idx_Cb, artifact_idx_Cb, Cb_data.pwr_out_Cb] = sleep_classification2(Cb_data.LFP,1,Cb_data.Fs_LFP);
            saveas(gcf, [rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_Classification_Cb.fig']);
            close all
            
            M1_data.artifact_idx = artifact_idx_M1 | artifact_idx_Cb;
            M1_data.sleep_idx = sleep_idx_M1;
            Cb_data.artifact_idx = M1_data.artifact_idx;
            Cb_data.sleep_idx = M1_data.sleep_idx;
            
            shared_data.sleep_dur(day, block) = sum(sleep_idx_M1)/M1_data.Fs_LFP;
            
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_data.mat'], 'M1_data', 'Cb_data')
                        
            M1_channel_data.artifact_idx = artifact_idx_M1 | artifact_idx_Cb;
            M1_channel_data.sleep_idx = sleep_idx_M1;
            Cb_channel_data.artifact_idx = M1_data.artifact_idx;
            Cb_channel_data.sleep_idx = M1_data.sleep_idx;
            
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Channel_Sleep_data.mat'], 'M1_channel_data', 'Cb_channel_data','-v7.3')
        end
    end
    
    save([rootpath,animal,'/Shared_Data.mat'], 'shared_data');
    rmpath(genpath('C:\Users\FleischerP\Documents\MATLAB\power-based sleep detection code'))
end

%% Use Jaekyung Code to find Delta and Slow Oscillations (use with 9 and 11 as an alternative to 3-8) (10) 

if enabled(10)
    disp('Block 10...')
    addpath(genpath('C:\Users\FleischerP\Documents\MATLAB\SO Delta detection code'))
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_data.mat']);
            
            M1_data.so_delta=detect_so_delta(M1_data.LFP,M1_data.Fs_LFP,...
                'sleep_idx',M1_data.sleep_idx,...
                'artifact_idx',M1_data.artifact_idx,...
                'PLOT',1,...
                'mnl_parm',[85 40 .15 .5]);
            figHandles = findobj('Type', 'figure');
            saveas(figHandles(1),[rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/M1_',get(figHandles(1),'Name'),'_Waves.fig'])
            saveas(figHandles(2),[rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/M1_',get(figHandles(2),'Name'),'_Waves.fig'])
            close all;
            
            Cb_data.so_delta=detect_so_delta(Cb_data.LFP,Cb_data.Fs_LFP,...
                'sleep_idx',Cb_data.sleep_idx,...
                'artifact_idx',Cb_data.artifact_idx,...
                'PLOT',1,...
                'mnl_parm',[85 40 .15 .5]);
            figHandles = findobj('Type', 'figure');
            saveas(figHandles(1),[rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Cb_',get(figHandles(1),'Name'),'_Waves.fig'])
            saveas(figHandles(2),[rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Cb_',get(figHandles(2),'Name'),'_Waves.fig'])
            close all;
            
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_data.mat'], 'M1_data', 'Cb_data')
            
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Channel_Sleep_data.mat'])
            for M1_chan = 1:size(M1_channel_data.LFP,2)
                M1_channel_data.so_delta(:,M1_chan)=detect_so_delta(M1_channel_data.LFP(:,M1_chan),M1_channel_data.Fs_LFP,...
                    'sleep_idx',M1_channel_data.sleep_idx,...
                    'artifact_idx',M1_channel_data.artifact_idx,...
                    'PLOT',1,...
                    'mnl_parm',[85 40 .15 .5]);
                close all;
            end
            for Cb_chan = 1:size(Cb_channel_data.LFP,2)
                Cb_channel_data.so_delta(:,Cb_chan)=detect_so_delta(Cb_channel_data.LFP(:,Cb_chan),Cb_channel_data.Fs_LFP,...
                    'sleep_idx',Cb_channel_data.sleep_idx,...
                    'artifact_idx',Cb_channel_data.artifact_idx,...
                    'PLOT',1,...
                    'mnl_parm',[85 40 .15 .5]);
                close all;
            end
            
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Channel_Sleep_data.mat'], 'M1_channel_data', 'Cb_channel_data','-v7.3')
        end
    end
    rmpath(genpath('C:\Users\FleischerP\Documents\MATLAB\SO Delta detection code'))
end

%% Use Jaekyung Code to find Sleep spindles (use with 9 and 10 as an alternative to 3-8) (11) 

if enabled(11)
    disp('Block 11...')
    addpath(genpath('C:\Users\FleischerP\Documents\MATLAB\Spindle detection code'))
    spindle_count_M1 = nan(param.days, param.sleep_blocks);
    spindle_count_Cb = nan(param.days, param.sleep_blocks);
    sleep_duration = nan(param.days, param.sleep_blocks);
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_data.mat']);
            
            M1_data.spindles = detect_spindles( mat2cell(M1_data.LFP, length(M1_data.LFP), [1]),...
                'Fs',M1_data.Fs_LFP,...
                'sleep_idx',mat2cell(M1_data.sleep_idx, length(M1_data.LFP), [1]),...
                'artifact_idx',mat2cell(M1_data.artifact_idx, length(M1_data.LFP), [1]),...
                'PLOT',1,...
                'sleep_classify',1);
            saveas(gcf,[rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/spindles_M1_fig.fig'])
            close all;
            
            Cb_data.spindles = detect_spindles( mat2cell(Cb_data.LFP, length(Cb_data.LFP), [1]),...
                'Fs',Cb_data.Fs_LFP,...
                'sleep_idx',mat2cell(Cb_data.sleep_idx, length(Cb_data.LFP), [1]),...
                'artifact_idx',mat2cell(Cb_data.artifact_idx, length(Cb_data.LFP), [1]),...
                'PLOT',1,...
                'sleep_classify',1);
            saveas(gcf,[rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/spindles_Cb_fig.fig'])
            close all;
            
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_data.mat'], 'M1_data', 'Cb_data')
            
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Channel_Sleep_data.mat'])
            for M1_chan = 1:size(M1_channel_data.LFP,2)
                M1_channel_data.spindles(:,M1_chan) = detect_spindles( mat2cell(M1_channel_data.LFP(:,M1_chan), size(M1_channel_data.LFP,1), [1]),...
                    'Fs',M1_channel_data.Fs_LFP,...
                    'sleep_idx',mat2cell(M1_channel_data.sleep_idx, size(M1_channel_data.LFP,1), [1]),...
                    'artifact_idx',mat2cell(M1_channel_data.artifact_idx, size(M1_channel_data.LFP,1), [1]),...
                    'PLOT',1,...
                    'sleep_classify',1);
                close all;
            end
            for Cb_chan = 1:size(Cb_channel_data.LFP,2)
                Cb_channel_data.spindles(:,Cb_chan) = detect_spindles( mat2cell(Cb_channel_data.LFP(:,Cb_chan), size(Cb_channel_data.LFP,1), [1]),...
                    'Fs',Cb_channel_data.Fs_LFP,...
                    'sleep_idx',mat2cell(Cb_channel_data.sleep_idx, size(Cb_channel_data.LFP,1), [1]),...
                    'artifact_idx',mat2cell(Cb_channel_data.artifact_idx, size(Cb_channel_data.LFP,1), [1]),...
                    'PLOT',1,...
                    'sleep_classify',1);
                close all;
            end
            
            spindle_count_M1(day, block) = length(M1_data.spindles{1}.pks);
            spindle_count_Cb(day, block) = length(Cb_data.spindles{1}.pks);
            sleep_duration(day, block) = sum(M1_data.sleep_idx)/M1_data.Fs_LFP;
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Channel_Sleep_data.mat'], 'M1_channel_data', 'Cb_channel_data','-v7.3')
        end
    end
    load([rootpath,animal,'/Shared_Data.mat'])
    shared_data.spindle_count_M1 = spindle_count_M1;
    shared_data.spindle_count_Cb = spindle_count_Cb;
    shared_data.sleep_duration = sleep_duration;
    save([rootpath,animal,'/Shared_Data.mat'], 'shared_data');
    rmpath(genpath('C:\Users\FleischerP\Documents\MATLAB\spindle detection code'))
end

%% Get Sleep/wake indexes (12)
%Classify all LFP into sleep/wake. After that it will be easy to implemnt whether or not wake should override NREM.

if enabled(12)
    disp('Block 12...')
    framerate = 1; %hz
    movement_threshold = 0.65;
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([origin_rootpath,animal,'/Day',num2str(day),'/Sleep_Vids/',param.sleep_block_vid_names{day, block},'.mat'])
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/WAV_Truncated.mat']);
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_data.mat']);
            start_idx = find(abs(diff(WAVE2_t))>0.05, 1);
            detected_sleep = motion <= movement_threshold;
            detected_sleep = repmat(detected_sleep,round(param.M1_Fs/framerate),1);
            detected_sleep = detected_sleep(:);
            sleep = zeros(1,length(WAVE2_t));
            sleep(start_idx:(length(detected_sleep)+start_idx-1)) = detected_sleep;
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_data.mat'], 'M1_data', 'Cb_data', 'sleep');
        end
    end
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Get spindles Frequency(13)

if enabled(13)
    disp('Block 13...')
    for day = 1:param.days
        spindle_freq_M1 = zeros(1,param.sleep_blocks);
        spindle_freq_Cb = zeros(1,param.sleep_blocks);
        spindle_count_M1 = zeros(1,param.sleep_blocks);
        spindle_count_Cb = zeros(1,param.sleep_blocks);
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/spindles_Sleep.mat']);
            
            % Unpack sindle cell arrays and count the spindles
            Epoch_idx = 1:size(spindle_ts_Cb,1);
            nEpochs = length(Epoch_idx);
            bfr = logical(~cellfun('isempty',spindle_ts_M1(Epoch_idx)).*~cellfun('isempty',spindle_ts_Cb(Epoch_idx)));
            spindle_ts_M1 = spindle_ts_M1(bfr);
            spindle_ts_Cb = spindle_ts_Cb(bfr);
            count_M1 = zeros(1,size(spindle_ts_M1,1));
            count_Cb = zeros(1,size(spindle_ts_M1,1));
            for j=1:size(spindle_ts_M1,1)
                count_M1(j) = size(spindle_ts_M1{j},1);
                count_Cb(j) = size(spindle_ts_Cb{j},1);
            end
            spindle_freq_M1(block) = sum(count_M1)/(nEpochs*param.epoch_length);
            spindle_freq_Cb(block) = sum(count_Cb)/(nEpochs*param.epoch_length);
            
            spindle_count_M1(block) = sum(count_M1);
            spindle_count_Cb(block) = sum(count_Cb);
        end
        subplot(2,2,1);
        bar(spindle_freq_M1);%ylim([0 800]);
        title('M1');
        ylabel('spindle Frequency (hz)');
        subplot(2,2,2);
        bar(spindle_freq_Cb);%ylim([0 800]);
        title('Cb');
        subplot(2,2,3);
        bar(spindle_count_M1);%ylim([0 800]);
        xlabel('Block');
        ylabel('spindle Count');
        subplot(2,2,4);
        bar(spindle_count_Cb);%ylim([0 800]);
        xlabel('Block');
        saveas(gcf,[rootpath,animal,'/Day',num2str((day)),'/spindle_stats.fig'])
        close all;
    end
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Detect Slow-Oscillations and Delta waves timestamps (14)

if enabled(14)
    disp('Block 14...')
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Events_Sleep.mat']);
  
            nEpochs = length(peaks_Cb);
            UP_State_M1 = zeros(1,nEpochs);
            UP_State_Cb = zeros(1,nEpochs);
            DOWN_State_M1 = zeros(1,nEpochs);
            DOWN_State_Cb = zeros(1,nEpochs);
            SO_timestamps_M1 = cell(1,nEpochs);
            SO_timestamps_Cb = cell(1,nEpochs);
            D_timestamps_M1 = cell(1,nEpochs);
            D_timestamps_Cb = cell(1,nEpochs);
            for j=1:nEpochs
                bfr1 = SO_d_M1_filtered(:,j);
                bfr2 = SO_d_Cb_filtered(:,j);
                
                % Define UP-state threshold
                UP_State_M1(j)   = prctile(bfr1(troughs_M1{j}),60);
                UP_State_Cb(j)   = prctile(bfr2(troughs_Cb{j}),60);
                
                % Define DOWN-state threshold
                DOWN_State_M1(j) = prctile(bfr1(peaks_M1{j}),80);
                DOWN_State_Cb(j) = prctile(bfr2(peaks_Cb{j}),80);
                
                % Amplitude of peaks and troughs
                tmp1_p = bfr1(peaks_M1{j});
                tmp1_t = bfr1(troughs_M1{j});
                tmp2_p = bfr2(peaks_Cb{j});
                tmp2_t = bfr2(troughs_Cb{j});
                
                % Detect SO timestamps
                bfr1_p = peaks_M1{j}(tmp1_p > DOWN_State_M1(j));
                bfr1_t = troughs_M1{j}(tmp1_t < UP_State_M1(j));
                bfr2_p = peaks_Cb{j}(tmp2_p > DOWN_State_Cb(j));
                bfr2_t = troughs_Cb{j}(tmp2_t < UP_State_Cb(j));
                
                ts = [];
                if ~isempty(bfr1_t) && ~isempty(bfr1_p) && ~isempty(bfr2_t) && ~isempty(bfr2_p)
                    bfr1_t = bfr1_t(bfr1_t>bfr1_p(1));
                    for k = 1:length(bfr1_t)
                        for l = 1:length(bfr1_p)
                            if((bfr1_t(k)-bfr1_p(l))/param.M1_Fs*1000 >= 150 && (bfr1_t(k)-bfr1_p(l))/param.M1_Fs*1000 <= 500)
                                z = zero_crossings_M1{j};
                                ts = [ts z(find(z<bfr1_p(l),1,'last')) bfr1_p(l) bfr1_t(k) z(find(z>bfr1_t(k),1,'first'))];
                                break;
                            end
                        end
                    end
                end
                SO_timestamps_M1{j} = sort(ts);
                if ~isempty(SO_timestamps_M1{j})
                    troughs_M1{j} = troughs_M1{j}(~ismember(troughs_M1{j},SO_timestamps_M1{j}));
                    tmp1_t = tmp1_t(~ismember(troughs_M1{j},SO_timestamps_M1{j}));
                end
                
                ts = [];
                if ~isempty(bfr1_t) && ~isempty(bfr1_p) && ~isempty(bfr2_t) && ~isempty(bfr2_p)
                    bfr2_t = bfr2_t(bfr2_t>bfr2_p(1));
                    for k = 1:length(bfr2_t)
                        for l = 1:length(bfr2_p)
                            if((bfr2_t(k)-bfr2_p(l))/param.Cb_Fs*1000 >= 150 && (bfr2_t(k)-bfr2_p(l))/param.Cb_Fs*1000 <= 500)
                                z = zero_crossings_Cb{j};
                                ts = [ts z(find(z<bfr2_p(l),1,'last')) bfr2_p(l) bfr2_t(k) z(find(z>bfr2_t(k),1,'first'))];
                                break;
                            end
                        end
                    end
                end
                SO_timestamps_Cb{j} = sort(ts);
                if ~isempty(SO_timestamps_Cb{j})
                    troughs_Cb{j} = troughs_Cb{j}(~ismember(troughs_Cb{j},SO_timestamps_Cb{j}));
                    tmp2_t = tmp2_t(~ismember(troughs_Cb{j},SO_timestamps_Cb{j}));
                end
                
                % Detect Delta timestamps
                bfr1_p = peaks_M1{j}(tmp1_p < DOWN_State_M1(j));
                bfr1_t = troughs_M1{j}(tmp1_t < UP_State_M1(j));
                bfr2_p = peaks_Cb{j}(tmp2_p < DOWN_State_Cb(j));
                bfr2_t = troughs_Cb{j}(tmp2_t < UP_State_Cb(j));
                
                ts = [];
                if ~isempty(bfr1_t) && ~isempty(bfr1_p) && ~isempty(bfr2_t) && ~isempty(bfr2_p)
                    bfr1_t = bfr1_t(bfr1_t>bfr1_p(1));
                    for k = 1:length(bfr1_t)
                        for l = 1:length(bfr1_p)
                            if((bfr1_t(k)-bfr1_p(l))/param.M1_Fs*1000 <= 500)
                                z = zero_crossings_M1{j};
                                ts = [ts z(find(z<bfr1_p(l),1,'last')) bfr1_p(l) bfr1_t(k) z(find(z>bfr1_t(k),1,'first'))];
                                break;
                            end
                        end
                    end
                end
                D_timestamps_M1{j} = sort(ts);
                
                ts = [];
                if ~isempty(bfr1_t) && ~isempty(bfr1_p) && ~isempty(bfr2_t) && ~isempty(bfr2_p)
                    bfr2_t = bfr2_t(bfr2_t>bfr2_p(1));
                    for k = 1:length(bfr2_t)
                        for l = 1:length(bfr2_p)
                            if((bfr2_t(k)-bfr2_p(l))/param.Cb_Fs*1000 <= 500)
                                z = zero_crossings_Cb{j};
                                ts = [ts z(find(z<bfr2_p(l),1,'last')) bfr2_p(l) bfr2_t(k) z(find(z>bfr2_t(k),1,'first'))];
                                break;
                            end
                        end
                    end
                end
                D_timestamps_Cb{j} = sort(ts);
            end
            % Save
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Events_SOD_Sleep.mat'],...
                'SO_timestamps_M1','SO_timestamps_Cb','D_timestamps_M1','D_timestamps_Cb',...
                'UP_State_M1','UP_State_Cb','DOWN_State_M1','DOWN_State_Cb');
        end
    end
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Classify LFP into Slow-Oscillations and Delta oscillations(15)

if enabled(15)
    disp('Block 15...')
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Events_SOD_Sleep.mat']);
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Epoch_Classification.mat']);
            
            nEpochs = length(UP_State_M1);
            SO_M1 = cell(1,nEpochs);
            SO_Cb = cell(1,nEpochs);
            D_M1 = cell(1,nEpochs);
            D_Cb = cell(1,nEpochs);
            for j=1:nEpochs
                
                % Filter M1 channels in 0.1-4Hz band
                lfps_M1 = SO_d_epochs_M1(:,:,j)';
                Fn = param.M1_Fs/2;
                CutOff_freqs = 0.1;
                Wn = CutOff_freqs./Fn;
                filterOrder = 2;
                [b,a] = butter(filterOrder,Wn,'high');
                bfr = filtfilt(b,a,double(lfps_M1));
                CutOff_freqs = 4;
                Wn = CutOff_freqs./Fn;
                filterOrder = 5;
                [b,a] = butter(filterOrder,Wn,'low');
                lfps_M1_filtered = filtfilt(b,a,bfr)';
                
                % Filter Cb channels in 0.1-4Hz band
                lfps_Cb = SO_d_epochs_Cb(:,:,j)';
                Fn = param.Cb_Fs/2;
                CutOff_freqs = 0.1;
                Wn = CutOff_freqs./Fn;
                filterOrder = 2;
                [b,a] = butter(filterOrder,Wn,'high');
                bfr = filtfilt(b,a,double(lfps_Cb));
                CutOff_freqs = 4;
                Wn = CutOff_freqs./Fn;
                filterOrder = 5;
                [b,a] = butter(filterOrder,Wn,'low');
                lfps_Cb_filtered = filtfilt(b,a,bfr)';
                
                % Extract Slow-oscillations
                ts = SO_timestamps_M1{j};
                if ~isempty(ts) && mod(length(ts),4)==0
                    ts = reshape(ts,4,length(ts)/4);
                    tmp = cell(1,size(ts,2));
                    for k=1:size(ts,2)
                        tmp{k} = lfps_M1_filtered(:,ts(1,k):ts(end,k));
                    end
                else
                    tmp = {};
                end
                SO_M1{j} = tmp;
                
                ts = SO_timestamps_Cb{j};
                if ~isempty(ts) && mod(length(ts),4)==0
                    ts = reshape(ts,4,length(ts)/4);
                    tmp = cell(1,size(ts,2));
                    for k=1:size(ts,2)
                        tmp{k} = lfps_Cb_filtered(:,ts(1,k):ts(end,k));
                    end
                else
                    tmp = {};
                end
                SO_Cb{j} = tmp;
                
                % Extract Delta-oscillations
                ts = D_timestamps_M1{j};
                if ~isempty(ts) && mod(length(ts),4)==0
                    ts = reshape(ts,4,length(ts)/4);
                    tmp = cell(1,size(ts,2));
                    for k=1:size(ts,2)
                        tmp{k} = lfps_M1_filtered(:,ts(1,k):ts(end,k));
                    end
                else
                    tmp = {};
                end
                D_M1{j} = tmp;
                
                ts = D_timestamps_Cb{j};
                if ~isempty(ts) && mod(length(ts),4)==0
                    ts = reshape(ts,4,length(ts)/4);
                    tmp = cell(1,size(ts,2));
                    for k=1:size(ts,2)
                        tmp{k} = lfps_Cb_filtered(:,ts(1,k):ts(end,k));
                    end
                else
                    tmp = {};
                end
                D_Cb{j} = tmp;
            end
            
            % Save
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Oscillations_trace.mat'],...
                'SO_M1','SO_Cb','D_M1','D_Cb');
        end
    end
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Detect Nested spindles(16)

if enabled(16)
    disp('Block 16...')
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/spindles_Sleep.mat']);
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Events_SOD_Sleep.mat']);
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_Epochs.mat']);
            
            nEpochs = size(spindle_ts_Cb,1);
            SO_nested_spindle_M1 = cell(1,nEpochs);
            SO_nested_spindle_Cb = cell(1,nEpochs);
            D_nested_spindle_M1 = cell(1,nEpochs);
            D_nested_spindle_Cb = cell(1,nEpochs);
            for j=1:nEpochs
                % Get epoch
                epoch_M1 = M1_epochs(:,:,j);
                epoch_Cb = Cb_epochs(:,:,j);
                len_epoch = 1:length(epoch_M1);
                
                % Detect SO nested spindles
                tmp = SO_timestamps_M1{j};
                if mod(length(tmp),4) == 0
                    bfr = [];
                    for k=3:4:length(tmp) %take only troughs
                        for l=1:size(spindle_ts_M1{j},1)
                            tmp2 = spindle_ts_M1{j}{l};
                            spindle_peak = tmp2(spindle_amp_M1{j}{l} == max(spindle_amp_M1{j}{l}));
                            diff_SO_spindle = tmp(k) - spindle_peak;
                            if diff_SO_spindle>-0.5*param.M1_Fs && diff_SO_spindle<1*param.M1_Fs
                                bfr = [bfr,spindle_peak];
                            end
                        end
                    end
                end
                if exist('bfr','var') == 1
                    SO_nested_spindle_M1{j} = bfr;
                    clear bfr;
                end
                
                tmp = SO_timestamps_Cb{j};
                if mod(length(tmp),4) == 0
                    bfr = [];
                    for k=3:4:length(tmp) %take only troughs
                        for l=1:size(spindle_ts_Cb{j},1)
                            tmp2 = spindle_ts_Cb{j}{l};
                            spindle_peak = tmp2(spindle_amp_Cb{j}{l} == max(spindle_amp_Cb{j}{l}));
                            diff_SO_spindle = tmp(k) - spindle_peak;
                            if diff_SO_spindle>-0.5*param.Cb_Fs && diff_SO_spindle<1*param.Cb_Fs
                                bfr = [bfr,spindle_peak];
                            end
                        end
                    end
                end
                if exist('bfr','var') == 1
                    SO_nested_spindle_Cb{j} = bfr;
                    clear bfr;
                end
                
                % Detect Delta nested spindles
                tmp = D_timestamps_M1{j};
                if mod(length(tmp),4) == 0
                    bfr = [];
                    for k=3:4:length(tmp) %take only troughs
                        for l=1:size(spindle_ts_M1{j},1)
                            tmp2 = spindle_ts_M1{j}{l};
                            spindle_peak = tmp2(spindle_amp_M1{j}{l} == max(spindle_amp_M1{j}{l}));
                            diff_D_spindle = tmp(k) - spindle_peak;
                            if diff_D_spindle>-0.5*param.M1_Fs && diff_D_spindle<1*param.M1_Fs
                                bfr = [bfr,spindle_peak];
                            end
                        end
                    end
                end
                if exist('bfr','var') == 1
                    D_nested_spindle_M1{j} = bfr;
                    clear bfr;
                end
                
                tmp = D_timestamps_Cb{j};
                if mod(length(tmp),4) == 0
                    bfr = [];
                    for k=3:4:length(tmp) %take only troughs
                        for l=1:size(spindle_ts_Cb{j},1)
                            tmp2 = spindle_ts_Cb{j}{l};
                            spindle_peak = tmp2(spindle_amp_Cb{j}{l} == max(spindle_amp_Cb{j}{l}));
                            diff_D_spindle = tmp(k) - spindle_peak;
                            if diff_D_spindle>-0.5*param.Cb_Fs && diff_D_spindle<1*param.Cb_Fs
                                bfr = [bfr,spindle_peak];
                            end
                        end
                    end
                end
                if exist('bfr','var') == 1
                    D_nested_spindle_Cb{j} = bfr;
                    clear bfr;
                end
            end
            
            % Save
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Nested_spindles_Sleep.mat'],...
                'SO_nested_spindle_M1','SO_nested_spindle_Cb','D_nested_spindle_M1','D_nested_spindle_Cb');
        end
    end
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Get frequency of Slow and Delta oscillations from block1 to block2 (17)

if enabled(17)
    disp('Block 17...')
    for day = 1:param.days
        SO_Freq_M1 = zeros(1,param.sleep_blocks);
        SO_Freq_Cb = zeros(1,param.sleep_blocks);
        D_Freq_M1 = zeros(1,param.sleep_blocks);
        D_Freq_Cb = zeros(1,param.sleep_blocks);
        for block = 1:param.sleep_blocks
            % Read data
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Events_SOD_Sleep.mat']);
            
            nEpochs = length(UP_State_M1);
            tmp1_SO = zeros(1,nEpochs);
            tmp2_SO = zeros(1,nEpochs);
            tmp1_D = zeros(1,nEpochs);
            tmp2_D = zeros(1,nEpochs);
            for j = 1:nEpochs
                % Count Slow-oscillations
                ts = SO_timestamps_M1{j};
                if ~isempty(ts) && mod(length(ts),4)==0
                    ts = reshape(ts,4,length(ts)/4);
                end
                tmp1_SO(j) = size(ts,1);
                
                ts = SO_timestamps_Cb{j};
                if ~isempty(ts) && mod(length(ts),4)==0
                    ts = reshape(ts,4,length(ts)/4);
                end
                tmp2_SO(j) = size(ts,1);
                
                % Count Delta-oscillations
                ts = D_timestamps_M1{j};
                if ~isempty(ts) && mod(length(ts),4)==0
                    ts = reshape(ts,4,length(ts)/4);
                end
                tmp1_D(j) = size(ts,1);
                
                ts = D_timestamps_Cb{j};
                if ~isempty(ts) && mod(length(ts),4)==0
                    ts = reshape(ts,4,length(ts)/4);
                end
                tmp2_D(j) = size(ts,1);
                
            end
            
            SO_Freq_M1(block) = sum(tmp1_SO)/(nEpochs*param.epoch_length);
            SO_Freq_Cb(block) = sum(tmp2_SO)/(nEpochs*param.epoch_length);
            
            D_Freq_M1(block)  = sum(tmp1_D)/(nEpochs*param.epoch_length);
            D_Freq_Cb(block)  = sum(tmp2_D)/(nEpochs*param.epoch_length);
        end
        subplot(2,2,1);
        bar(SO_Freq_M1);%ylim([0 800]);
        title('M1');
        ylabel('Slow Oscillation Freq (hz)');
        subplot(2,2,2);
        bar(SO_Freq_Cb);%ylim([0 800]);
        title('Cb');
        subplot(2,2,3);
        bar(D_Freq_M1);%ylim([0 800]);
        xlabel('Block');
        ylabel('Delta Oscillation Freq (hz)');
        subplot(2,2,4);
        bar(D_Freq_Cb);%ylim([0 800]);
        xlabel('Block');
        
        saveas(gcf,[rootpath,animal,'/Day',num2str((day)),'/Oscillation_Frequencies.fig']);
        close all;
    end
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Get Nested spindle Frequency (18)

if enabled(18)
    disp('Block 18...')
    for day = 1:param.days
        nested_spindle_freq_M1 = zeros(1,param.sleep_blocks);
        nested_spindle_freq_Cb = zeros(1,param.sleep_blocks);
        nested_spindle_count_M1 = zeros(1,param.sleep_blocks);
        nested_spindle_count_Cb = zeros(1,param.sleep_blocks);
        for block = 1:param.sleep_blocks
            % Read data
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Nested_spindles_Sleep.mat']);
            
            % Unpack sindle cell arrays and count the spindles
            nEpochs = 1:size(D_nested_spindle_Cb,2);
            
            bfr = logical(~cellfun('isempty',D_nested_spindle_M1(nEpochs))) & logical(~cellfun('isempty',D_nested_spindle_Cb(nEpochs)));
            SO_nested_spindle_M1 = D_nested_spindle_M1(bfr);
            SO_nested_spindle_Cb = D_nested_spindle_Cb(bfr);
            count_M1 = zeros(size(SO_nested_spindle_M1));
            for j=1:size(SO_nested_spindle_M1,2)
                count_M1(j) = size(SO_nested_spindle_M1{j},2);
            end
            count_Cb = zeros(size(SO_nested_spindle_Cb));
            for j=1:size(SO_nested_spindle_Cb,2)
                count_Cb(j) = size(SO_nested_spindle_Cb{j},2);
            end
            nested_spindle_freq_M1(block) = sum(count_M1)/(length(nEpochs)*param.epoch_length);
            nested_spindle_freq_Cb(block) = sum(count_Cb)/(length(nEpochs)*param.epoch_length);
            
            nested_spindle_count_M1(block) = sum(count_M1);
            nested_spindle_count_Cb(block) = sum(count_Cb);
        end
        subplot(2,2,1);
        bar(nested_spindle_freq_M1);%ylim([0 800]);
        title('M1');
        ylabel('spindle Frequency (hz)');
        subplot(2,2,2);
        bar(nested_spindle_freq_Cb);%ylim([0 800]);
        title('Cb');
        subplot(2,2,3);
        bar(nested_spindle_count_M1);%ylim([0 800]);
        xlabel('Block');
        ylabel('spindle Count');
        subplot(2,2,4);
        bar(nested_spindle_count_Cb);%ylim([0 800]);
        xlabel('Block');
        
        saveas(gcf,[rootpath,animal,'/Day',num2str((day)),'/SO_Nested_spindle_Frequencies.fig']);
        close all;
    end
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Plot Block LFP with 'bad', 'REM/Wake', and 'NREM' identifiers (19)

if enabled(19)
    disp('Block 19...')
    day = 1;
    block = 1;
    area = 'M1';
    classification_value = .0015;
    channel = 1;
    
    load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/LFP_',area,'_Truncated.mat']);
    load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/bad_Epochs.mat']);
    load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Epoch_Classification.mat']);
    if strcmp(area,'M1')
        M1_LFP_t(param.M1_sleep_bad_chans,:) = []; 
        epoch_classification = zeros(1,floor(size(M1_LFP_t,2)/(param.epoch_length*param.M1_Fs)));
        epoch_classification(bad_epochs_M1) = 4;
        epoch_classification(setdiff(bad_epochs_Cb,bad_epochs_M1)) = 5;
        epoch_classification(epoch_classification==0) = idx;
        y_val = repmat(classification_value,1,length(epoch_classification));
        
        if param.sleep_C(1,1) < param.sleep_C(2,1)
            label1 = 'REM/wake';
            label2 = 'NREM';
        else
            label2 = 'REM/wake';
            label1 = 'NREM';
        end
        hold on
        scatter((find(epoch_classification == 1)+0.5)*(param.epoch_length*param.M1_Fs),y_val(1:(sum(epoch_classification == 1))));
        scatter((find(epoch_classification == 2)+0.5)*(param.epoch_length*param.M1_Fs),y_val(1:(sum(epoch_classification == 2))));
        scatter((find(epoch_classification == 3)+0.5)*(param.epoch_length*param.M1_Fs),y_val(1:(sum(epoch_classification == 3))));
        scatter((find(epoch_classification == 4)+0.5)*(param.epoch_length*param.M1_Fs),y_val(1:(sum(epoch_classification == 4))));
        scatter((find(epoch_classification == 5)+0.5)*(param.epoch_length*param.M1_Fs),y_val(1:(sum(epoch_classification == 5))));
        plot(M1_LFP_t(channel,:));
    elseif strcmp(area,'Cb')
        Cb_LFP_t(param.Cb_sleep_bad_chans,:) = [];
        epoch_classification = zeros(1,floor(size(Cb_LFP_t,2)/(param.epoch_length*param.Cb_Fs)));
        epoch_classification(bad_epochs_Cb) = 5;
        epoch_classification(setdiff(bad_epochs_M1,bad_epochs_Cb)) = 4;
        epoch_classification(epoch_classification==0) = idx;
        y_val = repmat(classification_value,1,length(epoch_classification));
        
        if param.sleep_C(1,1) < param.sleep_C(2,1)
            label1 = 'REM/wake';
            label2 = 'NREM';
        else
            label2 = 'REM/wake';
            label1 = 'NREM';
        end
        
        hold on
        scatter((find(epoch_classification == 1)+0.5)*(param.epoch_length*param.Cb_Fs),y_val(1:(sum(epoch_classification == 1))));
        scatter((find(epoch_classification == 2)+0.5)*(param.epoch_length*param.Cb_Fs),y_val(1:(sum(epoch_classification == 1))));
        scatter((find(epoch_classification == 3)+0.5)*(param.epoch_length*param.Cb_Fs),y_val(1:(sum(epoch_classification == 1))));
        scatter((find(epoch_classification == 4)+0.5)*(param.epoch_length*param.Cb_Fs),y_val(1:(sum(epoch_classification == 1))));
        scatter((find(epoch_classification == 5)+0.5)*(param.epoch_length*param.Cb_Fs),y_val(1:(sum(epoch_classification == 1))));
        plot(Cb_LFP_t(channel,:));
    end
    
    
    legend(label1,label2,'short_NREM','bad M1','bad Cb','LFP');
    hold off
    input('Press enter when done.')
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Rename Spike Timestamp Data (19.5)

if false
    addpath(genpath('Z:\BMI_Analysis'))
    addpath(genpath('Z:\Matlab Offline Files SDK'))
    
    param.M1_spike_wave_Fs = 24414;
    param.M1_neurons = 6;
    param.Cb_spike_wave_Fs = 24414;
    param.Cb_neurons = 6;
    param.M1_mua_neurons = cell(5,1);
    param.Cb_mua_neurons = cell(5,1);
    
    for day=1:param.days
        if strcmp(animal,'I060')
            param.M1_mua_neurons{day} = zeros(param.M1_chans,param.M1_neurons);
            param.Cb_mua_neurons{day} = zeros(param.Cb_chans,param.Cb_neurons);
            %Super simple, everything is already done
            for block=1:param.blocks
                load([origin_rootpath,animal,'/',param.LFP_path,'/',param.sleep_block_names{day, block},'/Timestamps_S',num2str((block)),'.mat']);
                load([origin_rootpath,animal,'/',param.LFP_path,'/',param.sleep_block_names{day, block},'/Waves_S',num2str((block)),'.mat']);
            end
            TimeStamps2_S1(setdiff(1:param.Cb_chans,1:4:param.Cb_chans),:) = cell(24,6); %#ok<SAGROW>
            TimeStamps2_S2(setdiff(1:param.Cb_chans,1:4:param.Cb_chans),:) = cell(24,6); %#ok<SAGROW>
        %elseif strcmp(animal,'I061')
            
        else
            error('Unrecognized animal ID');
        end
        
        %Rename and save
        var_idx = 1;
        for block = 1:length(param.s_block_names)
            M1_spike_timestamps = eval(['TimeStamps1_S', num2str(var_idx)]);
            M1_spike_waves = eval(['Waves1_S', num2str(var_idx)]);
            Cb_spike_timestamps = eval(['TimeStamps2_S', num2str(var_idx)]);
            Cb_spike_waves = eval(['Waves2_S', num2str(var_idx)]);
            
            save([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_timestamps.mat'], 'M1_spike_timestamps', 'Cb_spike_timestamps', 'M1_spike_waves', 'Cb_spike_waves','-v7.3');
            var_idx = var_idx+1;
        end
        
        
        clearvars -except rootpath origin_rootpath animal param enabled day;
    end
    save([rootpath,animal,'/Parameters.mat'],'param');
    rmpath(genpath('Z:\BMI_Analysis'))
    rmpath(genpath('Z:\Matlab Offline Files SDK'))
    clear day;
end


%% Filter LFPs into spindle band (20)
if enabled(20)
    disp('Block 20...')
    addpath(genpath('Z:\Matlab for analysis\eeglab\functions'))
    rmpath(genpath('Z:\Matlab for analysis\eeglab\functions\octavefunc\signal'))
    M1_spindle_band = [10 16]; %in hz
    Cb_spindle_band = [10 16]; %in hz
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_data.mat'])
            M1_filt_data = M1_data;
            Cb_filt_data = Cb_data;
            M1_filt_data.LFP = eegfilt(M1_data.LFP', M1_data.Fs_LFP, M1_spindle_band(1), M1_spindle_band(2));
            Cb_filt_data.LFP = eegfilt(Cb_data.LFP', Cb_data.Fs_LFP, Cb_spindle_band(1), Cb_spindle_band(2));
            
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Filtered_Sleep_data.mat'], 'M1_filt_data', 'Cb_filt_data')
        end
    end
    rmpath(genpath('Z:\Matlab for analysis\eeglab\functions'))
    clearvars -except rootpath origin_rootpath animal param enabled;
end
            
            

%% Spike-spike Cross correlation around sleep events (21)

if enabled(21)
    disp('Block 21...')
    pre_margin = 0.5; %in seconds
    post_margin = 0.5; %in seconds
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_data.mat'])
            load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_timestamps.mat']);

            delta_idxs = round(M1_data.Fs_LFP * M1_data.so_delta.delta_up_states);
            so_idxs = round(M1_data.Fs_LFP * M1_data.so_delta.so_up_states);
            spindle_idxs = round(M1_data.Fs_LFP * M1_data.spindles{1}.pks);
            
            %Calculate spike-spike correlation at delta waves
            start_idxs = M1_data.so_delta.delta_up_states - pre_margin;
            end_idxs = M1_data.so_delta.delta_up_states + post_margin;
            cross_corrs = nan(sum(param.M1_task_related_neurons{day}(:)), sum(param.Cb_task_related_neurons{day}(:)), length(start_idxs), ((pre_margin + post_margin)*2000)-1);
            %cross_corrs dims:M1 neuron, Cb neuron, event index, time lag bin(ms)
            M1_hist = nan(sum(param.M1_task_related_neurons{day}(:)),length(start_idxs),(pre_margin + post_margin)*1000);
            Cb_hist = nan(sum(param.Cb_task_related_neurons{day}(:)),length(start_idxs),(pre_margin + post_margin)*1000);
            %Compile M1 histograms for each Task-related neuron and event
            n_idx = 0;
            for chan = 1:param.M1_chans
                for neuron = 1:param.M1_neurons
                    if param.M1_task_related_neurons{day}(chan,neuron)
                        n_idx = n_idx+1;
                        for event = 1:length(start_idxs)
                            M1_hist(n_idx,event,:) = histcounts(M1_spike_timestamps{chan,neuron}(M1_spike_timestamps{chan,neuron} >= start_idxs(event) & M1_spike_timestamps{chan,neuron} <= end_idxs(event)),start_idxs(event):.001:end_idxs(event));
                        end
                    end
                end
            end
            %Compile Cb histograms for each Task-related neuron and event
            n_idx = 0;
            for chan = 1:param.Cb_chans
                for neuron = 1:param.Cb_neurons
                    if param.Cb_task_related_neurons{day}(chan,neuron)
                        n_idx = n_idx+1;
                        for event = 1:length(start_idxs)
                            Cb_hist(n_idx,event,:) = histcounts(Cb_spike_timestamps{chan,neuron}(Cb_spike_timestamps{chan,neuron} >= start_idxs(event) & Cb_spike_timestamps{chan,neuron} <= end_idxs(event)),start_idxs(event):.001:end_idxs(event));
                        end
                    end
                end
            end
            
            %Calculate each cross correlation
            for M1_neuron = 1:size(M1_hist,1)
                for Cb_neuron = 1:size(Cb_hist,1)
                    for event = 1:length(start_idxs)
                        cross_corrs(M1_neuron,Cb_neuron,event,:) = xcorr(squeeze(M1_hist(M1_neuron,event,:)), squeeze(Cb_hist(Cb_neuron,event,:)));
                    end
                end
            end
            delta_cross_corrs = cross_corrs;
            delta_M1_hist = M1_hist;
            delta_Cb_hist = Cb_hist;
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/delta_data.mat'], 'delta_cross_corrs', 'delta_M1_hist', 'delta_Cb_hist', '-v7.3')
            
            %Calculate spike-spike correlation at SO waves
            start_idxs = M1_data.so_delta.so_up_states - pre_margin;
            end_idxs = M1_data.so_delta.so_up_states + post_margin;
            cross_corrs = nan(sum(param.M1_task_related_neurons{day}(:)), sum(param.Cb_task_related_neurons{day}(:)), length(start_idxs), ((pre_margin + post_margin)*2000)-1);
            %cross_corrs dims:M1 neuron, Cb neuron, event index, time lag bin(ms)
            M1_hist = nan(sum(param.M1_task_related_neurons{day}(:)),length(start_idxs),(pre_margin + post_margin)*1000);
            Cb_hist = nan(sum(param.Cb_task_related_neurons{day}(:)),length(start_idxs),(pre_margin + post_margin)*1000);
            %Compile M1 histograms for each Task-related neuron and event
            n_idx = 0;
            for chan = 1:param.M1_chans
                for neuron = 1:param.M1_neurons
                    if param.M1_task_related_neurons{day}(chan,neuron)
                        n_idx = n_idx+1;
                        for event = 1:length(start_idxs)
                            M1_hist(n_idx,event,:) = histcounts(M1_spike_timestamps{chan,neuron}(M1_spike_timestamps{chan,neuron} >= start_idxs(event) & M1_spike_timestamps{chan,neuron} <= end_idxs(event)),start_idxs(event):.001:end_idxs(event));
                        end
                    end
                end
            end
            %Compile Cb histograms for each Task-related neuron and event
            n_idx = 0;
            for chan = 1:param.Cb_chans
                for neuron = 1:param.Cb_neurons
                    if param.Cb_task_related_neurons{day}(chan,neuron)
                        n_idx = n_idx+1;
                        for event = 1:length(start_idxs)
                            Cb_hist(n_idx,event,:) = histcounts(Cb_spike_timestamps{chan,neuron}(Cb_spike_timestamps{chan,neuron} >= start_idxs(event) & Cb_spike_timestamps{chan,neuron} <= end_idxs(event)),start_idxs(event):.001:end_idxs(event));
                        end
                    end
                end
            end
            %Calculate each cross correlation
            for M1_neuron = 1:size(M1_hist,1)
                for Cb_neuron = 1:size(Cb_hist,1)
                    for event = 1:length(start_idxs)
                        cross_corrs(M1_neuron,Cb_neuron,event,:) = xcorr(squeeze(M1_hist(M1_neuron,event,:)), squeeze(Cb_hist(Cb_neuron,event,:)));
                    end
                end
            end
            so_cross_corrs = cross_corrs;
            so_M1_hist = M1_hist;
            so_Cb_hist = Cb_hist;
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/SO_data.mat'], 'so_cross_corrs', 'so_M1_hist', 'so_Cb_hist', '-v7.3')
            
            %Calculate spike-spike correlation at spindles
            start_idxs = M1_data.spindles{1}.pks - pre_margin;
            end_idxs = M1_data.spindles{1}.pks + post_margin;
            cross_corrs = nan(sum(param.M1_task_related_neurons{day}(:)), sum(param.Cb_task_related_neurons{day}(:)), length(start_idxs), ((pre_margin + post_margin)*2000)-1);
            %cross_corrs dims:M1 neuron, Cb neuron, event index, time lag bin(ms)
            M1_hist = nan(sum(param.M1_task_related_neurons{day}(:)),length(start_idxs),(pre_margin + post_margin)*1000);
            Cb_hist = nan(sum(param.Cb_task_related_neurons{day}(:)),length(start_idxs),(pre_margin + post_margin)*1000);
            %Compile M1 histograms for each Task-related neuron and event
            n_idx = 0;
            for chan = 1:param.M1_chans
                for neuron = 1:param.M1_neurons
                    if param.M1_task_related_neurons{day}(chan,neuron)
                        n_idx = n_idx+1;
                        for event = 1:length(start_idxs)
                            M1_hist(n_idx,event,:) = histcounts(M1_spike_timestamps{chan,neuron}(M1_spike_timestamps{chan,neuron} >= start_idxs(event) & M1_spike_timestamps{chan,neuron} <= end_idxs(event)),start_idxs(event):.001:end_idxs(event));
                        end
                    end
                end
            end
            %Compile Cb histograms for each Task-related neuron and event
            n_idx = 0;
            for chan = 1:param.Cb_chans
                for neuron = 1:param.Cb_neurons
                    if param.Cb_task_related_neurons{day}(chan,neuron)
                        n_idx = n_idx+1;
                        for event = 1:length(start_idxs)
                            Cb_hist(n_idx,event,:) = histcounts(Cb_spike_timestamps{chan,neuron}(Cb_spike_timestamps{chan,neuron} >= start_idxs(event) & Cb_spike_timestamps{chan,neuron} <= end_idxs(event)),start_idxs(event):.001:end_idxs(event));
                        end
                    end
                end
            end
            %Calculate each cross correlation
            for M1_neuron = 1:size(M1_hist,1)
                for Cb_neuron = 1:size(Cb_hist,1)
                    for event = 1:length(start_idxs)
                        cross_corrs(M1_neuron,Cb_neuron,event,:) = xcorr(squeeze(M1_hist(M1_neuron,event,:)), squeeze(Cb_hist(Cb_neuron,event,:)));
                    end
                end
            end
            spindle_cross_corrs = cross_corrs;
            spindle_M1_hist = M1_hist;
            spindle_Cb_hist = Cb_hist;
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/spindle_data.mat'], 'spindle_cross_corrs', 'spindle_M1_hist', 'spindle_Cb_hist', '-v7.3')
            
            mean_delta_corr = squeeze(mean(mean(mean(delta_cross_corrs,1),2),3));
            mean_so_corr = squeeze(mean(mean(mean(so_cross_corrs,1),2),3));
            mean_spindle_corr = squeeze(mean(mean(mean(spindle_cross_corrs,1),2),3));
            
            ss_fit = fit((1:length(mean_delta_corr))', mean_delta_corr, 'smoothingspline', 'SmoothingParam', param.smoothing_param);
            smooth_delta_corr = ss_fit(1:length(mean_delta_corr));
            line(1:length(mean_delta_corr), mean_delta_corr, 'Color', [.7 .7 .7])
            line(1:length(mean_delta_corr), smooth_delta_corr, 'Color', [.4 .4 .4], 'LineWidth', 2)
            saveas(gcf,[rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/delta_cross_corr.fig'])
            close all
            
            ss_fit = fit((1:length(mean_so_corr))', mean_so_corr, 'smoothingspline', 'SmoothingParam', param.smoothing_param);
            smooth_so_corr = ss_fit(1:length(mean_so_corr));
            line(1:length(mean_so_corr), mean_so_corr, 'Color', [.7 .7 .7])
            line(1:length(mean_so_corr), smooth_so_corr, 'Color', [.4 .4 .4], 'LineWidth', 2)
            saveas(gcf,[rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/SO_cross_corr.fig'])
            close all
            
            ss_fit = fit((1:length(mean_spindle_corr))', mean_spindle_corr, 'smoothingspline', 'SmoothingParam', param.smoothing_param);
            smooth_spindle_corr = ss_fit(1:length(mean_spindle_corr));
            line(1:length(mean_spindle_corr), mean_spindle_corr, 'Color', [.7 .7 .7])
            line(1:length(mean_spindle_corr), smooth_spindle_corr, 'Color', [.4 .4 .4], 'LineWidth', 2)
            saveas(gcf,[rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/spindle_cross_corr.fig'])
            close all
        end
    end
    rmpath(genpath('Z:\Matlab for analysis\eeglab\functions'))
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Spike-phase counts around spindles (22)

if enabled(22)
    disp('Block 22...')
    addpath(genpath('Z:\Matlab for analysis\eeglab\functions'))
    rmpath(genpath('Z:\Matlab for analysis\eeglab\functions\octavefunc\signal'))
    cross_area_bool = false;
    max_idx_err = inf;
    pre_cycles = 5;
    post_cycles = 5;
    Cb_center_window = [0.3 0]; %pre/post (in seconds)
    Cb_ref_to_use = 'Cb_spin'; %Options: 'Cb_spin', 'M1_spin', 'nearest_Cb_peak'  -  
    control_offsets = [-5, -10];
    
    if cross_area_bool
        save_name_add = '_cross';
    else
        save_name_add = '';
    end
    
    if ~any(strcmp(Cb_ref_to_use, {'Cb_spin', 'M1_spin', 'nearest_Cb_peak'}))
        error('Unrecognized Cb reference option')
    end
    
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Filtered_Sleep_data.mat'])
            load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_timestamps.mat']);
            
            if cross_area_bool
                data_temp = M1_filt_data;
                M1_filt_data = Cb_filt_data;
                Cb_filt_data = data_temp;
                clear data_temp
            end
            
            for offset = [0, control_offsets]
                
                M1_spindle_idxs = round(M1_filt_data.Fs_LFP * (M1_filt_data.spindles{1}.pks + offset));
                
                hilbert_LFP = hilbert(M1_filt_data.LFP);
                inst_phase = unwrap(angle(hilbert_LFP));%inst phase
                M1_radial_phase = mod(inst_phase,2*pi);
                M1_peaks = diff(M1_radial_phase) < 0;
                
                %align spindle_indexes with a peak
                M1_peak_idxs = find(M1_peaks);
                %should spindles that are too close to the begining/end be skipped for all offsets?
                M1_spindle_idxs = M1_spindle_idxs(M1_spindle_idxs > M1_peak_idxs(1+pre_cycles));
                M1_spindle_idxs = M1_spindle_idxs(M1_spindle_idxs < M1_peak_idxs(end-post_cycles));
                for i = 1:length(M1_spindle_idxs)
                    [idx_diff, M1_spindle_idxs(i)] = min(abs(M1_peak_idxs - M1_spindle_idxs(i)));
                    if idx_diff > max_idx_err
                        warning('A spindle index had greater offset from the cycle peak than is allowed by the threshold');
                    end
                end
                
                if strcmp(Cb_ref_to_use, 'Cb_spin')
                    
                    Cb_spindle_idxs = round(Cb_filt_data.Fs_LFP * (Cb_filt_data.spindles{1}.pks + offset));
                    
                    hilbert_LFP = hilbert(Cb_filt_data.LFP);
                    inst_phase = unwrap(angle(hilbert_LFP));%inst phase
                    Cb_radial_phase = mod(inst_phase,2*pi);
                    Cb_peaks = diff(Cb_radial_phase) < 0;
                    
                    %align spindle_indexes with a peak
                    Cb_peak_idxs = find(Cb_peaks);
                    %should spindles that are too close to the begining/end be skipped for all offsets?
                    Cb_spindle_idxs = Cb_spindle_idxs(Cb_spindle_idxs > Cb_peak_idxs(1+pre_cycles));
                    Cb_spindle_idxs = Cb_spindle_idxs(Cb_spindle_idxs < Cb_peak_idxs(end-post_cycles));
                    for i = 1:length(Cb_spindle_idxs)
                        [idx_diff, Cb_spindle_idxs(i)] = min(abs(Cb_peak_idxs - Cb_spindle_idxs(i)));
                        if idx_diff > max_idx_err
                            warning('A spindle index had greater offset from the cycle peak than is allowed by the threshold');
                        end
                    end
                
                elseif strcmp(Cb_ref_to_use, 'M1_spin')
                    Cb_spindle_idxs = M1_spindle_idxs;
                    Cb_peak_idxs = M1_peak_idxs;
                
                elseif strcmp(Cb_ref_to_use, 'nearest_Cb_peak')
                    Cb_spindle_idxs = nan(size(M1_spindle_idxs));
                    
                    hilbert_LFP = hilbert(Cb_filt_data.LFP);
                    inst_phase = unwrap(angle(hilbert_LFP));%inst phase
                    Cb_radial_phase = mod(inst_phase,2*pi);
                    Cb_peaks = diff(Cb_radial_phase) < 0;
                    
                    Cb_peak_idxs = find(Cb_peaks);
                    
                    for i = 1:length(M1_spindle_idxs)
                        possible_Cb_centers = Cb_peak_idxs(Cb_peak_idxs >= round(((M1_peak_idxs(M1_spindle_idxs(i))/param.M1_Fs) - Cb_center_window(1)) * param.Cb_Fs) & Cb_peak_idxs <= round(((M1_peak_idxs(M1_spindle_idxs(i))/param.M1_Fs) + Cb_center_window(2)) * param.Cb_Fs));
                        [~, center_idx] = max(Cb_filt_data.LFP(possible_Cb_centers));
                        Cb_spindle_idxs(i) = center_idx;
                    end
                    Cb_spindle_idxs(isnan(Cb_spindle_idxs)) = [];
                end
                
                
                prePeak_cycle_phases_M1_spindles_M1_spikes = cell(param.M1_chans,param.M1_neurons,pre_cycles);
                postPeak_cycle_phases_M1_spindles_M1_spikes = cell(param.M1_chans,param.M1_neurons,post_cycles);
                prePeak_cycle_phases_Cb_spindles_Cb_spikes = cell(param.Cb_chans,param.Cb_neurons,pre_cycles);
                postPeak_cycle_phases_Cb_spindles_Cb_spikes = cell(param.Cb_chans,param.Cb_neurons,post_cycles);
                
                for M1_spin_idx = M1_spindle_idxs'
                    for chan = 1:param.M1_chans
                        for neuron = 1:param.M1_neurons
                            for cycle = 1:pre_cycles
                                cycle_start = M1_peak_idxs(M1_spin_idx - cycle)/param.M1_Fs;
                                cycle_end = M1_peak_idxs(M1_spin_idx - (cycle-1))/param.M1_Fs;
                                cycle_spikes = M1_spike_timestamps{chan,neuron}(M1_spike_timestamps{chan,neuron} >= cycle_start & M1_spike_timestamps{chan,neuron} <= cycle_end);
                                prePeak_cycle_phases_M1_spindles_M1_spikes{chan,neuron,cycle} = cat(2,prePeak_cycle_phases_M1_spindles_M1_spikes{chan,neuron,cycle},M1_radial_phase(round(cycle_spikes*param.M1_Fs)));
                            end
                            for cycle = 1:post_cycles
                                cycle_start = M1_peak_idxs(M1_spin_idx + (cycle-1))/param.M1_Fs;
                                cycle_end = M1_peak_idxs(M1_spin_idx + cycle)/param.M1_Fs;
                                cycle_spikes = M1_spike_timestamps{chan,neuron}(M1_spike_timestamps{chan,neuron} >= cycle_start & M1_spike_timestamps{chan,neuron} <= cycle_end);
                                postPeak_cycle_phases_M1_spindles_M1_spikes{chan,neuron,cycle} = cat(2,postPeak_cycle_phases_M1_spindles_M1_spikes{chan,neuron,cycle},M1_radial_phase(round(cycle_spikes*param.M1_Fs)));
                            end
                        end
                    end
                end
                
                for Cb_spin_idx = Cb_spindle_idxs'                   
                    for chan = 1:param.Cb_chans
                        for neuron = 1:param.Cb_neurons
                            for cycle = 1:pre_cycles
                                cycle_start = Cb_peak_idxs(Cb_spin_idx - cycle)/param.Cb_Fs;
                                cycle_end = Cb_peak_idxs(Cb_spin_idx - (cycle-1))/param.Cb_Fs;
                                cycle_spikes = Cb_spike_timestamps{chan,neuron}(Cb_spike_timestamps{chan,neuron} >= cycle_start & Cb_spike_timestamps{chan,neuron} <= cycle_end);
                                prePeak_cycle_phases_Cb_spindles_Cb_spikes{chan,neuron,cycle} = cat(2,prePeak_cycle_phases_Cb_spindles_Cb_spikes{chan,neuron,cycle},Cb_radial_phase(round(cycle_spikes*param.M1_Fs)));
                            end
                            for cycle = 1:post_cycles
                                cycle_start = Cb_peak_idxs(Cb_spin_idx + (cycle-1))/param.M1_Fs;
                                cycle_end = Cb_peak_idxs(Cb_spin_idx + cycle)/param.M1_Fs;
                                cycle_spikes = Cb_spike_timestamps{chan,neuron}(Cb_spike_timestamps{chan,neuron} >= cycle_start & Cb_spike_timestamps{chan,neuron} <= cycle_end);
                                postPeak_cycle_phases_Cb_spindles_Cb_spikes{chan,neuron,cycle} = cat(2,postPeak_cycle_phases_Cb_spindles_Cb_spikes{chan,neuron,cycle},Cb_radial_phase(round(cycle_spikes*param.M1_Fs)));
                            end
                        end
                    end
                end
                
                if offset ==0
                    save([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_spike_phasestamps', save_name_add, '.mat'], 'prePeak_cycle_phases_M1_spindles_M1_spikes', 'postPeak_cycle_phases_M1_spindles_M1_spikes', 'prePeak_cycle_phases_Cb_spindles_Cb_spikes', 'postPeak_cycle_phases_Cb_spindles_Cb_spikes')
                else
                    save([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_spike_phasestamps_C', num2str(offset), save_name_add, '.mat'], 'prePeak_cycle_phases_M1_spindles_M1_spikes', 'postPeak_cycle_phases_M1_spindles_M1_spikes', 'prePeak_cycle_phases_Cb_spindles_Cb_spikes', 'postPeak_cycle_phases_Cb_spindles_Cb_spikes')
                end
            end
        end
    end
    rmpath(genpath('Z:\Matlab for analysis\eeglab\functions'))
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Plot Spike-phase counts around spindles (23)

if enabled(23)
    disp('Block 23...')
    % for each neuron and each cycle (different days have different neurons, different blocks of the same day have the same neurons)
    %    One each of mean spike count (spike count sum / spindle count), prefered phase (phase of the averaged vectors), and phase locking value(length of the averaged vectors)
    
    addpath(genpath('Z:\Matlab for analysis'))
    cross_area_bool = false;
    
    if cross_area_bool
        save_name_add = '_cross';
    else
        save_name_add = '';
    end
    
    load([rootpath,animal,'/Day1/',param.s_block_names{1},'/Spindle_spike_phasestamps', save_name_add, '.mat'])
    load([rootpath,animal,'/Shared_Data.mat'])
    pre_cycles = size(prePeak_cycle_phases_M1_spindles_M1_spikes,3);
    post_cycles = size(postPeak_cycle_phases_M1_spindles_M1_spikes,3);
    
    if cross_area_bool
        shared_data.M1_spindle_cycle_spike_counts_cross = cell(1,param.days);
        shared_data.M1_spindle_cycle_prefered_phases_cross = cell(1,param.days);
        shared_data.M1_spindle_cycle_phase_locking_values_cross = cell(1,param.days);
        shared_data.Cb_spindle_cycle_spike_counts_cross = cell(1,param.days);
        shared_data.Cb_spindle_cycle_prefered_phases_cross = cell(1,param.days);
        shared_data.Cb_spindle_cycle_phase_locking_values_cross = cell(1,param.days);
    else
        shared_data.M1_spindle_cycle_spike_counts = cell(1,param.days);
        shared_data.M1_spindle_cycle_prefered_phases = cell(1,param.days);
        shared_data.M1_spindle_cycle_phase_locking_values = cell(1,param.days);
        shared_data.Cb_spindle_cycle_spike_counts = cell(1,param.days);
        shared_data.Cb_spindle_cycle_prefered_phases = cell(1,param.days);
        shared_data.Cb_spindle_cycle_phase_locking_values = cell(1,param.days);
    end
    
    
    M1_day_spindle_count = zeros(1,param.days);
    Cb_day_spindle_count = zeros(1,param.days);
    
    control_offsets = [-5, -10];
    
    control_M1_spindle_cycle_spike_counts = cell(1,param.days);
    control_M1_spindle_cycle_prefered_phases = cell(1,param.days);
    control_M1_spindle_cycle_phase_locking_values = cell(1,param.days);
    control_Cb_spindle_cycle_spike_counts = cell(1,param.days);
    control_Cb_spindle_cycle_prefered_phases = cell(1,param.days);
    control_Cb_spindle_cycle_phase_locking_values = cell(1,param.days);
    
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Sleep_data.mat'])
            M1_day_spindle_count(day) = M1_day_spindle_count(day) + length(M1_data.spindles{1}.pks);
            Cb_day_spindle_count(day) = Cb_day_spindle_count(day) + length(Cb_data.spindles{1}.pks);
        end
    end
    
    for offset = [0, control_offsets]
        for day = 1:param.days
            M1_neuron_idx = 0;
            Cb_neuron_idx = 0;
            M1_day_neurons = sum(param.M1_task_related_neurons{day}(:));
            Cb_day_neurons = sum(param.Cb_task_related_neurons{day}(:));
            if cross_area_bool
                shared_data.M1_spindle_cycle_spike_counts_cross{day} = nan(pre_cycles + post_cycles,M1_day_neurons);
                shared_data.M1_spindle_cycle_prefered_phases_cross{day} = nan(pre_cycles + post_cycles,M1_day_neurons);
                shared_data.M1_spindle_cycle_phase_locking_values_cross{day} = nan(pre_cycles + post_cycles,M1_day_neurons);
                shared_data.Cb_spindle_cycle_spike_counts_cross{day} = nan(pre_cycles + post_cycles,Cb_day_neurons);
                shared_data.Cb_spindle_cycle_prefered_phases_cross{day} = nan(pre_cycles + post_cycles,Cb_day_neurons);
                shared_data.Cb_spindle_cycle_phase_locking_values_cross{day} = nan(pre_cycles + post_cycles,Cb_day_neurons);
            else
                shared_data.M1_spindle_cycle_spike_counts{day} = nan(pre_cycles + post_cycles,M1_day_neurons);
                shared_data.M1_spindle_cycle_prefered_phases{day} = nan(pre_cycles + post_cycles,M1_day_neurons);
                shared_data.M1_spindle_cycle_phase_locking_values{day} = nan(pre_cycles + post_cycles,M1_day_neurons);
                shared_data.Cb_spindle_cycle_spike_counts{day} = nan(pre_cycles + post_cycles,Cb_day_neurons);
                shared_data.Cb_spindle_cycle_prefered_phases{day} = nan(pre_cycles + post_cycles,Cb_day_neurons);
                shared_data.Cb_spindle_cycle_phase_locking_values{day} = nan(pre_cycles + post_cycles,Cb_day_neurons);
            end
            
            M1_prePeak_spike_phases = cell(param.M1_chans,param.M1_neurons,pre_cycles);
            M1_postPeak_spike_phases = cell(param.M1_chans,param.M1_neurons,post_cycles);
            Cb_prePeak_spike_phases = cell(param.Cb_chans,param.Cb_neurons,pre_cycles);
            Cb_postPeak_spike_phases = cell(param.Cb_chans,param.Cb_neurons,post_cycles);
            
            for block = 1:param.sleep_blocks
                if offset == 0
                    load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_spike_phasestamps', save_name_add, '.mat'])
                else
                    load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_spike_phasestamps_C', num2str(offset), save_name_add, '.mat'])
                end
                M1_prePeak_spike_phases = cellfun(@horzcat,M1_prePeak_spike_phases,prePeak_cycle_phases_M1_spindles_M1_spikes,'UniformOutput',false);
                M1_postPeak_spike_phases = cellfun(@horzcat,M1_postPeak_spike_phases,postPeak_cycle_phases_M1_spindles_M1_spikes,'UniformOutput',false);
                Cb_prePeak_spike_phases = cellfun(@horzcat,Cb_prePeak_spike_phases,prePeak_cycle_phases_Cb_spindles_Cb_spikes,'UniformOutput',false);
                Cb_postPeak_spike_phases = cellfun(@horzcat,Cb_postPeak_spike_phases,postPeak_cycle_phases_Cb_spindles_Cb_spikes,'UniformOutput',false);
            end
            
            if M1_day_spindle_count(day) == 0
                if cross_area_bool
                    shared_data.M1_spindle_cycle_spike_counts_cross{day} = shared_data.M1_spindle_cycle_spike_counts_cross{day}(:,[]);
                    shared_data.M1_spindle_cycle_prefered_phases_cross{day} = shared_data.M1_spindle_cycle_prefered_phases_cross{day}(:,[]);
                    shared_data.M1_spindle_cycle_phase_locking_values_cross{day} = shared_data.M1_spindle_cycle_phase_locking_values_cross{day}(:,[]);
                else
                    shared_data.M1_spindle_cycle_spike_counts{day} = shared_data.M1_spindle_cycle_spike_counts{day}(:,[]);
                    shared_data.M1_spindle_cycle_prefered_phases{day} = shared_data.M1_spindle_cycle_prefered_phases{day}(:,[]);
                    shared_data.M1_spindle_cycle_phase_locking_values{day} = shared_data.M1_spindle_cycle_phase_locking_values{day}(:,[]);
                end
            else
                for chan = 1:param.M1_chans
                    for neuron = 1:param.M1_neurons
                        if param.M1_task_related_neurons{day}(chan,neuron)
                            M1_neuron_idx = M1_neuron_idx+1;
                            
                            if cross_area_bool
                                shared_data.M1_spindle_cycle_spike_counts_cross{day}(:,M1_neuron_idx) = cat(3,cellfun(@length,M1_prePeak_spike_phases(chan,neuron,end:-1:1)), cellfun(@length,M1_postPeak_spike_phases(chan,neuron,:)));
                                shared_data.M1_spindle_cycle_spike_counts_cross{day}(:,M1_neuron_idx) = shared_data.M1_spindle_cycle_spike_counts_cross{day}(:,M1_neuron_idx)/M1_day_spindle_count(day);
                            else
                                shared_data.M1_spindle_cycle_spike_counts{day}(:,M1_neuron_idx) = cat(3,cellfun(@length,M1_prePeak_spike_phases(chan,neuron,end:-1:1)), cellfun(@length,M1_postPeak_spike_phases(chan,neuron,:)));
                                shared_data.M1_spindle_cycle_spike_counts{day}(:,M1_neuron_idx) = shared_data.M1_spindle_cycle_spike_counts{day}(:,M1_neuron_idx)/M1_day_spindle_count(day);
                            end
                            
                            pre_sum_polar = nan(pre_cycles, 2);
                            for cycle = 1:pre_cycles
                                %Sum unit vectors to find prefered phase and phase locking for the cycles
                                [vect_x, vect_y] = pol2cart(M1_prePeak_spike_phases{chan,neuron,cycle},1);
                                [pre_sum_polar(cycle,1), pre_sum_polar(cycle,2)] = cart2pol(mean(vect_x),mean(vect_y));
                            end
                            
                            post_sum_polar = nan(post_cycles, 2);
                            for cycle = 1:post_cycles
                                %Sum unit vectors to find prefered phase and phase locking for the cycles
                                [vect_x, vect_y] = pol2cart(M1_postPeak_spike_phases{chan,neuron,cycle},1);
                                [post_sum_polar(cycle,1), post_sum_polar(cycle,2)] = cart2pol(mean(vect_x),mean(vect_y));
                            end
                            
                            if cross_area_bool
                                shared_data.M1_spindle_cycle_prefered_phases_cross{day}(:,M1_neuron_idx) = [pre_sum_polar(end:-1:1,1); post_sum_polar(:,1)];
                                shared_data.M1_spindle_cycle_phase_locking_values_cross{day}(:,M1_neuron_idx) = [pre_sum_polar(end:-1:1,2); post_sum_polar(:,2)];
                            else
                                shared_data.M1_spindle_cycle_prefered_phases{day}(:,M1_neuron_idx) = [pre_sum_polar(end:-1:1,1); post_sum_polar(:,1)];
                                shared_data.M1_spindle_cycle_phase_locking_values{day}(:,M1_neuron_idx) = [pre_sum_polar(end:-1:1,2); post_sum_polar(:,2)];
                            end
                        end
                    end
                end
            end
            
            if Cb_day_spindle_count(day) == 0
                if cross_area_bool
                    shared_data.Cb_spindle_cycle_spike_counts_cross{day} = shared_data.Cb_spindle_cycle_spike_counts_cross{day}(:,[]);
                    shared_data.Cb_spindle_cycle_prefered_phases_cross{day} = shared_data.Cb_spindle_cycle_prefered_phases_cross{day}(:,[]);
                    shared_data.Cb_spindle_cycle_phase_locking_values_cross{day} = shared_data.Cb_spindle_cycle_phase_locking_values_cross{day}(:,[]);
                else
                    shared_data.Cb_spindle_cycle_spike_counts{day} = shared_data.Cb_spindle_cycle_spike_counts{day}(:,[]);
                    shared_data.Cb_spindle_cycle_prefered_phases{day} = shared_data.Cb_spindle_cycle_prefered_phases{day}(:,[]);
                    shared_data.Cb_spindle_cycle_phase_locking_values{day} = shared_data.Cb_spindle_cycle_phase_locking_values{day}(:,[]);
                end
            else
                for chan = 1:param.Cb_chans
                    for neuron = 1:param.Cb_neurons
                        if param.Cb_task_related_neurons{day}(chan,neuron)
                            Cb_neuron_idx = Cb_neuron_idx+1;
                            
                            if cross_area_bool
                                shared_data.Cb_spindle_cycle_spike_counts_cross{day}(:,Cb_neuron_idx) = cat(3,cellfun(@length,Cb_prePeak_spike_phases(chan,neuron,end:-1:1)), cellfun(@length,Cb_postPeak_spike_phases(chan,neuron,:)));
                                shared_data.Cb_spindle_cycle_spike_counts_cross{day}(:,Cb_neuron_idx) = shared_data.Cb_spindle_cycle_spike_counts_cross{day}(:,Cb_neuron_idx)/Cb_day_spindle_count(day);
                            else
                                shared_data.Cb_spindle_cycle_spike_counts{day}(:,Cb_neuron_idx) = cat(3,cellfun(@length,Cb_prePeak_spike_phases(chan,neuron,end:-1:1)), cellfun(@length,Cb_postPeak_spike_phases(chan,neuron,:)));
                                shared_data.Cb_spindle_cycle_spike_counts{day}(:,Cb_neuron_idx) = shared_data.Cb_spindle_cycle_spike_counts{day}(:,Cb_neuron_idx)/Cb_day_spindle_count(day);
                            end
                            
                            pre_sum_polar = nan(pre_cycles, 2);
                            for cycle = 1:pre_cycles
                                %Sum unit vectors to find prefered phase and phase locking for the cycles
                                [vect_x, vect_y] = pol2cart(Cb_prePeak_spike_phases{chan,neuron,cycle},1);
                                [pre_sum_polar(cycle,1), pre_sum_polar(cycle,2)] = cart2pol(mean(vect_x),mean(vect_y));
                            end
                            
                            post_sum_polar = nan(post_cycles, 2);
                            for cycle = 1:post_cycles
                                %Sum unit vectors to find prefered phase and phase locking for the cycles
                                [vect_x, vect_y] = pol2cart(Cb_postPeak_spike_phases{chan,neuron,cycle},1);
                                [post_sum_polar(cycle,1), post_sum_polar(cycle,2)] = cart2pol(mean(vect_x),mean(vect_y));
                            end
                            
                            if cross_area_bool
                                shared_data.Cb_spindle_cycle_prefered_phases_cross{day}(:,Cb_neuron_idx) = [pre_sum_polar(end:-1:1,1); post_sum_polar(:,1)];
                                shared_data.Cb_spindle_cycle_phase_locking_values_cross{day}(:,Cb_neuron_idx) = [pre_sum_polar(end:-1:1,2); post_sum_polar(:,2)];
                            else
                                shared_data.Cb_spindle_cycle_prefered_phases{day}(:,Cb_neuron_idx) = [pre_sum_polar(end:-1:1,1); post_sum_polar(:,1)];
                                shared_data.Cb_spindle_cycle_phase_locking_values{day}(:,Cb_neuron_idx) = [pre_sum_polar(end:-1:1,2); post_sum_polar(:,2)];
                            end
                        end
                    end
                end
            end
        end
        %To make a graph for each day count the task-related neurons for each day and use that to divide up the neurons
        if offset == 0
            save([rootpath,animal,'/Shared_Data.mat'], 'shared_data');
        else
            if cross_area_bool
                control_M1_spindle_cycle_spike_counts = cellfun(@horzcat,control_M1_spindle_cycle_spike_counts,shared_data.M1_spindle_cycle_spike_counts_cross,'UniformOutput',false);
                control_M1_spindle_cycle_prefered_phases = cellfun(@horzcat,control_M1_spindle_cycle_prefered_phases,shared_data.M1_spindle_cycle_prefered_phases_cross,'UniformOutput',false);
                control_M1_spindle_cycle_phase_locking_values = cellfun(@horzcat,control_M1_spindle_cycle_phase_locking_values,shared_data.M1_spindle_cycle_phase_locking_values_cross,'UniformOutput',false);
                control_Cb_spindle_cycle_spike_counts = cellfun(@horzcat,control_Cb_spindle_cycle_spike_counts,shared_data.Cb_spindle_cycle_spike_counts_cross,'UniformOutput',false);
                control_Cb_spindle_cycle_prefered_phases = cellfun(@horzcat,control_Cb_spindle_cycle_prefered_phases,shared_data.Cb_spindle_cycle_prefered_phases_cross,'UniformOutput',false);
                control_Cb_spindle_cycle_phase_locking_values = cellfun(@horzcat,control_Cb_spindle_cycle_phase_locking_values,shared_data.Cb_spindle_cycle_phase_locking_values_cross,'UniformOutput',false);
                
            else
                control_M1_spindle_cycle_spike_counts = cellfun(@horzcat,control_M1_spindle_cycle_spike_counts,shared_data.M1_spindle_cycle_spike_counts,'UniformOutput',false);
                control_M1_spindle_cycle_prefered_phases = cellfun(@horzcat,control_M1_spindle_cycle_prefered_phases,shared_data.M1_spindle_cycle_prefered_phases,'UniformOutput',false);
                control_M1_spindle_cycle_phase_locking_values = cellfun(@horzcat,control_M1_spindle_cycle_phase_locking_values,shared_data.M1_spindle_cycle_phase_locking_values,'UniformOutput',false);
                control_Cb_spindle_cycle_spike_counts = cellfun(@horzcat,control_Cb_spindle_cycle_spike_counts,shared_data.Cb_spindle_cycle_spike_counts,'UniformOutput',false);
                control_Cb_spindle_cycle_prefered_phases = cellfun(@horzcat,control_Cb_spindle_cycle_prefered_phases,shared_data.Cb_spindle_cycle_prefered_phases,'UniformOutput',false);
                control_Cb_spindle_cycle_phase_locking_values = cellfun(@horzcat,control_Cb_spindle_cycle_phase_locking_values,shared_data.Cb_spindle_cycle_phase_locking_values,'UniformOutput',false);
            end
        end
    end
    x_idx = 1:(pre_cycles+post_cycles);
    load([rootpath,animal,'/Shared_Data.mat']) %Yes, I know this is a wierd, awkward way of doing this. Blame poor early design choices and an unwillingness to go back and redo everything.
    
    if cross_area_bool
        M1_spindle_cycle_spike_counts = cat(2,shared_data.M1_spindle_cycle_spike_counts_cross{:});
        M1_spindle_cycle_prefered_phases = cat(2,shared_data.M1_spindle_cycle_prefered_phases_cross{:});
        M1_spindle_cycle_phase_locking_values = cat(2,shared_data.M1_spindle_cycle_phase_locking_values_cross{:});
        Cb_spindle_cycle_spike_counts = cat(2,shared_data.Cb_spindle_cycle_spike_counts_cross{:});
        Cb_spindle_cycle_prefered_phases = cat(2,shared_data.Cb_spindle_cycle_prefered_phases_cross{:});
        Cb_spindle_cycle_phase_locking_values = cat(2,shared_data.Cb_spindle_cycle_phase_locking_values_cross{:});
        
        shared_data.c_M1_spindle_cycle_spike_counts_cross = control_M1_spindle_cycle_spike_counts;
        shared_data.c_M1_spindle_cycle_prefered_phases_cross = control_M1_spindle_cycle_prefered_phases;
        shared_data.c_M1_spindle_cycle_phase_locking_values_cross = control_M1_spindle_cycle_phase_locking_values;
        shared_data.c_Cb_spindle_cycle_spike_counts_cross = control_Cb_spindle_cycle_spike_counts;
        shared_data.c_Cb_spindle_cycle_prefered_phases_cross = control_Cb_spindle_cycle_prefered_phases;
        shared_data.c_Cb_spindle_cycle_phase_locking_values_cross = control_Cb_spindle_cycle_phase_locking_values;
    else
        M1_spindle_cycle_spike_counts = cat(2,shared_data.M1_spindle_cycle_spike_counts{:});
        M1_spindle_cycle_prefered_phases = cat(2,shared_data.M1_spindle_cycle_prefered_phases{:});
        M1_spindle_cycle_phase_locking_values = cat(2,shared_data.M1_spindle_cycle_phase_locking_values{:});
        Cb_spindle_cycle_spike_counts = cat(2,shared_data.Cb_spindle_cycle_spike_counts{:});
        Cb_spindle_cycle_prefered_phases = cat(2,shared_data.Cb_spindle_cycle_prefered_phases{:});
        Cb_spindle_cycle_phase_locking_values = cat(2,shared_data.Cb_spindle_cycle_phase_locking_values{:});
        
        shared_data.c_M1_spindle_cycle_spike_counts = control_M1_spindle_cycle_spike_counts;
        shared_data.c_M1_spindle_cycle_prefered_phases = control_M1_spindle_cycle_prefered_phases;
        shared_data.c_M1_spindle_cycle_phase_locking_values = control_M1_spindle_cycle_phase_locking_values;
        shared_data.c_Cb_spindle_cycle_spike_counts = control_Cb_spindle_cycle_spike_counts;
        shared_data.c_Cb_spindle_cycle_prefered_phases = control_Cb_spindle_cycle_prefered_phases;
        shared_data.c_Cb_spindle_cycle_phase_locking_values = control_Cb_spindle_cycle_phase_locking_values;
    end
    
    
    save([rootpath,animal,'/Shared_Data.mat'], 'shared_data');
    
    c_M1_spindle_cycle_spike_counts = cat(2,control_M1_spindle_cycle_spike_counts{:});
    c_M1_spindle_cycle_prefered_phases = cat(2,control_M1_spindle_cycle_prefered_phases{:});
    c_M1_spindle_cycle_phase_locking_values = cat(2,control_M1_spindle_cycle_phase_locking_values{:});
    c_Cb_spindle_cycle_spike_counts = cat(2,control_Cb_spindle_cycle_spike_counts{:});
    c_Cb_spindle_cycle_prefered_phases = cat(2,control_Cb_spindle_cycle_prefered_phases{:});
    c_Cb_spindle_cycle_phase_locking_values = cat(2,control_Cb_spindle_cycle_phase_locking_values{:});
    
    c_M1_spindle_cycle_spike_counts = c_M1_spindle_cycle_spike_counts([pre_cycles, pre_cycles+1],:);
    c_M1_spindle_cycle_prefered_phases = c_M1_spindle_cycle_prefered_phases([pre_cycles, pre_cycles+1],:);
    c_M1_spindle_cycle_phase_locking_values = c_M1_spindle_cycle_phase_locking_values([pre_cycles, pre_cycles+1],:);
    c_Cb_spindle_cycle_spike_counts = c_Cb_spindle_cycle_spike_counts([pre_cycles, pre_cycles+1],:);
    c_Cb_spindle_cycle_prefered_phases = c_Cb_spindle_cycle_prefered_phases([pre_cycles, pre_cycles+1],:);
    c_Cb_spindle_cycle_phase_locking_values = c_Cb_spindle_cycle_phase_locking_values([pre_cycles, pre_cycles+1],:);
    
    %plot M1 figures
    mean_spikes = mean(M1_spindle_cycle_spike_counts,2);
    err_spikes = std(M1_spindle_cycle_spike_counts,0,2)/sqrt(size(M1_spindle_cycle_spike_counts,2));
    perfered_phases = mean(M1_spindle_cycle_prefered_phases,2,'omitnan');
    err_phases = std(M1_spindle_cycle_prefered_phases,0,2,'omitnan')./sqrt(sum(~isnan(M1_spindle_cycle_prefered_phases),2));
    phase_locking_vals = mean(M1_spindle_cycle_phase_locking_values,2,'omitnan');
    err_PLV = std(M1_spindle_cycle_phase_locking_values,0,2,'omitnan')./sqrt(sum(~isnan(M1_spindle_cycle_phase_locking_values),2));
    
    c_mean_spikes = mean(c_M1_spindle_cycle_spike_counts(:));
    c_err_spikes = std(c_M1_spindle_cycle_spike_counts(:))/sqrt(length(c_M1_spindle_cycle_spike_counts(:)));
    c_perfered_phases = mean(c_M1_spindle_cycle_prefered_phases(:),'omitnan');
    c_err_phases = std(c_M1_spindle_cycle_prefered_phases(:),'omitnan')/sqrt(sum(~isnan(c_M1_spindle_cycle_prefered_phases(:))));
    c_phase_locking_vals = mean(c_M1_spindle_cycle_phase_locking_values(:),'omitnan');
    c_err_PLV = std(c_M1_spindle_cycle_phase_locking_values(:),'omitnan')/sqrt(sum(~isnan(c_M1_spindle_cycle_phase_locking_values(:))));
    
    hold on
    shadedErrorBar(x_idx,repmat(c_mean_spikes,[size(x_idx)]),repmat([(c_mean_spikes+c_err_spikes);(c_mean_spikes-c_err_spikes)],[size(x_idx)]));
    scatter(x_idx,mean_spikes,'filled')
    errorbar(x_idx,mean_spikes,err_spikes)
    xlim([(x_idx(1)-0.5), (x_idx(end)+0.5)])
    saveas(gcf, [rootpath,animal,'/Mean_M1_spikes_per_M1_cycle', save_name_add, '.fig'])
    hold off
    close all
    
    hold on
    shadedErrorBar(x_idx,repmat(c_perfered_phases,[size(x_idx)]),repmat([(c_perfered_phases+c_err_phases);(c_perfered_phases-c_err_phases)],[size(x_idx)]));
    scatter(x_idx,perfered_phases,'filled')
    errorbar(x_idx,perfered_phases,err_phases)
    xlim([(x_idx(1)-0.5), (x_idx(end)+0.5)])
    ylim([-pi pi])
    saveas(gcf, [rootpath,animal,'/M1_spikes_M1_cycle_prefered_phase', save_name_add, '.fig'])
    hold off
    close all
    
    hold on
    shadedErrorBar(x_idx,repmat(c_phase_locking_vals,[size(x_idx)]),repmat([(c_phase_locking_vals+c_err_PLV);(c_phase_locking_vals-c_err_PLV)],[size(x_idx)]));
    scatter(x_idx,phase_locking_vals,'filled')
    errorbar(x_idx,phase_locking_vals,err_PLV)
    xlim([(x_idx(1)-0.5), (x_idx(end)+0.5)])
    saveas(gcf, [rootpath,animal,'/M1_cycle_phase_M1_spike_locking_value', save_name_add, '.fig'])
    hold off
    close all
    
    %plot Cb figures
    mean_spikes = mean(Cb_spindle_cycle_spike_counts,2);
    err_spikes = std(Cb_spindle_cycle_spike_counts,0,2)/sqrt(size(Cb_spindle_cycle_spike_counts,2));
    perfered_phases = mean(Cb_spindle_cycle_prefered_phases,2,'omitnan');
    err_phases = std(Cb_spindle_cycle_prefered_phases,0,2,'omitnan')./sqrt(sum(~isnan(Cb_spindle_cycle_prefered_phases),2));
    phase_locking_vals = mean(Cb_spindle_cycle_phase_locking_values,2,'omitnan');
    err_PLV = std(Cb_spindle_cycle_phase_locking_values,0,2,'omitnan')./sqrt(sum(~isnan(Cb_spindle_cycle_phase_locking_values),2));
    
    c_mean_spikes = mean(c_Cb_spindle_cycle_spike_counts(:));
    c_err_spikes = std(c_Cb_spindle_cycle_spike_counts(:))/sqrt(length(c_Cb_spindle_cycle_spike_counts(:)));
    c_perfered_phases = mean(c_Cb_spindle_cycle_prefered_phases(:),'omitnan');
    c_err_phases = std(c_Cb_spindle_cycle_prefered_phases(:),'omitnan')/sqrt(sum(~isnan(c_Cb_spindle_cycle_prefered_phases(:))));
    c_phase_locking_vals = mean(c_Cb_spindle_cycle_phase_locking_values(:),'omitnan');
    c_err_PLV = std(c_Cb_spindle_cycle_phase_locking_values(:),'omitnan')/sqrt(sum(~isnan(c_Cb_spindle_cycle_phase_locking_values(:))));
    
    hold on
    shadedErrorBar(x_idx,repmat(c_mean_spikes,[size(x_idx)]),repmat([(c_mean_spikes+c_err_spikes);(c_mean_spikes-c_err_spikes)],[size(x_idx)]));
    scatter(x_idx,mean_spikes,'filled')
    errorbar(x_idx,mean_spikes,err_spikes)
    xlim([(x_idx(1)-0.5), (x_idx(end)+0.5)])
    saveas(gcf, [rootpath,animal,'/Mean_Cb_spikes_per_Cb_cycle', save_name_add, '.fig'])
    hold off
    close all
    
    hold on
    shadedErrorBar(x_idx,repmat(c_perfered_phases,[size(x_idx)]),repmat([(c_perfered_phases+c_err_phases);(c_perfered_phases-c_err_phases)],[size(x_idx)]));
    scatter(x_idx,perfered_phases,'filled')
    errorbar(x_idx,perfered_phases,err_phases)
    xlim([(x_idx(1)-0.5), (x_idx(end)+0.5)])
    ylim([-pi pi])
    saveas(gcf, [rootpath,animal,'/Cb_spikes_Cb_cycle_prefered_phase', save_name_add, '.fig'])
    hold off
    close all
    
    hold on
    shadedErrorBar(x_idx,repmat(c_phase_locking_vals,[size(x_idx)]),repmat([(c_phase_locking_vals+c_err_PLV);(c_phase_locking_vals-c_err_PLV)],[size(x_idx)]));
    scatter(x_idx,phase_locking_vals,'filled')
    errorbar(x_idx,phase_locking_vals,err_PLV)
    xlim([(x_idx(1)-0.5), (x_idx(end)+0.5)])
    saveas(gcf, [rootpath,animal,'/Cb_cycle_phase_Cb_spike_locking_value', save_name_add, '.fig'])
    hold off
    close all
    
    rmpath(genpath('Z:\Matlab for analysis'))
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Plot Cross-area LFP comparisons over an event-rich time period (24)

if enabled(24)
    disp('Block 24...')
    day = 1;
    block = 1;
    delta_so_band = [.1,4];
    spindle_band = [10,15];
    %for i = 1:1000
    start_time = 454; %in seconds %454 %1518 %1581 %(i*10)+1;
    end_time = 465;   %in seconds %465 %1540 %1636 %start_time+10;
    
    load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Sleep_data.mat'])
    
    M1_idxs = round(start_time * M1_data.Fs_LFP): round(end_time * M1_data.Fs_LFP);
    Cb_idxs = round(start_time * Cb_data.Fs_LFP): round(end_time * Cb_data.Fs_LFP);
    M1_x_vals = (M1_idxs./M1_data.Fs_LFP);
    Cb_x_vals = (Cb_idxs./Cb_data.Fs_LFP);
    
    if sum(M1_data.artifact_idx(M1_idxs)) > 0 || sum(Cb_data.artifact_idx(Cb_idxs)) > 0
        warning('Artifacts in selection')
        %continue 
    end
    if sum(M1_data.sleep_idx(M1_idxs)) == 0
        warning('Contains no sleep')
        %continue 
    end
    
    M1_LFP = M1_data.LFP;
    Cb_LFP = Cb_data.LFP;
    M1_LFP(M1_data.artifact_idx) = normrnd(mean(M1_LFP(~M1_data.artifact_idx)),std(M1_LFP(~M1_data.artifact_idx)),[sum(M1_data.artifact_idx) 1]);
    Cb_LFP(M1_data.artifact_idx) = normrnd(mean(Cb_LFP(~Cb_data.artifact_idx)),std(Cb_LFP(~Cb_data.artifact_idx)),[sum(Cb_data.artifact_idx) 1]);
    
    M1_LFP_norm = zscore(M1_LFP);
    Cb_LFP_norm = zscore(Cb_LFP);
    
    fpass = delta_so_band;
    % filter M1 for delta/slow wave
    [b,a] = butter(2,fpass(1)/(M1_data.Fs_LFP/2),'high');
    lfp_delta = filtfilt(b,a,M1_LFP);
    % lowpass
    [b,a] = butter(4,fpass(2)/(M1_data.Fs_LFP/2),'low');
    lfp_delta = filtfilt(b,a,lfp_delta);
    M1_so_filt_norm = zscore(lfp_delta);
    
    % filter Cb for delta/slow wave
    [b,a] = butter(2,fpass(1)/(Cb_data.Fs_LFP/2),'high');
    lfp_delta = filtfilt(b,a,Cb_LFP);
    % lowpass
    [b,a] = butter(4,fpass(2)/(Cb_data.Fs_LFP/2),'low');
    lfp_delta = filtfilt(b,a,lfp_delta);
    Cb_so_filt_norm = zscore(lfp_delta);
    
    fpass = spindle_band;
    % filter M1 for spindles
    [bhigh,ahigh] = butter(6,fpass(1)/(M1_data.Fs_LFP/2),'high');
    [blow,alow] = butter(8,fpass(2)/(M1_data.Fs_LFP/2),'low');
    flfp = filtfilt(bhigh,ahigh,M1_LFP);
    flfp = filtfilt(blow,alow,flfp);
    M1_spind_filt_norm = zscore(flfp);
    
    % filter Cb for spindles
    [bhigh,ahigh] = butter(6,fpass(1)/(Cb_data.Fs_LFP/2),'high');
    [blow,alow] = butter(8,fpass(2)/(Cb_data.Fs_LFP/2),'low');
    flfp = filtfilt(bhigh,ahigh,Cb_LFP);
    flfp = filtfilt(blow,alow,flfp);
    Cb_spind_filt_norm = zscore(flfp);
    
    M1_delta = M1_data.so_delta.delta_up_states((M1_data.so_delta.delta_up_states > start_time) & (M1_data.so_delta.delta_up_states < end_time));
    M1_so = M1_data.so_delta.so_up_states((M1_data.so_delta.so_up_states > start_time) & (M1_data.so_delta.so_up_states < end_time));
    M1_spin = M1_data.spindles{1}.pks((M1_data.spindles{1}.pks > start_time) & (M1_data.spindles{1}.pks < end_time));
    
    Cb_delta = Cb_data.so_delta.delta_up_states((Cb_data.so_delta.delta_up_states > start_time) & (Cb_data.so_delta.delta_up_states < end_time));
    Cb_so = Cb_data.so_delta.so_up_states((Cb_data.so_delta.so_up_states > start_time) & (Cb_data.so_delta.so_up_states < end_time));
    Cb_spin = Cb_data.spindles{1}.pks((Cb_data.spindles{1}.pks > start_time) & (Cb_data.spindles{1}.pks < end_time));

    plot(M1_x_vals,M1_LFP_norm(M1_idxs) + 10,'-k')
    hold on
    plot(M1_x_vals,M1_so_filt_norm(M1_idxs),'-k')
    plot(M1_x_vals,M1_spind_filt_norm(M1_idxs) - 10,'-k')
    for timestamp = M1_delta'
        if ~isempty(timestamp)
            plot([timestamp timestamp],[-15 15],'LineWidth',2,'Color',[1 0.5 0.5])
        end
    end
    for timestamp = M1_so'
        if ~isempty(timestamp)
            plot([timestamp timestamp],[-15 15],'LineWidth',2,'Color',[0.6 0 0])
        end
    end
    for timestamp = M1_spin'
        if ~isempty(timestamp)
            plot([timestamp timestamp],[-15 15],'LineWidth',2,'Color',[1 0.1 0.1])
        end
    end
    for timestamp = Cb_delta'
        if ~isempty(timestamp)
            %plot([timestamp timestamp],[-15 15],'LineWidth',2,'Color',[0.5 0.5 1])
        end
    end
    for timestamp = Cb_so'
        if ~isempty(timestamp)
            %plot([timestamp timestamp],[-15 15],'LineWidth',2,'Color',[0 0 0.6])
        end
    end
    for timestamp = Cb_spin'
        if ~isempty(timestamp)
            %plot([timestamp timestamp],[-15 15],'LineWidth',2,'Color',[0.1 0.1 1])
        end
    end
    axis([start_time-1, end_time+1, -20, 20])
    saveas(gcf,[rootpath,animal,'/LFP_rhythms_M1.fig'])
    %close all
    
    figure
    plot(Cb_x_vals,Cb_LFP_norm(M1_idxs) + 10,'-k')
    hold on
    plot(Cb_x_vals,Cb_so_filt_norm(M1_idxs),'-k')
    plot(Cb_x_vals,Cb_spind_filt_norm(M1_idxs) - 10,'-k')
    for timestamp = M1_delta'
        if ~isempty(timestamp)
            plot([timestamp timestamp],[-15 15],'LineWidth',2,'Color',[1 0.5 0.5])
        end
    end
    for timestamp = M1_so'
        if ~isempty(timestamp)
            plot([timestamp timestamp],[-15 15],'LineWidth',2,'Color',[0.6 0 0])
        end
    end
    for timestamp = M1_spin'
        if ~isempty(timestamp)
            plot([timestamp timestamp],[-15 15],'LineWidth',2,'Color',[1 0.1 0.1])
        end
    end
    for timestamp = Cb_delta'
        if ~isempty(timestamp)
            %plot([timestamp timestamp],[-15 15],'LineWidth',2,'Color',[0.5 0.5 1])
        end
    end
    for timestamp = Cb_so'
        if ~isempty(timestamp)
            %plot([timestamp timestamp],[-15 15],'LineWidth',2,'Color',[0 0 0.6])
        end
    end
    for timestamp = Cb_spin'
        if ~isempty(timestamp)
            %plot([timestamp timestamp],[-15 15],'LineWidth',2,'Color',[0.1 0.1 1])
        end
    end
    axis([start_time-1, end_time+1, -20, 20])
    saveas(gcf,[rootpath,animal,'/LFP_rhythms_Cb.fig'])
    close all
    %end
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Online and offline changes in LFP coherence and reach correlation (25)

if enabled(25)
    disp('Block 25...')
    addpath(genpath('Z:/Matlab for analysis/chronux_2_10/chronux/spectral_analysis'))
    
    do_full_coh_hist = false;
    
    coh_params.Fs = param.M1_Fs;
    coh_params.fpass = [0 20];
    coh_params.tapers = [10 19];
    coh_params.trialave = 0;
    coh_params.pad = 1;
    coh_params.err = [2 0.05];
    M1_example_chan = 1;
    Cb_example_chan = 7;
    hz_min = 4;
    hz_max = 8;
    
    example_pair_coh = nan(param.days, param.sleep_blocks,41);
    all_pair_cohs = nan(param.days, param.sleep_blocks, length(param.M1_good_chans), length(param.M1_good_chans), size(example_pair_coh,3));
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/LFP_M1_Truncated.mat']);
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/LFP_Cb_Truncated.mat']);
            load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Sleep_data.mat'])
            
            M1_LFP = M1_LFP_t(param.M1_sleep_good_chans,M1_data.sleep_idx);
            Cb_LFP = Cb_LFP_t(param.Cb_sleep_good_chans,Cb_data.sleep_idx);
            
            %subtract median at each timepoint
            M1_med = median(M1_LFP,1);
            Cb_med = median(Cb_LFP,1);
            M1_LFP = M1_LFP - M1_med;
            Cb_LFP = Cb_LFP - Cb_med;
            
            for M1_chan = 1:size(M1_LFP,1)
                for Cb_chan = 1:size(Cb_LFP,1)
                    if param.M1_sleep_good_chans(M1_chan) == M1_example_chan && param.Cb_sleep_good_chans(Cb_chan) == Cb_example_chan
                        [coh,phi_cmr,~,~,~,coh_times,coh_freqs,~,~,~] = cohgramc(squeeze(M1_LFP(M1_chan,:))', squeeze(Cb_LFP(Cb_chan,:))', [1 .025], coh_params);
                        coh_score = mean(coh,1);
                        all_pair_cohs(day,block,M1_chan,Cb_chan,:) = coh_score;
                        example_pair_coh(day,block,:) = coh_score;
                    elseif do_full_coh_hist
                        [coh,phi_cmr,~,~,~,coh_times,coh_freqs,~,~,~] = cohgramc(squeeze(M1_LFP(M1_chan,:))', squeeze(Cb_LFP(Cb_chan,:))', [1 .025], coh_params);
                        coh_score = mean(coh,1);
                        all_pair_cohs(day,block,M1_chan,Cb_chan,:) = coh_score;
                    end
                end
            end
        end
    end
    
    %Make Fig 2b
    for day = 1:param.days-1
        coh_pre = squeeze(example_pair_coh(day,1,:));
        coh_post_train = squeeze(example_pair_coh(day,2,:));
        coh_post_offline = squeeze(example_pair_coh(day+1,1,:));
        
        hold on
        plot(coh_freqs, coh_pre, 'Color', [.6 .6 .6], 'LineWidth', 1, 'LineStyle', '-')
        plot(coh_freqs, coh_post_train, 'Color', [.6 .6 .6], 'LineWidth', 1, 'LineStyle', '--')
        plot(coh_freqs, coh_post_offline, 'Color', [0 0 0], 'LineWidth', 1, 'LineStyle', '-')
        hold off
        saveas(gcf, [rootpath,animal,'/LFP_coherence_spectrum_',num2str(day),'.fig'])
        close all
    end
    
    %Make Fig 2c
    hold on
    plot_len = param.days*param.blocks;
    coh_means = mean(example_pair_coh(:,:,(coh_freqs <= hz_max) & (coh_freqs >= hz_min)),3);
    coh_means = reshape(coh_means', [1 plot_len]);
    
    for idx = 1:plot_len-1
        if mod(idx,param.blocks)
            line([idx idx+1],coh_means([idx idx+1]), 'Color', [.5 .5 .5], 'LineWidth', 2)
        end
    end
    for idx = param.blocks+1:param.blocks:plot_len
        line([idx-1 idx],coh_means([idx-1 idx]),'Color', [.7 .6 0], 'LineWidth', 2)
    end
    scatter(1:plot_len,coh_means,'Marker', '.', 'MarkerEdgeColor', [0 0 0], 'SizeData', 300)
    
    %Load reach velocity corrilation and add to plot
    
    hold off
    saveas(gcf, [rootpath,animal,'/LFP_coherence_by_day.fig'])
    close all
    
    %Make Fig 2d
    if do_full_coh_hist
        all_pair_cohs = all_pair_cohs(:,:,:,:,(coh_freqs <= hz_max) & (coh_freqs >= hz_min));
        online_coh_changes = all_pair_cohs(:,end,:,:,:) - all_pair_cohs(:,1,:,:,:);
        offline_coh_changes = all_pair_cohs(2:param.days,1,:,:,:) - all_pair_cohs(1:(param.days-1),end,:,:,:);
        
        online_coh_changes = squeeze(mean(mean(online_coh_changes(:,:,:,:,:),5) ,1));
        offline_coh_changes = squeeze(mean(mean(offline_coh_changes(:,:,:,:,:),5) ,1));
        
        save([rootpath,animal,'/online_offline_coh_changes.mat'], 'online_coh_changes', 'offline_coh_changes');
    end
        
    rmpath(genpath('Z:/Matlab for analysis/chronux_2_10/chronux/spectral_analysis'))
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Parallel online and offline changes in LFP coherence and reach corrilation (25.5)

if false
    if isunix  %#ok<UNRCH>
        addpath(genpath('/common/fleischerp/HPC_code'))
        addpath(genpath('/common/fleischerp/3rd_party_code'))
    else
        addpath(genpath('Z:/Matlab for analysis/chronux_2_10/chronux/spectral_analysis'))
    end
    
    worker_num = 24;
    do_full_coh_hist = true;
    
    coh_params.Fs = param.M1_Fs;
    coh_params.fpass = [0 20];
    coh_params.tapers = [10 19];
    coh_params.trialave = 0;
    coh_params.pad = 1;
    coh_params.err = [2 0.05];
    M1_example_chan = 1;
    Cb_example_chan = 7;
    hz_min = 4;
    hz_max = 8;
    
    example_pair_coh = nan(param.days, param.sleep_blocks,41);
    all_pair_cohs = nan(param.days, param.sleep_blocks, length(param.M1_good_chans), length(param.Cb_good_chans), size(example_pair_coh,3));
    parpool('local',worker_num);
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/LFP_M1_Truncated.mat']);
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/LFP_Cb_Truncated.mat']);
            load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Sleep_data.mat'])
            
            M1_LFP = M1_LFP_t(param.M1_sleep_good_chans,M1_data.sleep_idx);
            Cb_LFP = Cb_LFP_t(param.Cb_sleep_good_chans,Cb_data.sleep_idx);
            
            %subtract median at each timepoint
            M1_med = median(M1_LFP,1);
            Cb_med = median(Cb_LFP,1);
            M1_LFP = M1_LFP - M1_med;
            Cb_LFP = Cb_LFP - Cb_med;
            
            [~,~,~,~,~,coh_times,coh_freqs,~,~,~] = cohgramc(squeeze(M1_LFP(1,:))', squeeze(Cb_LFP(1,:))', [1 .025], coh_params);

            
            for M1_chan = 1:size(M1_LFP,1)
                parfor Cb_chan = 1:size(Cb_LFP,1)
                    if do_full_coh_hist
                        [coh,phi_cmr,~,~,~,~,~,~,~,~] = cohgramc(squeeze(M1_LFP(M1_chan,:))', squeeze(Cb_LFP(Cb_chan,:))', [1 .025], coh_params);
                        coh_score = mean(coh,1);
                        all_pair_cohs(day,block,M1_chan,Cb_chan,:) = coh_score;
                    end
                end
            end
        end
    end
    delete(gcp('nocreate'))
    
    %Make Fig 2d
    if do_full_coh_hist
        all_pair_cohs = all_pair_cohs(:,:,:,:,(coh_freqs <= hz_max) & (coh_freqs >= hz_min));
        online_coh_changes = all_pair_cohs(:,end,:,:,:) - all_pair_cohs(:,1,:,:,:);
        offline_coh_changes = all_pair_cohs(2:param.days,1,:,:,:) - all_pair_cohs(1:(param.days-1),end,:,:,:);
        
        online_coh_changes = squeeze(mean(mean(online_coh_changes(:,:,:,:,:),5) ,1));
        offline_coh_changes = squeeze(mean(mean(offline_coh_changes(:,:,:,:,:),5) ,1));
        
        save([rootpath,animal,'/online_offline_coh_changes.mat'], 'online_coh_changes', 'offline_coh_changes');
        %create and plot histograms and/or save data for multi-animal analysis
    end
        
    if isunix
        rmpath(genpath('/common/fleischerp/HPC_code'))
        rmpath(genpath('/common/fleischerp/3rd_party_code'))
    else
        rmpath(genpath('Z:/Matlab for analysis/chronux_2_10/chronux/spectral_analysis'))
    end
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Spike-spike Cross correlation within spindles using shuffle subtraction (26)

if enabled(26)
    disp('Block 26...')
    addpath(genpath('Z:\Matlab for analysis\eeglab\functions'))
    rmpath(genpath('Z:\Matlab for analysis\eeglab\functions\octavefunc\signal'))
    max_idx_err = 10;
    pre_cycles = 5;
    post_cycles = 5;
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Filtered_Sleep_data.mat'])
            load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_timestamps.mat']);
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/LFP_M1_Truncated.mat']);
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/LFP_Cb_Truncated.mat']);

            spindle_idxs = round(M1_data.Fs_LFP * M1_data.spindles{1}.pks);
            
            %get spindle phase data
            hilbert_LFP = hilbert(M1_filt_data.LFP);
            inst_phase = unwrap(angle(hilbert_LFP));%inst phase
            radial_phase = mod(inst_phase,2*pi);
            peaks = diff(radial_phase) < 0;
            %align spindle_indexes with a peak
            peak_idxs = find(peaks);
            for i = 1:length(spindle_idxs)
                [idx_diff, spindle_idxs(i)] = min(abs(peak_idxs - spindle_idxs(i)));
                if idx_diff > max_idx_err
                    warning('A spindle index had greater offset from the cycle peak than is allowed by the threshold');
                end
            end
            
            %TODO for spin_idx = spindle_idxs'
                        
            %Calculate spike-spike correlation at spindles
            %TODO start_idxs = ;
            %TODO end_idxs = ;
            cross_corrs = nan(sum(param.M1_task_related_neurons{day}(:)), sum(param.Cb_task_related_neurons{day}(:)), length(start_idxs), ((pre_margin + post_margin)*2000)-1);
            %cross_corrs dims:M1 neuron, Cb neuron, event index, time lag bin(ms)
            M1_hist = nan(sum(param.M1_task_related_neurons{day}(:)),length(start_idxs),(pre_margin + post_margin)*1000);
            Cb_hist = nan(sum(param.Cb_task_related_neurons{day}(:)),length(start_idxs),(pre_margin + post_margin)*1000);
            
            %Compile M1 histograms for each Task-related neuron and event
            n_idx = 0;
            for chan = 1:param.M1_chans
                for neuron = 1:param.M1_neurons
                    if param.M1_task_related_neurons{day}(chan,neuron)
                        n_idx = n_idx+1;
                        for event = 1:length(start_idxs)
                            M1_hist(n_idx,event,:) = histcounts(M1_spike_timestamps{chan,neuron}(M1_spike_timestamps{chan,neuron} >= start_idxs(event) & M1_spike_timestamps{chan,neuron} <= end_idxs(event)),start_idxs(event):.001:end_idxs(event));
                        end
                    end
                end
            end
            %Compile Cb histograms for each Task-related neuron and event
            n_idx = 0;
            for chan = 1:param.Cb_chans
                for neuron = 1:param.Cb_neurons
                    if param.Cb_task_related_neurons{day}(chan,neuron)
                        n_idx = n_idx+1;
                        for event = 1:length(start_idxs)
                            Cb_hist(n_idx,event,:) = histcounts(Cb_spike_timestamps{chan,neuron}(Cb_spike_timestamps{chan,neuron} >= start_idxs(event) & Cb_spike_timestamps{chan,neuron} <= end_idxs(event)),start_idxs(event):.001:end_idxs(event));
                            snapshot = Cb_spike_timestamps{chan,neuron}(Cb_spike_timestamps{chan,neuron} >= start_idxs(event) & Cb_spike_timestamps{chan,neuron} <= end_idxs(event));
                            shuffled_snapshot = nan(size(snapshot));
                            
                            %Shuffle spike cycles while holding cycle phase constant
                            if length(snapshot) > 1
                                disp('Here')
                            end
                            for spiketime = snapshot
                                %find phase that it occurs at
                                %select random cycle (1-5, +/-)
                                %find timestamp in new cycle with the phase closest to the spike's phase
                                %record timestamp as shuffled spike's timestamp
                            end
                            
                            %add shuffled window to shuffled histogram
                            Cb_hist_shuff(n_idx,event,:) = histcounts(shuffled_snapshot,start_idxs(event):.001:end_idxs(event));
                        end
                    end
                end
            end
            %Calculate each cross correlation
            for M1_neuron = 1:size(M1_hist,1)
                for Cb_neuron = 1:size(Cb_hist,1)
                    for event = 1:length(start_idxs)
                        cross_corrs(M1_neuron,Cb_neuron,event,:) = xcorr(squeeze(M1_hist(M1_neuron,event,:)), squeeze(Cb_hist(Cb_neuron,event,:)));
                    end
                end
            end
            spindle_cross_corrs = cross_corrs;
            spindle_M1_hist = M1_hist;
            spindle_Cb_hist = Cb_hist;
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/spindle_data.mat'], 'spindle_cross_corrs', 'spindle_M1_hist', 'spindle_Cb_hist', '-v7.3')
            
            mean_delta_corr = squeeze(mean(mean(mean(delta_cross_corrs,1),2),3));
            mean_so_corr = squeeze(mean(mean(mean(so_cross_corrs,1),2),3));
            mean_spindle_corr = squeeze(mean(mean(mean(spindle_cross_corrs,1),2),3));
            
            ss_fit = fit((1:length(mean_delta_corr))', mean_delta_corr, 'smoothingspline', 'SmoothingParam', param.smoothing_param);
            smooth_delta_corr = ss_fit(1:length(mean_delta_corr));
            line(1:length(mean_delta_corr), mean_delta_corr, 'Color', [.7 .7 .7])
            line(1:length(mean_delta_corr), smooth_delta_corr, 'Color', [.4 .4 .4], 'LineWidth', 2)
            saveas(gcf,[rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/delta_cross_corr.fig'])
            close all
            
            ss_fit = fit((1:length(mean_so_corr))', mean_so_corr, 'smoothingspline', 'SmoothingParam', param.smoothing_param);
            smooth_so_corr = ss_fit(1:length(mean_so_corr));
            line(1:length(mean_so_corr), mean_so_corr, 'Color', [.7 .7 .7])
            line(1:length(mean_so_corr), smooth_so_corr, 'Color', [.4 .4 .4], 'LineWidth', 2)
            saveas(gcf,[rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/SO_cross_corr.fig'])
            close all
            
            ss_fit = fit((1:length(mean_spindle_corr))', mean_spindle_corr, 'smoothingspline', 'SmoothingParam', param.smoothing_param);
            smooth_spindle_corr = ss_fit(1:length(mean_spindle_corr));
            line(1:length(mean_spindle_corr), mean_spindle_corr, 'Color', [.7 .7 .7])
            line(1:length(mean_spindle_corr), smooth_spindle_corr, 'Color', [.4 .4 .4], 'LineWidth', 2)
            saveas(gcf,[rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/spindle_cross_corr.fig'])
            close all
        end
    end
    rmpath(genpath('Z:\Matlab for analysis\eeglab\functions'))
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Normalized Spike-spike Cross correlation around sleep events using shuffled rhythms (27)

if enabled(27)
    disp('Block 27...')
    pre_margin = 0.5; %in seconds
    post_margin = 0.5; %in seconds
    shuffle_reps = 25;
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/delta_data.mat'])
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/SO_data.mat'])
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/spindle_data.mat'])
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_data.mat'])
            load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_timestamps.mat']);
            
            delta_idxs = round(M1_data.Fs_LFP * M1_data.so_delta.delta_up_states);
            so_idxs = round(M1_data.Fs_LFP * M1_data.so_delta.so_up_states);
            spindle_idxs = round(M1_data.Fs_LFP * M1_data.spindles{1}.pks);
            
            shuff_delta_mean = nan(shuffle_reps, ((pre_margin + post_margin)*2000)-1);
            shuff_so_mean = nan(shuffle_reps, ((pre_margin + post_margin)*2000)-1);
            shuff_spindle_mean = nan(shuffle_reps, ((pre_margin + post_margin)*2000)-1);
            for rep = 1:shuffle_reps
                
                shuff_delta_cross_corrs = nan(sum(param.M1_task_related_neurons{day}(:)), sum(param.Cb_task_related_neurons{day}(:)), length(M1_data.so_delta.delta_up_states), ((pre_margin + post_margin)*2000)-1);
                shuff_so_cross_corrs = nan(sum(param.M1_task_related_neurons{day}(:)), sum(param.Cb_task_related_neurons{day}(:)), length(M1_data.so_delta.so_up_states), ((pre_margin + post_margin)*2000)-1);
                shuff_spindle_cross_corrs = nan(sum(param.M1_task_related_neurons{day}(:)), sum(param.Cb_task_related_neurons{day}(:)), length(M1_data.spindles{1}.pks), ((pre_margin + post_margin)*2000)-1);
                
                %Calculate spike-spike correlation at delta waves
                shuffled_idx = randperm(length(M1_data.so_delta.delta_up_states));
                M1_start_idxs = M1_data.so_delta.delta_up_states - pre_margin;
                M1_end_idxs = M1_data.so_delta.delta_up_states + post_margin;
                Cb_start_idxs = M1_data.so_delta.delta_up_states(shuffled_idx) - pre_margin;
                Cb_end_idxs = M1_data.so_delta.delta_up_states(shuffled_idx) + post_margin;
                M1_hist = nan(sum(param.M1_task_related_neurons{day}(:)),length(M1_start_idxs),(pre_margin + post_margin)*1000);
                Cb_hist = nan(sum(param.Cb_task_related_neurons{day}(:)),length(Cb_start_idxs),(pre_margin + post_margin)*1000);
                %Compile M1 histograms for each Task-related neuron and event
                n_idx = 0;
                for chan = 1:param.M1_chans
                    for neuron = 1:param.M1_neurons
                        if param.M1_task_related_neurons{day}(chan,neuron)
                            n_idx = n_idx+1;
                            for event = 1:length(M1_start_idxs)
                                M1_hist(n_idx,event,:) = histcounts(M1_spike_timestamps{chan,neuron}(M1_spike_timestamps{chan,neuron} >= M1_start_idxs(event) & M1_spike_timestamps{chan,neuron} <= M1_end_idxs(event)),M1_start_idxs(event):.001:M1_end_idxs(event));
                            end
                        end
                    end
                end
                %Compile Cb histograms for each Task-related neuron and event
                n_idx = 0;
                for chan = 1:param.Cb_chans
                    for neuron = 1:param.Cb_neurons
                        if param.Cb_task_related_neurons{day}(chan,neuron)
                            n_idx = n_idx+1;
                            for event = 1:length(M1_start_idxs)
                                Cb_hist(n_idx,event,:) = histcounts(Cb_spike_timestamps{chan,neuron}(Cb_spike_timestamps{chan,neuron} >= Cb_start_idxs(event) & Cb_spike_timestamps{chan,neuron} <= Cb_end_idxs(event)),Cb_start_idxs(event):.001:Cb_end_idxs(event));
                            end
                        end
                    end
                end
                
                %Calculate each cross correlation
                for M1_neuron = 1:size(M1_hist,1)
                    for Cb_neuron = 1:size(Cb_hist,1)
                        for event = 1:length(M1_start_idxs)
                            shuff_delta_cross_corrs(M1_neuron,Cb_neuron,event,:) = xcorr(squeeze(M1_hist(M1_neuron,event,:)), squeeze(Cb_hist(Cb_neuron,event,:)));
                        end
                    end
                end
                
                
                %Calculate spike-spike correlation at SO waves
                shuffled_idx = randperm(length(M1_data.so_delta.so_up_states));
                M1_start_idxs = M1_data.so_delta.so_up_states - pre_margin;
                M1_end_idxs = M1_data.so_delta.so_up_states + post_margin;
                Cb_start_idxs = M1_data.so_delta.so_up_states(shuffled_idx) - pre_margin;
                Cb_end_idxs = M1_data.so_delta.so_up_states(shuffled_idx) + post_margin;
                M1_hist = nan(sum(param.M1_task_related_neurons{day}(:)),length(M1_start_idxs),(pre_margin + post_margin)*1000);
                Cb_hist = nan(sum(param.Cb_task_related_neurons{day}(:)),length(Cb_start_idxs),(pre_margin + post_margin)*1000);
                %Compile M1 histograms for each Task-related neuron and event
                n_idx = 0;
                for chan = 1:param.M1_chans
                    for neuron = 1:param.M1_neurons
                        if param.M1_task_related_neurons{day}(chan,neuron)
                            n_idx = n_idx+1;
                            for event = 1:length(M1_start_idxs)
                                M1_hist(n_idx,event,:) = histcounts(M1_spike_timestamps{chan,neuron}(M1_spike_timestamps{chan,neuron} >= M1_start_idxs(event) & M1_spike_timestamps{chan,neuron} <= M1_end_idxs(event)),M1_start_idxs(event):.001:M1_end_idxs(event));
                            end
                        end
                    end
                end
                %Compile Cb histograms for each Task-related neuron and event
                n_idx = 0;
                for chan = 1:param.Cb_chans
                    for neuron = 1:param.Cb_neurons
                        if param.Cb_task_related_neurons{day}(chan,neuron)
                            n_idx = n_idx+1;
                            for event = 1:length(M1_start_idxs)
                                Cb_hist(n_idx,event,:) = histcounts(Cb_spike_timestamps{chan,neuron}(Cb_spike_timestamps{chan,neuron} >= Cb_start_idxs(event) & Cb_spike_timestamps{chan,neuron} <= Cb_end_idxs(event)),Cb_start_idxs(event):.001:Cb_end_idxs(event));
                            end
                        end
                    end
                end
                
                %Calculate each cross correlation
                for M1_neuron = 1:size(M1_hist,1)
                    for Cb_neuron = 1:size(Cb_hist,1)
                        for event = 1:length(M1_start_idxs)
                            shuff_so_cross_corrs(M1_neuron,Cb_neuron,event,:) = xcorr(squeeze(M1_hist(M1_neuron,event,:)), squeeze(Cb_hist(Cb_neuron,event,:)));
                        end
                    end
                end

                
                %Calculate spike-spike correlation at SO waves
                shuffled_idx = randperm(length(M1_data.spindles{1}.pks));
                M1_start_idxs = M1_data.spindles{1}.pks - pre_margin;
                M1_end_idxs = M1_data.spindles{1}.pks + post_margin;
                Cb_start_idxs = M1_data.spindles{1}.pks(shuffled_idx) - pre_margin;
                Cb_end_idxs = M1_data.spindles{1}.pks(shuffled_idx) + post_margin;
                M1_hist = nan(sum(param.M1_task_related_neurons{day}(:)),length(M1_start_idxs),(pre_margin + post_margin)*1000);
                Cb_hist = nan(sum(param.Cb_task_related_neurons{day}(:)),length(Cb_start_idxs),(pre_margin + post_margin)*1000);
                %Compile M1 histograms for each Task-related neuron and event
                n_idx = 0;
                for chan = 1:param.M1_chans
                    for neuron = 1:param.M1_neurons
                        if param.M1_task_related_neurons{day}(chan,neuron)
                            n_idx = n_idx+1;
                            for event = 1:length(M1_start_idxs)
                                M1_hist(n_idx,event,:) = histcounts(M1_spike_timestamps{chan,neuron}(M1_spike_timestamps{chan,neuron} >= M1_start_idxs(event) & M1_spike_timestamps{chan,neuron} <= M1_end_idxs(event)),M1_start_idxs(event):.001:M1_end_idxs(event));
                            end
                        end
                    end
                end
                %Compile Cb histograms for each Task-related neuron and event
                n_idx = 0;
                for chan = 1:param.Cb_chans
                    for neuron = 1:param.Cb_neurons
                        if param.Cb_task_related_neurons{day}(chan,neuron)
                            n_idx = n_idx+1;
                            for event = 1:length(M1_start_idxs)
                                Cb_hist(n_idx,event,:) = histcounts(Cb_spike_timestamps{chan,neuron}(Cb_spike_timestamps{chan,neuron} >= Cb_start_idxs(event) & Cb_spike_timestamps{chan,neuron} <= Cb_end_idxs(event)),Cb_start_idxs(event):.001:Cb_end_idxs(event));
                            end
                        end
                    end
                end
                
                %Calculate each cross correlation
                for M1_neuron = 1:size(M1_hist,1)
                    for Cb_neuron = 1:size(Cb_hist,1)
                        for event = 1:length(M1_start_idxs)
                            shuff_spindle_cross_corrs(M1_neuron,Cb_neuron,event,:) = xcorr(squeeze(M1_hist(M1_neuron,event,:)), squeeze(Cb_hist(Cb_neuron,event,:)));
                        end
                    end
                end
                
                
                shuff_delta_mean(rep,:) = mean(mean(mean(shuff_delta_cross_corrs,1),2),3);
                shuff_so_mean(rep,:) = mean(mean(mean(shuff_so_cross_corrs,1),2),3);
                shuff_spindle_mean(rep,:) = mean(mean(mean(shuff_spindle_cross_corrs,1),2),3);
            end
            
            
            mean_delta_corr = squeeze(mean(mean(mean(delta_cross_corrs,1),2),3));
            mean_so_corr = squeeze(mean(mean(mean(so_cross_corrs,1),2),3));
            mean_spindle_corr = squeeze(mean(mean(mean(spindle_cross_corrs,1),2),3));
            
            mean_shuff_delta_corr = mean(shuff_delta_mean,1)';
            mean_shuff_so_corr = mean(shuff_so_mean,1)';
            mean_shuff_spindle_corr = mean(shuff_spindle_mean,1)';
            
            mean_delta_corr = mean_delta_corr - mean_shuff_delta_corr;
            mean_so_corr = mean_so_corr - mean_shuff_so_corr;
            mean_spindle_corr = mean_spindle_corr - mean_shuff_spindle_corr;
            
            ss_fit = fit((1:length(mean_delta_corr))', mean_delta_corr, 'smoothingspline', 'SmoothingParam', param.smoothing_param);
            smooth_delta_corr = ss_fit(1:length(mean_delta_corr));
            line(1:length(mean_delta_corr), mean_delta_corr, 'Color', [.7 .7 .7])
            line(1:length(mean_delta_corr), smooth_delta_corr, 'Color', [.4 .4 .4], 'LineWidth', 2)
            saveas(gcf,[rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/delta_cross_corr_norm.fig'])
            close all
            
            ss_fit = fit((1:length(mean_so_corr))', mean_so_corr, 'smoothingspline', 'SmoothingParam', param.smoothing_param);
            smooth_so_corr = ss_fit(1:length(mean_so_corr));
            line(1:length(mean_so_corr), mean_so_corr, 'Color', [.7 .7 .7])
            line(1:length(mean_so_corr), smooth_so_corr, 'Color', [.4 .4 .4], 'LineWidth', 2)
            saveas(gcf,[rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/SO_cross_corr_norm.fig'])
            close all
            
            ss_fit = fit((1:length(mean_spindle_corr))', mean_spindle_corr, 'smoothingspline', 'SmoothingParam', param.smoothing_param);
            smooth_spindle_corr = ss_fit(1:length(mean_spindle_corr));
            line(1:length(mean_spindle_corr), mean_spindle_corr, 'Color', [.7 .7 .7])
            line(1:length(mean_spindle_corr), smooth_spindle_corr, 'Color', [.4 .4 .4], 'LineWidth', 2)
            saveas(gcf,[rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/spindle_cross_corr_norm.fig'])
            close all
        end
    end
    rmpath(genpath('Z:\Matlab for analysis\eeglab\functions'))
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Normalized Spike-spike cross area cofiring in M1 spindle cycles using shuffling (28)


if enabled(28)
    disp('Block 28...')
    addpath(genpath('Z:\Matlab for analysis\eeglab\functions'))
    rmpath(genpath('Z:\Matlab for analysis\eeglab\functions\octavefunc\signal'))
    pre_margin = 0.05; %in seconds
    post_margin = 0.05; %in seconds
    pre_cycles = 5;
    post_cycles = 5;
    shuffle_reps = 25;
    max_idx_err = 10;
    kernel_stdv = 5;
    control_offsets = [-5, -10];
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Filtered_Sleep_data.mat'])
            load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_timestamps.mat']);
            
            for offset = [0, control_offsets]
                spindle_idxs = round(M1_filt_data.Fs_LFP *(M1_filt_data.spindles{1}.pks + offset));
                
                %get spindle phase data
                hilbert_LFP = hilbert(M1_filt_data.LFP);
                inst_phase = unwrap(angle(hilbert_LFP));%inst phase
                radial_phase = mod(inst_phase,2*pi);
                peaks = diff(radial_phase) < 0;
                %align spindle_indexes with a peak
                peak_idxs = find(peaks);
                for i = 1:length(spindle_idxs)
                    [idx_diff, spindle_idxs(i)] = min(abs(peak_idxs - spindle_idxs(i)));
                    if idx_diff > max_idx_err
                        warning('A spindle index had greater offset from a cycle peak than is allowed by the threshold');
                    end
                end
                
                %full_raw_CCH = nan((param.M1_chans * param.M1_neurons),(param.Cb_chans * param.Cb_neurons),(pre_cycles+post_cycles),round((pre_margin+post_margin)*1000)+1);
                raw_CCH = nan(sum(param.Cb_task_related_neurons{day}(:)) * sum(param.M1_task_related_neurons{day}(:)),(pre_cycles+post_cycles),round((pre_margin+post_margin)*1000)+1);
                pair_idx = 0;
                for M1_neuron = 1:(param.M1_chans * param.M1_neurons)
                    if param.M1_task_related_neurons{day}(M1_neuron)
                        for Cb_neuron = 1:(param.Cb_chans * param.Cb_neurons)
                            if param.Cb_task_related_neurons{day}(Cb_neuron)
                                pair_idx = pair_idx+1;
                                for cycle = (1-pre_cycles):post_cycles %zero-th cycle is the cycle imeadiatly preceding the center
                                    lag_hists = nan(0,round((pre_margin+post_margin)*1000)+1);
                                    for spindle = 1:length(spindle_idxs)
                                        
                                        cycle_start_idx = peak_idxs(spindle_idxs(spindle)+(cycle-1));
                                        cycle_end_idx = peak_idxs(spindle_idxs(spindle)+(cycle));
                                        
                                        M1_cycle_spikes = M1_spike_timestamps{M1_neuron}(M1_spike_timestamps{M1_neuron} > (cycle_start_idx/param.M1_Fs)  &  M1_spike_timestamps{M1_neuron} < (cycle_end_idx/param.M1_Fs));
                                        for M1_spiketime = M1_cycle_spikes
                                            Cb_spikes = Cb_spike_timestamps{Cb_neuron}(Cb_spike_timestamps{Cb_neuron} >= (M1_spiketime - pre_margin) & Cb_spike_timestamps{Cb_neuron} <= (M1_spiketime + post_margin));
                                            Cb_spikes = Cb_spikes - M1_spiketime; %It should be Cb - M1 for both
                                            lag_hists = cat(1,lag_hists,histcounts(Cb_spikes,(-(round(pre_margin*1000)+0.5)):(round(post_margin*1000)+0.5)));
                                        end
                                        
                                        Cb_cycle_spikes = Cb_spike_timestamps{Cb_neuron}(Cb_spike_timestamps{Cb_neuron} > (cycle_start_idx/param.Cb_Fs)  &  Cb_spike_timestamps{Cb_neuron} < (cycle_end_idx/param.Cb_Fs));
                                        for Cb_spiketime = Cb_cycle_spikes
                                            M1_spikes = M1_spike_timestamps{M1_neuron}(M1_spike_timestamps{M1_neuron} >= (Cb_spiketime - pre_margin) & M1_spike_timestamps{M1_neuron} <= (Cb_spiketime + post_margin));
                                            M1_spikes = Cb_spiketime - M1_spikes; %It should be Cb - M1 for both
                                            lag_hists = cat(1,lag_hists,histcounts(M1_spikes,(-(round(pre_margin*1000)+0.5)):(round(post_margin*1000)+0.5)));
                                        end
                                    end
                                    %average hist arrays across spindles and store it
                                    if isempty(lag_hists)
                                        %full_raw_CCH(M1_neuron,Cb_neuron,(cycle + pre_cycles),:) = 0;
                                        raw_CCH(pair_idx,(cycle + pre_cycles),:) = 0;
                                    else
                                        %full_raw_CCH(M1_neuron,Cb_neuron,(cycle + pre_cycles),:) = mean(lag_hists,1);
                                        raw_CCH(pair_idx,(cycle + pre_cycles),:) = mean(lag_hists,1);
                                    end
                                    
                                end
                            end
                        end
                    end
                end
                
                
                %full_shuffled_CCH = nan((param.M1_chans * param.M1_neurons),(param.Cb_chans * param.Cb_neurons),(pre_cycles+post_cycles),shuffle_reps,round((pre_margin+post_margin)*1000)+1);
                shuffled_CCH = nan(sum(param.Cb_task_related_neurons{day}(:)) * sum(param.M1_task_related_neurons{day}(:)),(pre_cycles+post_cycles),shuffle_reps,round((pre_margin+post_margin)*1000)+1);
                shuffled_spindle_idxs = nan(shuffle_reps,length(spindle_idxs));
                for shuffle = 1:shuffle_reps
                    shuffled_spindle_idxs(shuffle,:) = randperm(length(spindle_idxs));
                end
                pair_idx = 0;
                for M1_neuron = 1:(param.M1_chans * param.M1_neurons)
                    if param.M1_task_related_neurons{day}(M1_neuron)
                        for Cb_neuron = 1:(param.Cb_chans * param.Cb_neurons)
                            if param.Cb_task_related_neurons{day}(Cb_neuron)
                                pair_idx = pair_idx+1;
                                for cycle = (1-pre_cycles):post_cycles %zero-th cycle is the cycle imeadiatly preceding the center
                                    for shuffle = 1:shuffle_reps
                                        lag_hists = nan(0,round((pre_margin+post_margin)*1000)+1);
                                        for spindle = 1:length(spindle_idxs)
                                            
                                            cycle_start_idx = peak_idxs(spindle_idxs(spindle)+(cycle-1));
                                            cycle_end_idx = peak_idxs(spindle_idxs(spindle)+(cycle));
                                            
                                            
                                            M1_cycle_spikes = M1_spike_timestamps{M1_neuron}(M1_spike_timestamps{M1_neuron} > (cycle_start_idx/param.M1_Fs)  &  M1_spike_timestamps{M1_neuron} < (cycle_end_idx/param.M1_Fs));
                                            M1_cycle_spikes = M1_cycle_spikes + (peak_idxs(spindle_idxs(shuffled_spindle_idxs(shuffle,spindle)))/param.M1_Fs) - (peak_idxs(spindle_idxs(spindle))/param.M1_Fs);
                                            for M1_spiketime = M1_cycle_spikes
                                                Cb_spikes = Cb_spike_timestamps{Cb_neuron}(Cb_spike_timestamps{Cb_neuron} >= (M1_spiketime - pre_margin) & Cb_spike_timestamps{Cb_neuron} <= (M1_spiketime + post_margin));
                                                Cb_spikes = Cb_spikes - M1_spiketime; %It should be Cb - M1 for both
                                                lag_hists = cat(1,lag_hists,histcounts(Cb_spikes,(-(round(pre_margin*1000)+0.5)):(round(post_margin*1000)+0.5)));
                                            end
                                            
                                            Cb_cycle_spikes = Cb_spike_timestamps{Cb_neuron}(Cb_spike_timestamps{Cb_neuron} > (cycle_start_idx/param.Cb_Fs)  &  Cb_spike_timestamps{Cb_neuron} < (cycle_end_idx/param.Cb_Fs));
                                            Cb_cycle_spikes = Cb_cycle_spikes + (peak_idxs(spindle_idxs(shuffled_spindle_idxs(shuffle,spindle)))/param.Cb_Fs) - (peak_idxs(spindle_idxs(spindle))/param.Cb_Fs);
                                            for Cb_spiketime = Cb_cycle_spikes
                                                M1_spikes = M1_spike_timestamps{M1_neuron}(M1_spike_timestamps{M1_neuron} >= (Cb_spiketime - pre_margin) & M1_spike_timestamps{M1_neuron} <= (Cb_spiketime + post_margin));
                                                M1_spikes = Cb_spiketime - M1_spikes; %It should be Cb - M1 for both
                                                lag_hists = cat(1,lag_hists,histcounts(M1_spikes,(-(round(pre_margin*1000)+0.5)):(round(post_margin*1000)+0.5)));
                                            end
                                        end
                                        %average hist array across spindles and store it
                                        if isempty(lag_hists)
                                            %full_shuffled_CCH(M1_neuron,Cb_neuron,(cycle + pre_cycles),shuffle,:) = 0;
                                            shuffled_CCH(pair_idx,(cycle + pre_cycles),shuffle,:) = 0;
                                        else
                                            %full_shuffled_CCH(M1_neuron,Cb_neuron,(cycle + pre_cycles),shuffle,:) = mean(lag_hists,1);
                                            shuffled_CCH(pair_idx,(cycle + pre_cycles),shuffle,:) = mean(lag_hists,1);
                                        end
                                        
                                    end
                                end
                            end
                        end
                    end
                end
                
                
                ave_shuffled_CCH = squeeze(mean(shuffled_CCH,3));
                corrected_CCH = raw_CCH-ave_shuffled_CCH;
                kernel = gausswin(size(corrected_CCH,3),size(corrected_CCH,3)/kernel_stdv);
                smoothed_CCH = nan(size(corrected_CCH));
                for i = 1:size(corrected_CCH,1)
                    for j = 1:size(corrected_CCH,2)
                        smoothed_CCH(i,j,:) = conv(squeeze(corrected_CCH(i,j,:)),kernel,'same');
                    end
                end
                [peak_val,peak_idx] = max(smoothed_CCH,[],3);
                peak_time = peak_idx - (round(pre_margin * 1000) + 1);
                pair_peaks_mean = mean(peak_val,1);
                pair_peaks_stdv = std(peak_val,0,1);
                pair_times_mean = mean(peak_time,1);
                pair_times_stdv = std(peak_time,0,1);
                pair_peaks_err = pair_peaks_stdv/sqrt(size(peak_time,1));
                pair_times_err = pair_times_stdv/sqrt(size(peak_time,1));
                
                x_axis = [1-pre_cycles:post_cycles];
                hold on
                errorbar(x_axis,pair_peaks_mean,pair_peaks_err,pair_peaks_err,zeros(1,pre_cycles+post_cycles),zeros(1,pre_cycles+post_cycles),'.')
                scatter(x_axis,pair_peaks_mean,'Marker', 's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 1], 'SizeData', 100)
                xlim([-4.5 5.5])
                hold off
                if offset == 0
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_prob.fig'])
                else
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_prob_C', num2str(offset), '.fig'])
                end
                close all
                
                hold on
                errorbar(x_axis,pair_times_mean,pair_times_err,pair_times_err,zeros(1,pre_cycles+post_cycles),zeros(1,pre_cycles+post_cycles),'.')
                scatter(x_axis,pair_times_mean,'Marker', 's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 1], 'SizeData', 100)
                xlim([-4.5 5.5])
                hold off
                if offset == 0
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_peak_times.fig'])
                else
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_peak_times_C', num2str(offset), '.fig'])
                end
                close all
            end
        end
    end
    rmpath(genpath('Z:\Matlab for analysis\eeglab\functions'))
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% LFP amplitude and neuron firing rate around rhythms (29)

if enabled(29)
    disp('Block 29...')
    delta_so_band = [.1,4];
    spindle_band = [10,15];
    pre_margin = 0.5; %in seconds
    post_margin = 0.5; %in seconds
    lfp_pre = round(pre_margin*param.M1_Fs);
    lfp_post = round(post_margin*param.M1_Fs);
    bin_width = .01; %in seconds
    task_related_neurons_only = false;
    max_to_plot = 200;
    
    control_offsets = [-5, -10];
    
    x_axis_LFP = round(-pre_margin*param.M1_Fs):round(post_margin*param.M1_Fs);
    x_axis_LFP = x_axis_LFP*1000/param.M1_Fs;
    x_axis_spike = round(-pre_margin*1000):round(bin_width*1000):round(post_margin*1000);
    x_axis_spike = x_axis_spike(1:(end-1)) + round(bin_width*500);
    
    M1_delta_snapshots = nan(0, length(x_axis_LFP));
    M1_SO_snapshots = nan(0, length(x_axis_LFP));
    M1_spind_snapshots = nan(0, length(x_axis_LFP));
    Cb_delta_snapshots = nan(0, length(x_axis_LFP));
    Cb_SO_snapshots = nan(0, length(x_axis_LFP));
    Cb_spind_snapshots = nan(0, length(x_axis_LFP));
    
    for offset = [0, control_offsets]
        for day = 1:param.days
            for block = 1:param.sleep_blocks
                load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_data.mat'])
                load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_timestamps.mat']);
                
                if task_related_neurons_only
                    M1_neurons = param.M1_task_related_neurons{day};
                    Cb_neurons = param.Cb_task_related_neurons{day};
                else
                    M1_neurons = ~cellfun(@isempty,M1_spike_timestamps);
                    Cb_neurons = ~cellfun(@isempty,Cb_spike_timestamps);
                end
                
                M1_LFP = M1_data.LFP;
                Cb_LFP = Cb_data.LFP;
                M1_LFP(M1_data.artifact_idx) = normrnd(mean(M1_LFP(~M1_data.artifact_idx)),std(M1_LFP(~M1_data.artifact_idx)),[sum(M1_data.artifact_idx) 1]);
                Cb_LFP(Cb_data.artifact_idx) = normrnd(mean(Cb_LFP(~Cb_data.artifact_idx)),std(Cb_LFP(~Cb_data.artifact_idx)),[sum(Cb_data.artifact_idx) 1]);
                
                fpass = delta_so_band;
                % filter M1 for delta/slow wave
                [b,a] = butter(2,fpass(1)/(M1_data.Fs_LFP/2),'high');
                lfp_delta = filtfilt(b,a,M1_LFP);
                % lowpass
                [b,a] = butter(4,fpass(2)/(M1_data.Fs_LFP/2),'low');
                lfp_delta = filtfilt(b,a,lfp_delta);
                M1_so_filt_norm = zscore(lfp_delta);
                
                % filter Cb for delta/slow wave
                [b,a] = butter(2,fpass(1)/(Cb_data.Fs_LFP/2),'high');
                lfp_delta = filtfilt(b,a,Cb_LFP);
                % lowpass
                [b,a] = butter(4,fpass(2)/(Cb_data.Fs_LFP/2),'low');
                lfp_delta = filtfilt(b,a,lfp_delta);
                Cb_so_filt_norm = zscore(lfp_delta);
                
                fpass = spindle_band;
                % filter M1 for spindles
                [bhigh,ahigh] = butter(6,fpass(1)/(M1_data.Fs_LFP/2),'high');
                [blow,alow] = butter(8,fpass(2)/(M1_data.Fs_LFP/2),'low');
                flfp = filtfilt(bhigh,ahigh,M1_LFP);
                flfp = filtfilt(blow,alow,flfp);
                M1_spind_filt_norm = zscore(flfp);
                
                % filter Cb for spindles
                [bhigh,ahigh] = butter(6,fpass(1)/(Cb_data.Fs_LFP/2),'high');
                [blow,alow] = butter(8,fpass(2)/(Cb_data.Fs_LFP/2),'low');
                flfp = filtfilt(bhigh,ahigh,Cb_LFP);
                flfp = filtfilt(blow,alow,flfp);
                Cb_spind_filt_norm = zscore(flfp);
                
                %Delta
                all_rhythms_M1 = M1_data.so_delta.delta_up_states + offset;
                all_rhythms_M1 = all_rhythms_M1(all_rhythms_M1 > (0+pre_margin));
                all_rhythms_M1 = all_rhythms_M1(all_rhythms_M1 < ((length(M1_LFP)/M1_data.Fs_LFP)-post_margin));
                all_rhythms_Cb = Cb_data.so_delta.delta_up_states + offset;
                all_rhythms_Cb = all_rhythms_Cb(all_rhythms_Cb > (0+pre_margin));
                all_rhythms_Cb = all_rhythms_Cb(all_rhythms_Cb < ((length(Cb_LFP)/Cb_data.Fs_LFP)-post_margin));
                M1_LFP = M1_so_filt_norm;
                Cb_LFP = Cb_so_filt_norm;
                M1_rhythm_snapshots = nan(length(all_rhythms_M1), round(pre_margin*param.M1_Fs)+round(post_margin*param.M1_Fs)+1);
                Cb_rhythm_snapshots = nan(length(all_rhythms_Cb), round(pre_margin*param.Cb_Fs)+round(post_margin*param.Cb_Fs)+1);
                M1_spike_hist_sum = zeros(1, round((pre_margin + post_margin)/bin_width));
                Cb_spike_hist_sum = zeros(1, round((pre_margin + post_margin)/bin_width));
                for rhythm = 1:length(all_rhythms_M1)
                    center = round(M1_data.Fs_LFP * all_rhythms_M1(rhythm));
                    M1_rhythm_snapshots(rhythm,:) = M1_LFP((center-round(pre_margin*param.M1_Fs)):(center+round(post_margin*param.M1_Fs)));
                    for spike_times_cell = M1_spike_timestamps(M1_neurons)'
                        all_spike_times = spike_times_cell{1};
                        spike_times = all_spike_times(all_spike_times >= ((center/param.M1_Fs) - pre_margin) & all_spike_times <= ((center/param.M1_Fs) + post_margin));
                        spike_times = spike_times - (center/param.M1_Fs);
                        cur_hist = histcounts(spike_times, -pre_margin:bin_width:post_margin);
                        M1_spike_hist_sum = M1_spike_hist_sum + cur_hist;
                    end
                end
                
                
                for rhythm = 1:length(all_rhythms_Cb)
                    center = round(Cb_data.Fs_LFP * all_rhythms_Cb(rhythm));
                    Cb_rhythm_snapshots(rhythm,:) = Cb_LFP((center-round(pre_margin*param.Cb_Fs)):(center+round(post_margin*param.Cb_Fs)));
                    for spike_times_cell =  Cb_spike_timestamps(Cb_neurons)' %Does this work?
                        all_spike_times = spike_times_cell{1};
                        spike_times = all_spike_times(all_spike_times >= ((center/param.Cb_Fs) - pre_margin) & all_spike_times <= ((center/param.Cb_Fs) + post_margin));
                        spike_times = spike_times - (center/param.Cb_Fs);
                        cur_hist = histcounts(spike_times, -pre_margin:bin_width:post_margin);
                        Cb_spike_hist_sum = Cb_spike_hist_sum + cur_hist;
                    end
                end
                %average LFP snapshots
                M1_mean_LFP = mean(M1_rhythm_snapshots,1);
                Cb_mean_LFP = mean(Cb_rhythm_snapshots,1);
                
                M1_delta_snapshots = cat(1, M1_delta_snapshots, M1_rhythm_snapshots);
                Cb_delta_snapshots = cat(1, Cb_delta_snapshots, Cb_rhythm_snapshots);
                
                %multiply summed histograms by 1000 (to get spikes per sec), then divide by sum(param.??_task_related_neurons{day}(:)) * number of rythms (to get mean)
                M1_mean_spike_hist = M1_spike_hist_sum / (length(all_rhythms_M1) * sum(M1_neurons(:)) * bin_width);
                Cb_mean_spike_hist = Cb_spike_hist_sum / (length(all_rhythms_Cb) * sum(Cb_neurons(:)) * bin_width);
                
                %plot and save
                plot(x_axis_LFP, M1_mean_LFP, 'Color', [0.6 0 0], 'LineWidth', 1);
                hold on
                plot(x_axis_LFP, Cb_mean_LFP, 'Color', [0 0 0.6], 'LineWidth', 1);
                hold off
                xlim([-500 500])
                if offset == 0
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/LFP_mean_amp_delta.fig'])
                else
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/LFP_mean_amp_delta_C', num2str(offset), '.fig'])
                end
                close all
                
                plot(x_axis_spike, M1_mean_spike_hist, 'Color', [0.6 0 0], 'LineWidth', 1);
                hold on
                plot(x_axis_spike, Cb_mean_spike_hist, 'Color', [0 0 0.6], 'LineWidth', 1);
                hold off
                xlim([-500 500])
                if offset == 0
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_rate_delta.fig'])
                else
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_rate_delta_C', num2str(offset), '.fig'])
                end
                close all
                
                %Slow Oscillations
                all_rhythms_M1 = M1_data.so_delta.so_up_states + offset;
                all_rhythms_M1 = all_rhythms_M1(all_rhythms_M1 > (0+pre_margin));
                all_rhythms_M1 = all_rhythms_M1(all_rhythms_M1 < ((length(M1_LFP)/M1_data.Fs_LFP)-post_margin));
                all_rhythms_Cb = Cb_data.so_delta.so_up_states + offset;
                all_rhythms_Cb = all_rhythms_Cb(all_rhythms_Cb > (0+pre_margin));
                all_rhythms_Cb = all_rhythms_Cb(all_rhythms_Cb < ((length(Cb_LFP)/Cb_data.Fs_LFP)-post_margin));
                M1_LFP = M1_so_filt_norm;
                Cb_LFP = Cb_so_filt_norm;
                M1_rhythm_snapshots = nan(length(all_rhythms_M1), round(pre_margin*param.M1_Fs)+round(post_margin*param.M1_Fs)+1);
                Cb_rhythm_snapshots = nan(length(all_rhythms_Cb), round(pre_margin*param.Cb_Fs)+round(post_margin*param.Cb_Fs)+1);
                M1_spike_hist_sum = zeros(1, round((pre_margin + post_margin)/bin_width));
                Cb_spike_hist_sum = zeros(1, round((pre_margin + post_margin)/bin_width));
                for rhythm = 1:length(all_rhythms_M1)
                    center = round(M1_data.Fs_LFP * all_rhythms_M1(rhythm));
                    M1_rhythm_snapshots(rhythm,:) = M1_LFP((center-round(pre_margin*param.M1_Fs)):(center+round(post_margin*param.M1_Fs)));
                    for spike_times_cell = M1_spike_timestamps(M1_neurons)'
                        all_spike_times = spike_times_cell{1};
                        spike_times = all_spike_times(all_spike_times >= ((center/param.M1_Fs) - pre_margin) & all_spike_times <= ((center/param.M1_Fs) + post_margin));
                        spike_times = spike_times - (center/param.M1_Fs);
                        cur_hist = histcounts(spike_times, -pre_margin:bin_width:post_margin);
                        M1_spike_hist_sum = M1_spike_hist_sum + cur_hist;
                    end
                end
                for rhythm = 1:length(all_rhythms_Cb)
                    center = round(Cb_data.Fs_LFP * all_rhythms_Cb(rhythm));
                    Cb_rhythm_snapshots(rhythm,:) = Cb_LFP((center-round(pre_margin*param.Cb_Fs)):(center+round(post_margin*param.Cb_Fs)));
                    for spike_times_cell =  Cb_spike_timestamps(Cb_neurons)' %Does this work?
                        all_spike_times = spike_times_cell{1};
                        spike_times = all_spike_times(all_spike_times >= ((center/param.Cb_Fs) - pre_margin) & all_spike_times <= ((center/param.Cb_Fs) + post_margin));
                        spike_times = spike_times - (center/param.Cb_Fs);
                        cur_hist = histcounts(spike_times, -pre_margin:bin_width:post_margin);
                        Cb_spike_hist_sum = Cb_spike_hist_sum + cur_hist;
                    end
                end
                %average LFP snapshots
                M1_mean_LFP = mean(M1_rhythm_snapshots,1);
                Cb_mean_LFP = mean(Cb_rhythm_snapshots,1);
                
                M1_SO_snapshots = cat(1, M1_SO_snapshots, M1_rhythm_snapshots);
                Cb_SO_snapshots = cat(1, Cb_SO_snapshots, Cb_rhythm_snapshots);
                
                %multiply summed histograms by 1000 (to get spikes per sec), then divide by sum(param.??_task_related_neurons{day}(:)) * number of rythms (to get mean)
                M1_mean_spike_hist = M1_spike_hist_sum / (length(all_rhythms_M1) * sum(M1_neurons(:)) * bin_width);
                Cb_mean_spike_hist = Cb_spike_hist_sum / (length(all_rhythms_Cb) * sum(Cb_neurons(:)) * bin_width);
                
                %plot and save
                plot(x_axis_LFP, M1_mean_LFP, 'Color', [0.6 0 0], 'LineWidth', 1);
                hold on
                plot(x_axis_LFP, Cb_mean_LFP, 'Color', [0 0 0.6], 'LineWidth', 1);
                hold off
                xlim([-500 500])
                if offset == 0
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/LFP_mean_amp_SO.fig'])
                else
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/LFP_mean_amp_SO_C', num2str(offset), '.fig'])
                end
                close all
                
                plot(x_axis_spike, M1_mean_spike_hist, 'Color', [0.6 0 0], 'LineWidth', 1);
                hold on
                plot(x_axis_spike, Cb_mean_spike_hist, 'Color', [0 0 0.6], 'LineWidth', 1);
                hold off
                xlim([-500 500])
                if offset == 0
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_rate_SO.fig'])
                else
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_rate_SO_C', num2str(offset), '.fig'])
                end
                close all
                
                %Spindles
                all_rhythms_M1 = M1_data.spindles{1}.pks + offset;
                all_rhythms_M1 = all_rhythms_M1(all_rhythms_M1 > (0+pre_margin));
                all_rhythms_M1 = all_rhythms_M1(all_rhythms_M1 < ((length(M1_LFP)/M1_data.Fs_LFP)-post_margin));
                all_rhythms_Cb = Cb_data.spindles{1}.pks + offset;
                all_rhythms_Cb = all_rhythms_Cb(all_rhythms_Cb > (0+pre_margin));
                all_rhythms_Cb = all_rhythms_Cb(all_rhythms_Cb < ((length(Cb_LFP)/Cb_data.Fs_LFP)-post_margin));
                M1_LFP = M1_spind_filt_norm;
                Cb_LFP = Cb_spind_filt_norm;
                M1_rhythm_snapshots = nan(length(all_rhythms_M1), round(pre_margin*param.M1_Fs)+round(post_margin*param.M1_Fs)+1);
                Cb_rhythm_snapshots = nan(length(all_rhythms_Cb), round(pre_margin*param.Cb_Fs)+round(post_margin*param.Cb_Fs)+1);
                M1_spike_hist_sum = zeros(1, round((pre_margin + post_margin)/bin_width));
                Cb_spike_hist_sum = zeros(1, round((pre_margin + post_margin)/bin_width));
                for rhythm = 1:length(all_rhythms_M1)
                    center = round(M1_data.Fs_LFP * all_rhythms_M1(rhythm));
                    M1_rhythm_snapshots(rhythm,:) = M1_LFP((center-round(pre_margin*param.M1_Fs)):(center+round(post_margin*param.M1_Fs)));
                    for spike_times_cell = M1_spike_timestamps(M1_neurons)'
                        all_spike_times = spike_times_cell{1};
                        spike_times = all_spike_times(all_spike_times >= ((center/param.M1_Fs) - pre_margin) & all_spike_times <= ((center/param.M1_Fs) + post_margin));
                        spike_times = spike_times - (center/param.M1_Fs);
                        cur_hist = histcounts(spike_times, -pre_margin:bin_width:post_margin);
                        M1_spike_hist_sum = M1_spike_hist_sum + cur_hist;
                    end
                end
                for rhythm = 1:length(all_rhythms_Cb)
                    center = round(Cb_data.Fs_LFP * all_rhythms_Cb(rhythm));
                    Cb_rhythm_snapshots(rhythm,:) = Cb_LFP((center-round(pre_margin*param.Cb_Fs)):(center+round(post_margin*param.Cb_Fs)));
                    for spike_times_cell =  Cb_spike_timestamps(Cb_neurons)' %Does this work?
                        all_spike_times = spike_times_cell{1};
                        spike_times = all_spike_times(all_spike_times >= ((center/param.Cb_Fs) - pre_margin) & all_spike_times <= ((center/param.Cb_Fs) + post_margin));
                        spike_times = spike_times - (center/param.Cb_Fs);
                        cur_hist = histcounts(spike_times, -pre_margin:bin_width:post_margin);
                        Cb_spike_hist_sum = Cb_spike_hist_sum + cur_hist;
                    end
                end
                %average LFP snapshots
                M1_mean_LFP = mean(M1_rhythm_snapshots,1);
                Cb_mean_LFP = mean(Cb_rhythm_snapshots,1);
                
                M1_spind_snapshots = cat(1, M1_spind_snapshots, M1_rhythm_snapshots);
                Cb_spind_snapshots = cat(1, Cb_spind_snapshots, Cb_rhythm_snapshots);
                
                %multiply summed histograms by 1000 (to get spikes per sec), then divide by sum(param.??_task_related_neurons{day}(:)) * number of rythms (to get mean)
                M1_mean_spike_hist = M1_spike_hist_sum / (length(all_rhythms_M1) * sum(M1_neurons(:)) * bin_width);
                Cb_mean_spike_hist = Cb_spike_hist_sum / (length(all_rhythms_Cb) * sum(Cb_neurons(:)) * bin_width);
                
                %plot and save
                plot(x_axis_LFP, M1_mean_LFP, 'Color', [0.6 0 0], 'LineWidth', 1);
                hold on
                plot(x_axis_LFP, Cb_mean_LFP, 'Color', [0 0 0.6], 'LineWidth', 1);
                hold off
                xlim([-500 500])
                if offset == 0
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/LFP_mean_amp_spindle.fig'])
                else
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/LFP_mean_amp_spindle_C', num2str(offset), '.fig'])
                end
                close all
                
                plot(x_axis_spike, M1_mean_spike_hist, 'Color', [0.6 0 0], 'LineWidth', 1);
                hold on
                plot(x_axis_spike, Cb_mean_spike_hist, 'Color', [0 0 0.6], 'LineWidth', 1);
                hold off
                xlim([-500 500])
                if offset == 0
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_rate_spindle.fig'])
                else
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_rate_spindle_C', num2str(offset), '.fig'])
                end
                close all
            end
        end
        
        rand_idxs = randperm(size(M1_delta_snapshots,1));
        if size(M1_delta_snapshots,1) > max_to_plot
            rand_idxs = rand_idxs(1:max_to_plot);
        end
        plot(x_axis_LFP, M1_delta_snapshots(rand_idxs,:), 'Color', [0.6 0 0], 'LineWidth', 1);
        hold on
        rand_idxs = randperm(size(Cb_delta_snapshots,1));
        if size(Cb_delta_snapshots,1) > max_to_plot
            rand_idxs = rand_idxs(1:max_to_plot);
        end
        plot(x_axis_LFP, Cb_delta_snapshots(rand_idxs,:), 'Color', [0 0 0.6], 'LineWidth', 1);
        hold off
        xlim([-500 500])
        if offset == 0
            saveas(gcf,[rootpath,animal,'/LFP_mean_amp_delta.fig'])
        else
            saveas(gcf,[rootpath,animal,'/LFP_mean_amp_delta_C', num2str(offset), '.fig'])
        end
        close all
        
        rand_idxs = randperm(size(M1_SO_snapshots,1));
        if size(M1_SO_snapshots,1) > max_to_plot
            rand_idxs = rand_idxs(1:max_to_plot);
        end
        plot(x_axis_LFP, M1_SO_snapshots(rand_idxs,:), 'Color', [0.6 0 0], 'LineWidth', 1);
        hold on
        rand_idxs = randperm(size(Cb_SO_snapshots,1));
        if size(Cb_SO_snapshots,1) > max_to_plot
            rand_idxs = rand_idxs(1:max_to_plot);
        end
        plot(x_axis_LFP, Cb_SO_snapshots(rand_idxs,:), 'Color', [0 0 0.6], 'LineWidth', 1);
        hold off
        xlim([-500 500])
        if offset == 0
            saveas(gcf,[rootpath,animal,'/LFP_mean_amp_SO.fig'])
        else
            saveas(gcf,[rootpath,animal,'/LFP_mean_amp_SO_C', num2str(offset), '.fig'])
        end
        close all
        
        rand_idxs = randperm(size(M1_spind_snapshots,1));
        if size(M1_spind_snapshots,1) > max_to_plot
            rand_idxs = rand_idxs(1:max_to_plot);
        end
        plot(x_axis_LFP, M1_spind_snapshots(rand_idxs,:), 'Color', [0.6 0 0], 'LineWidth', 1);
        hold on
        rand_idxs = randperm(size(Cb_spind_snapshots,1));
        if size(Cb_spind_snapshots,1) > max_to_plot
            rand_idxs = rand_idxs(1:max_to_plot);
        end
        plot(x_axis_LFP, Cb_spind_snapshots(rand_idxs,:), 'Color', [0 0 0.6], 'LineWidth', 1);
        hold off
        xlim([-500 500])
        if offset == 0
            saveas(gcf,[rootpath,animal,'/LFP_mean_amp_spindle.fig'])
        else
            saveas(gcf,[rootpath,animal,'/LFP_mean_amp_spindle_C', num2str(offset), '.fig'])
        end
        close all
    end
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Create Polar Histogram of Selected Neurons' Spindle Spike Phase data (30)

if enabled(30)
    disp('Block 30...')
    day = 1;
    M1_neuron_chan = 1;
    M1_neuron_idx = 1;
    M1_neuron_num = 1; %If there is no data at {neuron_chan, neuron_idx}, the code uses the [neuron_num]th neuron with data
    M1_cycle = 1; % from 1 to 10. 1 and 10 are the left and right tails, 5 and 6 are the peak
    Cb_neuron_chan = 1;
    Cb_neuron_idx = 1;
    Cb_neuron_num = 1; %If there is no data at {neuron_chan, neuron_idx}, the code uses the [neuron_num]th neuron with data
    Cb_cycle = 1; % from 1 to 10. 1 and 10 are the left and right tails, 5 and 6 are the peak
    
    M1_phase_data = nan(0);
    Cb_phase_data = nan(0);
    
    for block = 1:param.sleep_blocks
        load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_spike_phasestamps.mat'])
        
%         if M1_cycle > size(prePeak_cycle_phases_M1_spindles_M1_spikes,??)
%             M1_phase_data = cat(??,M1_phase_data, postPeak_cycle_phases_M1_spindles_M1_spikes(M1_cycle - size(prePeak_cycle_phases_M1_spindles_M1_spikes,??)));
%         else
%             M1_phase_data = cat(??,M1_phase_data, prePeak_cycle_phases_M1_spindles_M1_spikes(abs(M1_cycle - size(prePeak_cycle_phases_M1_spindles_M1_spikes,??))));
%         end
%         
%         if Cb_cycle > size(prePeak_cycle_phases_Cb_spindles_Cb_spikes,??)
%             Cb_phase_data = cat(??,Cb_phase_data, postPeak_cycle_phases_Cb_spindles_Cb_spikes(Cb_cycle - size(prePeak_cycle_phases_Cb_spindles_Cb_spikes,??)));
%         else
%            Cb_phase_data = cat(??,Cb_phase_data, prePeak_cycle_phases_Cb_spindles_Cb_spikes(abs(Cb_cycle - size(prePeak_cycle_phases_Cb_spindles_Cb_spikes,??))));
%         end
    end
    
    
end

%% Normalized M1 Spike-spike within area cofiring in M1 spindle cycles using shuffling (31)

if enabled(31)
    disp('Block 31...')
    addpath(genpath('Z:\Matlab for analysis\eeglab\functions'))
    rmpath(genpath('Z:\Matlab for analysis\eeglab\functions\octavefunc\signal'))
    pre_margin = 0.05; %in seconds
    post_margin = 0.05; %in seconds
    trace_margin = 1; %in seconds - Spindles closer than this to the edge of the recording will be skipped to avoid index out of bounds errors
    pre_cycles = 5;
    post_cycles = 5;
    shuffle_reps = 25;
    max_idx_err = 10;
    kernel_stdv = 5;
    control_offsets = [-5, -10];
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Filtered_Sleep_data.mat'])
            load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_timestamps.mat']);
            
            for offset = [0, control_offsets]
                spindle_idxs = round(M1_filt_data.Fs_LFP * (M1_filt_data.spindles{1}.pks + offset));
                spindle_idxs(spindle_idxs < trace_margin * M1_filt_data.Fs_LFP) = [];
                spindle_idxs(spindle_idxs > length(M1_filt_data.LFP) - (trace_margin * M1_filt_data.Fs_LFP)) = [];
                
                %get spindle phase data
                hilbert_LFP = hilbert(M1_filt_data.LFP);
                inst_phase = unwrap(angle(hilbert_LFP));%inst phase
                radial_phase = mod(inst_phase,2*pi);
                peaks = diff(radial_phase) < 0;
                %align spindle_indexes with a peak
                peak_idxs = find(peaks);
                for i = 1:length(spindle_idxs)
                    [idx_diff, spindle_idxs(i)] = min(abs(peak_idxs - spindle_idxs(i)));
                    if idx_diff > max_idx_err
                        warning('A spindle index had greater offset from a cycle peak than is allowed by the threshold');
                    end
                end
                
                raw_CCH = nan((sum(param.M1_task_related_neurons{day}(:)) * (sum(param.M1_task_related_neurons{day}(:))-1))/2,(pre_cycles+post_cycles),round((pre_margin+post_margin)*1000)+1);
                pair_idx = 0;
                for M1_neuron1 = 1:(param.M1_chans * param.M1_neurons)
                    if param.M1_task_related_neurons{day}(M1_neuron1)
                        for M1_neuron2 = (M1_neuron1+1):(param.M1_chans * param.M1_neurons)
                            if param.M1_task_related_neurons{day}(M1_neuron2)
                                pair_idx = pair_idx+1;
                                disp(['Day ' num2str(day) ', Block ' num2str(block) ', Offset ' num2str(offset) ': Running pair ' num2str(pair_idx) ' out of ' num2str(size(raw_CCH,1)) '.'])
                                for cycle = (1-pre_cycles):post_cycles %zero-th cycle is the cycle imeadiatly preceding the center
                                    lag_hists = nan(0,round((pre_margin+post_margin)*1000)+1);
                                    for spindle = 1:length(spindle_idxs)
                                        
                                        cycle_start_idx = peak_idxs(spindle_idxs(spindle)+(cycle-1));
                                        cycle_end_idx = peak_idxs(spindle_idxs(spindle)+(cycle));
                                        
                                        M1_cycle_spikes1 = M1_spike_timestamps{M1_neuron1}(M1_spike_timestamps{M1_neuron1} > (cycle_start_idx/param.M1_Fs)  &  M1_spike_timestamps{M1_neuron1} < (cycle_end_idx/param.M1_Fs));
                                        for M1_spiketime1 = M1_cycle_spikes1
                                            M1_spikes2 = M1_spike_timestamps{M1_neuron2}(M1_spike_timestamps{M1_neuron2} >= (M1_spiketime1 - pre_margin) & M1_spike_timestamps{M1_neuron2} <= (M1_spiketime1 + post_margin));
                                            M1_spikes2 = M1_spikes2 - M1_spiketime1; %It should be 2 - 1 for both
                                            lag_hists = cat(1,lag_hists,histcounts(M1_spikes2,(-(round(pre_margin*1000)+0.5)):(round(post_margin*1000)+0.5)));
                                        end
                                        
                                        M1_cycle_spikes2 = M1_spike_timestamps{M1_neuron2}(M1_spike_timestamps{M1_neuron2} > (cycle_start_idx/param.M1_Fs)  &  M1_spike_timestamps{M1_neuron2} < (cycle_end_idx/param.M1_Fs));
                                        for M1_spiketime2 = M1_cycle_spikes2
                                            M1_spikes1 = M1_spike_timestamps{M1_neuron1}(M1_spike_timestamps{M1_neuron1} >= (M1_spiketime2 - pre_margin) & M1_spike_timestamps{M1_neuron1} <= (M1_spiketime2 + post_margin));
                                            M1_spikes1 = M1_spiketime2 - M1_spikes1; %It should be 2 - 1 for both
                                            lag_hists = cat(1,lag_hists,histcounts(M1_spikes1,(-(round(pre_margin*1000)+0.5)):(round(post_margin*1000)+0.5)));
                                        end
                                    end
                                    %average hist arrays across spindles and store it
                                    if isempty(lag_hists)
                                        raw_CCH(pair_idx,(cycle + pre_cycles),:) = 0;
                                    else
                                        raw_CCH(pair_idx,(cycle + pre_cycles),:) = mean(lag_hists,1);
                                    end
                                    
                                end
                            end
                        end
                    end
                end
                
                
                shuffled_CCH = nan((sum(param.M1_task_related_neurons{day}(:)) * (sum(param.M1_task_related_neurons{day}(:))-1))/2,(pre_cycles+post_cycles),round((pre_margin+post_margin)*1000)+1,shuffle_reps);
                shuffled_spindle_idxs = nan(shuffle_reps,length(spindle_idxs));
                for shuffle = 1:shuffle_reps
                    shuffled_spindle_idxs(shuffle,:) = randperm(length(spindle_idxs));
                end
                pair_idx = 0;
                for M1_neuron1 = 1:(param.M1_chans * param.M1_neurons)
                    if param.M1_task_related_neurons{day}(M1_neuron1)
                        for M1_neuron2 = 1:(param.M1_chans * param.M1_neurons)
                            if param.M1_task_related_neurons{day}(M1_neuron2) && (M1_neuron1 ~= M1_neuron2)
                                pair_idx = pair_idx+1;
                                disp(['Day ' num2str(day) ', Block ' num2str(block) ', Offset ' num2str(offset) ': Running shuffled baseline for pair ' num2str(pair_idx) ' out of ' num2str(size(shuffled_CCH,1)) '.'])
                                for cycle = (1-pre_cycles):post_cycles %zero-th cycle is the cycle imeadiatly preceding the center
                                    for shuffle = 1:shuffle_reps
                                        lag_hists = nan(0,round((pre_margin+post_margin)*1000)+1);
                                        for spindle = 1:length(spindle_idxs)
                                            
                                            cycle_start_idx = peak_idxs(spindle_idxs(spindle)+(cycle-1));
                                            cycle_end_idx = peak_idxs(spindle_idxs(spindle)+(cycle));
                                            
                                            
                                            M1_cycle_spikes1 = M1_spike_timestamps{M1_neuron1}(M1_spike_timestamps{M1_neuron1} > (cycle_start_idx/param.M1_Fs)  &  M1_spike_timestamps{M1_neuron1} < (cycle_end_idx/param.M1_Fs));
                                            M1_cycle_spikes1 = M1_cycle_spikes1 + (peak_idxs(spindle_idxs(shuffled_spindle_idxs(shuffle,spindle)))/param.M1_Fs) - (peak_idxs(spindle_idxs(spindle))/param.M1_Fs);
                                            for M1_spiketime1 = M1_cycle_spikes1
                                                M1_spikes2 = M1_spike_timestamps{M1_neuron2}(M1_spike_timestamps{M1_neuron2} >= (M1_spiketime1 - pre_margin) & M1_spike_timestamps{M1_neuron2} <= (M1_spiketime1 + post_margin));
                                                M1_spikes2 = M1_spikes2 - M1_spiketime1; %It should be 2 - 1 for both
                                                lag_hists = cat(1,lag_hists,histcounts(M1_spikes2,(-(round(pre_margin*1000)+0.5)):(round(post_margin*1000)+0.5)));
                                            end
                                            
                                            M1_cycle_spikes2 = M1_spike_timestamps{M1_neuron2}(M1_spike_timestamps{M1_neuron2} > (cycle_start_idx/param.M1_Fs)  &  M1_spike_timestamps{M1_neuron2} < (cycle_end_idx/param.M1_Fs));
                                            M1_cycle_spikes2 = M1_cycle_spikes2 + (peak_idxs(spindle_idxs(shuffled_spindle_idxs(shuffle,spindle)))/param.M1_Fs) - (peak_idxs(spindle_idxs(spindle))/param.M1_Fs);
                                            for M1_spiketime2 = M1_cycle_spikes2
                                                M1_spikes1 = M1_spike_timestamps{M1_neuron1}(M1_spike_timestamps{M1_neuron1} >= (M1_spiketime2 - pre_margin) & M1_spike_timestamps{M1_neuron1} <= (M1_spiketime2 + post_margin));
                                                M1_spikes1 = M1_spiketime2 - M1_spikes1; %It should be 2 - 1 for both
                                                lag_hists = cat(1,lag_hists,histcounts(M1_spikes1,(-(round(pre_margin*1000)+0.5)):(round(post_margin*1000)+0.5)));
                                            end
                                        end
                                        %average hist array across spindles and store it
                                        if isempty(lag_hists)
                                            shuffled_CCH(pair_idx,(cycle + pre_cycles),:,shuffle) = 0;
                                        else
                                            shuffled_CCH(pair_idx,(cycle + pre_cycles),:,shuffle) = mean(lag_hists,1);
                                        end
                                        
                                    end
                                end
                            end
                        end
                    end
                end
                
                
                ave_shuffled_CCH = mean(shuffled_CCH,4);
                corrected_CCH = raw_CCH-ave_shuffled_CCH;
                kernel = gausswin(size(corrected_CCH,3),size(corrected_CCH,3)/kernel_stdv);
                smoothed_CCH = nan(size(corrected_CCH));
                for i = 1:size(corrected_CCH,1)
                    for j = 1:size(corrected_CCH,2)
                        smoothed_CCH(i,j,:) = conv(squeeze(corrected_CCH(i,j,:)),kernel,'same');
                    end
                end
                [peak_val,peak_idx] = max(smoothed_CCH,[],3);
                peak_time = peak_idx - (round(pre_margin * 1000) + 1);
                pair_peaks_mean = mean(peak_val,1);
                pair_peaks_stdv = std(peak_val,0,1);
                pair_times_mean = mean(peak_time,1);
                pair_times_stdv = std(peak_time,0,1);
                pair_peaks_err = pair_peaks_stdv/sqrt(size(peak_time,1));
                pair_times_err = pair_times_stdv/sqrt(size(peak_time,1));
                
                x_axis = [1-pre_cycles:post_cycles];
                hold on
                errorbar(x_axis,pair_peaks_mean,pair_peaks_err,pair_peaks_err,zeros(1,pre_cycles+post_cycles),zeros(1,pre_cycles+post_cycles),'.')
                scatter(x_axis,pair_peaks_mean,'Marker', 's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 1], 'SizeData', 100)
                xlim([-4.5 5.5])
                hold off
                if offset == 0
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_prob_M1.fig'])
                else
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_prob_M1_C', num2str(offset), '.fig'])
                end
                close all
                
                hold on
                errorbar(x_axis,pair_times_mean,pair_times_err,pair_times_err,zeros(1,pre_cycles+post_cycles),zeros(1,pre_cycles+post_cycles),'.')
                scatter(x_axis,pair_times_mean,'Marker', 's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 1], 'SizeData', 100)
                xlim([-4.5 5.5])
                hold off
                if offset == 0
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_peak_times_M1.fig'])
                else
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_peak_times_M1_C', num2str(offset), '.fig'])
                end
                close all
            end
        end
    end
    rmpath(genpath('Z:\Matlab for analysis\eeglab\functions'))
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% HPC Normalized M1 Spike-spike within area cofiring in M1 spindle cycles using shuffling (31.5)

if false
    if isunix  %#ok<UNRCH>
        addpath(genpath('/common/fleischerp/HPC_code'))
        addpath(genpath('/common/fleischerp/3rd_party_code'))
        %rmpath(genpath('/common/fleischerp/3rd_party_code/octavefunc/signal'))
    else
        addpath(genpath('Z:\Matlab for analysis\eeglab\functions'))
        rmpath(genpath('Z:\Matlab for analysis\eeglab\functions\octavefunc\signal'))
    end
    worker_num = 24;
    %parpool('local',worker_num);
	
    pre_margin = 0.05; %in seconds
    post_margin = 0.05; %in seconds
    trace_margin = 1; %in seconds - Spindles closer than this to the edge of the recording will be skipped to avoid index out of bounds errors
    pre_cycles = 5;
    post_cycles = 5;
    shuffle_reps = 25;
    max_idx_err = 10;
    kernel_stdv = 5;
    control_offsets = [-5, -10];
    
    bin_edges = ((-(round(pre_margin*1000)+0.5)):(round(post_margin*1000)+0.5))/1000;
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Filtered_Sleep_data.mat'])
            load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_timestamps.mat']);
            
            for offset = [0, control_offsets]
                spindle_idxs = round(M1_filt_data.Fs_LFP * (M1_filt_data.spindles{1}.pks + offset));
                spindle_idxs(spindle_idxs < trace_margin * M1_filt_data.Fs_LFP) = [];
                spindle_idxs(spindle_idxs > length(M1_filt_data.LFP) - (trace_margin * M1_filt_data.Fs_LFP)) = [];
                
                %get spindle phase data
                hilbert_LFP = hilbert(M1_filt_data.LFP);
                inst_phase = unwrap(angle(hilbert_LFP));%inst phase
                radial_phase = mod(inst_phase,2*pi);
                peaks = diff(radial_phase) < 0;
                %align spindle_indexes with a peak
                peak_idxs = find(peaks);
                for i = 1:length(spindle_idxs)
                    [idx_diff, spindle_idxs(i)] = min(abs(peak_idxs - spindle_idxs(i)));
                    if idx_diff > max_idx_err
                        warning('A spindle index had greater offset from a cycle peak than is allowed by the threshold');
                    end
                end
                
                raw_CCH = nan((sum(param.M1_task_related_neurons{day}(:)) * (sum(param.M1_task_related_neurons{day}(:))-1))/2,(pre_cycles+post_cycles),round((pre_margin+post_margin)*1000)+1);
                neuron_idxs = find(param.M1_task_related_neurons{day});
                good_pairs = (neuron_idxs - neuron_idxs') > 0;
                pair_idxs = find(good_pairs);
                for pair_idx = 1:length(pair_idxs) %parfor
                    [M1_neuron1, M1_neuron2] = ind2sub(size(good_pairs),pair_idxs(pair_idx));
                    M1_neuron1 = neuron_idxs(M1_neuron1);
                    M1_neuron2 = neuron_idxs(M1_neuron2);
                    pair_data = nan(pre_cycles+post_cycles,round((pre_margin+post_margin)*1000)+1);
                    for cycle = (1-pre_cycles):post_cycles %zero-th cycle is the cycle imeadiatly preceding the center
                        lag_hists = nan(0,round((pre_margin+post_margin)*1000)+1);
                        for spindle = 1:length(spindle_idxs)
                            
                            cycle_start_idx = peak_idxs(spindle_idxs(spindle)+(cycle-1));
                            cycle_end_idx = peak_idxs(spindle_idxs(spindle)+(cycle));
                            
                            M1_cycle_spikes1 = M1_spike_timestamps{M1_neuron1}(M1_spike_timestamps{M1_neuron1} > (cycle_start_idx/param.M1_Fs)  &  M1_spike_timestamps{M1_neuron1} < (cycle_end_idx/param.M1_Fs));
                            for M1_spiketime1 = M1_cycle_spikes1
                                M1_spikes2 = M1_spike_timestamps{M1_neuron2}(M1_spike_timestamps{M1_neuron2} >= (M1_spiketime1 - pre_margin) & M1_spike_timestamps{M1_neuron2} <= (M1_spiketime1 + post_margin));
                                M1_spikes2 = M1_spikes2 - M1_spiketime1; %It should be 2 - 1 for both
                                lag_hists = cat(1,lag_hists,histcounts(M1_spikes2,bin_edges));
                            end
                            
                            M1_cycle_spikes2 = M1_spike_timestamps{M1_neuron2}(M1_spike_timestamps{M1_neuron2} > (cycle_start_idx/param.M1_Fs)  &  M1_spike_timestamps{M1_neuron2} < (cycle_end_idx/param.M1_Fs));
                            for M1_spiketime2 = M1_cycle_spikes2
                                M1_spikes1 = M1_spike_timestamps{M1_neuron1}(M1_spike_timestamps{M1_neuron1} >= (M1_spiketime2 - pre_margin) & M1_spike_timestamps{M1_neuron1} <= (M1_spiketime2 + post_margin));
                                M1_spikes1 = M1_spiketime2 - M1_spikes1; %It should be 2 - 1 for both
                                lag_hists = cat(1,lag_hists,histcounts(M1_spikes1,bin_edges));
                            end
                        end
                        %average hist arrays across spindles and store it
                        if isempty(lag_hists)
                            pair_data(cycle + pre_cycles,:) = 0;
                        else
                            pair_data(cycle + pre_cycles,:) = mean(lag_hists,1);
                        end
                    end
                    raw_CCH(pair_idx,:,:) = pair_data;
                end
                
                
                shuffled_CCH = nan((sum(param.M1_task_related_neurons{day}(:)) * (sum(param.M1_task_related_neurons{day}(:))-1))/2,(pre_cycles+post_cycles),round((pre_margin+post_margin)*1000)+1,shuffle_reps);
                shuffled_spindle_idxs = nan(shuffle_reps,length(spindle_idxs));
                for shuffle = 1:shuffle_reps
                    shuffled_spindle_idxs(shuffle,:) = randperm(length(spindle_idxs));
                end
                neuron_idxs = find(param.M1_task_related_neurons{day});
                good_pairs = (neuron_idxs - neuron_idxs') > 0;
                pair_idxs = find(good_pairs);
                for pair_idx = 1:length(pair_idxs) %parfor
                    [M1_neuron1, M1_neuron2] = ind2sub(size(good_pairs),pair_idxs(pair_idx));
                    M1_neuron1 = neuron_idxs(M1_neuron1);
                    M1_neuron2 = neuron_idxs(M1_neuron2);
                    pair_data = nan(pre_cycles+post_cycles,round((pre_margin+post_margin)*1000)+1);
                    for cycle = (1-pre_cycles):post_cycles %zero-th cycle is the cycle imeadiatly preceding the center
                        for shuffle = 1:shuffle_reps
                            lag_hists = nan(0,round((pre_margin+post_margin)*1000)+1);
                            for spindle = 1:length(spindle_idxs)
                                
                                cycle_start_idx = peak_idxs(spindle_idxs(spindle)+(cycle-1));
                                cycle_end_idx = peak_idxs(spindle_idxs(spindle)+(cycle));
                                
                                
                                M1_cycle_spikes1 = M1_spike_timestamps{M1_neuron1}(M1_spike_timestamps{M1_neuron1} > (cycle_start_idx/param.M1_Fs)  &  M1_spike_timestamps{M1_neuron1} < (cycle_end_idx/param.M1_Fs));
                                M1_cycle_spikes1 = M1_cycle_spikes1 + (peak_idxs(spindle_idxs(shuffled_spindle_idxs(shuffle,spindle)))/param.M1_Fs) - (peak_idxs(spindle_idxs(spindle))/param.M1_Fs);
                                for M1_spiketime1 = M1_cycle_spikes1
                                    M1_spikes2 = M1_spike_timestamps{M1_neuron2}(M1_spike_timestamps{M1_neuron2} >= (M1_spiketime1 - pre_margin) & M1_spike_timestamps{M1_neuron2} <= (M1_spiketime1 + post_margin));
                                    M1_spikes2 = M1_spikes2 - M1_spiketime1; %It should be 2 - 1 for both
                                    lag_hists = cat(1,lag_hists,histcounts(M1_spikes2,bin_edges));
                                end
                                
                                M1_cycle_spikes2 = M1_spike_timestamps{M1_neuron2}(M1_spike_timestamps{M1_neuron2} > (cycle_start_idx/param.M1_Fs)  &  M1_spike_timestamps{M1_neuron2} < (cycle_end_idx/param.M1_Fs));
                                M1_cycle_spikes2 = M1_cycle_spikes2 + (peak_idxs(spindle_idxs(shuffled_spindle_idxs(shuffle,spindle)))/param.M1_Fs) - (peak_idxs(spindle_idxs(spindle))/param.M1_Fs);
                                for M1_spiketime2 = M1_cycle_spikes2
                                    M1_spikes1 = M1_spike_timestamps{M1_neuron1}(M1_spike_timestamps{M1_neuron1} >= (M1_spiketime2 - pre_margin) & M1_spike_timestamps{M1_neuron1} <= (M1_spiketime2 + post_margin));
                                    M1_spikes1 = M1_spiketime2 - M1_spikes1; %It should be 2 - 1 for both
                                    lag_hists = cat(1,lag_hists,histcounts(M1_spikes1,bin_edges));
                                end
                            end
                            %average hist arrays across spindles and store it
                            if isempty(lag_hists)
                                pair_data(cycle + pre_cycles,:) = 0;
                            else
                                pair_data(cycle + pre_cycles,:) = mean(lag_hists,1);
                            end
                            shuffled_CCH(pair_idx,:,:,shuffle) = pair_data;
                        end
                    end
                end
                
                
                ave_shuffled_CCH = squeeze(mean(shuffled_CCH,4));
                corrected_CCH = raw_CCH-ave_shuffled_CCH;
                kernel = gausswin(size(corrected_CCH,3),size(corrected_CCH,3)/kernel_stdv);
                smoothed_CCH = nan(size(corrected_CCH));
                for i = 1:size(corrected_CCH,1)
                    for j = 1:size(corrected_CCH,2)
                        smoothed_CCH(i,j,:) = conv(squeeze(corrected_CCH(i,j,:)),kernel,'same');
                    end
                end
                [peak_val,peak_idx] = max(smoothed_CCH,[],3);
                peak_time = peak_idx - (round(pre_margin * 1000) + 1);
                pair_peaks_mean = mean(peak_val,1);
                pair_peaks_stdv = std(peak_val,0,1);
                pair_times_mean = mean(peak_time,1);
                pair_times_stdv = std(peak_time,0,1);
                pair_peaks_err = pair_peaks_stdv/sqrt(size(peak_time,1));
                pair_times_err = pair_times_stdv/sqrt(size(peak_time,1));
                
                if offset == 0
                    save([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_data_M1.mat'], 'peak_val', 'peak_time')
                else
                    save([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_data_M1_C', num2str(offset), '.mat'], 'peak_val', 'peak_time')
                end
                
                x_axis = [1-pre_cycles:post_cycles];
                hold on
                errorbar(x_axis,pair_peaks_mean,pair_peaks_err,pair_peaks_err,zeros(1,pre_cycles+post_cycles),zeros(1,pre_cycles+post_cycles),'.')
                scatter(x_axis,pair_peaks_mean,'Marker', 's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 1], 'SizeData', 100)
                xlim([-4.5 5.5])
                hold off
                if offset == 0
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_prob_M1.fig'])
                else
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_prob_M1_C', num2str(offset), '.fig'])
                end
                close all
                
                hold on
                errorbar(x_axis,pair_times_mean,pair_times_err,pair_times_err,zeros(1,pre_cycles+post_cycles),zeros(1,pre_cycles+post_cycles),'.')
                scatter(x_axis,pair_times_mean,'Marker', 's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 1], 'SizeData', 100)
                xlim([-4.5 5.5])
                hold off
                if offset == 0
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_peak_times_M1.fig'])
                else
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_peak_times_M1_C', num2str(offset), '.fig'])
                end
                close all
            end
        end
    end
    if isunix
        rmpath(genpath('/common/fleischerp/HPC_code'))
        rmpath(genpath('/common/fleischerp/3rd_party_code'))
    else
        rmpath(genpath('Z:\Matlab for analysis\eeglab\functions'))
    end
	delete(gcp('nocreate'));
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Normalized Cb Spike-spike within area cofiring in Cb spindle cycles using shuffling (32)

if enabled(32)
    disp('Block 32...')
    addpath(genpath('Z:\Matlab for analysis\eeglab\functions'))
    rmpath(genpath('Z:\Matlab for analysis\eeglab\functions\octavefunc\signal'))
    pre_margin = 0.05; %in seconds
    post_margin = 0.05; %in seconds
    trace_margin = 1; %in seconds - Spindles closer than this to the edge of the recording will be skipped to avoid index out of bounds errors
    pre_cycles = 5;
    post_cycles = 5;
    shuffle_reps = 25;
    max_idx_err = 10;
    kernel_stdv = 5;
    control_offsets = [-5, -10];
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Filtered_Sleep_data.mat'])
            load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_timestamps.mat']);
            
            for offset = [0, control_offsets]
                spindle_idxs = round(Cb_filt_data.Fs_LFP * (Cb_filt_data.spindles{1}.pks + offset));
                spindle_idxs(spindle_idxs < trace_margin * Cb_filt_data.Fs_LFP) = [];
                spindle_idxs(spindle_idxs > length(Cb_filt_data.LFP) - (trace_margin * Cb_filt_data.Fs_LFP)) = [];
                
                %get spindle phase data
                hilbert_LFP = hilbert(Cb_filt_data.LFP);
                inst_phase = unwrap(angle(hilbert_LFP));%inst phase
                radial_phase = mod(inst_phase,2*pi);
                peaks = diff(radial_phase) < 0;
                %align spindle_indexes with a peak
                peak_idxs = find(peaks);
                for i = 1:length(spindle_idxs)
                    [idx_diff, spindle_idxs(i)] = min(abs(peak_idxs - spindle_idxs(i)));
                    if idx_diff > max_idx_err
                        warning('A spindle index had greater offset from a cycle peak than is allowed by the threshold');
                    end
                end
                
                raw_CCH = nan((sum(param.Cb_task_related_neurons{day}(:)) * (sum(param.Cb_task_related_neurons{day}(:))-1))/2,(pre_cycles+post_cycles),round((pre_margin+post_margin)*1000)+1);
                pair_idx = 0;
                for Cb_neuron1 = 1:(param.Cb_chans * param.Cb_neurons)
                    if param.Cb_task_related_neurons{day}(Cb_neuron1)
                        for Cb_neuron2 = (Cb_neuron1+1):(param.Cb_chans * param.Cb_neurons)
                            if param.Cb_task_related_neurons{day}(Cb_neuron2)
                                pair_idx = pair_idx+1;
                                disp(['Day ' num2str(day) ', Block ' num2str(block) ', Offset ' num2str(offset) ': Running pair ' num2str(pair_idx) ' out of ' num2str(size(raw_CCH,1)) '.'])
                                for cycle = (1-pre_cycles):post_cycles %zero-th cycle is the cycle imeadiatly preceding the center
                                    lag_hists = nan(0,round((pre_margin+post_margin)*1000)+1);
                                    for spindle = 1:length(spindle_idxs)
                                        
                                        cycle_start_idx = peak_idxs(spindle_idxs(spindle)+(cycle-1));
                                        cycle_end_idx = peak_idxs(spindle_idxs(spindle)+(cycle));
                                        
                                        Cb_cycle_spikes1 = Cb_spike_timestamps{Cb_neuron1}(Cb_spike_timestamps{Cb_neuron1} > (cycle_start_idx/param.Cb_Fs)  &  Cb_spike_timestamps{Cb_neuron1} < (cycle_end_idx/param.Cb_Fs));
                                        for Cb_spiketime1 = Cb_cycle_spikes1
                                            Cb_spikes2 = Cb_spike_timestamps{Cb_neuron2}(Cb_spike_timestamps{Cb_neuron2} >= (Cb_spiketime1 - pre_margin) & Cb_spike_timestamps{Cb_neuron2} <= (Cb_spiketime1 + post_margin));
                                            Cb_spikes2 = Cb_spikes2 - Cb_spiketime1; %It should be 2 - 1 for both
                                            lag_hists = cat(1,lag_hists,histcounts(Cb_spikes2,(-(round(pre_margin*1000)+0.5)):(round(post_margin*1000)+0.5)));
                                        end
                                        
                                        Cb_cycle_spikes2 = Cb_spike_timestamps{Cb_neuron2}(Cb_spike_timestamps{Cb_neuron2} > (cycle_start_idx/param.Cb_Fs)  &  Cb_spike_timestamps{Cb_neuron2} < (cycle_end_idx/param.Cb_Fs));
                                        for Cb_spiketime2 = Cb_cycle_spikes2
                                            Cb_spikes1 = Cb_spike_timestamps{Cb_neuron1}(Cb_spike_timestamps{Cb_neuron1} >= (Cb_spiketime2 - pre_margin) & Cb_spike_timestamps{Cb_neuron1} <= (Cb_spiketime2 + post_margin));
                                            Cb_spikes1 = Cb_spiketime2 - Cb_spikes1; %It should be 2 - 1 for both
                                            lag_hists = cat(1,lag_hists,histcounts(Cb_spikes1,(-(round(pre_margin*1000)+0.5)):(round(post_margin*1000)+0.5)));
                                        end
                                    end
                                    %average hist arrays across spindles and store it
                                    if isempty(lag_hists)
                                        raw_CCH(pair_idx,(cycle + pre_cycles),:) = 0;
                                    else
                                        raw_CCH(pair_idx,(cycle + pre_cycles),:) = mean(lag_hists,1);
                                    end
                                    
                                end
                            end
                        end
                    end
                end
                
                
                shuffled_CCH = nan((sum(param.Cb_task_related_neurons{day}(:)) * sum(param.Cb_task_related_neurons{day}(:)))/2,(pre_cycles+post_cycles),round((pre_margin+post_margin)*1000)+1,shuffle_reps);
                shuffled_spindle_idxs = nan(shuffle_reps,length(spindle_idxs));
                for shuffle = 1:shuffle_reps
                    shuffled_spindle_idxs(shuffle,:) = randperm(length(spindle_idxs));
                end
                pair_idx = 0;
                for Cb_neuron1 = 1:(param.Cb_chans * param.Cb_neurons)
                    if param.Cb_task_related_neurons{day}(Cb_neuron1)
                        for Cb_neuron2 = 1:(param.Cb_chans * param.Cb_neurons)
                            if param.Cb_task_related_neurons{day}(Cb_neuron2) && (Cb_neuron1 ~= Cb_neuron2)
                                pair_idx = pair_idx+1;
                                disp(['Day ' num2str(day) ', Block ' num2str(block) ', Offset ' num2str(offset) ': Running shuffled baseline for pair ' num2str(pair_idx) ' out of ' num2str(size(shuffled_CCH,1)) '.'])
                                for cycle = (1-pre_cycles):post_cycles %zero-th cycle is the cycle imeadiatly preceding the center
                                    for shuffle = 1:shuffle_reps
                                        lag_hists = nan(0,round((pre_margin+post_margin)*1000)+1);
                                        for spindle = 1:length(spindle_idxs)
                                            
                                            cycle_start_idx = peak_idxs(spindle_idxs(spindle)+(cycle-1));
                                            cycle_end_idx = peak_idxs(spindle_idxs(spindle)+(cycle));
                                            
                                            
                                            Cb_cycle_spikes1 = Cb_spike_timestamps{Cb_neuron1}(Cb_spike_timestamps{Cb_neuron1} > (cycle_start_idx/param.Cb_Fs)  &  Cb_spike_timestamps{Cb_neuron1} < (cycle_end_idx/param.Cb_Fs));
                                            Cb_cycle_spikes1 = Cb_cycle_spikes1 + (peak_idxs(spindle_idxs(shuffled_spindle_idxs(shuffle,spindle)))/param.Cb_Fs) - (peak_idxs(spindle_idxs(spindle))/param.Cb_Fs);
                                            for Cb_spiketime1 = Cb_cycle_spikes1
                                                Cb_spikes2 = Cb_spike_timestamps{Cb_neuron2}(Cb_spike_timestamps{Cb_neuron2} >= (Cb_spiketime1 - pre_margin) & Cb_spike_timestamps{Cb_neuron2} <= (Cb_spiketime1 + post_margin));
                                                Cb_spikes2 = Cb_spikes2 - Cb_spiketime1; %It should be 2 - 1 for both
                                                lag_hists = cat(1,lag_hists,histcounts(Cb_spikes2,(-(round(pre_margin*1000)+0.5)):(round(post_margin*1000)+0.5)));
                                            end
                                            
                                            Cb_cycle_spikes2 = Cb_spike_timestamps{Cb_neuron2}(Cb_spike_timestamps{Cb_neuron2} > (cycle_start_idx/param.Cb_Fs)  &  Cb_spike_timestamps{Cb_neuron2} < (cycle_end_idx/param.Cb_Fs));
                                            Cb_cycle_spikes2 = Cb_cycle_spikes2 + (peak_idxs(spindle_idxs(shuffled_spindle_idxs(shuffle,spindle)))/param.Cb_Fs) - (peak_idxs(spindle_idxs(spindle))/param.Cb_Fs);
                                            for Cb_spiketime2 = Cb_cycle_spikes2
                                                Cb_spikes1 = Cb_spike_timestamps{Cb_neuron1}(Cb_spike_timestamps{Cb_neuron1} >= (Cb_spiketime2 - pre_margin) & Cb_spike_timestamps{Cb_neuron1} <= (Cb_spiketime2 + post_margin));
                                                Cb_spikes1 = Cb_spiketime2 - Cb_spikes1; %It should be 2 - 1 for both
                                                lag_hists = cat(1,lag_hists,histcounts(Cb_spikes1,(-(round(pre_margin*1000)+0.5)):(round(post_margin*1000)+0.5)));
                                            end
                                        end
                                        %average hist array across spindles and store it
                                        if isempty(lag_hists)
                                            shuffled_CCH(pair_idx,(cycle + pre_cycles),:,shuffle) = 0;
                                        else
                                            shuffled_CCH(pair_idx,(cycle + pre_cycles),:,shuffle) = mean(lag_hists,1);
                                        end
                                        
                                    end
                                end
                            end
                        end
                    end
                end
                
                
                ave_shuffled_CCH = mean(shuffled_CCH,4);
                corrected_CCH = raw_CCH-ave_shuffled_CCH;
                kernel = gausswin(size(corrected_CCH,3),size(corrected_CCH,3)/kernel_stdv);
                smoothed_CCH = nan(size(corrected_CCH));
                for i = 1:size(corrected_CCH,1)
                    for j = 1:size(corrected_CCH,2)
                        smoothed_CCH(i,j,:) = conv(squeeze(corrected_CCH(i,j,:)),kernel,'same');
                    end
                end
                [peak_val,peak_idx] = max(smoothed_CCH,[],3);
                peak_time = peak_idx - (round(pre_margin * 1000) + 1);
                pair_peaks_mean = mean(peak_val,1);
                pair_peaks_stdv = std(peak_val,0,1);
                pair_times_mean = mean(peak_time,1);
                pair_times_stdv = std(peak_time,0,1);
                pair_peaks_err = pair_peaks_stdv/sqrt(size(peak_time,1));
                pair_times_err = pair_times_stdv/sqrt(size(peak_time,1));
                
                x_axis = [1-pre_cycles:post_cycles];
                hold on
                errorbar(x_axis,pair_peaks_mean,pair_peaks_err,pair_peaks_err,zeros(1,pre_cycles+post_cycles),zeros(1,pre_cycles+post_cycles),'.')
                scatter(x_axis,pair_peaks_mean,'Marker', 's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 1], 'SizeData', 100)
                xlim([-4.5 5.5])
                hold off
                if offset == 0
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_prob_Cb.fig'])
                else
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_prob_Cb_C', num2str(offset), '.fig'])
                end
                close all
                
                hold on
                errorbar(x_axis,pair_times_mean,pair_times_err,pair_times_err,zeros(1,pre_cycles+post_cycles),zeros(1,pre_cycles+post_cycles),'.')
                scatter(x_axis,pair_times_mean,'Marker', 's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 1], 'SizeData', 100)
                xlim([-4.5 5.5])
                hold off
                if offset == 0
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_peak_times_Cb.fig'])
                else
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_peak_times_Cb_C', num2str(offset), '.fig'])
                end
                close all
            end
        end
    end
    rmpath(genpath('Z:\Matlab for analysis\eeglab\functions'))
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% HPC Normalized Cb Spike-spike within area cofiring in Cb spindle cycles using shuffling (32.5)

if false
    if isunix  %#ok<UNRCH>
        addpath(genpath('/common/fleischerp/HPC_code'))
        addpath(genpath('/common/fleischerp/3rd_party_code'))
        %rmpath(genpath('/common/fleischerp/3rd_party_code/octavefunc/signal'))
    else
        addpath(genpath('Z:\Matlab for analysis\eeglab\functions'))
        rmpath(genpath('Z:\Matlab for analysis\eeglab\functions\octavefunc\signal'))
    end
    worker_num = 24;
    %parpool('local', worker_num);
    %pair_block = [1 worker_num];
    %pair_block = pair_block + (24 * 2);
    
    pre_margin = 0.05; %in seconds
    post_margin = 0.05; %in seconds
    trace_margin = 1; %in seconds - Spindles closer than this to the edge of the recording will be skipped to avoid index out of bounds errors
    pre_cycles = 5;
    post_cycles = 5;
    shuffle_reps = 25;
    max_idx_err = 10;
    kernel_stdv = 5;
    control_offsets = [-5, -10];
    
    bin_edges = ((-(round(pre_margin*1000)+0.5)):(round(post_margin*1000)+0.5))/1000;
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Filtered_Sleep_data.mat'])
            load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_timestamps.mat']);
            
            for offset = [0, control_offsets]
                spindle_idxs = round(Cb_filt_data.Fs_LFP * (Cb_filt_data.spindles{1}.pks + offset));
                spindle_idxs(spindle_idxs < trace_margin * Cb_filt_data.Fs_LFP) = [];
                spindle_idxs(spindle_idxs > length(Cb_filt_data.LFP) - (trace_margin * Cb_filt_data.Fs_LFP)) = [];
                
                %get spindle phase data
                hilbert_LFP = hilbert(Cb_filt_data.LFP);
                inst_phase = unwrap(angle(hilbert_LFP));%inst phase
                radial_phase = mod(inst_phase,2*pi);
                peaks = diff(radial_phase) < 0;
                %align spindle_indexes with a peak
                peak_idxs = find(peaks);
                for i = 1:length(spindle_idxs)
                    [idx_diff, spindle_idxs(i)] = min(abs(peak_idxs - spindle_idxs(i)));
                    if idx_diff > max_idx_err && offset == 0
                        warning('A spindle index had greater offset from a cycle peak than is allowed by the threshold');
                    end
                end
                
                raw_CCH = nan((sum(param.Cb_task_related_neurons{day}(:)) * (sum(param.Cb_task_related_neurons{day}(:))-1))/2,(pre_cycles+post_cycles),round((pre_margin+post_margin)*1000)+1);
                neuron_idxs = find(param.Cb_task_related_neurons{day});
                good_pairs = (neuron_idxs - neuron_idxs') > 0;
                pair_idxs = find(good_pairs);
                for pair_idx = 1:length(pair_idxs) % parfor %pair_block(1):min(length(pair_idxs),pair_block(2))
                    [Cb_neuron1, Cb_neuron2] = ind2sub(size(good_pairs),pair_idxs(pair_idx));
                    Cb_neuron1 = neuron_idxs(Cb_neuron1);
                    Cb_neuron2 = neuron_idxs(Cb_neuron2);
                    pair_data = nan(pre_cycles+post_cycles,round((pre_margin+post_margin)*1000)+1);
                    for cycle = (1-pre_cycles):post_cycles %zero-th cycle is the cycle imeadiatly preceding the center
                        lag_hists = nan(0,round((pre_margin+post_margin)*1000)+1);
                        for spindle = 1:length(spindle_idxs)
                            
                            cycle_start_idx = peak_idxs(spindle_idxs(spindle)+(cycle-1));
                            cycle_end_idx = peak_idxs(spindle_idxs(spindle)+(cycle));
                            
                            Cb_cycle_spikes1 = Cb_spike_timestamps{Cb_neuron1}(Cb_spike_timestamps{Cb_neuron1} > (cycle_start_idx/param.Cb_Fs)  &  Cb_spike_timestamps{Cb_neuron1} < (cycle_end_idx/param.Cb_Fs));
                            for Cb_spiketime1 = Cb_cycle_spikes1
                                Cb_spikes2 = Cb_spike_timestamps{Cb_neuron2}(Cb_spike_timestamps{Cb_neuron2} >= (Cb_spiketime1 - pre_margin) & Cb_spike_timestamps{Cb_neuron2} <= (Cb_spiketime1 + post_margin));
                                Cb_spikes2 = Cb_spikes2 - Cb_spiketime1; %It should be 2 - 1 for both
                                lag_hists = cat(1,lag_hists,histcounts(Cb_spikes2,bin_edges));
                            end
                            
                            Cb_cycle_spikes2 = Cb_spike_timestamps{Cb_neuron2}(Cb_spike_timestamps{Cb_neuron2} > (cycle_start_idx/param.Cb_Fs)  &  Cb_spike_timestamps{Cb_neuron2} < (cycle_end_idx/param.Cb_Fs));
                            for Cb_spiketime2 = Cb_cycle_spikes2
                                Cb_spikes1 = Cb_spike_timestamps{Cb_neuron1}(Cb_spike_timestamps{Cb_neuron1} >= (Cb_spiketime2 - pre_margin) & Cb_spike_timestamps{Cb_neuron1} <= (Cb_spiketime2 + post_margin));
                                Cb_spikes1 = Cb_spiketime2 - Cb_spikes1; %It should be 2 - 1 for both
                                lag_hists = cat(1,lag_hists,histcounts(Cb_spikes1,bin_edges));
                            end
                        end
                        %average hist arrays across spindles and store it
                        c_idx = cycle + pre_cycles;
                        if isempty(lag_hists)
                            pair_data(cycle + pre_cycles,:) = 0;
                        else
                            pair_data(cycle + pre_cycles,:) = mean(lag_hists,1);
                        end
                    end
                    raw_CCH(pair_idx,:,:) = pair_data;
                end
                
                
                shuffled_CCH = nan((sum(param.Cb_task_related_neurons{day}(:)) * (sum(param.Cb_task_related_neurons{day}(:))-1))/2,(pre_cycles+post_cycles),round((pre_margin+post_margin)*1000)+1,shuffle_reps);
                shuffled_spindle_idxs = nan(shuffle_reps,length(spindle_idxs));
                for shuffle = 1:shuffle_reps
                    shuffled_spindle_idxs(shuffle,:) = randperm(length(spindle_idxs));
                end
                neuron_idxs = find(param.Cb_task_related_neurons{day});
                good_pairs = (neuron_idxs - neuron_idxs') > 0;
                pair_idxs = find(good_pairs);
                for pair_idx = 1:length(pair_idxs) %parfor
                    [Cb_neuron1, Cb_neuron2] = ind2sub(size(good_pairs),pair_idxs(pair_idx));
                    Cb_neuron1 = neuron_idxs(Cb_neuron1);
                    Cb_neuron2 = neuron_idxs(Cb_neuron2);
                    pair_data = nan(pre_cycles+post_cycles,round((pre_margin+post_margin)*1000)+1);
                    for cycle = (1-pre_cycles):post_cycles %zero-th cycle is the cycle imeadiatly preceding the center
                        for shuffle = 1:shuffle_reps
                            lag_hists = nan(0,round((pre_margin+post_margin)*1000)+1);
                            for spindle = 1:length(spindle_idxs)
                                
                                cycle_start_idx = peak_idxs(spindle_idxs(spindle)+(cycle-1));
                                cycle_end_idx = peak_idxs(spindle_idxs(spindle)+(cycle));
                                
                                
                                Cb_cycle_spikes1 = Cb_spike_timestamps{Cb_neuron1}(Cb_spike_timestamps{Cb_neuron1} > (cycle_start_idx/param.Cb_Fs)  &  Cb_spike_timestamps{Cb_neuron1} < (cycle_end_idx/param.Cb_Fs));
                                Cb_cycle_spikes1 = Cb_cycle_spikes1 + (peak_idxs(spindle_idxs(shuffled_spindle_idxs(shuffle,spindle)))/param.Cb_Fs) - (peak_idxs(spindle_idxs(spindle))/param.Cb_Fs);
                                for Cb_spiketime1 = Cb_cycle_spikes1
                                    Cb_spikes2 = Cb_spike_timestamps{Cb_neuron2}(Cb_spike_timestamps{Cb_neuron2} >= (Cb_spiketime1 - pre_margin) & Cb_spike_timestamps{Cb_neuron2} <= (Cb_spiketime1 + post_margin));
                                    Cb_spikes2 = Cb_spikes2 - Cb_spiketime1; %It should be 2 - 1 for both
                                    lag_hists = cat(1,lag_hists,histcounts(Cb_spikes2,bin_edges));
                                end
                                
                                Cb_cycle_spikes2 = Cb_spike_timestamps{Cb_neuron2}(Cb_spike_timestamps{Cb_neuron2} > (cycle_start_idx/param.Cb_Fs)  &  Cb_spike_timestamps{Cb_neuron2} < (cycle_end_idx/param.Cb_Fs));
                                Cb_cycle_spikes2 = Cb_cycle_spikes2 + (peak_idxs(spindle_idxs(shuffled_spindle_idxs(shuffle,spindle)))/param.Cb_Fs) - (peak_idxs(spindle_idxs(spindle))/param.Cb_Fs);
                                for Cb_spiketime2 = Cb_cycle_spikes2
                                    Cb_spikes1 = Cb_spike_timestamps{Cb_neuron1}(Cb_spike_timestamps{Cb_neuron1} >= (Cb_spiketime2 - pre_margin) & Cb_spike_timestamps{Cb_neuron1} <= (Cb_spiketime2 + post_margin));
                                    Cb_spikes1 = Cb_spiketime2 - Cb_spikes1; %It should be 2 - 1 for both
                                    lag_hists = cat(1,lag_hists,histcounts(Cb_spikes1,bin_edges));
                                end
                            end
                            %average hist arrays across spindles and store it
                            if isempty(lag_hists)
                                pair_data(cycle + pre_cycles,:) = 0;
                            else
                                pair_data(cycle + pre_cycles,:) = mean(lag_hists,1);
                            end
                            shuffled_CCH(pair_idx,:,:,shuffle) = pair_data;
                        end
                    end
                end
                
                
                ave_shuffled_CCH = mean(shuffled_CCH,4);
                corrected_CCH = raw_CCH-ave_shuffled_CCH;
                kernel = gausswin(size(corrected_CCH,3),size(corrected_CCH,3)/kernel_stdv);
                smoothed_CCH = nan(size(corrected_CCH));
                for i = 1:size(corrected_CCH,1)
                    for j = 1:size(corrected_CCH,2)
                        smoothed_CCH(i,j,:) = conv(squeeze(corrected_CCH(i,j,:)),kernel,'same');
                    end
                end
                [peak_val,peak_idx] = max(smoothed_CCH,[],3);
                peak_time = peak_idx - (round(pre_margin * 1000) + 1);
                pair_peaks_mean = mean(peak_val,1);
                pair_peaks_stdv = std(peak_val,0,1);
                pair_times_mean = mean(peak_time,1);
                pair_times_stdv = std(peak_time,0,1);
                pair_peaks_err = pair_peaks_stdv/sqrt(size(peak_time,1));
                pair_times_err = pair_times_stdv/sqrt(size(peak_time,1));
                
                if offset == 0
                    save([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_data_Cb.mat'], 'peak_val', 'peak_time')
                else
                    save([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_data_Cb_C', num2str(offset), '.mat'], 'peak_val', 'peak_time')
                end
                
                x_axis = [1-pre_cycles:post_cycles];
                hold on
                errorbar(x_axis,pair_peaks_mean,pair_peaks_err,pair_peaks_err,zeros(1,pre_cycles+post_cycles),zeros(1,pre_cycles+post_cycles),'.')
                scatter(x_axis,pair_peaks_mean,'Marker', 's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 1], 'SizeData', 100)
                xlim([-4.5 5.5])
                hold off
                if offset == 0
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_prob_Cb.fig'])
                else
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_prob_Cb_C', num2str(offset), '.fig'])
                end
                close all
                
                hold on
                errorbar(x_axis,pair_times_mean,pair_times_err,pair_times_err,zeros(1,pre_cycles+post_cycles),zeros(1,pre_cycles+post_cycles),'.')
                scatter(x_axis,pair_times_mean,'Marker', 's', 'MarkerEdgeColor', [0 0 0], 'MarkerFaceColor', [0 0 1], 'SizeData', 100)
                xlim([-4.5 5.5])
                hold off
                if offset == 0
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_peak_times_Cb.fig'])
                else
                    saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_cofiring_peak_times_Cb_C', num2str(offset), '.fig'])
                end
                close all
            end
        end
    end
    if isunix
        rmpath(genpath('/common/fleischerp/HPC_code'))
        rmpath(genpath('/common/fleischerp/3rd_party_code'))
    else
        rmpath(genpath('Z:\Matlab for analysis\eeglab\functions'))
    end
	delete(gcp('nocreate'));
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Count Proportion of Spindles that are nested with Slow Oscillations (33)

if enabled(33)
    disp('Block 33...')
    M1_spindle_SO_dist = cell(param.days,param.sleep_blocks);
    Cb_spindle_SO_dist = cell(param.days,param.sleep_blocks);
    edges = 0:.25:10;
    bin_centers = edges(1:end-1) + ((edges(2) - edges(1))/2);
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Sleep_data.mat'])
            
            so_timestmps = M1_data.so_delta.so_up_states;
            spindle_timestmps = M1_data.spindles{1}.pks;
            diff_matrix =  spindle_timestmps - so_timestmps';
            leading_SO_dist = diff_matrix;
            leading_SO_dist(diff_matrix < 0) = inf;
            leading_SO_dist = min(leading_SO_dist,[],2);
            M1_spindle_SO_dist{day,block} = leading_SO_dist;
            
            dist_counts = histcounts(leading_SO_dist,edges);
            bar(bin_centers,dist_counts)
            saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/M1_Spindle_Nesting_Distance.fig'])
            close all
                
            so_timestmps = Cb_data.so_delta.so_up_states;
            spindle_timestmps = Cb_data.spindles{1}.pks;
            diff_matrix =  spindle_timestmps - so_timestmps';
            leading_SO_dist = diff_matrix;
            leading_SO_dist(diff_matrix < 0) = inf;
            leading_SO_dist = min(leading_SO_dist,[],2);
            Cb_spindle_SO_dist{day,block} = leading_SO_dist;
            
            dist_counts = histcounts(leading_SO_dist,edges);
            bar(bin_centers,dist_counts)
            saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Cb_Spindle_Nesting_Distance.fig'])
            close all
        end
    end
    load([rootpath,animal,'/Shared_Data.mat'])
    shared_data.M1_spindle_SO_dist = M1_spindle_SO_dist;
    shared_data.Cb_spindle_SO_dist = Cb_spindle_SO_dist;
    save([rootpath,animal,'/Shared_Data.mat'], 'shared_data')
    save([rootpath,animal,'/Nested_Spindles.mat'], 'M1_spindle_SO_dist', 'Cb_spindle_SO_dist')
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Spindle Ratios of Individual Channels (34)

if enabled(34)
    disp('Block 34...')
    day_range = 1:3;
    M1_all_days = nan(2,length(param.M1_good_chans),param.days);
    Cb_all_days = nan(2,length(param.Cb_good_chans),param.days);
    nesting_threshold = 1.5; %in seconds
    for day = 1:param.days
        load([rootpath,animal,'/Day',num2str(day),'/Sleep1/Channel_Sleep_data.mat'])
        M1_spindle_count = nan(2,size(M1_channel_data.spindles,2));
        Cb_spindle_count = nan(2,size(Cb_channel_data.spindles,2));
        for M1_chan = 1:size(M1_channel_data.spindles,2)
            nested_idxs = find_nested(M1_channel_data.so_delta(M1_chan).so_up_states,M1_channel_data.spindles{M1_chan}.pks,nesting_threshold);
            M1_spindle_count(1,M1_chan) = sum(nested_idxs);
        end
        for Cb_chan = 1:size(Cb_channel_data.spindles,2)
            nested_idxs = find_nested(Cb_channel_data.so_delta(Cb_chan).so_up_states,Cb_channel_data.spindles{Cb_chan}.pks,nesting_threshold);
            Cb_spindle_count(1,Cb_chan) = sum(nested_idxs);
        end
        sleep_minutes = sum(M1_channel_data.sleep_idx) / (M1_channel_data.Fs_LFP * 60);
        M1_spindle_count(1,:) = M1_spindle_count(1,:)/sleep_minutes;
        Cb_spindle_count(1,:) = Cb_spindle_count(1,:)/sleep_minutes;
            
        load([rootpath,animal,'/Day',num2str(day),'/Sleep2/Channel_Sleep_data.mat'])
        for M1_chan = 1:size(M1_channel_data.spindles,2)
            nested_idxs = find_nested(M1_channel_data.so_delta(M1_chan).so_up_states,M1_channel_data.spindles{M1_chan}.pks,nesting_threshold);
            M1_spindle_count(2,M1_chan) = sum(nested_idxs);
        end
        for Cb_chan = 1:size(Cb_channel_data.spindles,2)
            nested_idxs = find_nested(Cb_channel_data.so_delta(Cb_chan).so_up_states,Cb_channel_data.spindles{Cb_chan}.pks,nesting_threshold);
            Cb_spindle_count(2,Cb_chan) = sum(nested_idxs);
        end
        sleep_minutes = sum(M1_channel_data.sleep_idx) / (M1_channel_data.Fs_LFP * 60);
        M1_spindle_count(2,:) = M1_spindle_count(2,:)/sleep_minutes;
        Cb_spindle_count(2,:) = Cb_spindle_count(2,:)/sleep_minutes;
        M1_all_days(:,:,day) = M1_spindle_count;
        Cb_all_days(:,:,day) = Cb_spindle_count;
        
        hold on
        scatter(M1_spindle_count(1,:),M1_spindle_count(2,:));
        line([0,20],[0,20])
        saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/M1_Nested_Spindle_S1vS2_Rate.fig'])
        close all
        
        hold on
        scatter(Cb_spindle_count(1,:),Cb_spindle_count(2,:));
        line([0,20],[0,20])
        saveas(gcf,[rootpath,animal,'/Day',num2str(day),'/Cb_Nested_Spindle_S1vS2_Rate.fig'])
        close all
    end
    M1_comb_data = nan(2,0);
    Cb_comb_data = nan(2,0);
    for day = day_range
        M1_comb_data = cat(2, M1_comb_data, M1_all_days(:,:,day));
        Cb_comb_data = cat(2, Cb_comb_data, Cb_all_days(:,:,day));
    end
    
    hold on
    scatter(M1_comb_data(1,:),M1_comb_data(2,:));
    line([0,20],[0,20])
    saveas(gcf,[rootpath,animal,'/M1_Spindle_S1vS2_Rate.fig'])
    close all
    
    hold on
    scatter(Cb_comb_data(1,:),Cb_comb_data(2,:));
    line([0,20],[0,20])
    saveas(gcf,[rootpath,animal,'/Cb_Spindle_S1vS2_Rate.fig'])
    close all
    
    load([rootpath,animal,'/Shared_Data.mat'])
    shared_data.M1_spindle_rate = M1_all_days;
    shared_data.Cb_spindle_rate = Cb_all_days;
    save([rootpath,animal,'/Shared_Data.mat'], 'shared_data')
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Indentify Co-occuring Slow Oscilations and Delta Waves (35)

if enabled(35)
    disp('Block 35...')
    proximity_threshold = 0.1; % in seconds
    delta_ratios = nan(param.days, param.sleep_blocks);
    SO_ratios = nan(param.days, param.sleep_blocks);
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Sleep_data.mat'])
            dists = abs(M1_data.so_delta.so_up_states - Cb_data.so_delta.so_up_states');
            if length(M1_data.so_delta.so_up_states) > length(Cb_data.so_delta.so_up_states)
                minimums = min(dists,[],1)';
                cooccur_idxs = minimums < proximity_threshold;
                co_SO = Cb_data.so_delta.so_up_states(cooccur_idxs);
            else
                minimums = min(dists,[],2);
                cooccur_idxs = minimums < proximity_threshold;
                co_SO = M1_data.so_delta.so_up_states(cooccur_idxs);
            end
            SO_ratios(day,block) = length(co_SO)/length(minimums);
            
            dists = abs(M1_data.so_delta.delta_up_states - Cb_data.so_delta.delta_up_states');
            if length(M1_data.so_delta.delta_up_states) > length(Cb_data.so_delta.delta_up_states)
                minimums = min(dists,[],1)';
                cooccur_idxs = minimums < proximity_threshold;
                co_delta = Cb_data.so_delta.delta_up_states(cooccur_idxs);
            else
                minimums = min(dists,[],2);
                cooccur_idxs = minimums < proximity_threshold;
                co_delta = M1_data.so_delta.delta_up_states(cooccur_idxs);
            end
            delta_ratios(day,block) = length(co_delta)/length(minimums);
            
            save([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_data.mat'], 'M1_data', 'Cb_data', 'co_SO', 'co_delta')
        end
    end
    bar(1:param.days, delta_ratios)
    saveas(gcf,[rootpath,animal,'/cooccuring_delta_ratios.fig'])
    close all;
    
    bar(1:param.days, SO_ratios)
    saveas(gcf,[rootpath,animal,'/cooccuring_SO_ratios.fig'])
    close all;
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Spindle-centered Spike-field coherence (36)

if enabled(36)
    disp('Block 36...')
    addpath(genpath('Z:/Matlab for analysis/chronux_2_10/chronux/spectral_analysis'))
    addpath(genpath('Z:/Matlab for analysis/eeglab/functions'))
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Sleep_data.mat'])
            load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_timestamps.mat']);
            %Create spindle-centered LFP snapshots
            LFP_snapshots = create_LFP_snapshots(M1_data.LFP',M1_data.spindles{1}.pks,[-4, 4],M1_data.Fs_LFP);
            LFP_cat = LFP_snapshots(:)';
            %Create spindle-centered spiking snapshots
            spike_snapshots = nan(length(M1_spike_timestamps(:)),size(LFP_snapshots,2),size(LFP_snapshots,3));
            
            for n = 1:length(spike_snapshots)
                spike_snapshots(n,:,:) = create_spike_bin_snapshots(M1_spike_timestamps{n},M1_data.spindles{1}.pks,[-4, 4],M1_data.Fs_LFP);
                spikes_cat = spike_snapshots(n,:,:);
                spikes_cat = spikes_cat(:)';
                %Call newcrossf() on LFP and spike snapshots
                [cross_trial_coh,~,coh_times,coh_freqs,~,~,cross_spec,x_pspec,y_pspec] = ...
                    newcrossf(LFP_cat, spikes_cat, size(LFP_snapshots,2), [-4000 4000], M1_data.Fs_LFP, [0.01 0.1], 'type', 'coher', 'freqs', [0 60]);
            end
            
        end
    end
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Sleep Reactivation of Reach Principal Components (37)

if enabled(37)
    disp('Block 37...')
    addpath(genpath('Z:\Matlab for analysis\PCA_cellassembly'))
    M1_neurons_of_interest = cellfun(@not,cellfun(@isnan,param.M1_neuron_chans,'UniformOutput', false),'UniformOutput', false);
    Cb_neurons_of_interest = cellfun(@not,cellfun(@isnan,param.Cb_neuron_chans,'UniformOutput', false),'UniformOutput', false);
    bin_width = 1; %in Fs
    %load([rootpath,animal,'/Shared_Data.mat']);
    %shared_data.M1_reach_assembly_reactivation = cell(param.days,2);
    %shared_data.Cb_reach_assembly_reactivation = cell(param.days,2);
    %shared_data.M1_reach_assembly_eigen_ratios = cell(param.days,1);
    %shared_data.Cb_reach_assembly_eigen_ratios = cell(param.days,1);
    for day = 1:param.days
        day_neurons_of_interest_M1 = find(M1_neurons_of_interest{day}(:));
        day_neurons_of_interest_Cb = find(Cb_neurons_of_interest{day}(:));
        for block = 1:2
            load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Sleep_data.mat'])
            load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_timestamps.mat']);
            load([rootpath,animal,'/Day',num2str(day),'/',param.block_names{1},'/PC_reach_patterns.mat']);
            
            if ~isempty(day_neurons_of_interest_M1)
                edges = 0.5:(length(M1_data.LFP)+0.5);
                edges = edges/M1_data.Fs_LFP;
                neuron_idx = 0;
                sleep_cat_hist = nan(length(day_neurons_of_interest_M1), length(M1_data.LFP));
                for full_neuron_idx = day_neurons_of_interest_M1'
                    neuron_idx = neuron_idx+1;
                    [chan, neuron] = ind2sub(size(M1_neurons_of_interest{day}), full_neuron_idx);
                    sleep_cat_hist(neuron_idx,:) = histcounts(M1_spike_timestamps{chan,neuron}, edges);
                end
                sleep_cat_hist = sleep_cat_hist(:,M1_data.sleep_idx);
                eval(['sleep' num2str(block) '_cat_hist_M1 = sleep_cat_hist;']);
            end
                
            if ~isempty(day_neurons_of_interest_Cb)
                edges = 0.5:(length(Cb_data.LFP)+0.5);
                edges = edges/Cb_data.Fs_LFP;
                neuron_idx = 0;
                sleep_cat_hist = nan(length(day_neurons_of_interest_Cb), length(Cb_data.LFP));
                for full_neuron_idx = day_neurons_of_interest_Cb'
                    neuron_idx = neuron_idx+1;
                    [chan, neuron] = ind2sub(size(Cb_neurons_of_interest{day}), full_neuron_idx);
                    sleep_cat_hist(neuron_idx,:) = histcounts(Cb_spike_timestamps{chan,neuron}, edges);
                end
                sleep_cat_hist = sleep_cat_hist(:,Cb_data.sleep_idx);
                eval(['sleep' num2str(block) '_cat_hist_Cb = sleep_cat_hist;']);
            end
        end
        
        if ~isempty(day_neurons_of_interest_M1)
            Activities_pre_M1 = assembly_activity(M1_patterns,sleep1_cat_hist_M1);
            Activities_post_M1 = assembly_activity(M1_patterns,sleep2_cat_hist_M1);
            Activities_pre_M1 = abs(Activities_pre_M1);
            Activities_post_M1 = abs(Activities_post_M1);
            %shared_data.M1_reach_assembly_reactivation{day,1} = Activities_pre_M1;
            %shared_data.M1_reach_assembly_reactivation{day,2} = Activities_post_M1;
            %shared_data.M1_reach_assembly_eigen_ratios{day} = M1_eigen_lambda_ratio;
            
            for i = 1:(size(M1_patterns,2))
                figure(1); clf;
                subplot(311)
                plot(Activities_pre_M1(i,:))
                %     ylim([-50 50])
                title('Ensemble Strength Pre-Learning Sleep Block');
                
                subplot(312)
                plot(Activities_post_M1(i,:))
                %     ylim([-50 50])
                title('Ensemble Strength Post-Learning Sleep Block');
                
                subplot(313)
                
                Mn(i,1) = mean(Activities_pre_M1(i,:));
                Mn(i,2)= mean(Activities_post_M1(i,:));
                Err(i,1) = std(Activities_pre_M1(i,:))/sqrt(length(Activities_pre_M1(i,:)));
                Err(i,2) = std(Activities_post_M1(i,:))/sqrt(length(Activities_pre_M1(i,:)));
                [~,p]=ttest2(Activities_pre_M1(i,:),Activities_post_M1(i,:)); %SIG TEST
                
                errorbar(Mn(i,:),Err(i,:));
                title('Mean Change in Activation Strength');
                pval = strcat('p = ',num2str(p));
                text(1.8,Mn(1),pval);
                filename=[rootpath,animal,'/Day',num2str(day),'/Pre-Post_Activation_Ensemble',num2str(i) '_M1.tiff'];
                saveas(gcf,filename);
                close all
                %     figure;
                max_pre=max(Activities_pre_M1(i,:));
                max_post=max(Activities_post_M1(i,:));
                
                maxY=350;
                figure(2)
                subplot(221);plot(abs(Activities_pre_M1(i,:)));axis([0 length(Activities_pre_M1) 0 maxY]);title('Pre Sleep Activation')
                subplot(222);plot(abs(Activities_post_M1(i,:)));axis([0 length(Activities_post_M1) 0 maxY]);title('Post Sleep Activation')
                pre_hist=hist(abs(Activities_pre_M1(i,:)),[0:1:350]);
                post_hist=hist(abs(Activities_post_M1(i,:)),[0:1:350]);
                
                
                subplot(223);bar(hist(abs(Activities_pre_M1(i,:)),[0:1:300]));axis([-10 200 0 500]);
                subplot(224);bar(hist(abs(Activities_post_M1(i,:)),[0:1:300]));axis([-10 200 0 500]);
                filename1=[rootpath,animal,'/Day',num2str(day),'/Pre-Post_Activation_Pre-Post_Activation',num2str(i) '_M1.tiff'];
                saveas(gcf,filename1);
                close all
                
                figure(3)
                plot([cumsum(post_hist)-cumsum(pre_hist)],'b');
                filename2=[rootpath,animal,'/Day',num2str(day),'/Pre-Post_Activation_Bin_Pre-Post_Diff',num2str(i) '_M1.tiff'];
                saveas(gcf,filename2);
                close all
                
                % q = size(zSpikeCount,1)/size(zSpikeCount,2);
                % lambda_max = ((1+sqrt(1/q))^2);
                
                %     filesav = strcat('Pre-Post Activation Bin',num2str(bin),'Ensemble',num2str(i));
                %     saveas(gcf,filesav,'tiff');
                %     close all;
            end
        else
            %shared_data.M1_reach_assembly_reactivation{day,1} = [];
            %shared_data.M1_reach_assembly_reactivation{day,2} = [];
            %shared_data.M1_reach_assembly_eigen_ratios{day} = [];
        end
        
        if ~isempty(day_neurons_of_interest_Cb)
            Activities_pre_Cb = assembly_activity(Cb_patterns,sleep1_cat_hist_Cb);
            Activities_post_Cb = assembly_activity(Cb_patterns,sleep2_cat_hist_Cb);
            Activities_pre_Cb = abs(Activities_pre_Cb);
            Activities_post_Cb = abs(Activities_post_Cb);
            %shared_data.Cb_reach_assembly_reactivation{day,1} = Activities_pre_Cb;
            %shared_data.Cb_reach_assembly_reactivation{day,2} = Activities_post_Cb;
            %shared_data.Cb_reach_assembly_eigen_ratios{day} = Cb_eigen_lambda_ratio;
            
            for i = 1:(size(Cb_patterns,2))
                
                figure(1); clf;
                subplot(311)
                plot(Activities_pre_Cb(i,:))
                %     ylim([-50 50])
                title('Ensemble Strength Pre-Learning Sleep Block');
                
                subplot(312)
                plot(Activities_post_Cb(i,:))
                %     ylim([-50 50])
                title('Ensemble Strength Post-Learning Sleep Block');
                
                subplot(313)
                
                Mn(i,1) = mean(Activities_pre_Cb(i,:));
                Mn(i,2)= mean(Activities_post_Cb(i,:));
                Err(i,1) = std(Activities_pre_Cb(i,:))/sqrt(length(Activities_pre_Cb(i,:)));
                Err(i,2) = std(Activities_post_Cb(i,:))/sqrt(length(Activities_pre_Cb(i,:)));
                [~,p]=ttest2(Activities_pre_Cb(i,:),Activities_post_Cb(i,:)); %SIG TEST
                
                errorbar(Mn(i,:),Err(i,:));
                title('Mean Change in Activation Strength');
                pval = strcat('p = ',num2str(p));
                text(1.8,Mn(1),pval);
                filename=[rootpath,animal,'/Day',num2str(day),'/Pre-Post_Activation_Ensemble',num2str(i) '_Cb.tiff'];
                saveas(gcf,filename);
                close all
                %     figure;
                max_pre=max(Activities_pre_Cb(i,:));
                max_post=max(Activities_post_Cb(i,:));
                
                maxY=350;
                figure(2)
                subplot(221);plot(abs(Activities_pre_Cb(i,:)));axis([0 length(Activities_pre_Cb) 0 maxY]);title('Pre Sleep Activation')
                subplot(222);plot(abs(Activities_post_Cb(i,:)));axis([0 length(Activities_post_Cb) 0 maxY]);title('Post Sleep Activation')
                pre_hist=hist(abs(Activities_pre_Cb(i,:)),[0:1:350]);
                post_hist=hist(abs(Activities_post_Cb(i,:)),[0:1:350]);
                
                
                subplot(223);bar(hist(abs(Activities_pre_Cb(i,:)),[0:1:300]));axis([-10 200 0 500]);
                subplot(224);bar(hist(abs(Activities_post_Cb(i,:)),[0:1:300]));axis([-10 200 0 500]);
                filename1=[rootpath,animal,'/Day',num2str(day),'/Pre-Post_Activation_Pre-Post_Activation',num2str(i) '_Cb.tiff'];
                saveas(gcf,filename1);
                close all
                
                figure(3)
                plot([cumsum(post_hist)-cumsum(pre_hist)],'b');
                filename2=[rootpath,animal,'/Day',num2str(day),'/Pre-Post_Activation_Bin_Pre-Post_Diff',num2str(i) '_Cb.tiff'];
                saveas(gcf,filename2);
                close all
                
                % q = size(zSpikeCount,1)/size(zSpikeCount,2);
                % lambda_max = ((1+sqrt(1/q))^2);
                
                %     filesav = strcat('Pre-Post Activation Bin',num2str(bin),'Ensemble',num2str(i));
                %     saveas(gcf,filesav,'tiff');
                %     close all;
            end
        else
            %shared_data.Cb_reach_assembly_reactivation{day,1} = [];
            %shared_data.Cb_reach_assembly_reactivation{day,2} = [];
            %shared_data.Cb_reach_assembly_eigen_ratios{day} = [];
        end
        
    end
    %save([rootpath,animal,'/Shared_Data.mat'], 'shared_data')
    rmpath(genpath('Z:\Matlab for analysis\PCA_cellassembly'))
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Plot example co-firing histograms from single neuron-neuron-spindle triple in CB only (38)

if enabled(38)
    disp('Block 38...')
    pre_margin = 0.05; %in seconds
    post_margin = 0.05; %in seconds
    trace_margin = 1; %in seconds - Spindles closer than this to the edge of the recording will be skipped to avoid index out of bounds errors
    pre_cycles = 5;
    post_cycles = 5;
    shuffle_reps = 25;
    max_idx_err = 10;
    kernel_stdv = 5;
    control_offsets = [-5, -10];
    
    %area = 'Cb';
    day = 1;
    block = 1;
    %spindle_idx = 10;
    neuron1_idx = 1;
    neuron2_idx = 2;
    
    bin_edges = ((-(round(pre_margin*1000)+0.5)):(round(post_margin*1000)+0.5))/1000; %In seconds, for computation
    bin_centers = (-(round(pre_margin*1000))):(round(post_margin*1000)); %In miliseconds, for display
    
    load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Filtered_Sleep_data.mat'])
    load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_timestamps.mat']);
    
    hilbert_LFP = hilbert(Cb_filt_data.LFP);
    inst_phase = unwrap(angle(hilbert_LFP));%inst phase
    radial_phase = mod(inst_phase,2*pi);
    peaks = diff(radial_phase) < 0;
    %align spindle_indexes with a peak
    peak_idxs = find(peaks);
    
    for offset = [0, control_offsets]
        spindle_idxs = round(Cb_filt_data.Fs_LFP * (Cb_filt_data.spindles{1}.pks + offset));
        spindle_idxs(spindle_idxs < trace_margin * Cb_filt_data.Fs_LFP) = [];
        spindle_idxs(spindle_idxs > length(Cb_filt_data.LFP) - (trace_margin * Cb_filt_data.Fs_LFP)) = [];
        
        %get spindle phase data
        hilbert_LFP = hilbert(Cb_filt_data.LFP);
        inst_phase = unwrap(angle(hilbert_LFP));%inst phase
        radial_phase = mod(inst_phase,2*pi);
        peaks = diff(radial_phase) < 0;
        %align spindle_indexes with a peak
        peak_idxs = find(peaks);
        for i = 1:length(spindle_idxs)
            [idx_diff, spindle_idxs(i)] = min(abs(peak_idxs - spindle_idxs(i)));
            if idx_diff > max_idx_err && offset == 0
                warning('A spindle index had greater offset from a cycle peak than is allowed by the threshold');
            end
        end
        
        raw_CCH = nan((pre_cycles+post_cycles),round((pre_margin+post_margin)*1000)+1);
        neuron_idxs = find(param.Cb_task_related_neurons{day});
        Cb_neuron1 = neuron_idxs(neuron1_idx);
        Cb_neuron2 = neuron_idxs(neuron2_idx);
        pair_data = nan(pre_cycles+post_cycles,round((pre_margin+post_margin)*1000)+1);
        
        for cycle = (1-pre_cycles):post_cycles %zero-th cycle is the cycle imeadiatly preceding the center
            lag_hists = nan(0,round((pre_margin+post_margin)*1000)+1);
            for spindle = 1:length(spindle_idxs)
                
                cycle_start_idx = peak_idxs(spindle_idxs(spindle)+(cycle-1));
                cycle_end_idx = peak_idxs(spindle_idxs(spindle)+(cycle));
                
                Cb_cycle_spikes1 = Cb_spike_timestamps{Cb_neuron1}(Cb_spike_timestamps{Cb_neuron1} > (cycle_start_idx/param.Cb_Fs)  &  Cb_spike_timestamps{Cb_neuron1} < (cycle_end_idx/param.Cb_Fs));
                for Cb_spiketime1 = Cb_cycle_spikes1
                    Cb_spikes2 = Cb_spike_timestamps{Cb_neuron2}(Cb_spike_timestamps{Cb_neuron2} >= (Cb_spiketime1 - pre_margin) & Cb_spike_timestamps{Cb_neuron2} <= (Cb_spiketime1 + post_margin));
                    Cb_spikes2 = Cb_spikes2 - Cb_spiketime1; %It should be 2 - 1 for both
                    lag_hists = cat(1,lag_hists,histcounts(Cb_spikes2,bin_edges));
                end
                
                Cb_cycle_spikes2 = Cb_spike_timestamps{Cb_neuron2}(Cb_spike_timestamps{Cb_neuron2} > (cycle_start_idx/param.Cb_Fs)  &  Cb_spike_timestamps{Cb_neuron2} < (cycle_end_idx/param.Cb_Fs));
                for Cb_spiketime2 = Cb_cycle_spikes2
                    Cb_spikes1 = Cb_spike_timestamps{Cb_neuron1}(Cb_spike_timestamps{Cb_neuron1} >= (Cb_spiketime2 - pre_margin) & Cb_spike_timestamps{Cb_neuron1} <= (Cb_spiketime2 + post_margin));
                    Cb_spikes1 = Cb_spiketime2 - Cb_spikes1; %It should be 2 - 1 for both
                    lag_hists = cat(1,lag_hists,histcounts(Cb_spikes1,bin_edges));
                end
            end
            %average hist arrays across spindles and store it
            c_idx = cycle + pre_cycles;
            if isempty(lag_hists)
                pair_data(cycle + pre_cycles,:) = 0;
            else
                pair_data(cycle + pre_cycles,:) = mean(lag_hists,1);
            end
        end
        raw_CCH(:,:) = pair_data;
        
        shuffled_CCH = nan((pre_cycles+post_cycles),round((pre_margin+post_margin)*1000)+1,shuffle_reps);
        shuffled_spindle_idxs = nan(shuffle_reps,length(spindle_idxs));
        for shuffle = 1:shuffle_reps
            shuffled_spindle_idxs(shuffle,:) = randperm(length(spindle_idxs));
        end
        neuron_idxs = find(param.Cb_task_related_neurons{day});
        good_pairs = (neuron_idxs - neuron_idxs') > 0;
        pair_idxs = find(good_pairs);
        pair_data = nan(pre_cycles+post_cycles,round((pre_margin+post_margin)*1000)+1);
        
        for cycle = (1-pre_cycles):post_cycles %zero-th cycle is the cycle imeadiatly preceding the center
            for shuffle = 1:shuffle_reps
                lag_hists = nan(0,round((pre_margin+post_margin)*1000)+1);
                for spindle = 1:length(spindle_idxs)
                    
                    cycle_start_idx = peak_idxs(spindle_idxs(spindle)+(cycle-1));
                    cycle_end_idx = peak_idxs(spindle_idxs(spindle)+(cycle));
                    
                    
                    Cb_cycle_spikes1 = Cb_spike_timestamps{Cb_neuron1}(Cb_spike_timestamps{Cb_neuron1} > (cycle_start_idx/param.Cb_Fs)  &  Cb_spike_timestamps{Cb_neuron1} < (cycle_end_idx/param.Cb_Fs));
                    Cb_cycle_spikes1 = Cb_cycle_spikes1 + (peak_idxs(spindle_idxs(shuffled_spindle_idxs(shuffle,spindle)))/param.Cb_Fs) - (peak_idxs(spindle_idxs(spindle))/param.Cb_Fs);
                    for Cb_spiketime1 = Cb_cycle_spikes1
                        Cb_spikes2 = Cb_spike_timestamps{Cb_neuron2}(Cb_spike_timestamps{Cb_neuron2} >= (Cb_spiketime1 - pre_margin) & Cb_spike_timestamps{Cb_neuron2} <= (Cb_spiketime1 + post_margin));
                        Cb_spikes2 = Cb_spikes2 - Cb_spiketime1; %It should be 2 - 1 for both
                        lag_hists = cat(1,lag_hists,histcounts(Cb_spikes2,bin_edges));
                    end
                    
                    Cb_cycle_spikes2 = Cb_spike_timestamps{Cb_neuron2}(Cb_spike_timestamps{Cb_neuron2} > (cycle_start_idx/param.Cb_Fs)  &  Cb_spike_timestamps{Cb_neuron2} < (cycle_end_idx/param.Cb_Fs));
                    Cb_cycle_spikes2 = Cb_cycle_spikes2 + (peak_idxs(spindle_idxs(shuffled_spindle_idxs(shuffle,spindle)))/param.Cb_Fs) - (peak_idxs(spindle_idxs(spindle))/param.Cb_Fs);
                    for Cb_spiketime2 = Cb_cycle_spikes2
                        Cb_spikes1 = Cb_spike_timestamps{Cb_neuron1}(Cb_spike_timestamps{Cb_neuron1} >= (Cb_spiketime2 - pre_margin) & Cb_spike_timestamps{Cb_neuron1} <= (Cb_spiketime2 + post_margin));
                        Cb_spikes1 = Cb_spiketime2 - Cb_spikes1; %It should be 2 - 1 for both
                        lag_hists = cat(1,lag_hists,histcounts(Cb_spikes1,bin_edges));
                    end
                end
                %average hist arrays across spindles and store it
                if isempty(lag_hists)
                    pair_data(cycle + pre_cycles,:) = 0;
                else
                    pair_data(cycle + pre_cycles,:) = mean(lag_hists,1);
                end
                shuffled_CCH(:,:,shuffle) = pair_data;
            end
        end
        ave_shuffled_CCH = mean(shuffled_CCH,3);
        corrected_CCH = raw_CCH-ave_shuffled_CCH;
        
        raw_peak_CCH = sum(raw_CCH([5,6],:),1);
        shuf_peak_CCH = sum(ave_shuffled_CCH([5,6],:),1);
        corr_peak_CCH = sum(corrected_CCH([5,6],:),1);
        
        bar(bin_centers, raw_peak_CCH)
        filename=[rootpath,animal,'/Cofiring_Raw_Peak_CCH_Cb_offset_', num2str(offset), '.fig'];
        saveas(gcf,filename)
        close all
        bar(bin_centers, shuf_peak_CCH)
        filename=[rootpath,animal,'/Cofiring_Shuf_Peak_CCH_Cb_offset_', num2str(offset), '.fig'];
        saveas(gcf,filename)
        close all
        bar(bin_centers, corr_peak_CCH)
        filename=[rootpath,animal,'/Cofiring_Corr_Peak_CCH_Cb_offset_', num2str(offset), '.fig'];
        saveas(gcf,filename)
        close all
        
        raw_tail_CCH = sum(raw_CCH([1,10],:),1);
        shuf_tail_CCH = sum(ave_shuffled_CCH([1,10],:),1);
        corr_tail_CCH = sum(corrected_CCH([1,10],:),1);
        
        bar(bin_centers, raw_tail_CCH)
        filename=[rootpath,animal,'/Cofiring_Raw_Tail_CCH_Cb_offset_', num2str(offset), '.fig'];
        saveas(gcf,filename)
        close all
        bar(bin_centers, shuf_tail_CCH)
        filename=[rootpath,animal,'/Cofiring_Shuf_Tail_CCH_Cb_offset_', num2str(offset), '.fig'];
        saveas(gcf,filename)
        close all
        bar(bin_centers, corr_tail_CCH)
        filename=[rootpath,animal,'/Cofiring_Corr_Tail_CCH_Cb_offset_', num2str(offset), '.fig'];
        saveas(gcf,filename)
        close all
    end
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Compare percentage of neurons that are modulated at spindel tails vs. at spindel peaks (39)

if enabled(39)
    disp('Block 39...')
    addpath('Z:\Matlab for analysis\CircStat2012a');
    control_offsets = [-5, -10];
    group_names = {'Tail' 'Peak'};
    group_cycles = [[1;10] [5;6]];%[[1;2;9;10] [4;5;6;7]];
    c_group_cycles = [5;6];%[4;5;6;7];
    
    edge_exp = -6:0.1:0;
    edges = exp(edge_exp);
    bin_centers = (diff(edges)/2) + edges(1:end-1);

    M1_Pval = nan(1,0);
    M1_Tstat = nan(1,0);
    Cb_Pval = nan(1,0);
    Cb_Tstat = nan(1,0);
    
    for day = 1:param.days
        M1_phases = cell([size(param.M1_task_related_neurons{1},1), size(param.M1_task_related_neurons{day},2), size(group_cycles,2)]);
        Cb_phases = cell([size(param.Cb_task_related_neurons{1},1), size(param.Cb_task_related_neurons{day},2), size(group_cycles,2)]);
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_spike_phasestamps.mat'])
            for group = 1:size(group_cycles,2)
                for cycle = group_cycles(:,group)'
                    if cycle > size(prePeak_cycle_phases_M1_spindles_M1_spikes,3)
                        %postPeak
                        M1_phases(:,:,group) = cellfun(@horzcat,M1_phases(:,:,group),postPeak_cycle_phases_M1_spindles_M1_spikes(:,:,cycle - size(prePeak_cycle_phases_M1_spindles_M1_spikes,3)),'UniformOutput',false);
                        Cb_phases(:,:,group) = cellfun(@horzcat,Cb_phases(:,:,group),postPeak_cycle_phases_Cb_spindles_Cb_spikes(:,:,cycle - size(prePeak_cycle_phases_Cb_spindles_Cb_spikes,3)),'UniformOutput',false);
                    else
                        %prePeak
                        M1_phases(:,:,group) = cellfun(@horzcat,M1_phases(:,:,group),prePeak_cycle_phases_M1_spindles_M1_spikes(:,:,cycle),'UniformOutput',false);
                        Cb_phases(:,:,group) = cellfun(@horzcat,Cb_phases(:,:,group),prePeak_cycle_phases_Cb_spindles_Cb_spikes(:,:,cycle),'UniformOutput',false);
                    end
                end
            end
        end
        
        M1_phases_c = cell([size(param.M1_task_related_neurons{1},1), size(param.M1_task_related_neurons{day},2)]);
        Cb_phases_c = cell([size(param.Cb_task_related_neurons{1},1), size(param.Cb_task_related_neurons{day},2)]);
        for offset = control_offsets
            for block = 1:param.sleep_blocks
                load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spindle_spike_phasestamps_C', num2str(offset), '.mat'])
                for cycle = c_group_cycles'
                    if cycle > size(prePeak_cycle_phases_M1_spindles_M1_spikes,3)
                        %postPeak
                        M1_phases_c(:,:) = cellfun(@horzcat,M1_phases_c(:,:),postPeak_cycle_phases_M1_spindles_M1_spikes(:,:,cycle - size(prePeak_cycle_phases_M1_spindles_M1_spikes,3)),'UniformOutput',false);
                        Cb_phases_c(:,:) = cellfun(@horzcat,Cb_phases_c(:,:),postPeak_cycle_phases_Cb_spindles_Cb_spikes(:,:,cycle - size(prePeak_cycle_phases_Cb_spindles_Cb_spikes,3)),'UniformOutput',false);
                    else
                        %prePeak
                        M1_phases_c = cellfun(@horzcat,M1_phases_c,prePeak_cycle_phases_M1_spindles_M1_spikes(:,:,cycle),'UniformOutput',false);
                        Cb_phases_c = cellfun(@horzcat,Cb_phases_c,prePeak_cycle_phases_Cb_spindles_Cb_spikes(:,:,cycle),'UniformOutput',false);
                    end
                end
            end
        end
    
        M1_Pval_day = nan(size(group_cycles,2),sum(param.M1_task_related_neurons{day}(:)));
        M1_Tstat_day = nan(size(group_cycles,2),sum(param.M1_task_related_neurons{day}(:)));
        M1_neuron_idx = 0;
        for chan = 1:param.M1_chans
            for neuron = 1:param.M1_neurons
                if param.M1_task_related_neurons{day}(chan,neuron)
                    undersized = cellfun(@length,M1_phases(chan,neuron,:)) < 10;
                    if sum([undersized(:)', (length(M1_phases_c{chan,neuron}) < 10)]) < 1
                        M1_neuron_idx = M1_neuron_idx+1;
                        for group = 1:size(group_cycles,2)
                            test_success = false;
                            sub_set_ratio = 1.0;
                            while ~test_success
                                phases_1 = M1_phases{chan,neuron,group};
                                phases_2 = M1_phases_c{chan,neuron};
                                phases_1 = phases_1(randperm(length(phases_1),round(length(phases_1)*sub_set_ratio)));
                                phases_2 = phases_2(randperm(length(phases_2),round(length(phases_2)*sub_set_ratio)));
                                try
                                    [M1_Pval_day(group, M1_neuron_idx), ~, M1_Tstat_day(group, M1_neuron_idx)] = circ_cmtest(phases_1, phases_2);
                                    test_success = true;
                                catch
                                    sub_set_ratio = sub_set_ratio/2;
                                    test_success = false;
                                end
                            end
                        end
                    end
                end
            end
        end
        
        Cb_Pval_day = nan(size(group_cycles,2),sum(param.Cb_task_related_neurons{day}(:)));
        Cb_Tstat_day = nan(size(group_cycles,2),sum(param.Cb_task_related_neurons{day}(:)));
        Cb_neuron_idx = 0;
        for chan = 1:param.Cb_chans
            for neuron = 1:param.Cb_neurons
                if param.Cb_task_related_neurons{day}(chan,neuron)
                    undersized = cellfun(@length,Cb_phases(chan,neuron,:)) < 10;
                    if sum([undersized(:)', (length(Cb_phases_c{chan,neuron}) < 10)]) < 1
                        Cb_neuron_idx = Cb_neuron_idx+1;
                        for group = 1:size(group_cycles,2)
                            test_success = false;
                            sub_set_ratio = 1.0;
                            while ~test_success
                                phases_1 = Cb_phases{chan,neuron,group};
                                phases_2 = Cb_phases_c{chan,neuron};
                                phases_1 = phases_1(randperm(length(phases_1),round(length(phases_1)*sub_set_ratio)));
                                phases_2 = phases_2(randperm(length(phases_2),round(length(phases_2)*sub_set_ratio)));
                                try
                                    [Cb_Pval_day(group, Cb_neuron_idx), ~, Cb_Tstat_day(group, Cb_neuron_idx)] = circ_cmtest(phases_1, phases_2);
                                    test_success = true;
                                catch
                                    sub_set_ratio = sub_set_ratio/2;
                                    test_success = false;
                                end
                            end
                        end
                    end
                end
            end
        end
        
        M1_Pval = [M1_Pval M1_Pval_day];
        M1_Tstat = [M1_Tstat M1_Tstat_day];
        Cb_Pval = [Cb_Pval Cb_Pval_day];
        Cb_Tstat = [Cb_Tstat Cb_Tstat_day];
    end
    
    %M1
    %Tail Line
    [p_counts, ~] = histcounts(M1_Pval(1,:),edges);
    cum_count = cumsum(p_counts);
    p_ratio_cum = cum_count/size(M1_Pval,2);
    line(bin_centers,p_ratio_cum','Color',[0.5, 0.5, 1]);
    %Peak Line
    [p_counts, ~] = histcounts(M1_Pval(2,:),edges);
    cum_count = cumsum(p_counts);
    p_ratio_cum = cum_count/size(M1_Pval,2);
    line(bin_centers,p_ratio_cum','Color',[0, 0, 1]);
    
    line([.05 .05], [0 1], 'Color', [0 0 0], 'LineStyle','--');
    set(gca, 'XScale', 'log')
    xlim([bin_centers(1) bin_centers(end)])
    saveas(gcf,[rootpath,animal,'/Different_Median_Spindel_Spike-phase_M1', '.fig']);
    close all;
    
    %Cb
    %Tail Line
    [p_counts, ~] = histcounts(Cb_Pval(1,:),edges);
    cum_count = cumsum(p_counts);
    p_ratio_cum = cum_count/size(Cb_Pval,2);
    line(bin_centers,p_ratio_cum','Color',[0.5, 0.5, 1]);
    %Peak Line
    [p_counts, ~] = histcounts(Cb_Pval(2,:),edges);
    cum_count = cumsum(p_counts);
    p_ratio_cum = cum_count/size(Cb_Pval,2);
    line(bin_centers,p_ratio_cum','Color',[0, 0, 1]);
    
    line([.05 .05], [0 1], 'Color', [0 0 0], 'LineStyle','--');
    set(gca, 'XScale', 'log')
    xlim([bin_centers(1) bin_centers(end)])
    saveas(gcf,[rootpath,animal,'/Different_Median_Spindel_Spike-phase_Cb', '.fig']);
    close all;
    
    rmpath('Z:\Matlab for analysis\CircStat2012a');
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Make Spindel-centered Rastors (40)

if enabled(40)
    disp('Block 40...')
    addpath('Z:\Matlab for analysis\CircStat2012a');
    cross_area_bool = true;
    
    spindle_band = [10,15];
    pre_margin = 0.5; %in seconds - Spindles closer than this to the edge of the recording will be skipped to avoid index out of bounds errors
    post_margin = 0.5; %in seconds - Spindles closer than this to the edge of the recording will be skipped to avoid index out of bounds errors
    lfp_pre = round(pre_margin*param.M1_Fs);
    lfp_post = round(post_margin*param.M1_Fs);
    bin_width = .001; %in seconds
    M1_neuron = 1;
    Cb_neuron = 2;
    task_related_neurons_only = true;
    max_to_plot = 100;
    peak_range = [0 1];
    tail_range = [4 5];
    max_idx_err = 10;
    
    x_axis_LFP = round(-pre_margin*param.M1_Fs):round(post_margin*param.M1_Fs);
    x_axis_LFP = x_axis_LFP*1000/param.M1_Fs;
    x_axis_spike = round(-pre_margin*1000):round(bin_width*1000):round(post_margin*1000);
    x_axis_spike = x_axis_spike(1:(end-1)) + round(bin_width*500);
    
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Filtered_Sleep_data.mat'])
            load([rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Spike_timestamps.mat']);
            if cross_area_bool
                temp = M1_filt_data;
                M1_filt_data = Cb_filt_data;
                Cb_filt_data = temp;
                clear temp
            end
            
            spindle_peak_idxs_M1 = round(M1_filt_data.Fs_LFP * (M1_filt_data.spindles{1}.pks));
            spindle_peak_idxs_M1(spindle_peak_idxs_M1 < pre_margin * M1_filt_data.Fs_LFP) = [];
            spindle_peak_idxs_M1(spindle_peak_idxs_M1 > length(M1_filt_data.LFP) - (post_margin * M1_filt_data.Fs_LFP)) = [];
            
            %get spindle phase data
            hilbert_LFP = hilbert(M1_filt_data.LFP);
            inst_phase = unwrap(angle(hilbert_LFP));%inst phase
            radial_phase_M1 = mod(inst_phase,2*pi);
            peaks = diff(radial_phase_M1) < 0;
            %align spindle_indexes with a peak
            peak_idxs_M1 = find(peaks);
            for i = 1:length(spindle_peak_idxs_M1)
                [idx_diff, spindle_peak_idxs_M1(i)] = min(abs(peak_idxs_M1 - spindle_peak_idxs_M1(i)));
                if idx_diff > max_idx_err && offset == 0
                    warning('A spindle index had greater offset from a cycle peak than is allowed by the threshold');
                end
            end
            
            spindle_peak_idxs_Cb = round(Cb_filt_data.Fs_LFP * (Cb_filt_data.spindles{1}.pks));
            spindle_peak_idxs_Cb(spindle_peak_idxs_Cb < pre_margin * Cb_filt_data.Fs_LFP) = [];
            spindle_peak_idxs_Cb(spindle_peak_idxs_Cb > length(Cb_filt_data.LFP) - (post_margin * Cb_filt_data.Fs_LFP)) = [];
            
            %get spindle phase data
            hilbert_LFP = hilbert(Cb_filt_data.LFP);
            inst_phase = unwrap(angle(hilbert_LFP));%inst phase
            radial_phase_Cb = mod(inst_phase,2*pi);
            peaks = diff(radial_phase_Cb) < 0;
            %align spindle_indexes with a peak
            peak_idxs_Cb = find(peaks);
            for i = 1:length(spindle_peak_idxs_Cb)
                [idx_diff, spindle_peak_idxs_Cb(i)] = min(abs(peak_idxs_Cb - spindle_peak_idxs_Cb(i)));
                if idx_diff > max_idx_err && offset == 0
                    warning('A spindle index had greater offset from a cycle peak than is allowed by the threshold');
                end
            end
            
            
            if task_related_neurons_only
                M1_neurons = param.M1_task_related_neurons{day};
                Cb_neurons = param.Cb_task_related_neurons{day};
            else
                M1_neurons = ~cellfun(@isempty,M1_spike_timestamps);
                Cb_neurons = ~cellfun(@isempty,Cb_spike_timestamps);
            end
            M1_neuron = find(M1_neurons, M1_neuron);
            if ~isempty(M1_neuron)
                M1_neuron = M1_neuron(end);
            end
            Cb_neuron = find(Cb_neurons, Cb_neuron);
            if ~isempty(Cb_neuron)
                Cb_neuron = Cb_neuron(end);
            end
            
            all_rhythms_M1 = M1_filt_data.spindles{1}.pks;
            all_rhythms_M1 = all_rhythms_M1(all_rhythms_M1 > (0+pre_margin));
            all_rhythms_M1 = all_rhythms_M1(all_rhythms_M1 < ((length(M1_filt_data.LFP)/M1_filt_data.Fs_LFP)-post_margin));
            all_rhythms_Cb = Cb_filt_data.spindles{1}.pks;
            all_rhythms_Cb = all_rhythms_Cb(all_rhythms_Cb > (0+pre_margin));
            all_rhythms_Cb = all_rhythms_Cb(all_rhythms_Cb < ((length(Cb_filt_data.LFP)/Cb_filt_data.Fs_LFP)-post_margin));
            M1_spike_bools = nan(0, round((pre_margin + post_margin)/bin_width));
            Cb_spike_bools = nan(0, round((pre_margin + post_margin)/bin_width));
            
            M1_peak_spike_phases = zeros(1,0);
            M1_tail_spike_phases = zeros(1,0);
            Cb_peak_spike_phases = zeros(1,0);
            Cb_tail_spike_phases = zeros(1,0);
            for rhythm = 1:length(all_rhythms_M1)
                center = round(M1_filt_data.Fs_LFP * all_rhythms_M1(rhythm));
                
                if isempty(M1_neuron)
                    all_spike_times = [];
                else
                    all_spike_times = M1_spike_timestamps{M1_neuron};
                end
                spike_times = all_spike_times(all_spike_times >= ((center/param.M1_Fs) - pre_margin) & all_spike_times <= ((center/param.M1_Fs) + post_margin));
                spike_times = spike_times - (center/param.M1_Fs);
                cur_hist = histcounts(spike_times, -pre_margin:bin_width:post_margin);
                M1_spike_bools = cat(1, M1_spike_bools, cur_hist);
                
                peak_idxs_pre = all_spike_times >= (peak_idxs_M1(spindle_peak_idxs_M1(rhythm) - peak_range(2))/param.M1_Fs) & all_spike_times <= (peak_idxs_M1(spindle_peak_idxs_M1(rhythm) - peak_range(1))/param.M1_Fs);
                peak_idxs_post = all_spike_times >= (peak_idxs_M1(spindle_peak_idxs_M1(rhythm) + peak_range(1))/param.M1_Fs) & all_spike_times <= (peak_idxs_M1(spindle_peak_idxs_M1(rhythm) + peak_range(2))/param.M1_Fs);
                peak_spike_times = all_spike_times(peak_idxs_pre | peak_idxs_post);
                peak_spike_phases = round(peak_spike_times * param.M1_Fs);
                peak_spike_phases = radial_phase_M1(peak_spike_phases);
                M1_peak_spike_phases = [M1_peak_spike_phases peak_spike_phases];
                
                tail_idxs_pre = all_spike_times >= (peak_idxs_M1(spindle_peak_idxs_M1(rhythm) - tail_range(2))/param.M1_Fs) & all_spike_times <= (peak_idxs_M1(spindle_peak_idxs_M1(rhythm) - tail_range(1))/param.M1_Fs);
                tail_idxs_post = all_spike_times >= (peak_idxs_M1(spindle_peak_idxs_M1(rhythm) + tail_range(1))/param.M1_Fs) & all_spike_times <= (peak_idxs_M1(spindle_peak_idxs_M1(rhythm) + tail_range(2))/param.M1_Fs);
                tail_spike_times = all_spike_times(tail_idxs_pre | tail_idxs_post);
                tail_spike_phases = round(tail_spike_times * param.M1_Fs);
                tail_spike_phases = radial_phase_M1(tail_spike_phases);
                M1_tail_spike_phases = [M1_tail_spike_phases tail_spike_phases];
            end
            
            for rhythm = 1:length(all_rhythms_Cb)
                center = round(Cb_filt_data.Fs_LFP * all_rhythms_Cb(rhythm));
                
                if isempty(Cb_neuron)
                    all_spike_times = [];
                else
                    all_spike_times = Cb_spike_timestamps{Cb_neuron};
                end
                spike_times = all_spike_times(all_spike_times >= ((center/param.Cb_Fs) - pre_margin) & all_spike_times <= ((center/param.Cb_Fs) + post_margin));
                spike_times = spike_times - (center/param.Cb_Fs);
                cur_hist = histcounts(spike_times, -pre_margin:bin_width:post_margin);
                Cb_spike_bools = cat(1, Cb_spike_bools, cur_hist);
                
                peak_idxs_pre = all_spike_times >= (peak_idxs_Cb(spindle_peak_idxs_Cb(rhythm) - peak_range(2))/param.Cb_Fs) & all_spike_times <= (peak_idxs_Cb(spindle_peak_idxs_Cb(rhythm) - peak_range(1))/param.Cb_Fs);
                peak_idxs_post = all_spike_times >= (peak_idxs_Cb(spindle_peak_idxs_Cb(rhythm) + peak_range(1))/param.Cb_Fs) & all_spike_times <= (peak_idxs_Cb(spindle_peak_idxs_Cb(rhythm) + peak_range(2))/param.Cb_Fs);
                peak_spike_times = all_spike_times(peak_idxs_pre | peak_idxs_post);
                peak_spike_phases = round(peak_spike_times * param.Cb_Fs);
                peak_spike_phases = radial_phase_Cb(peak_spike_phases);
                Cb_peak_spike_phases = [Cb_peak_spike_phases peak_spike_phases];
                
                tail_idxs_pre = all_spike_times >= (peak_idxs_Cb(spindle_peak_idxs_Cb(rhythm) - tail_range(2))/param.Cb_Fs) & all_spike_times <= (peak_idxs_Cb(spindle_peak_idxs_Cb(rhythm) - tail_range(1))/param.Cb_Fs);
                tail_idxs_post = all_spike_times >= (peak_idxs_Cb(spindle_peak_idxs_Cb(rhythm) + tail_range(1))/param.Cb_Fs) & all_spike_times <= (peak_idxs_Cb(spindle_peak_idxs_Cb(rhythm) + tail_range(2))/param.Cb_Fs);
                tail_spike_times = all_spike_times(tail_idxs_pre | tail_idxs_post);
                tail_spike_phases = round(tail_spike_times * param.Cb_Fs);
                tail_spike_phases = radial_phase_Cb(tail_spike_phases);
                Cb_tail_spike_phases = [Cb_tail_spike_phases tail_spike_phases];
            end
            
            %plot([raster_x; raster_x],[0; 1]+raster_y,'LineWidth',1,'Color',rgb('Black'));
            hold on
            for x = 1:size(M1_spike_bools,2)
                for y = 1:min([size(M1_spike_bools,1), max_to_plot])
                    if M1_spike_bools(y,x)
                        x_ms = ((x-0.5)*bin_width)-pre_margin;
                        plot([x_ms;x_ms], [y;y+1], 'LineWidth', 1, 'Color', [0,0,0]);
                    end
                end
            end
            if cross_area_bool
                title('CROSS!');
            end
            saveas(gcf, [rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Single_Neuron_Spindle_Rastor_M1.fig'])
            close all
                        
            t = circ_mean(M1_peak_spike_phases');
            r = circ_r(M1_peak_spike_phases');
            
            polarhistogram(M1_peak_spike_phases,20,'FaceColor', [0.4 0.4 0.8], 'FaceAlpha', 1, 'Normalization', 'pdf');
            hold on
            polarplot([t t], [0 (r)]) %*length(M1_peak_spike_phases)
            if cross_area_bool
                title('CROSS!');
            end
            close all
            
            t = circ_mean(M1_tail_spike_phases');
            r = circ_r(M1_tail_spike_phases');
            polarhistogram(M1_tail_spike_phases,20,'FaceColor', [0.3 0.7 1.0], 'FaceAlpha', 1, 'Normalization', 'pdf');
            hold on
            polarplot([t t], [0 (r)])
            if cross_area_bool
                title('CROSS!');
            end
            close all
            
            hold on
            for x = 1:size(Cb_spike_bools,2)
                for y = 1:min([size(Cb_spike_bools,1), max_to_plot])
                    if Cb_spike_bools(y,x)
                        x_ms = ((x-0.5)*bin_width)-pre_margin;
                        plot([x_ms;x_ms], [y;y+1], 'LineWidth', 1, 'Color', [0,0,0]);
                    end
                end
            end
            if cross_area_bool
                title('CROSS!');
            end
            saveas(gcf, [rootpath,animal,'/Day',num2str(day),'/',param.s_block_names{block},'/Single_Neuron_Spindle_Rastor_Cb.fig'])
            close all
            
            t = circ_mean(Cb_peak_spike_phases');
            r = circ_r(Cb_peak_spike_phases');
            polarhistogram(Cb_peak_spike_phases,20,'FaceColor', [0.4 0.4 0.8], 'FaceAlpha', 1, 'Normalization', 'pdf');
            hold on
            polarplot([t t], [0 (r)])
            if cross_area_bool
                title('CROSS!');
            end
            close all
            
            t = circ_mean(Cb_tail_spike_phases');
            r = circ_r(Cb_tail_spike_phases');
            polarhistogram(Cb_tail_spike_phases,20,'FaceColor', [0.3 0.7 1.0], 'FaceAlpha', 1, 'Normalization', 'pdf');
            hold on
            polarplot([t t], [0 (r)])
            if cross_area_bool
                title('CROSS!');
            end
            close all
        end
    end
    rmpath('Z:\Matlab for analysis\CircStat2012a');
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Cross area Inter-Spindle Intervals (41)

if enabled(41)
    disp('Block 41...')
    load([rootpath,animal,'/Shared_Data.mat'])
    shared_data.M12Cb_intervals = cell(param.days, param.sleep_blocks);
    shared_data.Cb2M1_intervals = cell(param.days, param.sleep_blocks);
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_data.mat'])
            M1_spin = M1_data.spindles{1}.pks;
            Cb_spin = Cb_data.spindles{1}.pks;
            M12Cb_intervals = nan(2,length(M1_spin)); %Row 1 is backward (negative intervals), 2 is forward
            Cb2M1_intervals = nan(2,length(Cb_spin));
            M1_idx = 1;
            Cb_idx = 1;
            if isempty(M1_spin) || isempty(Cb_spin)
                continue
            end
            if M1_spin(M1_idx) < Cb_spin(Cb_idx)
                while M1_spin(M1_idx) < Cb_spin(Cb_idx)
                    M1_idx = M1_idx +1;
                end
            else
                while Cb_spin(Cb_idx) < M1_spin(M1_idx)
                    Cb_idx = Cb_idx +1;
                end
            end
            
            %Spindles at the ends that are not bracketed are left as NaNs. NaNs are not removed to retain pairity with *data.spindles{1}.pks
            while M1_idx <= length(M1_spin) && Cb_idx <= length(Cb_spin)
                if M1_spin(M1_idx) == Cb_spin(Cb_idx)
                    %Set forward and back dist to 0 for both
                    M12Cb_intervals(:,M1_idx) = [0;0];
                    Cb2M1_intervals(:,Cb_idx) = [0;0];
                    %increment both idxs
                    M1_idx = M1_idx+1;
                    Cb_idx = Cb_idx+1;
                elseif M1_spin(M1_idx) < Cb_spin(Cb_idx)
                    M12Cb_intervals(1,M1_idx) = Cb_spin(Cb_idx-1) - M1_spin(M1_idx);
                    M12Cb_intervals(2,M1_idx) = Cb_spin(Cb_idx) - M1_spin(M1_idx);
                    M1_idx = M1_idx+1;
                else
                    Cb2M1_intervals(1,Cb_idx) = M1_spin(M1_idx-1) - Cb_spin(Cb_idx);
                    Cb2M1_intervals(2,Cb_idx) = M1_spin(M1_idx) - Cb_spin(Cb_idx);
                    Cb_idx = Cb_idx+1;
                end
            end
            shared_data.M12Cb_intervals{day,block} = M12Cb_intervals;
            shared_data.Cb2M1_intervals{day,block} = Cb2M1_intervals;
        end
    end
    save([rootpath,animal,'/Shared_Data.mat'], 'shared_data')
    clearvars -except rootpath origin_rootpath animal param enabled;
end

%% Plot sample spindle and corisponding control (42)

if enabled(42)
    disp('Block 42...')
    area = 'Cb';
    day = 1;
    block = 1;
    spindle_idx = 5;
    control_offset = -5; %in seconds
    
    load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Filtered_Sleep_data.mat'])
    
    if strcmp(area, 'M1')
        filt_data = M1_filt_data;
    elseif strcmp(area, 'Cb')
        filt_data = Cb_filt_data;
    else
        error('Unrecognized area.')
    end
    win_start_offset = -0.5 * filt_data.Fs_LFP;
    win_end_offset = 0.5 * filt_data.Fs_LFP;
    
    spin_timestamp = filt_data.spindles{1}.pks(spindle_idx) * filt_data.Fs_LFP;
    LFP_data = filt_data.LFP(round(spin_timestamp + win_start_offset) : round(spin_timestamp + win_end_offset));
    time_data = 0:(length(LFP_data)-1);
    time_data = ((time_data * 1000)/length(LFP_data)) - 500;
    plot(time_data, LFP_data)
    close all
    
    spin_timestamp = (filt_data.spindles{1}.pks(spindle_idx) + control_offset) * filt_data.Fs_LFP;
    LFP_data = filt_data.LFP(round(spin_timestamp + win_start_offset) : round(spin_timestamp + win_end_offset));
    time_data = 0:(length(LFP_data)-1);
    time_data = ((time_data * 1000)/length(LFP_data)) - 500;
    plot(time_data, LFP_data)    
    close all
end

%% Sum non-REM time and spindle count (43)

if enabled(43)
    disp('Block 43...')
    sleep_time = 0; %in sec
    M1_spindles = 0;
    Cb_spindles = 0;
    for day = 1:param.days
        for block = 1:param.sleep_blocks
            load([rootpath,animal,'/Day',num2str((day)),'/',param.s_block_names{block},'/Sleep_data.mat'])
            sleep = sum(M1_data.sleep_idx);
            sleep = sleep/M1_data.Fs_LFP;
            sleep_time = sleep_time + sleep;
            M1_spindles = M1_spindles + length(M1_data.spindles{1}.pks);
            Cb_spindles = Cb_spindles + length(Cb_data.spindles{1}.pks);
            disp(['Block Sleep dur: ', num2str(sleep)])
        end
    end
    disp(['Sleep dur: ', num2str(sleep_time), '.'])
    disp(['M1 spindles: ', num2str(M1_spindles), '.'])
    disp(['Cb spindles: ', num2str(Cb_spindles), '.'])
end

beep
disp 'Analysis Complete.'