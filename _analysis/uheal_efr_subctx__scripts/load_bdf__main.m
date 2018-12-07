function load_bdf__main

%path_bdf_file = ['/Volumes/EXT_DISK/encina/1_recordings' filesep...
%    '6_deaff_mod_depth_4m_pilot' filesep 'ando90' filesep 'ando90_raw_data'...
%    filesep 'subj-ando90_exp-deaff_mod_depth_4m_pilot_ear-Right_lvl-81dB_mod-85.bdf'];
%path_bdf_file   = ['/Users/Gerard/Desktop' filesep 'uheal_tmp' filesep 'JM_DG.bdf'];
path_bdf_file = '/Users/jmarcher/Documents/git/SSR/_data/raw_bdf/pilot3_2xsamtone+oldx10min.bdf';
epoch_time      = 3;

%% Read data using FieldTrip

fprintf('#################################### \n');
fprintf('# Reading bdf file using FIELDTRIP # \n');
fprintf('#################################### \n\n');

% Prepare config structure pointing to the bdf file
cfg_read                     = [];
cfg_read.dataset             = path_bdf_file;

% Define trial length (in second) and trigger
cfg_read.trialdef.prestim    = 0; % 0 s prestimulus window
cfg_read.trialdef.poststim   = epoch_time; % 1 s poststimulus window
cfg_read.trialdef.eventtype  = 'STATUS';
cfg_read.trialdef.eventvalue = 255;

%% Pre-processing and Re-referencing

cfg_read.trialfun            = 'ft_trialfun_general';

% Fitering options (Band-pass from 0.5 to 400 Hz)
cfg_read.bpfilter            = 'yes';       %Band-pass Butterworth
cfg_read.bpfiltord           = 2;           %2nd Order
cfg_read.bpfreq              = [1 400];    %Hz
cfg_read.bpfiltdir           = 'twopass';
cfg_read.bpfilttype          = 'but';

% Re-refererring all channels to Cz
cfg_read.reref               = 'yes';
cfg_read.channel             = 'all';
cfg_read.refchannel          = 'Cz';     % Cz is used as the new reference

%Define trials
cfg_in_trials               = ft_definetrial(cfg_read);

%Pre-process data
data_in_epochs              = ft_preprocessing(cfg_in_trials);

 %% Remove unused channels
% 
% %Find the target channels
% chan_idx_ref        = find(strcmp(data_in_epochs.label, 'Cz'));
% chan_idx_status     = find(strcmp(data_in_epochs.label, 'Status'));
% chan_idx_fc5        = find(strcmp(data_in_epochs.label, 'FC5'));
% chan_idx_fc6        = find(strcmp(data_in_epochs.label, 'FC6'));
% chan_idx_t7         = find(strcmp(data_in_epochs.label, 'T7'));
% chan_idx_t8         = find(strcmp(data_in_epochs.label, 'T8'));
% chan_idx_p7         = find(strcmp(data_in_epochs.label, 'P7'));
% chan_idx_p8         = find(strcmp(data_in_epochs.label, 'P8'));
% 
% %Channels to remove from the general list
% chan_unused_rmv             = [chan_idx_ref, chan_idx_status];
% chan_no_apply_thres_artif   = [chan_idx_ref, chan_idx_status, chan_idx_fc5,...
%     chan_idx_fc6, chan_idx_t7, chan_idx_t8, chan_idx_p7, chan_idx_p8];
% 
% %Channel list with the unused channels
% chan_unused_list                    =  data_in_epochs.label;
% chan_unused_list(chan_unused_rmv)   = [];
% 
% %Channel list without threshold artifact channels
% chan_list_thres_artif                               = data_in_epochs.label;
% chan_list_thres_artif(chan_no_apply_thres_artif)    = [];

%% Identify bad channels

% %Compute absolute value to all epochs
% data_in_epochs_abs  = cellfun(@abs, data_in_epochs.trial,...
%     'UniformOutput', false);
% data_in_epochs_abs_max  = cellfun(@(in) max(in, [], 2), ...
%     data_in_epochs_abs, 'UniformOutput', false);
% data_in_epochs_abs_max_mean  = mean(cell2mat(data_in_epochs_abs_max), 2);
% 
% %Find bad channels
% thres_bad_channel   = 80;    %80 muV
% idx_bad_channel     = data_in_epochs_abs_max_mean > thres_bad_channel;
% name_bad_channel    = data_in_epochs.label(idx_bad_channel);
% 
% %Verify that bad channels are also in the threshold artifact channel list
% idx_bad_channel_in_list_thres_artif = ismember(chan_list_thres_artif,...
%     name_bad_channel);
% name_bad_channel_verif = chan_list_thres_artif(idx_bad_channel_in_list_thres_artif);

%% Re-sample the data to a fixed sample rate of 2048 Hz

fs_raw_data_homogen     = 2048;

if data_in_epochs.fsample > fs_raw_data_homogen
    
    cfg_decimate                 = [];
    cfg_decimate.resamplefs = fs_raw_data_homogen;
    cfg_decimate.detrend    = 'no';
    cfg_decimate.demean     = 'no';
    cfg_decimate.feedback   = 'text';
    cfg_decimate.trials     = 'all';
    
    data_in_epochs_resample = ft_resampledata(cfg_decimate,...
        data_in_epochs);
    
elseif data_in_epochs.fsample == fs_raw_data_homogen
    
    data_in_epochs_resample = data_in_epochs;
    
else
    
    data_in_epochs_resample = data_in_epochs;
    fprintf('\n############################ WARNING ############################\n')
    fprintf('The sample frequency in the original bdf file is %i Hz,\nwhich is lower than the recommended sample frequency %i Hz!\n\n',...
        data_in_epochs.fsample, fs_raw_data_homogen)
end

%% Remove bad channels

% %Prepare electrode layout
% cfg_neig = [];
% cfg_neig.method = 'template';
% cfg_neig.layout = 'biosemi32.lay';
% 
% %Select channel repair method
% cfg_repair              = [];
% cfg_repair.method       = 'nearest';
% cfg_repair.trials       = 'all';
% cfg_repair.neighbours   = ft_prepare_neighbours(cfg_neig);
% cfg_repair.elec         = ft_read_sens('standard_1020.elc', 'senstype', 'eeg');
% 
% cfg_repair.badchannel   = name_bad_channel_verif;
% 
% data_in_epochs_clean    = ft_channelrepair(cfg_repair, data_in_epochs_resample);

%% Remove non-used channels (Cz, Status)

% cfg_rmv_chan                    = [];
% cfg_rmv_chan.channel            = chan_unused_list;
% data_in_epochs_preproc_resample = ft_preprocessing(cfg_rmv_chan,...
%     data_in_epochs_clean);

%% Re-organize data to custom structure (epochs x channels x epochs)

%Pre-allocate matrix
samples_epoch           = size(data_in_epochs_resample.trial{1,1}, 2);
num_chan                = size(data_in_epochs_resample.trial{1,1}, 1);
num_epochs              = size(data_in_epochs_resample.trialinfo, 1);
bdf_data                = zeros(samples_epoch, num_chan, num_epochs);

for idx = 1:num_epochs
    bdf_data(:, :, idx)   = data_in_epochs_resample.trial{1, idx}';
end

% Remove eeg data from data_in_epochs to save the rest as bdf_info
bdf_info = data_in_epochs_resample;
bdf_info = rmfield(bdf_info, 'trial');

%% Select one channel or the mean of several channels

chan_idx = strcmp(bdf_info.label, 'EXG1');
%chan_idx = logical(sum([strcmp(bdf_info.label, 'P10'), strcmp(bdf_info.label, 'EXG2')], 2));


data_in__weighavg = squeeze(mean(bdf_data(:, chan_idx, :), 2));

%% Weighted averaging

fs = bdf_info.fsample;
epochs_x_trial = 10;

%Weighted Averaging
[data__weighavg, data__weighavg_mtx, data__weighavg_weights] = ...
    art_rej_weighted_avg__uheal( data_in__weighavg, fs, epochs_x_trial); %#ok<ASGLU>

%% Plot

aux_time_spectrum_viewer( data__weighavg, fs, 'sig_lgth', 0, 1000, -80, 0, false, 'JM, fm=195, m=85%' )

end

%Eof