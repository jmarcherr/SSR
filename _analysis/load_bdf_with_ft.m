function [ bdf_data, bdf_info ] = load_bdf_with_ft( path_bdf_file,...
    epoch_time )
%load_bdf_with_ft - Read a single .bdf (Biosemi) file using FiedTrip.
%Perform a simple pre-processing to the imported data which consist on:
% - Filter all the trials with a BandPass 4th order Butterworth zero-phase
%   filter from 60 to 300 Hz
% - Re-referencing all the recorded channels to the vertex 'Cz'
% - Remove the 'Cz' channel (self-referred) and a dummy channel named
%   'Status' that Biosemi always include.
%To finalize, the function converts the FieldTrip data structure to the
%data structure used in the rest of custom-made functions.
%
%
% Syntax:  [ bdf_data, bdf_info ] = load_bdf_with_ft( path_bdf_file,
%                                   epoch_time )
%
% Inputs: path_bdf_file -> Path where the .bdf file is stored
%         epoch_time    -> Epoch duration in seconds (usually 1 sec)
%
% Outputs: bdf_data -> Imported and pre-processed eeg data
%          bdf_info -> Information of the bdf file and the processing done
%                      using FieldTrip
%
% Example:
%    Function load_bdf_with_ft.m is a stand alone function.
%    Call function: Get_ASSR_Compr_MultiFreq.m
%
% Other m-files required: FieldTrip toolbox,

% Subfunctions:
% MAT-files required: No
%
% See also: IO_function, Read_Bdf_IOfunct, Ftest_IO_function
%   Plot_IO_function, Plot_Lvl_growth

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
cfg_read.trialdef.eventvalue = 192;

%% Pre-processing and Re-referencing

cfg_read.trialfun            = 'ft_trialfun_general';

% Fitering options (Band-pass from 60 to 300 Hz)
cfg_read.bpfilter            = 'yes';       %Band-pass Butterworth
cfg_read.bpfiltord           = 2;           %2nd Order
cfg_read.bpfreq              = [60 300];    %Hz
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

%Find the target channels
chan_idx_ref            = find(strcmp(data_in_epochs.label, 'Cz'));
chan_idx_status         = find(strcmp(data_in_epochs.label, 'Status'));
chan_idx_p9             = find(strcmp(data_in_epochs.label, 'P9'));
chan_idx_p10            = find(strcmp(data_in_epochs.label, 'P10'));
chan_idx_exg1           = find(strcmp(data_in_epochs.label, 'EXG1'));
chan_idx_exg2           = find(strcmp(data_in_epochs.label, 'EXG2'));

%Channels to remove from the general list
chan_unused_rmv             = [chan_idx_ref, chan_idx_status];
chan_no_apply_thres_artif   = [chan_idx_ref, chan_idx_status, chan_idx_p9,...
    chan_idx_p10, chan_idx_exg1, chan_idx_exg2];

%Channel list with the unused channels
chan_unused_list                    =  data_in_epochs.label;
chan_unused_list(chan_unused_rmv)   = [];

%Channel list without threshold artifact channels
chan_list_thres_artif                               = data_in_epochs.label;
chan_list_thres_artif(chan_no_apply_thres_artif)    = [];

%% Identify bad channels

%Compute absolute value to all epochs
data_in_epochs_abs  = cellfun(@abs, data_in_epochs.trial,...
    'UniformOutput', false);
data_in_epochs_abs_max  = cellfun(@(in) max(in, [], 2), ...
    data_in_epochs_abs, 'UniformOutput', false);
data_in_epochs_abs_max_mean  = mean(cell2mat(data_in_epochs_abs_max), 2);

%Find bad channels
thres_bad_channel   = 100;    %100 muV
idx_bad_channel     = data_in_epochs_abs_max_mean > thres_bad_channel;
name_bad_channel    = data_in_epochs.label(idx_bad_channel);

%Verify that bad channels are also in the threshold artifact channel list
idx_bad_channel_in_list_thres_artif = ismember(chan_list_thres_artif,...
    name_bad_channel);
name_bad_channel_verif = chan_list_thres_artif(idx_bad_channel_in_list_thres_artif);

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

%Prepare electrode layout
cfg_neig = [];
cfg_neig.method = 'template';
cfg_neig.layout = 'biosemi64.lay';

%Select channel repair method
cfg_repair              = [];
cfg_repair.method       = 'nearest';
cfg_repair.trials       = 'all';
cfg_repair.neighbours   = ft_prepare_neighbours(cfg_neig);
cfg_repair.elec         = ft_read_sens('standard_1020.elc', 'senstype', 'eeg');

cfg_repair.badchannel   = name_bad_channel_verif;

data_in_epochs_clean    = ft_channelrepair(cfg_repair, data_in_epochs_resample);

%% Remove non-used channels (Cz, Status)

cfg_rmv_chan                    = [];
cfg_rmv_chan.channel            = chan_unused_list;
data_in_epochs_preproc_resample = ft_preprocessing(cfg_rmv_chan,...
    data_in_epochs_clean);

%% Re-organize data to custom structure (epochs x channels x epochs)

%Pre-allocate matrix
samples_epoch           = size(data_in_epochs_preproc_resample.trial{1,1}, 2);
num_chan                = size(data_in_epochs_preproc_resample.trial{1,1}, 1);
num_epochs              = size(data_in_epochs_preproc_resample.trialinfo, 1);
bdf_data                = zeros(samples_epoch, num_chan, num_epochs);

for idx = 1:num_epochs
    bdf_data(:, :, idx)   = data_in_epochs_preproc_resample.trial{1, idx}';
end

% Remove eeg data from data_in_epochs to save the rest as bdf_info
bdf_info = data_in_epochs_preproc_resample;
bdf_info = rmfield(bdf_info, 'trial');


end

%Eof