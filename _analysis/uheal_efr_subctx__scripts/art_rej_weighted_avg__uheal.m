function [data_weighavg, data_weighavg_mtx, data_weighavg_weights] = ...
    art_rej_weighted_avg__uheal(data_2process, fs, epochs_x_trial)

%art_rej_weighted_avg - Artifact rejection for ASSR data based on
%weighted averaging.
%The function averages the input data weighting each epoch by the inverse
%of its variance. It is assumed that the more noisy an epoch is the more
%intrinsict variance have. Thus, a smaller weight is given to the epochs
%with large variance and the less variable epochs (less noisy) will get
%bigger weight. Prior to computing the epoch variance, it is filtered using
%a zero-phase second order bandpass butterworth filter from 70 to 110 Hz.
%Prior to that, an absolute threshold is applied to all epochs so that any
%epoch with absolute amplitude over 40 muV is discarded.
%Trials are formed by concatenating 16 epochs together. Trials are weighted
%averaged which results in the final weighted averaged trial.
%
%
% Syntax:  [data_weighavg, data_weighavg_mtx, data_weighavg_weights] = ...
%    art_rej_weighted_avg(data_2process, fs, fm1, fm2, fm3, fm4 ,...
%    epochs_x_trial)
%
% Inputs:   data_2process   --> Input data
%           fs              --> Sample frequency
%           epochs_x_trial  --> Num. of epochs to concatenate (Default: 16)
%
% Outputs:  data_weighavg           --> Weighted averaged output trial
%           data_weighavg_mtx       --> Weighted data prior to average
%           data_weighavg_weights   --> Weights
%
% Example:
%    Function art_rej_weighted_avg.m is usually called by another function
%    from the "Get_data_func" group
%
%
% See also: Get_ASSR_Compr_MultiFreq Get_ASSR_ModDepth

%% Initial variables and Filter design

% Pre-allocate variables space
num_epochs              = size(data_2process, 2);
sample_x_epoch          = size(data_2process, 1);
% Matrix of weighted epochs
epochs_weighted_mtx     = zeros(size(data_2process));
% Epoch variance vector
epoch_var               = zeros(num_epochs, 1);
% Count the number of rejected epochs
rejt_epochs             = 0;
% Index to compensate for the rejected epochs
i_rejt_comp             = 0;

% Artifact rejection Filter (2nd Butterworth filter from 170-210 Hz)
hd_art_rej_filt = designfilt('bandpassiir','FilterOrder',2, ...
    'HalfPowerFrequency1',170,'HalfPowerFrequency2',210, ...
    'SampleRate',fs,'designmethod','butter');

%% Absolute amplitude thresholding (+-40 muV) and weights computing

for idx_var = 1:num_epochs
    epoch_currt_anlys = data_2process(:,idx_var);
    
    if max(abs(epoch_currt_anlys)) > 80  %Absolute threshold above +-80 muV
        rejt_epochs = rejt_epochs+1;
        i_rejt_comp = i_rejt_comp+1;
    else
        %Filter the current epoch from 70-110 Hz
        epoch_currt_anlys_filt = filtfilt(hd_art_rej_filt,...
            epoch_currt_anlys);
        
        %Compute filtered variance of the current epoch
        epoch_var(idx_var-i_rejt_comp) = var(epoch_currt_anlys_filt);
        
        %Weight the current epoch by its inverse variance (1/var)
        epochs_weighted_mtx(:,idx_var-i_rejt_comp) = ...
            epoch_currt_anlys./epoch_var(idx_var-i_rejt_comp);
    end
end

%Cut the generated matrices to remove space of rejected epochs
epoch_var = epoch_var(1:end-rejt_epochs);
epoch_var = epoch_var.^(-1);
epochs_weighted_mtx = epochs_weighted_mtx(:, 1:end-rejt_epochs);

%Number of valid epochs
valid_epochs = num_epochs - rejt_epochs;
fprintf('Number of valid epochs: %d out of %d \n',valid_epochs, num_epochs)
%Rejected epochs ratio (in percentage)
rejt_ratio = rejt_epochs*100/num_epochs;
fprintf('Percentage of rejected epochs %.2f %% \n \n', rejt_ratio)

if rejt_ratio > 25
    warning('Epochs rejected ratio over 25% !!!')
end

%% Concatenate epochs to form Trials (Default: 16 epochs)

% Create a Trial linking epochs_x_trial (default: 16) epochs
num_full_trials = floor(valid_epochs/epochs_x_trial);

data_weighavg_mtx = reshape(epochs_weighted_mtx(:,...
    1:epochs_x_trial*num_full_trials), sample_x_epoch*epochs_x_trial,...
    num_full_trials);

% Vector of the summed weights for each trial
data_weighavg_weights= reshape(epoch_var(1:num_full_trials*...
    epochs_x_trial, 1), epochs_x_trial, num_full_trials);

summed_weigths = sum(data_weighavg_weights, 2);

% Summed trials
summed_trial = sum(data_weighavg_mtx, 2);

%Compensate each epoch of the final summed trial by the sum of the weigths
%of all the summed epochs
idx_sum_weight = 1;
for idx_sum = 1:sample_x_epoch:length(summed_trial)
    data_weighavg(idx_sum:idx_sum + sample_x_epoch-1, 1) = ...
        summed_trial(idx_sum:idx_sum+sample_x_epoch-1)./...
        summed_weigths(idx_sum_weight);
    
    idx_sum_weight = idx_sum_weight + 1;
end

end

%Eof