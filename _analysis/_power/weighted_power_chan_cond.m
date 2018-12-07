clear all
assr_startup
ft_defaults

%% Select preprocessed .mat file to process
cd(datadir)
load(['pilot3_2xsamtone+oldx10min.mat']);
cd(rootdir)
close all
for kk=1:3 % trial loop
    chan_idx = 0;
    for chan=[find(strcmp(data_DG.label,'EXG1')) find(strcmp(data_DG.label,'EXG2'))];
        %% rename data
        data_all= data_DG;fs = data_all.fsample;
        cfg = [];
        cfg.channel = chan;%{'EXG1'};%{'all','-Status','-Cz'} % select channel
        data_all = ft_selectdata(cfg,data_all)
        chan_idx = chan_idx +1;
        
        %% divide into trials
        trial_id = kk;
        epoched_data = epoch_data(data_all,trial_id,100);
        
        %% artifact rejection and weighting
        reject_epoch = 0
        reject_idx = 0;
        filt_coef = [180 220];
        filt_def = designfilt('bandpassiir','FilterOrder',2, ...
            'HalfPowerFrequency1',filt_coef(1),'HalfPowerFrequency2',filt_coef(2), ...
            'SampleRate',fs,'designmethod','butter');
        
        for ii=1:size(epoched_data,1)
            data_tmp = epoched_data(ii,:);
            reject_idx = reject_idx+1;
            if max(abs(data_tmp)) > 80 % thresholding 80mV
                reject_epoch = reject_epoch +1;
                reject_idx = reject_idx-1;
            else
                %filter between filt_coeff
                filtered_data = filtfilt(filt_def, data_tmp);
                %compute variance of filtered epoch
                epoch_var(reject_idx) = var(filtered_data);
                %Weight the current epoch by its inverse variance (1/var)
                epoch_weighted(reject_idx,:) = data_tmp./epoch_var(reject_idx);
            end
        end
        
        epoch_var = epoch_var.^(-1); % weights
        
        %Reporting
        clc
        %Number of valid epochs
        valid_epochs = size(epoched_data,1) - reject_epoch;
        fprintf('Number of valid epochs: %d out of %d \n',valid_epochs, size(epoched_data,1))
        %Rejected epochs ratio (in percentage)
        rejt_ratio = reject_epoch*100/size(epoched_data,1);
        fprintf('Percentage of rejected epochs %.2f %% \n \n', rejt_ratio)
        
        if rejt_ratio > 25
            warning('Epochs rejected ratio over 25% !!!')
        end
        
        %%
        %Concatenating into epochs
        
        %number of epochs per trial
        epochs_x_trial = 5;
        %number of epochs
        idx_tmp = size(epoch_weighted,1);
        %number of valid trials
        nr_trials = floor(idx_tmp/epochs_x_trial);
        %length of each epoch in samples
        sample_x_epoch = length(epoch_weighted);
        
        % generate epochs
        cc=0;
        epoched_data = [];
        
        trial_data = reshape(epoch_weighted(1:nr_trials*epochs_x_trial,:)',sample_x_epoch*epochs_x_trial, nr_trials);
        trial_weights = reshape(epoch_var(1:nr_trials*epochs_x_trial),epochs_x_trial, nr_trials);
        
        % summing over trials
        summed_data = sum(trial_data,2);
        summed_weights = sum(trial_weights,2);
        
        
        %Compensate each epoch of the final summed trial by the sum of the weigths
        %of all the summed epochs
        idx_sum_weight = 1;
        for idx_sum = 1:sample_x_epoch:length(summed_data)
            data_weighavg(idx_sum:idx_sum + sample_x_epoch-1, 1) = ...
                summed_data(idx_sum:idx_sum+sample_x_epoch-1)./...
                summed_weights(idx_sum_weight);
            
            idx_sum_weight = idx_sum_weight + 1;
        end
        M=[];
        M = data_weighavg;
        
        %% %% perform FFT on epoched data
        f_fft = [];
        f_fft = fft(M)/(length(M)/2); % length of vector
        
        %f_fft_mean = squeeze(nanmean(f_fft,1));
        f_fft_pow = abs(f_fft.^2);
        f_fft_pow = (f_fft_pow(1:end/2+1));
        
        f = fs/2*linspace(0,1,length(f_fft_pow));
        %%
        figure(1)
        subplot(3,2,chan_idx+(kk-1)*2)
        plot(f,10*log10(f_fft_pow));
        title(data_DG.label(chan))
        xlim([163 257])
        ylim([-40 -20])
    end
end