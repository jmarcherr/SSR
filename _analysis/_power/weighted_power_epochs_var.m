clear all
assr_startup
ft_defaults

%% Select preprocessed .mat file to process
subs = ['pilot_17_1_19_SAM_subj_4.mat';'pilot_18_1_19_SAM_subj_5.mat';
cd(datadir)
for ss=1:3
load(['pilot_17_1_19_SAM_subj_4.mat']);

cd(rootdir)
close all
ff=0;
fms=[4 207];
% trial loop
chan_idx = 0;
chan=[find(strcmp(data_DG.label,'Cz'))]% find(strcmp(data_DG.label,'EXG2'))];
%% rename data
data_all= data_DG;fs = data_all.fsample;
cfg = [];
cfg.channel = [chan];%{'EXG1'};%{'all','-Status','-Cz'} % select channel
data_all = ft_selectdata(cfg,data_all);
chan_idx = chan_idx +1;

%% divide into trials
for ff=1:2
    fm_freq_oi = fms(ff);
    %init
    epoched_data = [];epoch_var=[];epoch_weighted=[];trial_data = [];trial_weights=[];data_tmp=[];
    % trigger
    for kk=1:4
        trial_id = kk;
        %init
        epoched_data = [];epoch_var=[];epoch_weighted=[];trial_data = [];
        trial_weights=[];data_tmp=[];
        %steps
        e_step = 10;e_idx = 0;
        for ee = e_step:e_step:length(find(data_all.trialinfo==trial_id))
            e_idx = e_idx+1;
            trials_oi = ee-e_step+1:ee;
            
            epoched_data = epoch_data_2(data_all,trial_id,trials_oi);
            
            %% artifact rejection and weighting
            reject_epoch = 0;
            reject_idx = 0;
            if fm_freq_oi ==4
            filt_coef = [fm_freq_oi-fm_freq_oi/5 fm_freq_oi+fm_freq_oi/5]; %4 Hz
            else 
            filt_coef = [fm_freq_oi-fm_freq_oi/15 fm_freq_oi+fm_freq_oi/15]; % 207 Hz
            end
            
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
            
            
            %%Concatenating into epochs
            
            %number of epochs per trial
            if fm_freq_oi < 207
                epochs_x_trial = 5;
            else
                epochs_x_trial = 5;
            end
            %number of epochs
            idx_tmp = size(epoch_weighted,1);
            %number of valid trials
            nr_trials = floor(idx_tmp/epochs_x_trial);
            %length of each epoch in samples
            sample_x_epoch = length(epoch_weighted);
            
            % generate epochs
            cc=0;
            epoched_data = [];
            
            trial_data{e_idx} = reshape(epoch_weighted(1:nr_trials*epochs_x_trial,:)',sample_x_epoch*epochs_x_trial, nr_trials);
            trial_weights{e_idx} = reshape(epoch_var(1:nr_trials*epochs_x_trial),epochs_x_trial, nr_trials);
            
        end
        
        trial_data_cat{kk} = cell2mat(trial_data);
        trial_weights_cat{kk} = cell2mat(trial_weights);
        
    end
    
    %% Summing over trials
    
    for mm=1:4
         % common number of valid 15s blocks
        m(mm) = length(trial_data_cat{mm}(1,:));
    end
    for cond = 1:4 %looping for test-retest

        N=min(m)
        % random selection
        day_idx = randperm(N);
        
        D_data = trial_data_cat{cond};
        D_weights = trial_weights_cat{cond};
        
        for n=1:N
            summed_data = sum(D_data(:,day_idx(1:n)),2); % picking random segments
            summed_weights = sum(D_weights(:,day_idx(1:n)),2); % finding weights
            
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
            %FFT
            f_fft = fft(M)/(length(M)/2);
            %Convert to power
            f_fft_pow = abs(f_fft.^2); %
            %Truncate negative freqencies
            f_fft_pow = (f_fft_pow(1:end/2+1));
            %Frequency vector
            f = fs/2*linspace(0,1,length(f_fft_pow));
            
            %% finding peak
            %fc's
            fc_idx = [fm_freq_oi];
            fc = find(f==fc_idx);
            
            bg_freq = [];
            % background (fc +/-3 Hz)
            if fm_freq_oi == 4
            bg_freq = [find(f>=fc_idx-.5 & f<fc_idx), ...
                find(f<=fc_idx+.5 & f>fc_idx)];
            else
            bg_freq = [find(f>=fc_idx-3 & f<fc_idx), ...
                find(f<=fc_idx+3 & f>fc_idx)];
            end
          
            
            %Power at fm_freq_oi
            Psn{ff,cond,n} = f_fft_pow(fc);
            %Averaged power in surrounding noise bands
            Pn{ff,cond,n} = mean(f_fft_pow(bg_freq));
            %SNR(dB)
            SNR{ff,cond,n} = 10*log10(f_fft_pow(fc))-10*log10(mean(f_fft_pow(bg_freq)));
            %F-statistic
            F{ff,cond,n}=f_fft_pow(fc)/mean(f_fft_pow(bg_freq));
            %Signal power
            Ps{ff,cond,n} = f_fft_pow(fc)-mean(f_fft_pow(bg_freq));
            F_crit{ff} = finv(0.95,2,2*length(bg_freq))
        end
        
        
    end
end

%%
close all
colorcodes = {'k','g','b','r'};
for f=1:2 %frequency
    figure(f)
for cc=1:4 %condition
    for ii=1:N %number of segments
        p(cc)=plot(ii,SNR{f,cc,ii},[colorcodes{cc}, 'o'])
        hold on
        %plot(ii,Pn{f,cc,ii},[colorcodes{cc},'x'])
    end
    
    plot(1:N,ones(1,N)*db(sqrt(F_crit{ff}-1)))
end

legend([p(1),p(2),p(3),p(4)],'25% mod, 75dB','85% mod, 40dB','25% mod, 40dB','85% mod, 75dB')

end