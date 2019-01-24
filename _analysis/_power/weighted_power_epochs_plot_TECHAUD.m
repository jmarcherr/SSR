clear all
assr_startup
ft_defaults

%% Select preprocessed .mat file to process
cd(datadir)
load(['pilot3_2xsamtone+oldx10min.mat']);
cd(rootdir)
close all
ff=0;
for fm_freq_oi=[4 207];
    ff=ff+1;
    for kk=1:3 % trial loop
        chan_idx = 0;
        for chan=[find(strcmp(data_DG.label,'EXG1'))]% find(strcmp(data_DG.label,'EXG2'))];
            %% rename data
            data_all= data_DG;fs = data_all.fsample;
            cfg = [];
            cfg.channel = [chan];%{'EXG1'};%{'all','-Status','-Cz'} % select channel
            data_all = ft_selectdata(cfg,data_all);
            chan_idx = chan_idx +1;
            
            %% divide into trials
            
            trial_id = kk;
            pow_207 = [];
            e_step = 10;e_idx = 0;
            for ee = e_step:e_step:length(find(data_all.trialinfo==kk))
                e_idx = e_idx+1;
                epoched_data = [];epoch_var=[];epoch_weighted=[];trial_data = [];trial_weights=[];data_tmp=[];
                epoched_data = epoch_data(data_all,trial_id,ee);
                
                %% artifact rejection and weighting
                reject_epoch = 0;
                reject_idx = 0;
                filt_coef = [fm_freq_oi-fm_freq_oi/3 fm_freq_oi+fm_freq_oi/3];
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
                bg_freq = [find(f>=fc_idx-1.5 & f<fc_idx), ...
                    find(f<=fc_idx+1.5 & f>fc_idx)];
                
                %Power at fm_freq_oi
                Psn{ff,kk,e_idx} = f_fft_pow(fc);
                %Averaged power in surrounding noise bands
                Pn{ff,kk,e_idx} = mean(f_fft_pow(bg_freq));
                %SNR(dB)
                SNR{ff,kk,e_idx} = 10*log10(f_fft_pow(fc))-10*log10(mean(f_fft_pow(bg_freq)));
                %F-statistic
                F{ff,kk,e_idx}=f_fft_pow(fc)/mean(f_fft_pow(bg_freq));
                F_crit(ff) = finv(0.99,2,2*length(bg_freq)); % critial value at a<0.01
            end
        end
            pow_tmp{kk}= f_fft_pow;
    end

end

%% SNR plot
close all
pow = SNR;%cellfun(@minus,Psn,Pn,'Un',0);
cp = cbrewer('qual','Set1',3);
sym = {'o','sq','v'};
%figure(99)
fid = [4 207];
for ff=1:2
    figure(ff)
    for k=1:3  
        if k==3
            index = 10;
        else
            index = 20;
        end
        for i=1:index
            Fdb = db(sqrt(F_crit(ff)-1));
            if squeeze(cell2mat(pow(ff,k,i)))>Fdb
            p(k)=plot(i,squeeze(cell2mat(pow(ff,k,i))),['k-',sym{k}],'MarkerFacecolor',cp(k,:),'markersize',14,'linewidth',.5);
            else
            p(k)=plot(i,squeeze(cell2mat(pow(ff,k,i))),['k-',sym{k}],'MarkerFacecolor',[1 1 1],'markersize',14,'linewidth',.5); 
            end
            hold on
            set(gca,'xtick',0:3:length(pow(ff,k,:)),'xticklabel', (0:e_step*3:length(pow(ff,k,:))*e_step)*3,'fontsize',20);
            if ff>1
                %xlabel('Stim time (s)')
            end
            ylabel('SNR (dB)')
            ylim([-10 15])
            xlim([0 21])
            
        end
    end
    
    title(['fm@ ', num2str(fid(ff)), ' Hz', ', fc@2005']);
    p(k+1)=plot(ones(size(squeeze(cell2mat(pow(ff,2,:)))))*db(sqrt(F_crit(ff)-1)),'-.','color',[0.5 0.5 0.5]);
    
    
    set(gcf,'Position',[680 351 514 454]);
    box off;
    
    hleg = legend(p(:),'Burst AM','Double SAM','Continuous double SAM','F_{sig}','location','Best');
    
    %hleg.Box = 'off'
end

%% signal + noise plot

%close all
pow_signal = cellfun(@minus,Psn,Pn,'Un',0);
pow_noise = Pn;
pow = pow_noise;
cp = cbrewer('qual','Set1',3);
sym = {'o','sq','v'};
%Fstat test here
%sqrt(F_crit(ff)-1)

for ff = 1:2
    figure(ff+2)
    for k=1:3
        if k==3
            index = 10;
        else
            index = 20;
        end
        for i=1:index
            % fstat vector goes here
            Fstat = sqrt(F_crit(ff)-1)< cell2mat(pow_signal(ff,k,i))/cell2mat(pow_noise(ff,k,i));
            
            if Fstat
                plot(i,db(squeeze(cell2mat(pow_signal(ff,k,i)))),['k-',sym{k}],'MarkerFacecolor',cp(k,:),'markersize',14);
            else
                plot(i,db(squeeze(cell2mat(pow_signal(ff,k,i)))),['k-',sym{k}],'MarkerFacecolor',[1 1 1],'markersize',14);
            end
            hold on
        end
        plot(1:index,db(squeeze(cell2mat(pow_noise(ff,k,:)))),'-','color',cp(k,:),'MarkerFacecolor',cp(k,:),'markersize',14);
        
    end
    set(gca,'xtick',0:3:length(pow(ff,k,:)),'xticklabel', (0:e_step*3:length(pow(ff,k,:))*e_step)*3,'fontsize',20);
    
    xlabel('Stim time (s)');
    ylabel('Amplitude (dB)');
    title(['fm@ ', num2str(fid(ff)), ' Hz', ', fc@2005']);
    box off;
    set(gcf,'Position',[680 351 514 454]);
end



