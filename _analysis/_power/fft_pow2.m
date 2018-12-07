clear all
assr_startup
ft_defaults

%% Select preprocessed .mat file to process
cd(datadir)
load(['pilot3_2xsamtone+oldx10min.mat']);
cd(rootdir)
%% rename data
data_all= data_DG;
cfg = [];
cfg.channel = {'all','-Status','-Cz'}
data_all = ft_selectdata(cfg,data_all)

%% Select trials based on triggers (see preprocessing)
close all
%%
for kk=1:2
idx=find(data_all.trialinfo==kk)
idx_tmp=idx%(1:200);
% how many epochs to group together
epoch = 100;
% generate epochs
cc=0;
epoched_data = [];
for ii = 1:epoch:length(idx_tmp)
    cc=cc+1;
    cfg = [];
    cfg.trials = idx_tmp(ii:ii+epoch-1);
    data_tmp = ft_selectdata(cfg,data_all);
    epoched_data(cc,:,:) = cell2mat(data_tmp.trial(:)'); % epoch x chan x time

end
 

    
%% perform FFT on epoched data
f_fft = [];fs = data_all.fsample;
M = [];
for kk=1:size(epoched_data,1) % run through epochs
    M=squeeze(epoched_data(kk,:,:));
    f_fft(kk,:,:) = fft(M,length(M),2)/(size(M,2)/2); % length of vector
    %f_fft(kk,:,:) = fft(M,2.^nextpow2(length(M)),2); %nextpow2
end
f_fft_mean = squeeze(nanmean(f_fft,1));
f_fft_pow = abs(f_fft_mean.^2);
f_fft_pow = (f_fft_pow(:,1:end/2+1));

f = fs/2*linspace(0,1,length(f_fft_mean)/2+1);

%%
figure(kk+1)
for ii=1:length(data_tmp.label)
subplot(3,3,ii)
%plot(f,10*log10(f_fft_pow(find(strcmp(data_tmp.label,'EXG1')),:)));
plot(f,10*log10(f_fft_pow(find(strcmp(data_tmp.label,data_tmp.label{ii})),:)));
xlim([1 5])
title(data_tmp.label{ii})
end

%% finding peak
pow_207 = [];
pow_4 = [];
%fc's
fc_idx = [4 207];
fc(1)=find(f==fc_idx(1));
fc(2) = find(f==fc_idx(2));

bg_freq = [];
for ii=1:length(fc)
% background (fc +/-1.5-0.5Hz)
bg_freq{ii} = [find(f>=fc_idx(ii)-1.5 & f<=fc_idx(ii)-0.5), ...
    find(f<=fc_idx(ii)+1.5 & f>=fc_idx(ii)+0.5)];
end


for i=1:length(data_tmp.label)
    pow_207(i) = 10*log10(f_fft_pow(i,fc(2)))/10*log10(mean(f_fft_pow(i,bg_freq{2}),2));
    pow_4(i)   = 10*log10(f_fft_pow(i,fc(1)))/10*log10(mean(f_fft_pow(i,bg_freq{2}),2));
end
figure(99)
subplot(2,1,1)
plot(pow_207)
hold on
set(gca,'xticklabels',data_tmp.label)
subplot(2,1,2)
plot(pow_4)
hold on
set(gca,'xticklabels',data_tmp.label)
hold off
end
legend('old','new','new continuous')


%%

