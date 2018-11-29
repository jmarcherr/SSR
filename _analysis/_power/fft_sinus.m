clear all

fs = 1024;
t = 0:1/fs:3-(1/fs);

test_sig = sin(2*pi*93*t);
plot(t,test_sig)

test_s = repmat(test_sig,1,10);
M=test_s;
    f_fft = fft(M)%./(size(M,2)./2); % length of vector
    %f_fft(kk,:,:) = fft(M,2.^nextpow2(length(M)),2); %nextpow2

f_fft_mean = squeeze(nanmean(f_fft,1));
f_fft_pow = abs(f_fft_mean.^2);
f = linspace(0,fs-(1/fs),length(f_fft_pow)); 


plot(f(1:end/2),f_fft_pow(1:end/2))



%%

% Frequency defition
t = 0:1/fs:(30-(1/fs));
df = 1/(30);
fn = length(M)/30;