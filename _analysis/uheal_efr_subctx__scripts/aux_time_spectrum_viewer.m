function aux_time_spectrum_viewer( data_in, fs, nfft_flag, freq_ax_min,...
    freq_ax_max, y_lim_min, y_lim_max, play_sound_flag, fig_title )

% time_spectrum_viewer
%
%   Function that plots temporal waveform and Spectrum of a signal
%   aux_time_spectrum_viewer( data_in, fs, nfft_flag, freq_ax_min,...
%    freq_ax_max, y_lim_min, y_lim_max, play_sound_flag, fig_title )
%
% EXAMPLE:
%
%       aux_time_spectrum_viewer( sig, 48000, 'sig_lgth', 100,...
%    15000, -20, 50, false, 'Pure tone in noise' )

%% Initial paremeters
len_data_in = numel(data_in);
time_vect = (0:1/fs:len_data_in/fs-1/fs)';

switch nfft_flag
    case 'nxt_pwr_2'
        nfft = 2^nextpow2(len_data_in);
    case 'sig_lgth'
        nfft = len_data_in;
    otherwise
        error('ERROR! Invalid string statement. Choose between "nxt_pwr_2" or "sig_lgth".')
        return;
end

f_axis = fs/2*linspace(0,1,nfft/2+1);

%% Compute FFT
data_in_fft         = fft(data_in, nfft)/(len_data_in/2);
data_in_fft_magn    = 10 * log10( abs(data_in_fft).^2 );
data_in_fft_angle   = angle(data_in_fft);

%% Plot waveform and Spectrum

% Waveform
figure('Name', fig_title, 'NumberTitle', 'off'),
subplot(311);
plot(time_vect, data_in)
xlabel('Time [seconds]')
ylabel('Amplitude')
title('Temporal waveform')
ylim([min(data_in)-0.1*min(data_in) max(data_in)+0.1*max(data_in)])
xlim([time_vect(1) time_vect(end)])

subplot(312)
semilogx(f_axis, data_in_fft_magn(1:floor(nfft/2)+1))
xlim([freq_ax_min freq_ax_max])
ylim([y_lim_min y_lim_max])
xlabel('Frequency [Hz]')
ylabel('Magnitude [dB]')
title('Spectrum')

subplot(313)
semilogx(f_axis, unwrap( data_in_fft_angle(1:floor(nfft/2)+1))/(2*pi) )
xlim([freq_ax_min freq_ax_max])
phase_vector = unwrap( data_in_fft_angle(1:floor(nfft/2)+1))/(2*pi);
ylim([min(phase_vector)-0.1*min(phase_vector) max(phase_vector)+0.1*max(phase_vector)])
xlabel('Frequency [Hz]')
ylabel('Phase [cycles]')
title('Phase')

if play_sound_flag
    sound(data_in, fs);
end

%EoF
end
