cfg = struct;
cfg.fc = fc;
%cfg.fmod = fmod;
cfg.moddepth = moddepth;
cfg.stimrate = stimrate;
fm_ctx = cfg.stimrate*2;
cfg.fs = fs;
cfg.Lstim = Lstim;
cfg.L = L;
cfg.stimtype = stimtype;
cfg.gain_factor = gain_factor;
t = 0:1/fs:Lstim-1/fs;


sam_double = ( ( 1 + moddepth*sin(2*pi* fm_ctx*t - pi/2) )/2 )...
    .* ( ( 1 + moddepth*sin(2*pi* fmod*t- pi/2) )/2 ) .* ...
    sin(2*pi*fc*t);

amp=1;
fs = 48000;
m_ctx = .85;
fm_ctx = 4;
m_subctx = .85;
fm_subctx= 207;
t_vect = 0:1/fs:3;
fc = 2005;



plot(sam_double)
soundsc(sam_double,fs)
