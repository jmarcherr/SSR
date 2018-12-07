% startup
rootdir = fileparts(which('assr_startup.m'));
datadir = ([rootdir,'/_data']);
bdfdir = ([rootdir,'/_data/raw_bdf'])


try %#ok
    rng(1); 
end  


fprintf('\n project directory now added to the current path \n')

if ~exist(fileparts(which('ft_defaults.m')))
    fprintf('remember to add fieldtrip to you path! \n')
end


addpath(fullfile('_analysis/'))
addpath _analysis/_preprocessing
addpath _analysis/_power
addpath _analysis/_itpc
addpath _analysis/uheal_efr_subctx__scripts
addpath _genstim
addpath _genstim/_func
addpath _data




fprintf('\n directory addded to the path')

SUBJECTS = {'JM2'};% subjects = SUBJECTS;
CP       = [0];


