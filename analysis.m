
%% ------     SETUP		
% To do:
% - add fieldtrip to the path (only main folder)
% - navigate into the folder containing the contents of the zip file
% - please note: this code will create a number of subfolders
ft_defaults
addpath(fullfile(pwd, 'helper functions'))
paths                       = [];
paths.root					= pwd;
paths.home                  = enpath(fullfile(paths.root, 'analysis'));
paths.meta                  = enpath(fullfile(paths.home, 'meta'));
channel						= {'F3_1', 'F4_1', 'C3_1', 'C4_1', 'P3_1', 'P4_1'}; 
sdata                       = ic_subjectdata;

%% ------     CREATE CFG for trial definition
% This will create on cfg with all trial definitions. Those include the cue
% offset (column 4) and the sleep stage at the tone time point (column 5).
% This structure will only used once to derive artifact structures from
% in the next step. All further steps are based on that artifact structure.
paths.trialdef	= fullfile(paths.meta, 'trialdef.mat');

% Cut data around tone
trialdef = {};
for iSj = 1:numel(sdata)
	trialdef{iSj}.dataformat                = 'BrainVision';
	trialdef{iSj}.headerformat              = 'BrainVision';
	trialdef{iSj}.dataset                   = sdata(iSj).eeg;
	trialdef{iSj}.continuous                = 'no';
	trialdef{iSj}.trialdef.pre				= -4;
	trialdef{iSj}.trialdef.post	    		= 8;
	trialdef{iSj}.hypnogram					= sdata(iSj).hypno;
	trialdef{iSj}.cond_id					= sdata(iSj).condition;
	trialdef{iSj}.trialfun                  = 'ic_trialfun';
	trialdef{iSj}.id                        = sdata(iSj).id; % unique recording ID for future reference
    trialdef{iSj}.id_name                   = sdata(iSj).id_name;
	trialdef{iSj}							= ft_definetrial(trialdef{iSj});
    temp_numtrl = size(trialdef{iSj}.trl, 1);
	
 	% Note down which trials we're rejecting and which trials we're keeping
    trialdef{iSj}.trl                       = trialdef{iSj}.trl(any(trialdef{iSj}.trl(:, 5) == [2 3 4], 2),:);
    if temp_numtrl ~= size(trialdef{iSj}.trl, 1)
        warning(['Trials were excluded in subject ' trialdef{iSj}.id_name ' (id: ' num2str(trialdef{iSj}.id) ').']);
    end
end

% Get rid of empty entries in trialdef
trialdef = trialdef(~cellfun('isempty',trialdef));

%% ------     ARTIFACT REJECTION	
% For each subject, artifact information will accumulate in this cfg until
% in the very end, all artifacts will be rejected. Steps can be repeated
% and refined, in all cases the existing cfg will be loaded, altered, and
% saved again.
%
% Takes the trialdefinition cfg from above and reads in unpreprocessed data
% from disk. This allows padding of the trials (e.g. for filtering).
% Creates a file for each entry in the preprocessing cfg that contains all
% artifact information. Every artifact procedure will save its data in this
% file. In the end they can all be rejected together.
%
% The end goal is to have pairs/triplets of preprocessed data:
% . raw data before artifact rejection ('raw')
% . an artifact definition for each data set ('artifacts')
% . artifact-free data  ('clean')

paths.artifacts				= enpath(fullfile(paths.home, 'artifacts'));
paths.clean					= enpath(fullfile(paths.home, 'clean'));

for iEntry = 1:numel(trialdef)
	% --------------- Set up artifact structure ---------------
	arts                 = [];
	arts.id              = trialdef{iEntry}.id; % brings up an error
	arts.id_name         = trialdef{iEntry}.id_name;
	arts.continuous      = 'no';
	arts.trl             = trialdef{iEntry}.trl;
	arts.dataset         = trialdef{iEntry}.dataset;
	arts.artfctdef       = [];
	
	% --------------- Muscle (automatic, z-value-based) ---------------
	if ~isfield(arts.artfctdef, 'zvalue') % only if none has been set already
		arts.artfctdef.zvalue.cutoff      = 8;
	end
	arts.artfctdef.zvalue.channel     = channel; % add EMG to find muscle artifacts less obvious in the EEG
	arts.artfctdef.zvalue.trlpadding  = 0;
	arts.artfctdef.zvalue.fltpadding  = 0.2;   % only used for filtering before artifact detection (tutorial: .1)
	arts.artfctdef.zvalue.artpadding  = 0.2;   % window around artifacts still rejected
	arts.artfctdef.zvalue.detrend     = 'yes';
	arts.artfctdef.zvalue.bpfilter    = 'yes';
	arts.artfctdef.zvalue.bpfreq      = [70 90]; % cannot look higher than sampling freq / 2
	arts.artfctdef.zvalue.bpfiltord   = 12;
	arts.artfctdef.zvalue.bpfilttype  = 'but';
	arts.artfctdef.zvalue.hilbert     = 'yes';         % ?
	arts.artfctdef.zvalue.boxcar      = 0.2;           % ?
	arts.artfctdef.zvalue.interactive = 'yes';
	disp(['Showing you data with ID: ' trialdef{iEntry}.id_name '.'])
	temp                           = ft_artifact_zvalue(arts); % the second output equals cfg.artfctdef.zvalue.artifact
	arts.artfctdef                    = temp.artfctdef; % add interesting parts of current detection
	
	% --------------- Visual inspection ---------------
	% Mark remaining artifacts
	cfg_pp			= [];
	cfg_pp.trl		= arts.trl;
	cfg_pp.dataset	= arts.dataset;
	temp			= ft_preprocessing(cfg_pp);
	
	cfg_db					= arts;
	cfg_db.continuous		= 'no';
	cfg_db.ylim				= [-100 100];
	cfg_db.blocksize		= 30;
	cfg_db.selectmode		= 'markartifact';
	cfg_db.viewmode       	= 'vertical';
	cfg_db.preproc.demean  	= 'yes';
	cfg_db.preproc.detrend	= 'yes';
	cfg_db					= ft_databrowser(cfg_db, temp);
	clear temp
	arts.artfctdef  = cfg_db.artfctdef;
	
	% Ask for  channels to reject
	suspchans				= inputdlg('Enter space-separated numbers (1 = F3, 2 = F4, 3 = C3, 4 = C4, 5 = P3, 6 = P4):', 'Which channels to reject?', [1 50]);
	suspchans				= strsplit(suspchans{:});
	
	% We will use all discovererd artifacts, combine them, and cut them out.
	paths.clean					= enpath(fullfile(paths.home, 'clean'));
	path_arts					= paths.artifacts;
	path_result					= paths.clean;
	
	cfg						= trialdef{iEntry};
	cfg.headerformat		= 'brainvision_vhdr';
	cfg.dataformat			= 'brainvision_eeg';
	
	% Do actual artifact rejection
% 	cfg.dataset				= abpath(cfg.dataset);
	cfg.artfctdef.reject    = 'complete';
	cfg.artfctdef.feedback  = 'no';         % yes gives you a plot
	cfg                     = ft_rejectartifact(cfg);
	
	chans					= channel;
	if ~isempty(suspchans) && ~isempty(suspchans{1})
		chans(cellfun(@str2num,suspchans)) = [];
	end
	
	% Extend the trial length by a second to each side
	cfg.trl(:,1)			= cfg.trl(:,1) - 200;
	cfg.trl(:,2)			= cfg.trl(:,2) + 200;
	cfg.trl(:,3)			= cfg.trl(:,3) - 200;
	
	cfg.channel             = chans;
	data                    = ft_preprocessing(cfg);
	data.label				= cellfun(@(x) x(1:2), data.label, 'UniformOutput', false);
	data.id					= cfg.id;
	data.id_name			= cfg.id_name;
	realsave(fullfile(path_result, [data.id_name '_clean.mat']), data);
end

%% ------	  ERPs
paths.erp			= enpath(fullfile(paths.home, 'erp'));
path_origin			= paths.clean;
path_result			= paths.erp;
for iSj = 1:numel(sdata)
	data							= load_file(path_origin, sdata(iSj).id_name);
	
	cfg_pp							= [];
	cfg_pp.lpfilter					= 'yes';
	cfg_pp.lpfreq					= 10;
	data_pp							= ft_preprocessing(cfg_pp, data);
	
	cfg								= [];
	cfg.offset						= -1 .* data.trialinfo(:,1);
	data_re							= ft_redefinetrial(cfg, data_pp);
	
	cfg_sel							= [];
	cfg_sel.channel					= {'C3', 'C4'};
	cfg_sel.avgoverchan				= 'yes';
	cfg_sel.latency					= [-7 5];
	cfg_sel.nanmean					= 'yes';
	data_sel						= ft_selectdata(cfg_sel, data_re);
	data_sel.label					= {'AVG (C3-C4)'};
	
	cfg_tl							= [];
	cfg_tl.removemean				= 'yes';
	data_erp						= ft_timelockanalysis(cfg_tl, data_sel);
	data_erp.id						= data.id;
	
	realsave(fullfile(path_result, [data.id_name '_erp_lpfilter' num2str(cfg_pp.lpfreq) '.mat']), data_erp);
end

%% ------	  TIME-FREQUENCY TRANSFORMATION			CUE
% Takes the output of the preprocessing pipeline
tfrlimits						= [3 22]; 
paths.tfr						= enpath(fullfile(paths.home, ['tfr cue-aligned ' num2str(tfrlimits(1)) '-' num2str(tfrlimits(2))]));
paths.tfr_coll					= enpath(fullfile(paths.tfr, 'collected dataset'));

plotpathsuffix					= 'v1';
plotfilesuffix					= 'v1';
path_result						= paths.tfr;

for iSj = 1:numel(sdata)
	data							= load_file(paths.clean, sdata(iSj).id_name);
	name							= [data.id_name '_' num2str(tfrlimits(1)) '-' num2str(tfrlimits(2)) 'Hz']; %data.id
	cycles							= 10; 
	
	cfg								= [];
	cfg.offset						= -1 .* data.trialinfo(:,1);
	data_re							= ft_redefinetrial(cfg, data);
	
	% Time-Frequency Transformation
	cfg                             = [];
	cfg.method                      = 'wavelet';                        
	cfg.output                      = 'pow';
	cfg.width                       = cycles;                             
	cfg.toi                         = -7:0.005:5;
	cfg.foi                         = tfrlimits(1):0.005:tfrlimits(2);    
	data_freq						= ft_freqanalysis(cfg, data_re);
	data_freq.id					= data.id;
	data_freq.id_name				= data.id_name;
	realsave(fullfile(path_result, [data.id_name '.mat']), data_freq); 
	
	% Plotting	
	cfg						= [];
	cfg.baseline			= [-4 -3]; % the tone should start at about -2.8
	cfg.baselinetype		= 'relative';
	
	cfg.zlim				= [0 2.5];
	cfg.ylim				= [3 20];
	cfg.xlim				= [-3.5 2.5];
	cfg.showlabels			= 'no';
	cfg.layout				= 'elec1020.lay';
	
	if ismember('C3', data_freq.label)
		cfg.channel				= 'C3';
	else
		cfg.channel				= 'C4';
	end
	cfg.title				= [data_freq.id_name ' ' cfg.channel ' cue-aligned'];
	ft_singleplotTFR(cfg, data_freq) 
		
	cfg.title				= [data.id_name ' cue-aligned'];
	cfg.channel				= 'all';
	ft_multiplotTFR(cfg, data_freq)
end

% Collect datasets for each condition separately and put them into one
% array
RC = cell(1,1); RW = cell(1,1); NR = cell(1,1);
RC_cnt = 1; RW_cnt = 1; NR_cnt = 1;
files = get_filenames(path_result, 'full');
for iFile = 1:numel(files)
	[path, name, ext]				= fileparts(files{iFile});
	data	= load_file(files{iFile});
	s		= sdata(cellfun(@(x) isequal(x, data.id), {sdata.id})); % get the correct sdata enty
	
	switch   s.condition
		case 1
			RC{RC_cnt} = data; RC_cnt = RC_cnt + 1;
		case 2
			RW{RW_cnt} = data; RW_cnt = RW_cnt + 1;
		case 3
			NR{NR_cnt} = data; NR_cnt = NR_cnt + 1;
	end
	
end
realsave(fullfile(paths.tfr_coll, ['all_datasets_RC.mat']), RC)
% realsave(fullfile(paths.tfr_coll, ['all_datasets_NR.mat']), NR) 
% realsave(fullfile(paths.tfr_coll, ['all_datasets_RW.mat']), RW)

% Grand Average plots - ONLY MAKES SENSE FOR MORE THAN ONE DATASET
% files = get_filenames(paths.tfr_coll, 'full');
% for iFile = 1:numel(files)
% 	[path, name, ext]		= fileparts(files{iFile});
% 	data					= load_file(files{iFile});
% 	if ~isempty(data{1})
% 		% Combine contra-lateral channels (nanmean)
% 		for iData = 1:numel(data)
% 			labs1  = {'C3','C4'};
% 			labs2  = {'F3','F4'};
% 			labs3  = {'P3','P4'};
% 			data{iData}.powspctrm_merged = cat(1,...
% 				nanmean(data{iData}.powspctrm(ismember(data{iData}.label,labs1),:,:),1),...
% 				nanmean(data{iData}.powspctrm(ismember(data{iData}.label,labs2),:,:),1),...
% 				nanmean(data{iData}.powspctrm(ismember(data{iData}.label,labs3),:,:),1));
% 			data{iData}.label = {'C', 'F', 'P'};
% 			data{iData}.powspctrm = data{iData}.powspctrm_merged;
% 			data{iData} = rmfield(data{iData}, 'powspctrm_merged');
% 		end
% 		
% 		cfg                     = [];
% 		data_ga                 = ft_freqgrandaverage(cfg, data{:});
% 		data_ga.id				= name;
% 		
% 		cfg						= [];
% 		cfg.baseline			= [-4 -3];
% 		cfg.baselinetype		= 'relative';
% 		cfg.channel				= 'C';
% 		
% 		cfg.ylim				= [3 20];
% 		cfg.xlim				= [-3.5 2.5];
% 		cfg.showlabels			= 'no';
% 		cfg.layout				= 'elec1020.lay';
% 		cfg.title				= ['grand average ' name(end-1:end) ' ' cfg.channel  ' cue-aligned'];
% 		figure, ft_singleplotTFR(cfg, data_ga)
% 		
% 		path_plot = enpath(fullfile(paths.tfr_coll, 'plots GA'));
% 		
% 		cfg.title				= ['grand average ' name(end-1:end) ' cue-aligned'];
% 		cfg.channel				= 'all';
% 		figure, ft_multiplotTFR(cfg, data_freq)
% 	end
% end

%% ------     PREPARE TFR DATA FOR FINAL PLOTTING AND STATS (incl. averaging contralateral channel and baselining)
% Baselining here is relative (will be changed to relative change later)
path_origin						= paths.tfr_coll;
paths.tfr_collprep 				= enpath(fullfile(path_origin, 'prep'));
files							= get_filenames(path_origin, 'full');

for iFile = 1:numel(files)
	temp			= load_file(files{iFile});
	[~, name, ~]   = fileparts(files{iFile});
	% For each dataset inside the file, let's average contralateral
	% channels and baseline
	for iDs = 1:numel(temp)
		% Average the channels (with code above)
		labs1  = {'C3','C4'};
		labs2  = {'F3','F4'};
		labs3  = {'P3','P4'};
		
		temp{iDs}.powspctrm_merged = cat(1,...
			nanmean(temp{iDs}.powspctrm(ismember(temp{iDs}.label,labs1),:,:),1),...
			nanmean(temp{iDs}.powspctrm(ismember(temp{iDs}.label,labs2),:,:),1),...
			nanmean(temp{iDs}.powspctrm(ismember(temp{iDs}.label,labs3),:,:),1));
		temp{iDs}.label = {'C', 'F', 'P'};
		temp{iDs}.powspctrm = temp{iDs}.powspctrm_merged;
		temp{iDs} = rmfield(temp{iDs}, 'powspctrm_merged');
		
		% Do a relative baselining (we'll subtract -1 later on to make it
		% vary around 0)
		cfg_bl				= [];
		cfg_bl.baseline		= [-4 -3];	% [-4 -3]
		cfg_bl.baselinetype = 'relative';
		temp{iDs}			= ft_freqbaseline(cfg_bl, temp{iDs});
	end
	realsave(fullfile(paths.tfr_collprep, [name '_prep.mat']), temp);
end

%% ------     TIME-FREQUENCY STATISTICS       -  cannot be tested with only one dataset        
% Cannot be tested with only one dataset. Furthermore, this requires a lot
% of RAM. Stats were run on a cluster. For running it locally, set number
% of permutations to 1000 to reduce RAM usage.
paths.tfr_stats		= enpath(fullfile(paths.tfr_collprep, 'stats5000 sample-level .01'));
path_origin         = paths.tfr_collprep;

% Fish out those files belonging to the current experiment
files           = get_filenames(path_origin, 'full');

% Original code
[~,t,~]         = fileparts(files{1}); t = tokenize(t, '_'); if ~strcmp(t(3), 'NR'), error('Not the right dataset.'), end
[~,t,~]         = fileparts(files{2}); t = tokenize(t, '_'); if ~strcmp(t(3), 'RC'), error('Not the right dataset.'), end
[~,t,~]         = fileparts(files{3}); t = tokenize(t, '_'); if ~strcmp(t(3), 'RW'), error('Not the right dataset.'), end
NR              = load_file(files{1});
RC              = load_file(files{2});
RW              = load_file(files{3});

% Do the statistics
cfg                     = [];
cfg.latency             = [-3.5 2.5];
cfg.frequency           = 'all';
cfg.channel             = 'C';
%     cfg.avgoverchan			= 'yes';
cfg.correctm            = 'cluster';
cfg.method              = 'montecarlo';
cfg.statistic           = 'indepsamplesT';     % use actvsblT for activation against baseline
cfg.clusterstatistic    = 'maxsum';            % statistic used to decide cluster significance (sum of t-values within a cluster)
cfg.clustertail         = 0;
cfg.tail                = 0;
cfg.alpha               = 0.025;
cfg.numrandomization    = 5000;             
cfg.clusteralpha        = .01;              % threshold over which a triplet is chosen, e.g. .01 / .02 / .05
cfg.ivar                = 1;                 % condition (uvar would be the subjects)

% Contrast 1: RC vs. NR
% Design the statistical contrast
design                  = [];
design(1,:)             = [ones(1,length(RC)) ones(1,length(NR))*2];        % conditions, eg:   1 1 1 1 2 2 2 2
cfg.design              = design;
stats1                  = ft_freqstatistics(cfg, RC{:}, NR{:});
realsave(fullfile(paths.tfr_stats, ['Statistics RC vs NR' '_calpha' num2str(cfg.clusteralpha) '_alpha' num2str(cfg.alpha) '_' cfg.statistic '_' cfg.correctm '_C_cuealigned_pretone_baseline' '.mat']), stats1);
clear stats1

% Contrast 2: RW vs. NR
design                  = [];
design(1,:)             = [ones(1,length(RW)) ones(1,length(NR))*2];        % conditions, eg:   1 1 1 1 2 2 2 2
cfg.design              = design;
stats2                   = ft_freqstatistics(cfg, RW{:}, NR{:});
realsave(fullfile(paths.tfr_stats, ['Statistics RW vs NR' '_calpha' num2str(cfg.clusteralpha) '_alpha' num2str(cfg.alpha) '_' cfg.statistic '_' cfg.correctm '_C_cuealigned_pretone_baseline' '.mat']), stats2);
clear stats2

% Contrast 3: RW vs. RC
design                  = [];
design(1,:)             = [ones(1,length(RW)) ones(1,length(RC))*2];        % conditions, eg:   1 1 1 1 2 2 2 2
cfg.design              = design;
stats3                  = ft_freqstatistics(cfg, RW{:}, RC{:});
realsave(fullfile(paths.tfr_stats, ['Statistics RW vs RC' '_calpha' num2str(cfg.clusteralpha) '_alpha' num2str(cfg.alpha) '_' cfg.statistic '_' cfg.correctm '_C_cuealigned_pretone_baseline' '.mat']), stats3);
clear stats3

%% ------     COLLECT ERP DATA
paths.erp_coll      = enpath(fullfile(paths.erp, 'collected dataset'));
path_origin         = paths.erp;
path_result         = paths.erp_coll;

% Collect datasets for each condition separately and put them into one
% array
RC = cell(1,1); RW = cell(1,1); NR = cell(1,1);
RC_cnt = 1; RW_cnt = 1; NR_cnt = 1;
files = get_filenames(path_origin, 'full');
for iFile = 1:numel(files)
	[path, name, ext]				= fileparts(files{iFile});
	
	data	= load_file(files{iFile});
	s		= sdata(cellfun(@(x) isequal(x, data.id), {sdata.id})); % get the correct sdata enty
	
	switch   s.condition
		case 1
			RC{RC_cnt} = data; RC_cnt = RC_cnt + 1;
		case 2
			RW{RW_cnt} = data; RW_cnt = RW_cnt + 1;
		case 3
			NR{NR_cnt} = data; NR_cnt = NR_cnt + 1;
	end
end

if ~isempty(RC{1}), realsave(fullfile(path_result, ['all_datasets_RC_5000p.mat']), RC), end
if ~isempty(NR{1}), realsave(fullfile(path_result, ['all_datasets_NR_5000p.mat']), NR), end
if ~isempty(RW{1}), realsave(fullfile(path_result, ['all_datasets_RW_5000p.mat']), RW), end

%% ------     FINAL PLOTTING
path_stats = paths.tfr_stats; % can  either paths.tfr_stats or paths.tfr_stats_bl

% For each contrast, collect respective datasets (erp, tfr, stats)
data_erp	= cell(1,3);
files_erp	= get_filenames(paths.erp_coll, 'full');
for iFile = 1:numel(files_erp)          
	[path, name, ext]	= fileparts(files_erp{iFile});
	name				= tokenize(name, '_');
	data_erp{iFile}		= load_file(files_erp{iFile});
end

data_tfr	= cell(1,3);
files_tfr	= get_filenames(paths.tfr_collprep, 'full');
for iFile = 1:numel(files_tfr)
	[path, name, ext]	= fileparts(files_tfr{iFile});
	name				= tokenize(name, '_');
	data_tfr{iFile}		= load_file(files_tfr{iFile});
end

data_stats	= cell(1,3);
files_stats	= get_filenames(path_stats, 'full');
for iFile = 1:numel(files_stats)
	[path, name, ext]	= fileparts(files_stats{iFile});
	data_stats{iFile}	= load_file(files_stats{iFile});
end

% Now all datasets are in the same order (NR - RC - RW)

% In case there is only one dataset
for iCond = 1:numel(data_erp)
	if isstruct(data_erp{iCond}) % in case of only 1 dataset
		tmp = data_erp{iCond};
		data_erp{iCond} = [];
		data_erp{iCond}{1} = tmp;
	end
end
for iCond = 1:numel(data_tfr)
	if isstruct(data_tfr{iCond}) % in case of only 1 dataset
		tmp = data_tfr{iCond};
		data_tfr{iCond} = [];
		data_tfr{iCond}{1} = tmp;
	end
end

% ------------- Baseline ERP & TFR
% Average ERP of experimental groups incl. baseline
for iCond = 1:numel(data_erp)
	if ~isempty(data_erp{iCond})
		for iFile = 1:numel(data_erp{iCond})
			cfg						= [];
			cfg.baseline			= [-4 -3];
			data_erp{iCond}{iFile}	= ft_timelockbaseline(cfg, data_erp{iCond}{iFile});
		end
		cfg				= [];
		cfg.latency     = [-3.5 2.5];
		data_erp{iCond} = ft_timelockgrandaverage(cfg, data_erp{iCond}{:});
	end
end

% Average TFR of experimental groups, subtract 1 so it varies around 0
for iCond = 1:numel(data_tfr)
	if ~isempty(data_tfr{iCond})
		for iFile = 1:numel(data_tfr{iCond})
			cfg						= [];
			cfg.operation			= 'subtract';
			cfg.parameter			= 'powspctrm';
			cfg.scalar				= 1;
			data_tfr{iCond}{iFile}	= ft_math(cfg, data_tfr{iCond}{iFile});
		end
		cfg				= [];
		cfg.toilim      = [-3.5 2.5];
		cfg.channel		= 'C';
		data_tfr{iCond} = ft_freqgrandaverage(cfg, data_tfr{iCond}{:});
	end
end

% Plot three contrasts
% 1		RC vs NR
% 2		RW vs NR
% 3		RW vs RC

% ----------------------
% Plots are not possible with test dataset ant thus commented out.
% Please find code for plotting the test data set below.
% ----------------------

% % ------------- PLOT -- RC vs NR     ---------------------
% cfg_pl						= [];
% cfg_pl.ylim					= [3 22];
% cfg_pl.zlim					= [-.8 .8];
% cfg_pl.xlim					= [-3.5 2.5];
% cfg_pl.showlabels			= 'no';
% cfg_pl.layout				= 'elec1020.lay';
% cfg_pl.maskstyle			= 'opacity';
% cfg_pl.maskalpha			= .4;
% 
% % data_tfr{2}.mask			= data_stats{1}.mask; % no statistics for test dataset
% cfg_pl.maskparameter		= 'mask';
% cfg_pl.title				= ['RC vs NR, C, plotted cuealigned against baseline'];
% f = figure;
% ft_singleplotTFR(cfg_pl, data_tfr{2})
% 
% 
% % Adjust the plot and add ERP
% a1 = gca;
% set(f,'Position', [50 100 1800 900])
% set(a1,'tickdir','out');
% set(a1,'ticklength', [.01 .01]);
% set(a1, 'xlim', [-3.5 2.5]);
% set(a1, 'ylim', [3 22]);
% set(a1, 'fontsize',12);
% set(a1, 'YTick',2:2:22);
% % set(a1, 'YTickLabel',0:2:20);
% temp_pos	= get(a1,'pos');
% tempXLim	= get(a1,'XLim');
% tempYlim	= get(a1,'Ylim');
% 
% hold on;                                  % hold to superimpose the SO average
% a2 = axes;
% plot(data_erp{2}.time, data_erp{2}.avg, 'Color', 'k', 'LineWidth', 2);  % HERE (2x)
% set(a2, 'YAxisLocation', 'Right');         % bring the y axis to the right
% set(a2, 'color', 'none');                  % hide it
% set(a2, 'XLim', tempXLim);				  % scale the new x axis to fit the old one
% set(a2, 'XTick', []);                      % and hide it
% set(a2, 'YLim', [-25 25]);                 % scale the new y axis
% set(a2, 'YTick',-25:5:25);                  % define ticks on the y axis
% set(a2, 'Position',temp_pos)
% set(a2, 'fontsize',12);
% 
% % Move the colorbar a bit to the right
% cHandle			= findobj('tag','ft-colorbar');       % find all the colorbar handles
% cPosition		= get(cHandle,'position');
% cPosition(1)	= cPosition(1) + .02;
% set(cHandle,'position',cPosition);
% set(cHandle, 'Ticks', -1:.5:1)
% set(a1, 'Position',temp_pos)
% 
% %path_plot		= enpath(fullfile(path_stats, 'plots_5000per'));
% 
% path_plot = 'Y:\Julia\Incomplete cueing\analysis\freqanalysis 1.2\tfr cue-aligned 3-22\collected dataset\prep\plots stats5000 sample-level .01'
% alpha			= data_stats{1}.cfg.alpha;
% clusteralpha	= data_stats{1}.cfg.clusteralpha;
% statistic		= data_stats{1}.cfg.statistic;
% correctm		= data_stats{1}.cfg.correctm;
% export_fig(fullfile(path_plot, ['Statistics__5000perm_sample-level .01_RCvsNR' '_calpha' num2str(clusteralpha) '_alpha' num2str(alpha) '_' statistic '_' correctm '_C4_cuealigned_pretone_baseline'  '.png']),  '-nocrop', '-a4', '-CMYK', '-q90', '-m2');
% close all

% % ------------- PLOT -- RW vs NR     ---------------------
% cfg_pl						= [];
% cfg_pl.ylim					= [3 22];
% cfg_pl.zlim					= [-.8 .8];
% cfg_pl.xlim					= [-3.5 2.5];
% cfg_pl.showlabels			= 'no';
% cfg_pl.layout				= 'elec1020.lay';
% cfg_pl.maskstyle			= 'opacity';
% cfg_pl.maskalpha			= .4;
% 
% % newmask = or(stats1.posclusterslabelmat == 1, stats1.posclusterslabelmat == 2);
% data_tfr{3}.mask			= data_stats{2}.mask;
% cfg_pl.maskparameter		= 'mask';
% cfg_pl.title				= ['RW vs NR, C, plotted cuealigned against baseline'];
% f = figure;
% ft_singleplotTFR(cfg_pl, data_tfr{3})
% % close all
% 
% % Adjust the plot and add ERP
% a1 = gca;
% set(f,'Position', [50 100 1800 900])
% set(a1,'tickdir','out');
% set(a1,'ticklength', [.01 .01]);
% set(a1, 'xlim', [-3.5 2.5]);
% set(a1, 'ylim', [3 22]);
% set(a1, 'fontsize',12);
% set(a1, 'YTick',2:2:22);
% % set(a1, 'YTickLabel',0:2:20);
% temp_pos	= get(a1,'pos');
% tempXLim	= get(a1,'XLim');
% tempYlim	= get(a1,'Ylim');
% 
% hold on;                                  % hold to superimpose the SO average
% a2 = axes;
% plot(data_erp{3}.time, data_erp{3}.avg, 'Color', 'k', 'LineWidth', 2);  % HERE (2x)
% set(a2, 'YAxisLocation', 'Right');         % bring the y axis to the right
% set(a2, 'color', 'none');                  % hide it
% set(a2, 'XLim', tempXLim);				  % scale the new x axis to fit the old one
% set(a2, 'XTick', []);                      % and hide it
% set(a2, 'YLim', [-25 25]);                 % scale the new y axis
% set(a2, 'YTick',-25:5:25);                  % define ticks on the y axis
% set(a2, 'Position',temp_pos)
% set(a2, 'fontsize',12);
% 
% % Move the colorbar a bit to the right
% cHandle			= findobj('tag','ft-colorbar');       % find all the colorbar handles
% cPosition		= get(cHandle,'position');
% cPosition(1)	= cPosition(1) + .02;
% set(cHandle,'position',cPosition);
% set(cHandle, 'Ticks', -1:.5:1)
% set(a1, 'Position',temp_pos)
% 
% path_plot = 'Y:\Julia\Incomplete cueing\analysis\freqanalysis 1.2\tfr cue-aligned 3-22\collected dataset\prep\plots stats5000 sample-level .01'
% alpha			= data_stats{2}.cfg.alpha;
% clusteralpha	= data_stats{2}.cfg.clusteralpha;
% statistic		= data_stats{2}.cfg.statistic;
% correctm		= data_stats{2}.cfg.correctm;
% export_fig(fullfile(path_plot, ['Statistics__5000perm_sample-level .01_RWvsNR' '_calpha' num2str(clusteralpha) '_alpha' num2str(alpha) '_' statistic '_' correctm '_C4_cuealigned_pretone_baseline'  '.png']),  '-nocrop', '-a4', '-CMYK', '-q90', '-m2');
% close all

% % ------------- PLOT -- RC vs RW -------------
% cfg_pl						= [];
% cfg_pl.ylim					= [3 22];
% cfg_pl.zlim					= [-.8 .8];
% cfg_pl.xlim					= [-3.5 2.5];
% cfg_pl.showlabels			= 'no';
% cfg_pl.layout				= 'elec1020.lay';
% cfg_pl.maskstyle			= 'opacity';
% cfg_pl.maskalpha			= .4;
% 
% % Lets create a new TFR RW - RC
% temp_tfr = data_tfr{3};
% temp_tfr.powspctrm = temp_tfr.powspctrm - data_tfr{2}.powspctrm;
% 
% % Lets create a mask that also plots trends
% temp_mask = data_stats{3}.mask;
% %temp_mask(data_stats{3}.prob > .025 & data_stats{3}.prob <= .05) = 1; 
% temp_tfr.mask			= temp_mask;
% 
% cfg_pl.maskparameter		= 'mask';
% cfg_pl.title				= ['RW vs RC, C, plot is RW - RC, cuealigned against baseline'];
% f = figure;
% ft_singleplotTFR(cfg_pl, temp_tfr)
% % close all
% 
% % Adjust the plot and add ERP
% a1 = gca;
% set(f,'Position', [50 100 1800 900])
% set(a1,'tickdir','out');
% set(a1,'ticklength', [.01 .01]);
% set(a1, 'xlim', [-3.5 2.5]);
% set(a1, 'ylim', [3 22]);
% set(a1, 'fontsize',12);
% set(a1, 'YTick',2:2:22);
% % set(a1, 'YTickLabel',0:2:20);
% temp_pos	= get(a1,'pos');
% tempXLim	= get(a1,'XLim');
% tempYlim	= get(a1,'Ylim');
% 
% hold on;                                  % hold to superimpose the SO average
% a2 = axes;
% plot(data_erp{3}.time, data_erp{3}.avg - data_erp{2}.avg, 'Color', 'k', 'LineWidth', 2);  % HERE (2x)
% set(a2, 'YAxisLocation', 'Right');         % bring the y axis to the right
% set(a2, 'color', 'none');                  % hide it
% set(a2, 'XLim', tempXLim);				  % scale the new x axis to fit the old one
% set(a2, 'XTick', []);                      % and hide it
% set(a2, 'YLim', [-25 25]);                 % scale the new y axis
% set(a2, 'YTick',-25:5:25);                  % define ticks on the y axis
% set(a2, 'Position',temp_pos)
% set(a2, 'fontsize',12);
% 
% % Move the colorbar a bit to the right
% cHandle			= findobj('tag','ft-colorbar');       % find all the colorbar handles
% cPosition		= get(cHandle,'position');
% cPosition(1)	= cPosition(1) + .02;
% set(cHandle,'position',cPosition);
% set(cHandle, 'Ticks', -1:.5:1)
% set(a1, 'Position',temp_pos)
% 
% path_plot = 'Y:\Julia\Incomplete cueing\analysis\freqanalysis 1.2\tfr cue-aligned 3-22\collected dataset\prep\plots stats5000 sample-level .01'
% 
% alpha			= data_stats{3}.cfg.alpha;
% clusteralpha	= data_stats{3}.cfg.clusteralpha;
% statistic		= data_stats{3}.cfg.statistic;
% correctm		= data_stats{3}.cfg.correctm;
% export_fig(fullfile(path_plot, ['Statistics__5000perm_sample-level .01_RWvsRC' '_calpha' num2str(clusteralpha) '_alpha' num2str(alpha) '_' statistic '_' correctm '_C4_cuealigned_pretone_baseline'  '.png']),  '-nocrop', '-a4', '-CMYK', '-q90', '-m2');
% close all

% ------------- PLOT -- Test dataset     ---------------------
cfg_pl						= [];
cfg_pl.ylim					= [3 22];
cfg_pl.zlim					= [-.8 .8];
cfg_pl.xlim					= [-3.5 2.5];
cfg_pl.showlabels			= 'no';
cfg_pl.layout				= 'elec1020.lay';
cfg_pl.maskstyle			= 'opacity';
cfg_pl.maskalpha			= .4;

% data_tfr{2}.mask			= data_stats{1}.mask; % no statistics for test dataset
% cfg_pl.maskparameter		= 'mask';
f = figure;
ft_singleplotTFR(cfg_pl, data_tfr{1})


% Adjust the plot and add ERP
a1 = gca;
set(f,'Position', [50 100 1800 900])
set(a1,'tickdir','out');
set(a1,'ticklength', [.01 .01]);
set(a1, 'xlim', [-3.5 2.5]);
set(a1, 'ylim', [3 22]);
set(a1, 'fontsize',12);
set(a1, 'YTick',2:2:22);
% set(a1, 'YTickLabel',0:2:20);
temp_pos	= get(a1,'pos');
tempXLim	= get(a1,'XLim');
tempYlim	= get(a1,'Ylim');

hold on;                                  % hold to superimpose the SO average
a2 = axes;
plot(data_erp{1}.time, data_erp{1}.avg, 'Color', 'k', 'LineWidth', 2);  % HERE (2x)
set(a2, 'YAxisLocation', 'Right');         % bring the y axis to the right
set(a2, 'color', 'none');                  % hide it
set(a2, 'XLim', tempXLim);				  % scale the new x axis to fit the old one
set(a2, 'XTick', []);                      % and hide it
set(a2, 'YLim', [-25 25]);                 % scale the new y axis
set(a2, 'YTick',-25:5:25);                  % define ticks on the y axis
set(a2, 'Position',temp_pos)
set(a2, 'fontsize',12);

% Move the colorbar a bit to the right
cHandle			= findobj('tag','ft-colorbar');       % find all the colorbar handles
cPosition		= get(cHandle,'position');
cPosition(1)	= cPosition(1) + .02;
set(cHandle,'position',cPosition);
set(cHandle, 'Ticks', -1:.5:1)
set(a1, 'Position',temp_pos)
