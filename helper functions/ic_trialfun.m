function [trl, events] = ic_trialfun(cfg)
% Creates a fieldtrip trial structure based on the EEG recordings and sleep
% scoring
% cfg.dataset                   string; path to dataset
% cfg.hypnogram					string; path to hypnogram
%
% additional trl fields:
% column 4			offset of cue start
% column 5			sleep stage

%% Setup
% Check for required fields in the configuration data
requiredFields = {'dataset', 'hypnogram', 'cond_id'};
for i = requiredFields
	if ~isfield(cfg,i)
		error(['Required field missing in cfg: ' i{1} ' (%s).'], cfg.id);
	end
end

epoch_length_sec		= 30;       % length of epochs in hypnogram in s

%% Load and check data
hdr                 = ft_read_header(cfg.dataset);
events              = ft_read_event(cfg.dataset);
hyp					= load_hypnogram(cfg.hypnogram);
epoch_length_smpl	= epoch_length_sec * hdr.Fs;

trigger_tone    = events(strcmp('S  1', {events.value})); 
trigger_cue		= events(strcmp('S  2', {events.value}));

%% Sanity checks
% Data shouldn't be more than one (+1) epoch longer than the hypnogram
% (incomplete epochs are dropped by SchlafAUS) and never be shorter
if hdr.nSamples > (length(hyp)+1) * epoch_length_smpl || hdr.nSamples < length(hyp) * epoch_length_smpl
	error('Data header and hypnogram do not match.')
end

% Sanity check on the triggers
if numel(trigger_tone) ~=  numel(trigger_cue)
	warning('Unequal number of cue and tone triggers(S1/2). (%s, counter %s).\n', cfg.id, num2str(cfg.counter))
end

%% GO
trl = [];
for iTrigger = 1:numel(trigger_tone)
	trlbegin	= trigger_tone(iTrigger).sample + cfg.trialdef.pre * hdr.Fs;
	trlend		= trigger_tone(iTrigger).sample + cfg.trialdef.post * hdr.Fs;
	offset		= cfg.trialdef.pre * hdr.Fs;
	offset_cue	= trigger_cue(iTrigger).sample - trigger_tone(iTrigger).sample;
	epoch		= ceil(trigger_tone(iTrigger).sample / epoch_length_smpl);
	
	newtrl		= [trlbegin trlend offset offset_cue hyp(epoch, 1) iTrigger];
	trl			= [trl; newtrl];
end
if length(trl) ~= 30
	warning(['Subject only has ' num2str(length(trl)) ' trials (Subject ID: ' cfg.id_name ')']);
end