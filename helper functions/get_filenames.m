function filename = get_filenames(path, varargin)
% Returns a list of files in a given folder or the filename of the nth file
% in that folder or a filename starting with a given string or combinations
% of all those.
%
% Use as
% filename = get_filenames(path, n)
% eg.: filename = get_filenames('C:\files\')
%      filelist = get_filenames('C:\files\', 3)
%      filelist = get_filenames('C:\files\', 3, 'full')
%	   filelist = get_filenames('C:\files\', 'RC_051', 'full')
%
% INPUT VARIABLES:
% path              String; input folder
% full              string (optional)
%                   If 'full' the full path(s) will be returned
% match             string or int (optional)
%                   Every string other than 'full' will be matched with the
%                   beginning of all file names. If 'all' is not set (see below) only one matching filename
%                   will be returned. If there is more than one match an
%                   error will be thrown. If there is no match a warning
%                   will be thrown and the output will be empty. Instead of
%                   a matching string, the index of the desired file can
%                   be provided (ordered alphabetically). If a match string
%                   and an index n is provided, the nth file among those
%                   matching the string is returned.
%
%              !!   The matching feature is intended to allow the search
%                   for files belonging to a specific subject or recording
%                   ID. Because of that, the string is only considered a
%                   match IF IN THE FILENAME IT IS FOLLOWED BY '.' OR '_'
%                   (eg. 's42_rec23' -> 's42_rec23_filtered.mat').
% all				string (optional)
%					If 'all' all matches will be returned
%
% MFF-Files are treated as files, although windows treats them as folders.
% The DS_Store nonsense files that Macs leave behind are ignored.
%
% OUTPUT VARIABLES:
% filename          names or list of names (in a cell array) of requested
%                   files IN THE ONLY PROPER WAY TO SORT FILENAMES, even
%                   if matlab is of a different opinion!
%                   (a1.mat, a5.mat, a12.mat...)
%
%
% AUTHOR:
% Jens Klinzing, jens.klinzing@uni-tuebingen.de

%% SETUP
all             = dir(path);               % Get the data for the current directory
idx_dirs        = [all.isdir];             % Find the index for directories

% Make sure .mff files are treated as files not folders and Mac nonsense is
% treated as such and not even mentioned.
names = {all.name};
for i = 1:numel(names)
	if length(names{i}) > 4 && strcmpi(names{i}(end-3:end),'.mff')
		idx_dirs(i) = 0;
	elseif strcmpi(names{i},'.DS_Store') || strcmpi(names{i},'._.DS_Store')
		idx_dirs(i) = 1;
	end
end

fileList        = natsort({all(~idx_dirs).name}');  % ...delete them out of the list and sort them like real men

match           = [];                      % default: dont try to match a string
returnall		= false;				   % default: dont return all if there is more than one match
n               = [];                      % default: return all files in the directory
full            = false;                   % default: dont add path to filename

if nargin < 1 || nargin > 4
	error('Unexpected number of input arguments.');
elseif nargin > 1 % handle potential varargin arguments
	if any(strcmp(varargin, 'full'))       % look for 'full' argument
		full = true;
		varargin(strcmp(varargin, 'full')) = []; % makes life easier looking for another string
	end
	if any(strcmp(varargin, 'all'))       % look for 'full' argument
		returnall = true;
		varargin(strcmp(varargin, 'all')) = []; % makes life easier looking for another string
	end
	if any(cellfun(@ischar, varargin))     % look for match argument
		if sum(cellfun(@ischar, varargin)) == 1
			match = varargin{cellfun(@ischar, varargin)};
		else
			error('More than one matching string provided. Cannot compute.')
		end
	end
	if any(cellfun(@isnumeric, varargin))   % look for n argument
		if numel(varargin{cellfun(@isnumeric, varargin)}) ~= 1
			error('More than one numerical value provided. Cannot compute.')
		else
			n = varargin{cellfun(@isnumeric, varargin)};
			if returnall
				error('You cannot ask for all matches but also provide a file index.')
			end
		end
	end
end

% Get rid of underscores or dots in the end, those are obligatory anyways
if ~isempty(match) && (strcmp(match(end), '_') || strcmp(match(end), '.'))
	match = match(1:end-1);
end

%% START
% Are any selectors given?
m = [];		% file indices of all matches to the string
if ~isempty(match)
	% Find all matching files
	for iFile = 1:numel(fileList)
		if length(fileList{iFile}) >= length(match)+1 && (strcmp(fileList{iFile}(1:length(match)+1), [match '_']) || strcmp(fileList{iFile}(1:length(match)+1), [match '.']))
			m = [m iFile];
		end
	end
	if isempty(m)
		warning(['No file beginning with ''' match ''' was found.'])
		filename = [];
		return
	end
	
	% From those matches lets take the one selected by its index
	if ~isempty(n)
		if n <= numel(m)
			m = m(n);
		else
			error('A file index was provided that is higher than the number of files matched to the provided string')
		end
	end
	% .. or return a single match or all matches
	if numel(m) == 1
		fileList = fileList{m};
		disp(['Selecting file ' fileList '.'])
	elseif numel(m) > 1 
		if returnall
			fileList = fileList(m);
			disp(['Selecting ' num2str(numel(m)) ' files.'])
		else
			error('More than one matching string and no disambiguating file index provided. Cannot compute.')
		end
	else
		fileList = fileList{m};
		disp(['Selecting file ' fileList '.'])
	end
elseif ~isempty(n) % if no match but a file index is provided
	fileList = fileList{n};
	disp(['Selecting file ' fileList '.'])
end

% Filename or full path?
if full
	filename = fullfile(path, fileList);
else
	filename = fileList;
end


