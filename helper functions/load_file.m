function data = load_file(path, varargin)
% Returns a structure containing the fieldtrip-processed dataset found at
% the path given in path. Call this function either with a direct path to a
% file or with the path to a folder and the desired file number.
%
% Use as
% data = load_file(path, subjNo)
% eg. data = load_file('C:\data\file.mat')
%     data = load_file('C:\data\', 9)
%
% INPUT VARIABLES:
% path              String; path of a) input file or b) input folder
%                   containing files with fieldtrip data,
%                   though it can contain any number of folders. If it is a
%                   input folder, the second input argument is needed.
% fileNo            int; optional, number of desired file (folders are
%                   ignored in the count)
%
% OUTPUT VARIABLES:
% data              data taken from specified file/folder
%
% AUTHOR:
% Jens Klinzing, jens.klinzing@uni-tuebingen.de

switch nargin
    case 0
        error('At least one input argument needed.');
    case 1
        filename = path;
        [~,name,ext]   = fileparts(filename);
        disp(['Loading file ' name ext '.'])
    otherwise
        filename = get_filenames(path, varargin{:}, 'full');
end

if ~isempty(filename)
	temp = load(filename);
else
	error('No file found.')
end

names = fieldnames(temp);
if length(names) ~= 1
    error('Unexpected content in preprocessed data file. File should contain one single structure.')
end
data = temp.(names{1});
clear temp

end