function realsave(filename, data, ask)
% Proper save function
% Adds the obligatory '-v7.3' statement.
% Can be used within parfor loops.
% Input argument 'ask' can be set to 1, in which case you are asked first
% before a file is overwritten.
%
% realsave('/root/yourfile.mat', data)
% realsave('/root/yourfile.mat', data, 1)
%
% AUTHOR:
% Jens Klinzing, jens.klinzing@uni-tuebingen.de

%%  SETUP
if nargin < 2
	error('Unexpected number of input arguments.')
end
if nargin < 3
	ask = 0;
end
if ~ischar(filename)
	error('Unexpected first unput argument')
end

%% START
if exist(filename, 'file') || exist([filename '.mat'], 'file')
	[~, name, ext] = fileparts(filename);
	if ask == 0
		disp([name, ext, ' does already exist and will be overwritten.'])
		save(filename, 'data', '-v7.3')
	else
		a=input(['File ''' name, ext, ''' already exist. Do you want to overwrite? Y or N '], 's');
		switch lower(a)
			case 'y'
				disp('File will be overwritten.')
				save(filename, 'data', '-v7.3')
			case 'n'
				disp('File will NOT be overwritten. No data saved!')
		end
	end
else
	save(filename, 'data', '-v7.3')
end

end