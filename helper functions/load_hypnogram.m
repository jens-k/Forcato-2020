function output = load_hypnogram(filename)
%IMPORTFILE Import numeric data from a text file as column vectors.
%   [VARNAME1,VARNAME2] = IMPORTFILE(FILENAME) Reads data from text file
%   FILENAME for the default selection.
%
%   [VARNAME1,VARNAME2] = IMPORTFILE(FILENAME, STARTROW, ENDROW) Reads data
%   from rows STARTROW through ENDROW of text file FILENAME.
%
% Example:
%   [VarName1,VarName2] = importfile('vp01.txt',1, 161);
%
%    See also TEXTSCAN.
%
% MATLAB-generated, slightly changed by Jens Klinzing

%% Initialize variables.
startRow = 1;
endRow = inf;


%% Format string for each line of text:
%   column1: double (%f)
%	column2: double (%f)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%1f%f%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, endRow(1)-startRow(1)+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow(1)-1, 'ReturnOnError', false);
for block=2:length(startRow)
    frewind(fileID);
    dataArrayBlock = textscan(fileID, formatSpec, endRow(block)-startRow(block)+1, 'Delimiter', '', 'WhiteSpace', '', 'HeaderLines', startRow(block)-1, 'ReturnOnError', false);
    for col=1:length(dataArray)
        dataArray{col} = [dataArray{col};dataArrayBlock{col}];
    end
end

output(:,1)         = dataArray{1,1};
output(:,2)         = dataArray{1,2};

%% Close the text file.
fclose(fileID);



