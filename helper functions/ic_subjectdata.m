function subjectdata = ic_subjectdata()
% Returns subjectdata in consecutive notation.
%
% OUTPUT VARIABLE:
% subjectdata       Structure array containing the subject data in
%                   consecutive order. This order can change, e.g. if
%                   subjects are excluded. Save all filenames based on the
%                   field .id which is unique. Thats more important than
%                   forwarding the .id field at every analysis step.
%                   See also get_filenames.

% RC	Incomplete Reactivation (Sound + Syllable)
% RW	Complete Reactivation (Sound + Full Word)
% NR	No Reactivation

i = 1;
subjectdata(i).id = 101;
subjectdata(i).id_name = 'RC101';
subjectdata(i).eeg = 'ExpE01-S11.eeg';
subjectdata(i).hypno = 'E01S11.txt';
subjectdata(i).condition = 1;
subjectdata(i).experiment = 1;

% .... list further subjects here