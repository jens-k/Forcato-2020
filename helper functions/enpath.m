function path = enpath(path)
% Ensures that given path exists. If it doesnt exist its created.

if strcmp(path(1), '$')
    warning('You tried to use a relative path. I corrected that for you, but don''t do that.')
    path = abpath(path);
end

if ~exist(path,'dir')
    mkdir(path);
    fprintf('Created path %s. \n', path)
end