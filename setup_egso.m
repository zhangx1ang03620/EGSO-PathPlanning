function root = setup_egso()
% SETUP_EGSO  Add the EGSO code folders to the MATLAB path.
%
%   setup_egso            % run once per MATLAB session, from any folder
%   root = setup_egso();  % also returns the repository root
%
% Experiment scripts write to results/... relative to the current folder,
% so run them from the repository root (cd(setup_egso())).
%
% experiments/equal_nfe is deliberately NOT added: it contains a shadowing
% fitness.m that counts objective evaluations, and its scripts cd into that
% folder themselves.

root = fileparts(mfilename('fullpath'));
folders = {'egso', 'egso/variants', 'egso3d', ...
           'baselines/metaheuristics', 'baselines/planners', ...
           'fields', 'maps', 'stats', 'experiments', 'examples', 'tests'};
for i = 1:numel(folders)
    addpath(fullfile(root, folders{i}));
end
if nargout == 0
    fprintf('EGSO: %d folders added to the MATLAB path (root: %s)\n', numel(folders), root);
    clear root
end
end
