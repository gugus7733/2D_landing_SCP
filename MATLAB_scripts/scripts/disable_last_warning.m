% Get the last warning's ID
[~, id] = lastwarn;

if ~isempty(id)
    % Try to load the existing list
    filename = 'warnings_to_disable.mat';
    if exist(filename, 'file')
        S = load(filename);
        if isfield(S, 'warning_ids')
            warning_ids = S.warning_ids;
        else
            warning_ids = {};
        end
    else
        warning_ids = {};
    end
    
    % Add new id if not already in the list
    if ~any(strcmp(id, warning_ids))
        warning_ids{end+1} = id; %#ok<AGROW>
        save(filename, 'warning_ids');
        disp(['Warning ID "' id '" disabled and added to ' filename '.']);
    else
        disp(['Warning ID "' id '" was already in ' filename '.']);
    end

    % Disable the warning now
    warning('off', id);
else
    disp('No warning ID to disable.');
end