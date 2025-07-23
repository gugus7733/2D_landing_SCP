function varargout = control_replay_utils(action, varargin)
%CONTROL_REPLAY_UTILS Utility functions for the control replay system
%
% USAGE:
%   control_log = control_replay_utils('initialize')
%   control_replay_utils('log_control', control_log, sim_time, T_applied, delta_applied)
%   replay_log = control_replay_utils('load_replay', filename)
%   [T_applied, delta_applied] = control_replay_utils('get_replay_control', replay_log, sim_time)
%   control_replay_utils('save_log', control_log, filename)

switch lower(action)
    case 'initialize'
        % Initialize empty control log structure
        varargout{1} = struct('t', [], 'U', []);
        
    case 'log_control'
        % Log a control input at a given time
        % INPUTS: control_log, sim_time, T_applied, delta_applied
        control_log = varargin{1};
        sim_time = varargin{2};
        T_applied = varargin{3};
        delta_applied = varargin{4};
        
        if isempty(control_log.t)
            control_log.t = sim_time;
            control_log.U = [T_applied; delta_applied];
        else
            control_log.t(end+1) = sim_time;
            control_log.U(:,end+1) = [T_applied; delta_applied];
        end
        varargout{1} = control_log;
        
    case 'load_replay'
        % Load control log from file for replay
        % INPUTS: filename
        filename = varargin{1};
        
        if exist(filename, 'file')
            replay_data = load(filename);
            varargout{1} = replay_data.control_log;
        else
            varargout{1} = [];
        end
        
    case 'get_replay_control'
        % Get interpolated control from replay log at given time
        % INPUTS: replay_log, sim_time
        % OUTPUTS: T_applied, delta_applied
        replay_log = varargin{1};
        sim_time = varargin{2};
        
        if ~isempty(replay_log) && length(replay_log.t) > 1
            T_applied = interp1(replay_log.t, replay_log.U(1,:), sim_time, 'linear', replay_log.U(1,end));
            delta_applied = interp1(replay_log.t, replay_log.U(2,:), sim_time, 'linear', replay_log.U(2,end));
        else
            T_applied = 0;
            delta_applied = 0;
        end
        varargout{1} = T_applied;
        varargout{2} = delta_applied;
        
    case 'save_log'
        % Save control log to file
        % INPUTS: control_log, filename
        control_log = varargin{1};
        filename = varargin{2};
        
        save(filename, 'control_log');
        
    case 'validate_replay_coverage'
        % Check if replay log covers the requested time range
        % INPUTS: replay_log, t_replay_control
        % OUTPUTS: is_valid, max_time_available
        replay_log = varargin{1};
        t_replay_control = varargin{2};
        
        if isempty(replay_log) || isempty(replay_log.t)
            is_valid = false;
            max_time_available = 0;
        else
            max_time_available = max(replay_log.t);
            is_valid = max_time_available >= t_replay_control;
        end
        varargout{1} = is_valid;
        varargout{2} = max_time_available;
        
    otherwise
        error('Unknown action: %s', action);
end
end
