function varargout = slack_management_utils(action, varargin)
%SLACK_MANAGEMENT_UTILS Utility functions for enhanced slack management in SCP
%
% USAGE:
%   w_slack = slack_management_utils('compute_slack_weight', t_remaining, T_horizon, P)
%   tol_slack = slack_management_utils('compute_slack_tolerance', t_remaining, T_horizon, P)
%   s_max = slack_management_utils('compute_slack_bounds', t_remaining, T_horizon, P, retry_level)

switch lower(action)
    case 'compute_slack_weight'
        % Compute time-dependent slack weight
        % INPUTS: t_remaining, T_horizon, P
        % OUTPUT: w_slack
        t_remaining = varargin{1};
        T_horizon = varargin{2};
        P = varargin{3};
        
        % Time progress ratio (0 = start of horizon, 1 = end of horizon)
        time_progress = max(0, min(1, (T_horizon - t_remaining) / T_horizon));
        
        % Progressive slack weight increase using power law
        % w_slack(t) = w_initial + (w_terminal - w_initial) * (progress^alpha)
        weight_range = P.w_slack_terminal - P.w_slack_initial;
        progress_factor = time_progress ^ P.slack_decay_alpha;
        w_slack = P.w_slack_initial + weight_range * progress_factor;
        
        varargout{1} = w_slack;
        
    case 'compute_slack_tolerance'
        % Compute time-dependent slack tolerance for convergence
        % INPUTS: t_remaining, T_horizon, P
        % OUTPUT: tol_slack
        t_remaining = varargin{1};
        T_horizon = varargin{2};
        P = varargin{3};
        
        % Time progress ratio
        time_progress = max(0, min(1, (T_horizon - t_remaining) / T_horizon));
        
        % Progressive slack tolerance tightening
        tol_range = P.tol_slack_initial - P.tol_slack_terminal;
        progress_factor = time_progress ^ P.slack_decay_beta;
        tol_slack = P.tol_slack_initial - tol_range * progress_factor;
        
        varargout{1} = tol_slack;
        
    case 'compute_slack_bounds'
        % Compute time-dependent maximum slack bounds
        % INPUTS: t_remaining, T_horizon, P, retry_level (optional)
        % OUTPUT: s_max_vector [s_vx_max; s_vy_max; s_omega_max]
        t_remaining = varargin{1};
        T_horizon = varargin{2};
        P = varargin{3};
        retry_level = 0;
        if nargin >= 5
            retry_level = varargin{4};
        end
        
        % Time progress ratio
        time_progress = max(0, min(1, (T_horizon - t_remaining) / T_horizon));
        
        % Base slack bounds decrease with time
        bound_range = P.slack_max_initial - P.slack_max_terminal;
        progress_factor = time_progress ^ P.slack_decay_beta;
        s_max_base = P.slack_max_initial - bound_range * progress_factor;
        
        % Apply retry multiplier if needed
        s_max = s_max_base * (P.slack_retry_multiplier ^ retry_level);
        
        % Return as vector [s_vx_max; s_vy_max; s_omega_max]
        % Scale angular slack by reasonable factor (velocity units vs angular units)
        omega_scale = 10; % rad/s vs m/s scaling
        s_max_vector = [s_max; s_max; s_max * omega_scale];
        
        varargout{1} = s_max_vector;
        
    case 'log_slack_status'
        % Log current slack management status
        % INPUTS: t_remaining, T_horizon, P, slack_norm, retry_level (optional)
        t_remaining = varargin{1};
        T_horizon = varargin{2};
        P = varargin{3};
        slack_norm = varargin{4};
        retry_level = 0;
        if nargin >= 6
            retry_level = varargin{5};
        end
        
        % Compute current parameters
        w_slack = slack_management_utils('compute_slack_weight', t_remaining, T_horizon, P);
        tol_slack = slack_management_utils('compute_slack_tolerance', t_remaining, T_horizon, P);
        s_max = slack_management_utils('compute_slack_bounds', t_remaining, T_horizon, P, retry_level);
        
        % Time progress
        time_progress = max(0, min(1, (T_horizon - t_remaining) / T_horizon));
        
        if retry_level > 0
            fprintf('    SLACK RETRY %d: t_rem=%.2fs, prog=%.1f%%, w_slack=%.1e, tol=%.1e, s_max=%.1f, slack=%.1e\n', ...
                retry_level, t_remaining, time_progress*100, w_slack, tol_slack, s_max(1), slack_norm);
        else
            fprintf('    Slack status: t_rem=%.2fs, prog=%.1f%%, w_slack=%.1e, tol=%.1e, s_max=%.1f, slack=%.1e\n', ...
                t_remaining, time_progress*100, w_slack, tol_slack, s_max(1), slack_norm);
        end
        
    case 'validate_parameters'
        % Validate slack management parameters
        % INPUT: P
        % OUTPUT: is_valid, error_msg
        P = varargin{1};
        is_valid = true;
        error_msg = '';
        
        % Check required fields
        required_fields = {'w_slack_initial', 'w_slack_terminal', 'slack_decay_alpha', ...
                          'slack_decay_beta', 'slack_max_initial', 'slack_max_terminal', ...
                          'slack_retry_multiplier', 'max_slack_retries', ...
                          'tol_slack_initial', 'tol_slack_terminal'};
        
        for i = 1:length(required_fields)
            if ~isfield(P, required_fields{i})
                is_valid = false;
                error_msg = sprintf('Missing required field: %s', required_fields{i});
                break;
            end
        end
        
        % Check parameter ranges
        if is_valid
            if P.w_slack_initial >= P.w_slack_terminal
                is_valid = false;
                error_msg = 'w_slack_initial must be < w_slack_terminal';
            elseif P.slack_max_initial <= P.slack_max_terminal
                is_valid = false;
                error_msg = 'slack_max_initial must be > slack_max_terminal';
            elseif P.tol_slack_initial <= P.tol_slack_terminal
                is_valid = false;
                error_msg = 'tol_slack_initial must be > tol_slack_terminal';
            elseif P.slack_decay_alpha <= 0 || P.slack_decay_beta <= 0
                is_valid = false;
                error_msg = 'slack_decay_alpha and slack_decay_beta must be positive';
            end
        end
        
        varargout{1} = is_valid;
        varargout{2} = error_msg;
        
    otherwise
        error('Unknown action: %s', action);
end
end