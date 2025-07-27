
eps_set = [-0.15 -0.1 -0.05 0 +0.05 0.1 0.15];
all_costs = [];

% Generate the structures : dirty
if (isempty(last_scp_sol))
    [all_sols, all_logs] = run_scp_2d(current_state, P.time_to_touchdown, round(P.time_to_touchdown / 0.5), P, 0.5, 1, last_scp_sol, sim_time);
else
    all_sols = last_scp_sol;
    all_logs = last_scp_log;
end

all_tf = P.time_to_touchdown*(1+eps_set);
used_tf = [];
parfor i_eps = 1:numel(eps_set)
    tf_i = all_tf(i_eps);
    N_scp = round(tf_i / P.dt_scp);
    [sol_i, log_i] = run_scp_2d(current_state, tf_i, N_scp, P, P.dt_scp, 3, last_scp_sol, sim_time);

    if (~isempty(sol_i))
        used_tf(i_eps) = tf_i;
        all_costs(i_eps) = log_i.cost(end);
        all_sols(i_eps) = sol_i;
        all_logs(i_eps) = log_i;
    end
end

best_cost = inf;
for i_cost = 1:numel(all_costs)
    if  (all_costs(i_cost) < best_cost) && (~isempty(all_sols(i_cost)))
        best_cost = all_costs(i_cost);
        best_sol  = all_sols(i_cost);
        best_log  = all_logs(i_cost);
        best_tf    = used_tf(i_cost);
    end
end