function error = error_build_optimization(simulation, observation, mcmc_state, error)


if mcmc_state
    error = -0.5 * sum( rmmissing( ( (simulation - observation)...
        ./ (error / 100 * observation ) ).^2 ), 'all' ) / numel(observation);
else
    error_temp = (simulation - observation) ./ observation;
    error = rms(rmmissing( error_temp(not(isinf(error_temp))) ),'all');
end