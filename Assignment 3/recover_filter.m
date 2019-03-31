function filter=recover_filter()
% A function that feeds a white noise stimulus to a model (oned, written by Dr Christopher Pack);
% based on the model response, this function recovers the linear filter and static nonlinearity 
% described in Ch 2 of Dayan & Abbott
%
% stimulus -> linear filter -> static nonlinearity - > response
%
% Written by Thomas Christinck 01/2019

% Create a white noise vector (range (-1,1), length 1000) and give to the one 
% dimensional filter
white_noise = -1 + 2 * rand(1000,1);
response = oned(white_noise);

% Find cross-correlation; normalize by standard deviation of input squared and
% by experiment length
[x_corr,lags_trials]=xcorr(response, white_noise, 200,'none');
norm_x_corr = x_corr *  (1 / std(white_noise .^2)) * 0.001;

% Plot this normalized filter
figure(1)
plot(lags_trials, norm_x_corr)
title('Linear Filter Mapping Stimulus to Response')
xlabel('\tau')
ylabel('Filter Magnitude')
xlim([1 51])

% Now, get predicted output from this filter
lower_indices = find(lags_trials <= 52);
upper_indices = find(lags_trials > 1);
indices = intersect(lower_indices, upper_indices);

norm_x_corr = norm_x_corr(indices);
c=conv(norm_x_corr,white_noise);

figure(2)
plot(c, 'LineWidth', 1.5); hold on;
plot(response, 'LineWidth', 1.2); 
title('Response to White Noise Stimulus')
legend('Predicted Response (Linear Filter)', 'Observed Response')
ylabel('Response Magnitude')
xlabel('Time')
xlim([0 1050])
ylim([-175 150])

% Now, we want to optimize the static nonlinearity function to best fit the model
% static nonlinearity (in Figure 3). Visually from the plot, we see that the function
% should take the shape of a steep sigmoid with range (0,100).
figure(3)
mae_fun = @(x) sum(abs(response - (100./(1.0 + exp(x(1)*(x(2)-c))))));
mse_fun = @(x) sum((response - (100./(1.0 + exp(x(1)*(x(2)-c))))).^2);
initial_params = [10 1];
[params_mae,fminres_mae] = fminsearch(mae_fun,initial_params);
[params_mse,fminres_mse] = fminsearch(mse_fun,initial_params);
scatter(c, response); hold on;
fplot(@(x) 100./(1.0 + exp(params_mae(1)*(params_mae(2)-x))), [-150 150], 'LineWidth', 1.5); 
fplot(@(x) 100./(1.0 + exp(params_mse(1)*(params_mse(2)-x))), [-150 150], 'LineWidth', 1.5); 
legend('Model Static Nonlinearity', 'Mean Absolute Error Best Fit', 'Mean Square Error Best Fit', 'Location', 'E')
title('Visualization of Static Nonlinearity')
xlabel('Linear Filter Output')
ylabel('Model Output')

% Now we compute mse applying the static non-linearity to the linear filter output
lin_filter_mses(1:10, 1:5) = 0;
mse_nonlin_mses(1:10, 1:5) = 0;
mae_nonlin_mses(1:10, 1:5) = 0;
input_lengths = [100 500 1000 2000 5000];
nb_vector_lengths = length(input_lengths);
for trial = 1:10
	for i = 1:nb_vector_lengths
		white_noise_ips{i} = -1 + 2 * rand(input_lengths(i),1);
		responses{trial, i} = oned(white_noise_ips{i});
		lin_filter_response{trial, i} = conv(norm_x_corr,white_noise_ips{i});
		mae_responses{trial, i} = 100./(1.0 + exp(params_mae(1)*(params_mae(2)-lin_filter_response{trial, i})));
		mse_responses{trial, i} = 100./(1.0 + exp(params_mse(1)*(params_mse(2)-lin_filter_response{trial, i})));
		lin_filter_mses(trial, i) = immse(lin_filter_response{trial, i}, responses{trial, i});
		mse_nonlin_mses(trial, i) = immse(mse_responses{trial, i}, responses{trial, i});
		mae_nonlin_mses(trial, i) = immse(mae_responses{trial, i}, responses{trial, i});
	end
end

% Plot the average mean square error with standard deviation across trials
figure(4)
errorbar(input_lengths, mean(lin_filter_mses), std(lin_filter_mses), '.k'); hold on;
errorbar(input_lengths, mean(mse_nonlin_mses), std(mse_nonlin_mses), '.b'); 
errorbar(input_lengths, mean(mae_nonlin_mses), std(mae_nonlin_mses), '.r');
legend('Linear Filter', 'Mean Squared Error Optimized Nonlinearity', 'Mean Absolute Error Optimized Nonlinearity')
title('Mean Squared Error Between Estimated and Observed Responses')
xlabel('Length of Input Vector')
ylabel('Mean Squared Error')

% We'll also plot the best predicted response (MSE optimized nonlinearity applied to linear filter output)
figure(5)
r = 100./(1.0 + exp(params_mse(1)*(params_mse(2)-c)));
plot(r, 'LineWidth', 1.5); hold on;
plot(response, 'LineWidth', 1.2); 
title('Response to White Noise Stimulus')
legend('Predicted Response (with MSE optimized nonlinearity)', 'Observed Response')
ylabel('Response Magnitude')
xlabel('Time')
xlim([0 1050])
ylim([-5 125])

function r=oned(a)
% One dimensional filtering operation for simulating a reverse correlation
% experiment.  The input a is filtered by the simulated neuron and the
% response is returned in r.
% Christopher Pack, 11/07

k=25;  %Temporal scale factor 
slow_t=temp_imp_resp(5,k,0:.02:1);
fast_t=temp_imp_resp(3,k,0:.02:1);


b=slow_t + fast_t; % linear filter
c=conv(b,a); % convolve filter with stimulus
c=c./max(c); % normalize
c=25.0*c + 0.2;
r=100./(1.0 + exp(10*(0.5-c))); % static nonlinearity
%r=(c>0.0);

return;

function time_response=temp_imp_resp(n,k,t)
%time_response=temp_imp_resp(n,k,t)
%
%Produces a temporal impulse response function using the from from
%figure 1 in Adelson & Bergen (1985)
%
%It's pretty much a difference of Poisson functions with different
%time constants.

time_response=(k*t).^n .* exp(-k*t).*(1/factorial(n)-(k*t).^2/factorial(n+2));
return;