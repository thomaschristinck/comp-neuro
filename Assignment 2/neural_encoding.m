function neural_encoding
    % Load the data
    load('spiketrain1.mat');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PART 1 A)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf("------------------------- Part 1 A -------------------------\n");
    
    % Get the Vm sampling period and frequency; other relevant Vm stats (the number of timespoints)
    sample_period = timeaxis(2) - timeaxis(1);
    sample_freq = 1 / sample_period;
    fprintf("Sampling frequency is %.2f Hz \n", sample_freq);
    Vm_size = size(Vm);
    nb_timepoints = Vm_size(2);

    % Set a threshold for AP detection
    thresh = -60;

    % Get indices where an AP is firing; then with 'find' we return an array that is 1 at indices 
    % where an AP is firing and 0 elsewhere
    for i = 1:nb_timepoints
    	ap_timepoints(i) = Vm(i) > thresh && Vm(i -1) < thresh;
    end
    ap_indices = find(ap_timepoints);

    % Now, the spiketimes are given by the sample_period times the action potential indices.
    spiketimes = sample_period * ap_indices;

    % Get the number of spiketimes (we'll need to look through them now)
    size_spiketimes = size(spiketimes);
    nb_spiketimes = size_spiketimes(2);


    % For each time step (1 msec, 0.5 msec, 0.1 msec), we'll make a binary representation of the
    % spike train

    % First, for each timestep we get the number of bins associated, and then 
    timesteps = [0.001 0.0005 0.0001];
    nb_bins(1:3) = 0;

    for i = 1:3
        nb_bins(i) = floor(max(timeaxis) / timesteps(i));
    	binned_spiketimes{i} = round(spiketimes/timesteps(i));
    	bins{i}(1:nb_bins(i)) = 0;
    	for j = 1:nb_spiketimes
    		bins{i}(binned_spiketimes{i}(j)) = bins{i}(binned_spiketimes{i}(j)) + 1;
    	end
    end

    % Now, each bin is divided by its binwidth
    for i = 1:3
    	bins{i} = bins{i} / timesteps(i);
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PART 1 B)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf("\n------------------------- Part 1 B -------------------------\n");

    % Compose a histogram of interarrival times with bin size 1 ms on t E (0, 200 ms)
    interarrival_times = diff(spiketimes);
    hist_bins = linspace(0.0005,0.30005,300);
    fprintf("Figure 1 ...\n");
    figure(1)
    hist(interarrival_times, hist_bins)
    xlim([0 0.2])
    title('Interspike Interval Histogram')
    xlabel('Interspike Arrival Bin (ms)') 
    ylabel('Frequency of Interspike Interval Occurance (count per bin)')

    % Compute coefficient of variation (CV) for the neuron
    S = std(interarrival_times);
    M = mean(interarrival_times);
    CV = S / M;

    fprintf("Coefficient of variation is %.3f \n", CV);

    % Plot the interspike interval correlation coefficients
    fprintf("Figure 2 ...\n");
    figure(2)
    [ii_corr,ii_lags]=xcorr(interarrival_times - mean(interarrival_times), 200,'coeff');
    plot(ii_lags, ii_corr,'LineWidth', 1.5)
    title("Autocorrelation of Interspike Intervals")
    xlim([0 15])
    ylim([-0.4 1])
    ylabel('Autocorrelation')
    xlabel('j (Lag)')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PART 1 C)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf("\n------------------------- Part 1 C -------------------------\n");

    % Compute the autocorrelation of the binarized representations from A.
    max_lags = [35 70 350];
    colors = [[0, 0.4470, 0.7410]; [0.8500, 0.3250, 0.0980]; [0.9290, 0.6940, 0.1250]];
    fprintf("Figure 3 ...\n");
    figure(3)
    for i = 1:3
    	[c{i},lags{i}]=xcorr(bins{i} - mean(bins{i}), max_lags(i),'coeff');
    	 % Plot time course
        subplot(3,1,i);

        plot(lags{i} * timesteps(i), c{i}, 'LineWidth', 1.5, 'Color', colors(i, :));
    	xlim([0 0.035])

    	if i == 1
        	title("Autocorrelation of Binary Spike-Train Representations")
        	leg = legend('1 ms', 'Location', 'E');
    		title(leg,'Bin Size (dt)')
            hold on
    	end
        if i == 2
        	ylabel('Autocorrelation')
        	legend('0.5 ms', 'Location', 'E');
        end
        if i == 3
        	legend('0.1 ms', 'Location', 'E');
        	xlabel('\tau (ms)')
        end
    end
    drawnow();

    % Compute the power spectra of the binarized representations from A
    fprintf("Figure 4 ...\n");
    figure(4)
    sampling_freqs = [1000 2000 10000];
    for i = 1:3
    	[pxx{i},f{i}]=pwelch(bins{i},bartlett(2048),1024,2048,sampling_freqs(i));
    	 % Plot time course
        subplot(3,1,i);
        plot(f{i}, pxx{i}, 'LineWidth', 1.5, 'Color', colors(i, :));

    	if i == 1
        	title("Power Spectra of Binary Spike-Train Representations")
        	leg = legend('1 ms', 'Location', 'E');
    		title(leg,'Bin Size (dt)')
            hold on
    	end
        if i == 2
        	ylabel('Power Spectral Density')
        	legend('0.5 ms', 'Location', 'E');
        end
        if i == 3
        	legend('0.1 ms', 'Location', 'E');
        	xlabel('f (Hz)')
        end
    end
    drawnow();

    % Clear old data and load data for part 2
    clear 
    load('spiketrain2.mat');

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PART 2 A)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf("\n------------------------- Part 2 A -------------------------\n");

    % Plot the data in raster format
    fprintf("Figure 5 ...\n");
    figure(5)
    scatter(data(:,2), data(:,1), 's', 'filled', 'MarkerFaceColor',[0 .7 .7])
    xlim([0 500])
    ylim([0.5 40.5])
    title('Spike Times in First 0.5 s for Trials 1-40')
    xlabel('Time (s)')
    ylabel('Trial Number')
    ax = gca;
    ax.YGrid = 'off';
    ax.XGrid = 'on';

    % Now, there are 100 trials with 50 seconds of data. The 50 seconds has to be divided into 1 msec
    % bins, and we count the number of spikes in each bin.
    bin_width = 0.001; %ms
    nb_trials = max(data(:,1));
    trial_length = max(data(:,2));
    nb_bins = trial_length; % in ms
    edges = linspace(0,trial_length,nb_bins + 1);

    % Now, we discretize the spiketimes into bins and then the occurance of each index in bin_indices denotes
    % the occurance of a spike
    bin_indices = discretize(data(:,2), edges);
    spike_bins(1:50000) = 0;
    for i = 1:size(bin_indices)
    	spike_bins(bin_indices(i)) = spike_bins(bin_indices(i)) + 1;
    end

    % Get firing rate by dividing the # of spikes per bin by the number of epochs and by the bin width
    spike_rate = spike_bins * 1/nb_trials * 1/bin_width;
    fprintf("Figure 6 ...\n");
    figure(6)
    bar(spike_rate)
    xlim([0 500])
    title('Peristimulus Time Histogram')
    xlabel('Time (ms)')
    ylabel('Spikes per second')


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PART 2 B)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf("\n------------------------- Part 2 B -------------------------\n");

    % Make a binary representation of each trial at 2 kHz.
    sample_freq = 2000;
    sample_period = 1/ sample_freq * 1000; % ms 

    % Get the spike times for each trial
    spike_times = data(:, 2);
    trial_nbs = 100;
    for i = 1:trial_nbs
        trial_spike_times{i} = spike_times(data(:,1) == i);
    end

    % Put spike times into bins (number of bins will be same for each trial)
    fprintf('Binning spiketimes for 100 trials....\n');
    nb_bins = round(max(data(:,2)) / sample_period);
    for i = 1:100
        binned_spiketimes{i} = round(trial_spike_times{i}/sample_period);
        bins{i}(1:nb_bins) = 0;
        size_spiketimes = size(trial_spike_times{i});
        nb_spiketimes = size_spiketimes(1);
        for j = 1:nb_spiketimes
            bins{i}(binned_spiketimes{i}(j)) = bins{i}(binned_spiketimes{i}(j)) + 1;
        end
    end

    % Compute cross-correlations for all trials
    tau_range = 200;
    x_corr_totals(1:nb_trials, 1:tau_range * 2 + 1) = 0;
    for i = 1:nb_trials
        [x_corr_temp,lags_trials]=xcorr(bins{i}, stim, tau_range,'unbiased');
        x_corr_totals(i, :) = x_corr_temp(1, :);
        clear x_corr_temp
    end

    % Get average cross-correlation function and plot trial 1, 20 and average cross-corr
    avg_xcorr = mean(x_corr_totals);
    fprintf("Figure 7 ...\n");
    figure(7)
    plot(lags_trials * sample_period, x_corr_totals(1, :), 'LineWidth', 1.5)
    hold on
    plot(lags_trials * sample_period, x_corr_totals(20, :),'LineWidth', 1.5)
    hold on 
    plot(lags_trials * sample_period, avg_xcorr, 'LineWidth', 1.5)
    legend('Trial 1', 'Trial 20', 'Average across all trials', 'Location', 'NE');
    title('Stimulus-Response Cross-Correlation')
    xlabel('\tau (ms)')
    ylabel('Cross-Correlation')

    % Now we want to compute the stimulus-response cross-spectrum for trial/epoch 1.
    [pxy,f]=cpsd(bins{1}, stim,bartlett(2048),1024,2048,2000);
    [p_stim,f]=pwelch(stim,bartlett(2048),1024,2048,2000);
    
    % Take complex conjugate of cross-spectrum and divide by stimulus power to obtain
    % transfer function. Then normalize by gain at f = 1
    c_conjugate = conj(pxy);
    transfer_fn = c_conjugate ./ p_stim;
    
    % Normalize by gain at f = 1
    f_target = 1;
    [epsilon,index] = min(abs(f(:)-f_target));
    transfer_fn = transfer_fn / abs(transfer_fn(index));

    % Plot transfer function gain and phase
    fprintf("Figure 8 ...\n");
    figure(8)
    loglog(f, abs(transfer_fn),'LineWidth', 1.5)
    title('Transfer Function Normalized Gain')
    xlabel('Frequency (Hz)')
    ylabel('Gain')
    xlim([0 30])

    fprintf("Figure 9 ...\n");
    figure(9)
    semilogx(f, atan(imag(transfer_fn)/real(transfer_fn)),'LineWidth', 1.5)
    title('Transfer Function Phase')
    xlabel('Frequency (Hz)')
    ylabel('Phase (rad)')
    xlim([0 30])

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % PART 2 C)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    fprintf("\n------------------------- Part 2 C -------------------------\n");

    % Set up a matrix of the time bins. 100 trials with 100000 time points per trial
    bins_matrix(1:100, 1:100000) = 0;
    for i = 1:nb_trials
        bins_matrix(i, :) = bins{i};
    end
    % Take the mean response across all 100 trials
    mean_response = mean(bins_matrix);

    % Subtract the mean response from each trial; then get the power spectrum density of 
    % this response; sum each of these PSDs. We'll take the average (p_avg_noise) to get the
    % average noise power.
    p_tot(1:1025, 1) = 0;
    for i = 1:nb_trials
        norm_bins{i} = bins{i} - mean_response;
        [p_temp,f]=pwelch(norm_bins{i},bartlett(2048),1024,2048,2000);
        p_tot = p_tot + p_temp;
    end
    
    % Get power of average response; divide this by average noise power to get SNR
    [p_avg_response,f]=pwelch(mean_response,bartlett(2048),1024,2048,2000);
    p_avg_noise = p_tot / nb_trials;
    snr = p_avg_response ./ p_avg_noise;
    fprintf("Figure 10 ...\n");
    figure(10)
    loglog(f, snr,'LineWidth', 1.5)
    title('Signal-to-Noise Ratio')
    xlabel('Frequency (Hz)')
    ylabel('SNR')
end
