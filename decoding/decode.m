
function decode()
% Decoding lab 
% Loads data frm random dot experiment (monkey is trained to release a lever when random dot stimuli move 
% coherently) in 'decodingLabData.mat', which gives two matrices of neural responses where each row is a trial
% and each column is 1 ms (1 if there was aspike at that ms, 0 if not). The variable responseTime is a vector
% of the lever release time relative to the onset of the motion stimulus (NaN if there was no response).

% This code performs an ROC analysis of neurometric performance and neuronal activity correlation with motion 
% stimulus detection.

% Wtitten by Thomas Christinck (independently)

% Load data and get number of trials
load('decodingLabData.mat')
nb_trials = length(responseTime);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine how well each neuron informs an ideal observer that the motion stimulus has occurred on a trial-by-trial
% basis by comparing responses before stimulus onset with responses 40-140 ms after (40 ms because this is typical
% response time for MT neuron).

% Do this for both neurons 1 & 2
nonstim_neuron1 = neuron1(:, 1:100);
stim_neuron1 = neuron1(:, 541:640);
nonstim_neuron2 = neuron2(:, 1:100);
stim_neuron2 = neuron2(:, 541:640);

% Gives the number of spikes per second per trial 
trial_total_n1_stim = round(sum(stim_neuron1, 2)) .* 10 ;
trial_total_n1_nonstim = round(sum(nonstim_neuron1, 2)) .* 10 ;
trial_total_n2_stim = round(sum(stim_neuron2, 2)) .* 10 ;
trial_total_n2_nonstim = round(sum(nonstim_neuron2, 2)) .* 10 ;

% We want to plot a histogram of the number of spikes durng each 100 ms window
figure(1)
histogram(trial_total_n1_nonstim, 'BinWidth', 10); hold on; histogram(trial_total_n1_stim, 'BinWidth', 10);
xlabel('Frequency (Hz)')
ylabel('Number of Trials')
title('Neuron 1 Spike Rate During 100 ms Window Before and After Stimulus Onset')
legend('Pre-stimulus', 'Post-stimulus')

figure(2)
histogram(trial_total_n2_nonstim, 'BinWidth', 10); hold on; histogram(trial_total_n2_stim, 'BinWidth', 10);
xlabel('Frequency (Hz)')
ylabel('Number of Trials')
title('Neuron 2 Spike Rate During 100 ms Window Before and After Stimulus Onset')
legend('Pre-stimulus', 'Post-stimulus')

% Now generate ROC curves
thresh = [0 10 20 30 40 50 60];
[auc_neuron1, tp_scores_n1, fp_scores_n1] = roc(trial_total_n1_stim, trial_total_n1_nonstim, thresh, nb_trials, nb_trials);
[auc_neuron2, tp_scores_n2, fp_scores_n2] = roc(trial_total_n2_stim, trial_total_n2_nonstim, thresh, nb_trials, nb_trials);

fprintf("Neuron 1 neurometric area under ROC curve is %.2f \n", auc_neuron1);
fprintf("Neuron 2 neurometric area under ROC curve is %.2f \n", auc_neuron2);

figure(3)
plot(fp_scores_n1, tp_scores_n1, '-o', 'Linewidth', 1.5);
hold on;
plot(fp_scores_n2, tp_scores_n2, '-o', 'Linewidth', 1.5);
plot([0 1], [0 1], '--')
xlabel('False Positives')
ylabel('True Positives')
title('Neurometric ROC')
legend('Neuron 1', 'Neuron 2')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Part 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Separate spike counts during 40-140 ms after stimulus presentation into cases where the monkey responded and cases where they
% did not respond.
n1_responded = trial_total_n1_stim(~isnan(responseTime));
n1_noresponded = trial_total_n1_stim(isnan(responseTime));
n2_responded = trial_total_n2_stim(~isnan(responseTime));
n2_noresponded = trial_total_n2_stim(isnan(responseTime));

% Plot two histograms of the number of spikes occurring during the 100 ms window corresponding to correct and failed trials for each
% neuron
figure(4)
histogram(n1_responded, 'BinWidth', 10); hold on; histogram(n1_noresponded, 'BinWidth', 10);
xlabel('Frequency (Hz)')
ylabel('Number of Trials')
title('Neuron 1 Spike Rate During 100 ms Window After Stimulus Onset')
legend('Responded', 'Did not respond')

figure(5)
histogram(n2_responded, 'BinWidth', 10); hold on; histogram(n2_noresponded, 'BinWidth', 10);
xlabel('Frequency (Hz)')
ylabel('Number of Trials')
title('Neuron 2 Spike Rate During 100 ms Window After Stimulus Onset')
legend('Responded', 'Did not respond')

% Get the number of trials where we saw response and no response, and then generate ROC curves
nb_response_trials = length(n1_responded);
nb_noresponse_trials = length(n1_noresponded);
[auc_neuron1res, tp_scores_n1, fp_scores_n1] = roc(n1_responded, n1_noresponded, thresh, nb_response_trials, nb_noresponse_trials);
[auc_neuron2res, tp_scores_n2, fp_scores_n2] = roc(n2_responded, n2_noresponded, thresh, nb_response_trials, nb_noresponse_trials);

fprintf("Neuron 1 activity-motion stimulus detection area under ROC curve is %.2f \n", auc_neuron1res);
fprintf("Neuron 2 activity-motion stimulus detection area under ROC curve is %.2f \n", auc_neuron2res);

figure(6)
plot(fp_scores_n1, tp_scores_n1,'-o', 'Linewidth', 1.5); 
hold on;
plot(fp_scores_n2, tp_scores_n2,'-o', 'Linewidth', 1.5); 
plot([0 1], [0 1], '--')
xlabel('False Positives')
ylabel('True Positives')
title('ROC Detect Probability')
legend('Neuron 1', 'Neuron 2')

function [auc, tp_scores, fp_scores]=roc(spikecount1, spikecount2, thresh, nb_trials1, nb_trials2)
% Function for computing roc curve and auc. Iterates through thresholds in thresh.

nb_thresh = length(thresh);
tp_scores(1:nb_thresh) = 0;
fp_scores(1:nb_thresh) = 0;
for i = 1:nb_thresh
	threshed_sc1 = spikecount1 >= thresh(i);
	threshed_sc2 = spikecount2 >= thresh(i);
	total_sc1 = sum(threshed_sc1);
	total_sc2 = sum(threshed_sc2);
	tp_scores(i) = total_sc1 / nb_trials1;
	fp_scores(i) = total_sc2 / nb_trials2;
end

% Area under the curve will be negative because tp/fp scores are ordered backwards; multiply by -1 to make positive
auc = -1 * trapz(fp_scores, tp_scores);

return;
