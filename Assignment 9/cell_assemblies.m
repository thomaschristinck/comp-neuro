% Code written for NEUR 603 Assignment 9 for Dr Peyrache


% Part 1 - correlation matrix of PFC neurons during running, and pre- and post-sleep.
load DataPFC.mat % load data

%Compute correlation matrices
Hpre_corr = corr(Qpre);
Hrun_corr  = corr(Qrun);
Hpost_corr  = corr(Qpost);

% Zero-out diagonals for visualization
Hpre = Hpre_corr  - diag(diag(Hpre_corr ));
Hrun = Hrun_corr  - diag(diag(Hrun_corr ));
Hpost = Hpost_corr  - diag(diag(Hpost_corr ));

% Display
figure(1); imagesc(Hpre); colorbar;
title("Pre-sleep Correlation Matrix for PFC Neurons")

figure(2); imagesc(Hrun); colorbar;
title("Running Correlation Matrix for PFC Neurons")

figure(3); imagesc(Hpost); colorbar;
title("Post-sleep Correlation Matrix for PFC Neurons")

% Part 2 - correlation of neurons 20 and 26 during running, and pre- and post-sleep.
neuron_1 = 20;
neuron_2 = 26;
[Hpre, b] = xcorr (Qpre(:,neuron_1) , Qpre(:,neuron_2) , 100 , 'coeff');
Hrun = xcorr (Qrun(:,neuron_1) , Qrun(:,neuron_2) ,  100 , 'coeff');
Hpost = xcorr (Qpost(:,neuron_1) , Qpost(:,neuron_2) ,  100 , 'coeff');
t = b * 0.100; %Time bins are 100 ms

figure(4),clf
plot(t,Hpre,'b')
hold on, plot(t,Hpost,'r')
hold on,plot(t,Hrun,'k')
legend('Sleep PRE','Sleep POST','Running')
xlabel('Time-lag (s)')
ylabel('Correlation')
title('Cross-correlation between neurons 20 and 26 across states')

% Part 3 - PCA during running. Compute the PC weights, and the score associated with the
% first principle component. Spike trains are whitened (0 mean, stdev of 1) before they are 
% projected onto PCs.
PCweights_all = pcacov(Hrun_corr);
PCweights = PCweights_all(:,1);
PCscores = zscore(Qrun)*PCweights;

% Part 4 - Plot mean firing rate against the 1st principle component weights
runFR = mean(Qrun);
figure(5),clf
hold on, plot(runFR,PCweights,'o')
xlabel('Firing Rates (Hz)')
ylabel('Weights in PC 1')
title('Neuron firing rates and first principle component weights')

% Part 5  - Compute scores of first principle components during sleep
scorePRE = zscore(Qpre)*PCweights;
scorePOST = zscore(Qpost)*PCweights;

singleNpre = (zscore(Qpre).^2)*(PCweights.^2);
singleNpost = (zscore(Qpost).^2)*(PCweights.^2);

% Part 6 - Compute reactivation scores; plot reactivation strength during sleep-pre and post
reactPRE = scorePRE.^2 - singleNpre;
reactPOST = scorePOST.^2 - singleNpost;

figure(6); bar(reactPRE,1.0); title('Pre');
xlabel('Time (s)')
ylabel('Reactivation Strength')

figure(7); bar(reactPOST, 1.0); title('Post');
xlabel('Time (s)')
ylabel('Reactivation Strength')