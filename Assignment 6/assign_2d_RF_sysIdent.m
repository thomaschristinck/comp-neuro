% assign_2d_RF_sysIdent.m
%
% simulated estimation of 2d spatial receptive field of visual cortex simple cell
% uses 'scg' (scaled conjugate gradient) optimization, so no need to specify step size
% simple regression with netlab, using a glm, to estimate Gabor-RF (simulated neuron/filter response waveform)
% glm-fit uses ridge regression, with regularization parameter "alpha"
% use training dataset to get estimate of RF map 
%   using this estimated model, evaluates VAF on holdback validationdataset, for a given alpha

%   figure 1:  example stimulus images
%   figure 2:  model filter, and estimate of rfMap
 

% needs to run:
%   netlab - folder of functions from netlab toolbox
%   getStimulusMovies.m - construct movies, or read from files
%   makeModelRF.m - create linear gabor filter, for model
%   imagescZadj.m - imagesc, but force mid-range of color-map to be at zero
%   hwr.m - half wave rectification

clear all;  fprintf(1,'\n\n\n\n\n\n');   close all;  

% Add netlab
addpath('netlab3_3/')
addpath('nethelp3_3/')

rng('default');  % "standard" random number seed -> reproducible simulations 


%  choose stimulus   %%%%%%%%%%%%%%%
commandwindow
stimType = input('stimulus:  white noise (1) or natural images (2):  ');
if stimType==1
    option.stimulus = 'white';
else
    option.stimulus = 'McGill_clips';
end
clear stimType;

algoType = input('stimulus:  cross-correlation (1) or scg-regression (2):  ');
if algoType==1
    option.algorithm = 'crossCorr';
else
    option.algorithm = 'scg';
end
clear algoType;

% partion data:  most for training, some for validation, some for test 
nMovies.train  = 4;
nMovies.valid = 1; 
nMovies.test = 1;
nMovies.total = nMovies.train + nMovies.valid + nMovies.test; 
if (nMovies.total > 6)
    error('sorry, we only have 6 movie response datasets');
end

imgSiz  = 32;     % width/height of stimulus/filter/rfMap  
nPixels = imgSiz^2; 

durSec = 5;   refreshHz = 75;  % simulate 5 seconds at 75 hz frame rate
nFrames = durSec*refreshHz;  % e.g. 5 sec at 75 hz  (="ndata" in netlab)
if nFrames>375
    error('too many frames for these movie files !');
end

% specify model receptive field (Gabor function followed by half-power law)
model.lambda = 8;  
model.ori    = -45;  
model.pwrExp = 2; % input('power law exponent:   ');

%  graph specs
fig.stim.pos          = [400 200 300 300];     %[xOff yOff xSize ySize];
fig.stim.handle       = figure('position',fig.stim.pos,'name','stimulus');
fig.model.pos         = [50 600 300 400];  
fig.model.handle      = figure('position',fig.model.pos,'name','model');
fig.lcurve.pos         = [50 600 300 400];  
fig.lcurve.handle      = figure('position',fig.lcurve.pos,'name','learning curve');

%  create model filter, and plot in Figure 1 
rfModel = makeModelRF(model,imgSiz);    % creates model filter (Gabor function)
rfModelVec  = reshape(rfModel,1,nPixels);      % make a 1d version, for later use            

%  partition full dataset into 2 subsets, for training and validation 
stimMovie      = zeros(nPixels,nFrames);
stimMovieTrain   = [];
stimMovieValid  = [];
stimMovieTest = [];
respTrain      = [];
respValid     = [];
respTest    = [];

for iMovie=1:nMovies.total
    getStimulusMovies;            % -> stimMovie = nPixels x nFrames, range -1 to +1
    output = rfModelVec*stimMovie;    % linear filter response to the stimulus    
    output = hwr(output);  % half-wave rectify (set negative values to zero)
    output = output.^model.pwrExp;  % power-low for positive values
                   
    % accumulate results in dataset partitions: 
    if iMovie<=nMovies.train
        stimMovieTrain = [stimMovieTrain stimMovie]; % nPixels x nFrames*nMovies.train
        respTrain    = [respTrain    output];        %     1   x nFrames*nMovies.train  
    elseif iMovie > nMovies.train && (iMovie <= nMovies.valid + nMovies.train)
        stimMovieValid  = [stimMovieValid stimMovie];  % nPixels x nFrames*nMovies.valid
        respValid     = [respValid    output];          %     1   x nFrames*nMovies.valid 
    else
        stimMovieTest  = [stimMovieTest stimMovie];  % nPixels x nFrames*nMovies.valid
        respTest     = [respTest    output];          %     1   x nFrames*nMovies.valid   
    end                      
end  % end of iMovie-loop


if strcmp(option.algorithm,'scg')
    % initialize options for optimization
    nin  = imgSiz^2; % number of inputs
    nout = 1;        % number of outputs:  one neuron
    netOptions     = zeros (1,18); 
    netOptions(1)  = 0;
    netOptions(2)  = .0001;      % termination criterion: distance moved
    netOptions(3)  = netOptions(2); % for scg, use VERY small value, eg 10^-9
    netOptions(14) = 200;    % max no of iterations - should be >= no of dim.s ?
    
    commandwindow
    % regularization -> loop through alphas (alpha-loop should begin here)
    alphas = [10 90 100 175 200 225 350 400 450 500 600 700 900 1000 2500 5000 7500 10000];
    vafs_train(1:length(alphas)) = 0;
    vafs_valid(1:length(alphas)) = 0;
    loss_valid(1:length(alphas)) = 0;
    for i = 1:length(alphas)
        % estimate rfMap, for this alpha
        net = glm(nin, nout,'linear',alphas(i));       % initialize structure  
        net.w1 = 0*net.w1;    net.b1 = 0*net.b1;   % sparse prior
        [net, netOptions] = netopt(net,netOptions,stimMovieTrain',respTrain','scg');
        rfMap2d   = reshape(net.w1,imgSiz,imgSiz);  % reshape to 2d 
        rfMapVec    = reshape(rfMap2d,nPixels,1);    % make a 1-d version

        % use rfMap estimate to generate prediction of the training and validation responses)
        predRespTrain = rfMapVec'*stimMovieTrain; 
        predRespValid = rfMapVec'*stimMovieValid; 

        residTrainNew = respTrain - predRespTrain;  % Training error - might not use this
        residValidNew = respValid - predRespValid;  % residual - error in prediction of validation response
        loss_valid(i) = mean(residValidNew.^2);

        % could also calculate VAFs in addition to loss for validation dataset
    end;

    % Get the index of the alpha with the minimum loss
    [opt_loss, opt_alpha] = min(loss_valid);
    sample_alphas = [1 100000 alphas(opt_alpha) ];

    % Plot the RF model filter (Figure 4 will plot various alpha values)
    figure(4)
    subplot(2,2,1);
    imagescZadj(rfModel); axis image; axis off; colorbar; title('RF model filter');
    lim = caxis; hold on;

    for i = 1:length(sample_alphas)
        % Estimate rfMap, for this alpha
        net = glm(nin, nout,'linear',sample_alphas(i));       % initialize structure  
        net.w1 = 0*net.w1;    net.b1 = 0*net.b1;   % sparse prior
        [net, netOptions] = netopt(net,netOptions,stimMovieTrain',respTrain','scg');
        rfMap2d   = reshape(net.w1,imgSiz,imgSiz);  % reshape to 2d 
        rfMapVec    = reshape(rfMap2d,nPixels,1);
        % Graph estimated rfMap, below "actual" (model) receptive field:
        figure(4); 
        if i == 1
            subplot(2,2,2);  
            imagescZadj(rfMap2d); hold on; axis image; axis off; colorbar; caxis(lim);
            title(sprintf('SCG RF estimate for alpha = %d',sample_alphas(i)));  
        elseif i == 2  
            subplot(2,2,3);  
            imagescZadj(rfMap2d); hold on; axis image; axis off; colorbar; caxis(lim);
            title(sprintf('SCG RF estimate for alpha = %d',sample_alphas(i))); 
        else
            subplot(2,2,4);  
            imagescZadj(rfMap2d); hold on; axis image; axis off; colorbar; caxis(lim);
            title(sprintf('SCG RF estimate for alpha = %d',sample_alphas(i)));  
        end
    end;

    figure(fig.lcurve.handle);
    semilogx(alphas, loss_valid, 'b'); grid;
    if strcmp(option.stimulus, 'white')
        title(sprintf('Mean Square Error Loss for Varying Regularization \n Hyperparameter for White Noise Stimuli'))
    else
        title(sprintf('Mean Square Error Loss for Varying Regularization \n Hyperparameter for Natural Image Stimuli'))
    end;
    xlabel('Regularization Parameter \alpha')
    ylabel('MSE')


elseif strcmp(option.algorithm,'crossCorr')
    nLags=1; maxLag=0;  % (settings to make xcorr.m give us what we want)
    crossCorrAll = zeros(nPixels,nMovies.train*nFrames);
    for iPix=1:nPixels
        crossCorrAll(iPix,:) = xcorr(respTrain,stimMovieTrain(iPix,:),maxLag,'unbiased');  % cross-correlation
    end
    rfMap = crossCorrAll(:,end-nLags+1:end);  % only take positive lagged values
    rfMap2d = reshape(rfMap,imgSiz,imgSiz);  % reshape into 2d
    rfMapVec  = reshape(rfMap2d,nPixels,1);    % make a 1-d version
    clear iPix crossCorrAll rfMap nLags maxLag;
else
    error('unrecognized algorithm');
end

% Show some example stimulus images
figure(fig.stim.handle);
for ix=1:4
    stimImgVec = stimMovie(:,ix+20);
    stimImg = reshape(stimImgVec,imgSiz,imgSiz);
    subplot(2,2,ix);
    imagescZadj(stimImg); axis image; axis off; colormap('gray');
end

% Graph estimated rfMap, below "actual" (model) receptive field:
figure(fig.model.handle); 
subplot(2,1,1);
imagescZadj(rfModel);  hold on; axis image; axis off; colorbar; title('RF model filter');
lim = caxis;
subplot(2,1,2);  
imagescZadj(rfMap2d); colorbar; caxis(lim);
if strcmp(option.algorithm,'scg')
    title(sprintf('scg RF estimate for alpha = %d',alphas(5)));  
else
    if strcmp(option.stimulus, 'white')
        title(sprintf('RF estimate for cross-correlation \n for white noise stimuli')); 
    else
        title(sprintf('RF estimate for cross-correlation \n for natural image stimuli')); 
    end     
end
axis image;  axis off;  colorbar;  

% Use rfMap estimate to generate prediction of the training and validation responses)
predRespTrain = rfMapVec'*stimMovieTrain; 
predRespValid = rfMapVec'*stimMovieValid; 
predRespTest = rfMapVec'*stimMovieTest; 

residTestNew = respTest - predRespTest;  % residual - error in prediction of test response

% calculate VAF for test dataset
vaf.R_matrix = corrcoef(respTest,predRespTest);  % -> 2x2 matrix, ones on diagonal
vaf.offDiag = vaf.R_matrix(1,2);
vaf.vaf = vaf.offDiag^2.;
fprintf(1,'\nVAF for test dataset = %5.1f percent\n', 100*vaf.vaf);   
