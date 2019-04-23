% demo_1d_RF_sysIdent_overfit.m 
%
% demonstrate overfitting, for simple regression estimation of 1-d (spatial) receptive field
% use hold-back dataset, compare learning curves for training vs validation sets
%   use to illustrate how "early stopping" can work

% note: error = mean-squares not sum-squares

% for early stopping:  save RF estimate when prediction of validation set starts to get worse

clear all;  close all;  fprintf(1,'\n\n\n\n\n\n');

rng('default');   % "standard" random number seed -> reproducible simulations

nRFpts = 32;    % number of points in receptive field (== number of parameters to be estimated)
nMeasTrain = 70;    % number of measurements to use for receptive field estimation
nMeasValid = 30;

eta = 0.1;  %  learning rate 
early_stop = 1;

num_iterations = 50; % number of batch-mode iterations 

% define a model receptive field (Gabor function), and plot it
xPtsK = 1:1:nRFpts;
mu = nRFpts/2;   lambda = nRFpts/5;   sig = lambda*0.5;
env = exp(-(xPtsK-mu).^2/(2*sig^2));  % Gaussian envelope
receptiveField = env.*sin(2*pi*xPtsK/lambda);
figure(1);    
plot(xPtsK,receptiveField,'b-');      grid;

% create input signal (stimulus set):   white noise, range from -1 to +1
stimTrain = (rand(nRFpts,nMeasTrain) - 0.5);   % nMeasTrain measurements, for nRFpts pixels
stimValid = (rand(nRFpts, nMeasValid) - 0.5);

% simulate response of the model system (receptive field) to input signal:
respTrain = receptiveField*stimTrain + 0.3*randn(1,nMeasTrain);  % (with some added noise) 
respValid = receptiveField*stimValid + 0.3*randn(1,nMeasValid);
% stim   = nRFpts x nMeas
% resp   = 1 x nMeas           % note stim and resp are ~ zero-mean (as they need to be)
% w      = 1 x nRFpts

% $$ - replicate above, to produce stimValid and respValid 


w = zeros(1,nRFpts);  % initialize weights (receptive field estimate) - "sparse prior"

errTrain = zeros(num_iterations,1);     % initialize histories

for iteration = 1:num_iterations    % loop over iterations

   respCalc = w*stimTrain;     % predicted response for estimation dataset

   % gradient descent
   dw = (respCalc - respTrain)*stimTrain'; %  gradient
   w = w - eta*dw;   % learning rule:  update weights   
   
   errTrain(iteration) = mean((respTrain - respCalc).^2);   % record error-squared for history

   % $$  test how well we can now predict the validation dataset -> errValid
   respCalc = w*stimValid;

   errValid(iteration) = mean((respValid - respCalc).^2);

   % Check if validation error has increased
   if early_stop == 1 && iteration > 1 && ((errValid(iteration) - errValid(iteration - 1)) > 0.00001)
      num_iterations = iteration;
      errTrain = errTrain(1:num_iterations);
      errValid = errValid(1:num_iterations);
      fprintf('\nNumber of iterations: %d', iteration)
      fprintf('\nLearning rate: %.1f', eta)
      fprintf('\nValidation Loss: %.3f', errValid(iteration))
      fprintf('\nTraining Loss: %.3f \n', errTrain(iteration))

      break;
   end;

   % redraw plot of receptive field estimate 
   plot(xPtsK,receptiveField,'b-',xPtsK,w,'ro');   grid;  
   xMin = min(xPtsK);    xMax = max(xPtsK);   % set axis limits, to keep things stable
   yMin = 1.5*min(receptiveField);  yMax = 1.5*max(receptiveField);
   axis ([xMin xMax  yMin yMax]);
   legend('actual receptive field','estimated receptive field');
   if early_stop == 1
      title(sprintf('Actual and Learned Receptive Field Models with \n Early Stopping at Iteration %d', iteration + 1))
   else
      title(sprintf('Actual and Learned Receptive Field Models \n without Early Stopping'))
   end;
   drawnow  
end

figure(2); 
plot(1:1:num_iterations,errTrain,'b-'); hold on;
plot(1:1:num_iterations,errValid,'r-');
han.legend = legend('Training Error','Validation Error','Location','NorthEast');  % legend in upper right
set(han.legend,'FontSize',10);
set(han.legend,'EdgeColor','black');
grid on;  xlabel('Iterations');  ylabel('MSE');
if early_stop == 1
   title(sprintf('Learning Curve with Early Stopping and Learning Rate %.2f', eta))
else
   title(sprintf('Learning Curve without Early Stopping and Learning Rate %.2f', eta))
end;
drawnow;
% $$ modifiy above to also plot "errValid"

