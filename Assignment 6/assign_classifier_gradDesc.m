% assign_classifer_gradDesc.m 
%
% demo of single-layer binary classifer, using gradient descent
%   sum-square error, sigmoid activation
%   batch mode
%   derived from lab exercise in Bruno Olshausen's VS298 / lab2s.m - single neuron learning
%   uses Bruno's data files, "apples" and "oranges"
% finds line, w0 + w1*x1 + w2*x2 = 0, to best separate 2 categories of data pts

clear all;   close all;     
fprintf(1,'\n\n\n\n');

rng('default');   % "standard" random number seed -> reproducible simulations

%load data files (from Bruno's VS298), and initialize data array
load apples;    load oranges;   % -> pair of 2x10 matrices

data=[apples oranges];          % -> 2x20 matrix
[N nDataPts]=size(data);               % N=2, K=20
 

% initialize teacher (1x20 array)
teacher = [0*ones(1,nDataPts/2) +1*ones(1,nDataPts/2)]; 
lambda = 1;    % sigmoid parameter, in activation function

% learning rate
eta = input('learning rate:   ');
if isempty(eta)  error('sorry you MUST specify a learning rate !');  end

nTrials = input('number of batch-mode trials (e.g. 500):   ');
if isempty(nTrials)  ;  nTrials=500;  end

% Regularization parameter alpha
alpha = input('alpha:   ');
if isempty(alpha)  error('sorry you MUST specify a regularization coefficient !');  end

% initialize weights
w  = randn(2,1);          % 2x1 array
w0 = randn(1);            % scalar

% initialize data plot, and display hyperplane implied by initial-guess weights
if isunix
    figHanMain = figure('position',[60 1000 300 600]);
elseif ispc
    figHanMain = figure('position',[60   60 300 600]);
else
    error('unrecognized operating system');
end

subplot(2,1,1);
plot(apples(1,:),apples(2,:),'b+',oranges(1,:),oranges(2,:),'ro');  hold on
x1=0:4;
x2=-(w(1)*x1+w0)/w(2);
axis([0 4 -1 3])
h=plot(x1,x2);   grid on;   drawnow;

xlabel('x_1')
ylabel('x_2')
if alpha ~= 0
  title(sprintf('Function Performed by Neuron with \n Learning Rate %.1f and \n Regularization Hyperparameter %.3f', eta, alpha))
else
  title(sprintf('Function Performed by Neuron with \n Learning Rate %.1f', eta))
end

% Array to hold losses
loss_array = zeros(1, nTrials);
weight_array0 = zeros(1, nTrials);
weight_array1 = zeros(1, nTrials);
weight_array2 = zeros(1, nTrials);

for iTrial=1:nTrials    % loop over trials (iterations)
   % initialize dw's and loss_sum
   dw(1)  = 0;
   dw(2)  = 0;
   dw0    = 0;
   loss_sum = 0;
   
   % loop over training set
   for iDataPt=1:nDataPts
        % compute neuron output
        u = w0 + w(1)*data(1,iDataPt) + w(2)*data(2,iDataPt); 
        
        y = 1 / (1 + exp(-lambda*u));  % activation
            
        % compute error
        E = teacher(iDataPt) - y;     % scalar value of raw error
        
        % accumulate dw and dw0
        dw(1) = dw(1) + E*(exp(-lambda*u)/((1+exp(-lambda*u))^2))*data(1,iDataPt);
        dw(2) = dw(2) + E*(exp(-lambda*u)/((1+exp(-lambda*u))^2))*data(2,iDataPt);
        dw0   = dw0   + E*(exp(-lambda*u)/((1+exp(-lambda*u))^2));            
        
        loss = E^2;     % loss = error-squared
        
        % accumulate error
        loss_sum = loss_sum + loss;     % accumulated loss  
   end
   
   loss_array(iTrial) = loss_sum;

   % Plot the loss
   if iTrial == 1
      subplot(2,1,2)
      e=plot(1:nTrials, loss_array); grid on; drawnow;
      xlabel('Number of Trials')
      ylabel('Loss (Squared Error)')
      if alpha ~= 0
        title(sprintf('Learning Curve with Learning Rate %.1f and \n Regularization Hyperparameter %.3f', eta, alpha))
      else
        title(sprintf('Learning Curve with Learning Rate %.1f', eta))
      end

      figure(2)
      subplot(3,1,1)
      plot_w0=plot(1:nTrials, weight_array0); grid on; drawnow;
      if alpha ~= 0
        title(sprintf('Weight Evolution with Learning Rate %.1f and \n Regularization Hyperparameter %.3f', eta, alpha))
      else
        title(sprintf('Weight Evolution with Learning Rate %.1f', eta))
      end
      ylabel('w_0')
      subplot(3,1,2)
      plot_w1=plot(1:nTrials, weight_array1); grid on; drawnow;
      ylabel('w_1')
      subplot(3,1,3)
      plot_w2=plot(1:nTrials, weight_array2); grid on; drawnow;
      ylabel('w_2')
      xlabel('Number of Trials')
   end

   % update weights
   w(1) = w(1) + eta * (dw(1) - alpha * w(1)); 
   w(2) = w(2) + eta * (dw(2) - alpha * w(2));     
   w0   = w0   + eta * (dw0 - alpha * w0);       
   weight_array0(iTrial) = w0;
   weight_array1(iTrial) = w(1);
   weight_array2(iTrial) = w(2);

  
   % update display of separating hyperplane
   x2=-(w(1)*x1+w0)/w(2);
   set(h,'YData',x2);       % "h" from above, "h=plot(...)"
   drawnow
   set(e,'YData',loss_array);  
   drawnow  

   % Update 
   set(plot_w0,'YData',weight_array0);  
   drawnow  
   set(plot_w1,'YData',weight_array1);  
   drawnow 
   set(plot_w2,'YData',weight_array2);  
   drawnow 

end
hold off
fprintf(1,'\nfinal loss = %f\n', loss_sum);





