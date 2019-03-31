function filter=recover_receptive_field()
% A function that feeds a white noise stimulus to a model of a simple V1 neuron with two-
% dimensional spatial receptive field and temporal response that is the same as that in 
% recover_filter.m (threed.m, written by Dr Christopher Pack); based on the model response, 
% this function computes the spike-triggered average of the spatial receptive field and determines
% the neuron's preferred orientation
%
% Written by Thomas Christinck 01/2019

% Create a white noise image matrix that varies temporally (range (-1,1), 20 x 20,
% length 12000).and give to the three dimensional filter
white_noise = -1 + 2 * rand(20,20,12000,1);
response = threed(white_noise);

% Get spiketimes. Now, for each of these times, we look at the corresponding time of the stimulus
% train shifted by some amount tau. The average stimulus that preceded spiketimes by amount tau
% is the spike triggered avg
thresh = 0.1;
nb_timepoints = length(response);
for i = 1:nb_timepoints
    ap_timepoints(i) = response(i) > thresh && response(i -1) < thresh;
end
spiketimes = find(ap_timepoints);

% Each of the i's is a tau 
test_times = 8;
avg_stim(1:20, 1:20, 1:test_times) = 0;
for i = 1:test_times
	indices = spiketimes - i;
	if min(indices) <= 0 || max(indices) > length(white_noise)
		indices = indices(indices > 0);
		indices = indices(indices < length(white_noise));
	end
	avg_stim(:,:,i) = mean(white_noise(:,:,indices), 3);
end

% Now we plot the spike-triggered averages for tau 1-8
newmap = summer; 
f = figure('Position', [10 10 900 600]);
f.Name = 'Spike-Triggered Averages During Peak Temporal Response';

for i = 1:8
	adjust_stim = imadjust(avg_stim(:,:,i));
	subplot(3,4,i);
	imshow(adjust_stim, 'Colormap', newmap);
	title(strcat('\t',sprintf('au = %d', i)))
end

% We're also interested in the spike-triggered average in temporal response trough
% i.e. for tau ~ 15.
% For tau ~ 15 the stimulus and response are negatively correlated - in this case white
% noise is uniform and so the negative correlation indicates a decrease in firing rate. At this
% time lag, we should see what the stimulus looks like in order to quiet neuron activity

avg_stim(1:20, 1:20, 1:test_times) = 0;
for i = 1:test_times
	indices = spiketimes - (i + 12);
	if min(indices) <= 0 || max(indices) > length(white_noise)
		indices = indices(indices > 0);
		indices = indices(indices < length(white_noise));
	end
	avg_stim(:,:,i) = mean(white_noise(:,:,indices), 3);
end

% Now we plot the spike-triggered averages for tau 13-22
newmap = summer; 
f = figure('Position', [10 10 900 600]);
f.Name = 'Spike-Triggered Averages During Temporal Response Trough';

for i = 1:8
	adjust_stim = imadjust(avg_stim(:,:,i));
	subplot(3,4,i);
	imshow(adjust_stim, 'Colormap', newmap);
	title(strcat('\t',sprintf('au = %d', i + 12)))
end


% Three-dimensional filtering operation for simulating a reverse correlation
% experiment.  The three-dimensional input a is filtered by the simulated
% neuron and the response is returned as a one-dimensional vector r.
% ccp 1/21/08

function r=threed(a)
[x,y,z]=size(a);
if (x~=y)
    error('Use square input images.');
    return;
end;
SIZE = x;
SF = 0.15;
SIG = 7;
OR = 34*pi/180;
AR = 3;
PH = 0;

k=25;  %Temporal scale factor 
slow_t=temp_imp_resp(5,k,0:.02:1);
fast_t=temp_imp_resp(3,k,0:.02:1);

b=slow_t + fast_t; % linear temporal filter

xdata = meshgrid(1:SIZE,1:SIZE); 
temp1=(xdata-SIZE/2).*cos(OR)+(xdata'-SIZE/2).*sin(OR);
temp2=(-xdata+SIZE/2).*sin(OR)+(xdata'-SIZE/2).*cos(OR); % create oriented sinusoid

f1 = exp(-(temp1.*temp1+AR*AR*temp2.*temp2) / (2*SIG^2)); % window in a gaussian
f2=  cos(2*pi*SF*temp2+PH);

rf_image = f1.*f2;
rf_image = rf_image - mean(mean(rf_image)); % subtract off the mean

for i=1:z
    dp(i)=sum(sum(squeeze(a(:,:,i)).*rf_image));
end;

c=conv(b,dp); % convolve filter with stimulus
c=c./max(c); % normalize
c(c<0)=0;

for j=1:length(c) % for spiking output
    if (c(j)>.4)
        r(j)=1;
    else
        r(j)=0;
    end;
end;

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
