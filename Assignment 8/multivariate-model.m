load('MLM_demo.mat')
n = 630; % numbe rof nodes
k = n * (n - 1) /2; % number of connections
t = 73; % tps
mask = triu(ones(n), 1)> 0;
train_split = round(t * 0.75);
x = zeros(t, k);
for i = 1:t
	temp = squeeze(fc(:,:,i)); % slice of functional connectivity matrix, one for each tp
	x(i, :) = temp(mask);
end

% Split data into training and test sets
x_train = x(1:train_split - 1, :);
x_test = x(train_split:end, :);
y = panasx;
y_train = y(1:train_split - 1, :);
y_test = y(train_split:end, :);
ny = size(panasx, 2); % number of y variables

R_train = corr(x_train, y_train);
% Correlations between connections and mood scores
[u, s, v] = svd(R_train,0);


% What proportion of covariance is acounted for by singular values
s = diag(s); 
cb = s.^2 / sum(s.^2);

% So the first sv is accounting for the majority of our variance
wx = zeros(n); 
wx(mask) = u(:,1);
wx = wx + wx';

% This is how all the connections are weighted in our latent variable
[~,idx] = sortrows([ci mean(abs(wx))'], [1 -2]); % sort by membership in ci, then by weight magnitude
figure; imagesc(wx(idx,idx))

% scores 
vsc_train = y_train * v(:,1);
vlo = corr(y_train, vsc_train);
figure;
barh(vlo); set(gca, 'YTick', 1:ny, 'YTickLabel', labels)
% positive mood seems to be negatively correlated with this pattern of connectivity
usc_train = x_train * u(:,1);
figure; plot(usc_train, 'r'); hold on;
plot(vsc_train, 'b')
in_corr = corr(usc_train, vsc_train)
fprintf('In sample correlation of conncetivity and mood based scores (75/25 training/test split) is %.2f\n', in_corr)

% scores 
vsc_test = y_test * v(:,1);
vlo = corr(y_test, vsc_test);
figure;
barh(vlo); set(gca, 'YTick', 1:ny, 'YTickLabel', labels)
% positive mood seems to be negatively correlated with this pattern of connectivity
usc_test = x_test * u(:,1);
figure; plot(usc_test, 'r'); hold on;
plot(vsc_test, 'b')
in_corr = corr(usc_test, vsc_test)
fprintf('Out sample correlation of conncetivity and mood based scores (75/25 training/test split) is %.2f\n', in_corr)


% Permutation order (randomly order the 73 timepoints)
perm_order = randperm(size(x,1));
x_perm = x(perm_order, :);

% Correlate permuted x with y
R_perm = corr(x_perm, y);
% Correlations between connections and mood scores
[u_perm, s_perm, v_perm] = svd(R_perm,0);