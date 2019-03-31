% Changes made to the code provided by Prof Sjostrom in order to generate the figures etc. discussed in
% the report (choose number of distributions and k)

% Thomas Christinck (March 2019)


% Number of distributions
n_dist = input('Number of distributions:   ');

if isempty(n_dist) 
    disp("Using default 3 distributions") 
    n_dist = 3;
elseif n_dist ~= 2 && n_dist ~= 3
    error('The number of distributions must be 2 or 3!'); 
end;

k = input('Enter a k-value (2,3,4, or 5):   ');
if isempty(k) 
    disp("Using default k = 3 distributions") 
    k = 3;
elseif k ~= 2 && k ~= 3 && k ~= 4 && k ~= 5
    error('The k-value must be 2, 3, 4 or 5!'); 
end;


if n_dist == 3
    X = [randn(100,2)+1.5*ones(100,2);...
        randn(100,2)-2*ones(100,2);...
        randn(100,2)+[-3*ones(100,1) 2*ones(100,1)];];
else
    X = [randn(100,2)+1.5*ones(100,2);...
        randn(100,2)-2*ones(100,2)];
end;

opts = statset('Display','final');

% Get k-means clusters
[idx,ctrs] = kmeans(X,k,...
                    'Distance','sqeuclidean',...
                    'Replicates',5,...
                    'Options',opts);
figure(1);
clf;

colors = {'r.','b.','g.','y.','c.','m.'};

for i = 1:k
	plot(X(idx==i,1),X(idx==i,2),string(colors(i)),'MarkerSize',12)
	hold on;
end

plot(ctrs(:,1),ctrs(:,2),'kx',...
     'MarkerSize',12,'LineWidth',2)
plot(ctrs(:,1),ctrs(:,2),'ko',...
     'MarkerSize',12,'LineWidth',2)

if k == 2
    legend('Cluster 1','Cluster 2','Centroids',...
       'Location','NW')
elseif k == 3
    legend('Cluster 1','Cluster 2','Cluster 3','Centroids',...
       'Location','NW')
elseif k == 4
    legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Centroids',...
       'Location','NW')
elseif k == 5
    legend('Cluster 1','Cluster 2','Cluster 3','Cluster 4','Cluster 5','Centroids',...
       'Location','NW')
end

figure(2);
clf;
plot(X(:,1),X(:,2),'k.','MarkerSize',12)

figure(3);
[silh3,h] = silhouette(X,idx);
h = gca;
h.Children.EdgeColor = [0.8, 0.6, 0.7];
h.Children.FaceColor = [0.8, 0.6, 0.7];
xlabel('Silhouette Value');
ylabel('Cluster');