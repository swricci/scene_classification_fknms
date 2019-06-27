% Cluster feature vectors for LKSPA continuous test data

% LKP_08_featvec_results.mat includes the following:
    %ACI
    %ACT
    %BGN
    %EAS
    %ECV
    %ENT
    %EPS
    %EVN
    %HFC
    %LFC
    %SNR
    %pkfreq
    %tframe
% These variables are a summary index for each minute of recording time


%load feature vectors and concatenate
load LKP_08_featvec_results.mat;

features = cat(1,LKP_08_BGN_long,LKP_08_SNR_long,LKP_08_ACT_long,LKP_08_EVN_long,LKP_08_LFC_long,...
    LKP_08_HFC_long,LKP_08_ENT_long,LKP_08_EPS_long,LKP_08_EAS_long,LKP_08_ECV_long,...
    LKP_08_ACI_long,LKP_08_pkfreq_long)';

norm_features = normalize(features);

%normalization????
%calculate distance matrix (euclidean distance)
% rows = observations, columns = variables
feat_dist = pdist(features,'seuclidean'); % 1 x num of distance calculations in matrix



%PCA
[coeff,score,latent,~,explained] = pca(features);

% explained =
% 
%    97.4490
%     1.9575
%     0.5791
%     0.0070
%     0.0051
%     0.0012
%     0.0005
%     0.0004
%     0.0001
%     0.0000
%     0.0000


figure; plot(score(:,1),score(:,2),'.');
figure; scatter3(score(:,1),score(:,2),score(:,3)); axis equal;

[coeff,score,latent,~,explained] = pca(norm_features);


figure; plot(score(:,1),score(:,2),'.');
figure; scatter3(score(:,1),score(:,2),score(:,3)); axis equal;



corrfeat = corrcoef(features);
[COEFF,latent,explained] = pcacov(corrfeat); %same as pca with norm_features


featnames = {'BGN','SNR','ACT','EVN','LFC','HFC','ENT','EPS','EAS','ECV','ACI','pkfreq'};

figure; biplot(coeff(:,1:2),'Scores',score(:,1:2),'Varlabels',featnames);
figure; scatter3(score(:,1),score(:,2),score(:,3)); axis equal;

%dbscan

kD = pdist2(norm_features,norm_features,'euclidean','Smallest',50);
figure; plot(sort(kD(end,:)));
title('k-distance graph')
xlabel('Points sorted with 50th nearest distances')
ylabel('50th nearest distances')
grid

labels = dbscan(norm_features,2,50,'Distance','euclidean');




%tsne
Y = tsne(norm_features);
figure; gscatter(Y(:,1),Y(:,2)); %boat cluster, machine gun like cluster, perch/croaker type cluster

Y = tsne(norm_features);
figure; gscatter(Y(:,1),Y(:,2));

tic
Y2 = tsne(features,'Distance','seuclidean','Standardize',true);
figure; gscatter(Y2(:,1),Y2(:,2));
toc

Y3 = tsne(features, 'Standardize',true);
figure; gscatter(Y3(:,1),Y3(:,2));

Y4 = tsne(norm_features, 'Standardize',true);
figure; gscatter(Y4(:,1),Y4(:,2));


Y5 = tsne(norm_features,'Exaggeration',10);
figure; gscatter(Y5(:,1),Y5(:,2));

Y6 = tsne(norm_features,'Exaggeration',15);
figure; gscatter(Y6(:,1),Y6(:,2));

Y7 = tsne(norm_features,'Exaggeration',2);
figure; gscatter(Y7(:,1),Y7(:,2));

Y8 = tsne(norm_features,'NumDimensions',3);
figure; scatter3(Y8(:,1),Y8(:,2),Y8(:,3));


Y9 = tsne(norm_features,'Perplexity',50);
figure; gscatter(Y9(:,1),Y9(:,2));

Y10 = tsne(norm_features,'Perplexity',20);
figure; gscatter(Y10(:,1),Y10(:,2));

%Hierarchical clustering

%calculate distance matrix (euclidean distance)
% rows = observations, columns = variables
feat_dist = pdist(features,'seuclidean'); % 1 x num of distance calculations in matrix

feat_links = linkage(feat_dist,'average'); %other options here?

figure; dendrogram(feat_links);

feat_links = linkage(feat_dist); %other options here?

figure; dendrogram(feat_links);



% cutoff = median([feat_links(end-2,3) feat_links(end-1,3)]);
% figure; dendrogram(feat_links,'ColorThreshold',cutoff)

feat_clust = cluster(feat_links,'Maxclust',60);



dissim = squareform(feat_dist);
[Y] = cmdscale(dissim);

figure; plot(Y(:,1),Y(:,2),'.');