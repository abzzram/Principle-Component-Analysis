%Apply k-means clustering to Principle Component analysis outputs
%will compute optimal number of clusters and plot them against PC1 and PC2
%(3d axis) 

function [ fh ] = PCAkclust(Z, Names)

scr(:,1) = Z.score(:,1);scr(:,2) = Z.score(:,2); scr(:,3) = Z.score(:,3); 
%Compute optimal kmeans 
eva = evalclusters(Z.score(:,1:3),'kmeans','CalinskiHarabasz','KList',[1:10]);%test 10 clusters to find optimal k


fh = figure('Units', 'Normalized', 'OuterPosition', [.40 ,0.15, 0.35,.8]);
j = 0; figs = size(Names,2);
for i = 1:figs; j = j+1; c1 = subplot(figs,1,j);
    un = length(unique(Names(:,i)));
    if un >1;
    [clusters] = kmeans(scr, eva.OptimalK); %NOT sure what number of clusters hould be
gscatter3(scr(:,1), scr(:,2),clusters,Names(:,i)');
xlabel('PC1'); ylabel('PC2'); zlabel('cluster#');
title([num2str(eva.OptimalK), ' kmeans clusters']);
    else title('This suplot DNE, move on'); end
end 
suptitle(['K means clusts on 3 PCs: Total%Var ',num2str(sum(Z.explained(1:3)))]);


% fh = figure('Units', 'Normalized', 'OuterPosition', [.00 ,0.15, 0.35,.8]); 
% j = 0; figs = size(Names,2);
% for i = 1:figs; j = j+1; c1 = subplot(figs,1,j);
% gscatter3(scr(:,1), scr(:,2),scr(:,3),clusters)
% xlabel('PC1'); ylabel('PC2'); zlabel('PC3');
% end
 %more stuff from mathworks.com/products/demos/machine-learning/cluster-genes.html
% net = newsom(scr',[4 4]);
% net = train(net, scr');
% distances = dist(scr,net.IW{1}');
% [dnn ,center] = min(distances,[],2);
% 
% figure
% gscatter(scr(:,1),scr(:,2),center); %legend off;
% hold on
% plotsom(net.iw{1,1},net.layers{1}.distances);
% hold off
% 
% T = clusterdata(scr,'linkage','ward','SaveMemory','on','Maxclust',nclust);
% figure,
% scatter3(scr(:,1),scr(:,2),scr(:,3),10,T);
% xlabel('PC1'); ylabel('PC2'); zlabel('PC3');


end

