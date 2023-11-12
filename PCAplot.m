function [fh1, fh2, fh3 ] = PCAplot(Z,Names,PARAMS)
%Takes Paramterize data output to arrange data into cell by paramter matrix
%Then runs PCA on the matrix to generate plots. Do not use use function for
%window comparison analysis 

%INPUTS 
    %Z-is the output of PCApl. Structure array  with fields as outputs of
    %pca funciton. (coeff, score, latent, tsquared, explained). 
     %[Z.coeff, Z.score, Z.latent, Z.tsquared, Z.explained] = pca(data_matrix)
    %Names - other output of PCApl containing labels for each row of
    %Z.score
    %PARAMS - Cell array of strings containing names of
    %paramtemeters/dimensions of the data. 
addpath('\\albecklab.mcb.ucdavis.edu\data\Code\ARCode\Downloaded');% gscatter3 function    

%component weights
fh1 = figure('Units', 'Normalized', 'OuterPosition', [0, 0.15, 1, 0.75]); 
j=0; PC =1:3;
for i = PC; j =j+1; subplot(3,1,j);Y = zeros(1,size(Z.coeff,2));
    for o = 1:size(Z.coeff,2); Y(o) = Z.coeff(o,i); end
    bar(Y);
    ax= get(gcf,'children'); set(ax,'XTick',1:size(PARAMS,2)); set(ax, 'XTickLabel', PARAMS);
    xlabel(['PC', num2str(i)]); ylabel('weight')  
end
% suptitle('Parameter weights by component')
 sgtitle('Parameter weights by component')

%pareto plot
fh2=figure('Units', 'Normalized', 'OuterPosition', [.35 0.6, 0.4,.4]);
pareto(Z.explained); 
xlabel('Principal Component');ylabel('Variance Explained (%)');
title('Cummulative Variance Expalined')

%%asemble color map
%look through 
jj = 0;
fh3 = figure('Units', 'Normalized', 'OuterPosition', [.50 ,0.1, 0.4,.9]);
for iname =[1:size(Names,2)]
    jj = jj + 1;
U = unique(Names(:,iname));
% C = linspecer(size(U,1));
C = tab10(size(U,1));
C2 = zeros(size(Names(:,iname)));
for i = 1:size(U)
U2 = strcmp(U(i),Names(:,iname));
C2(U2) = i;
end
Colormap = C(C2,:);
%   

% fh{jj} = figure('Units', 'Normalized', 'OuterPosition', [.00 ,0.15, 0.75,.8]);
% fh{jj} = figure,
subplot(size(Names,2),1,jj)
scatter3(Z.score(:,1),Z.score(:,2),Z.score(:,3),[],Colormap,'.')
xlabel(['PC1 %Var=',num2str(Z.explained(1))]); ylabel(['PC2 %Var=',num2str(Z.explained(2))]); zlabel(['PC3 %Var=',num2str(Z.explained(3))]);

hold on %add color legend
for K = 1 : size(C, 1)
  H(K) = scatter3(nan, nan, nan, [], C(K,:));
end
legend(H, U,'Location','East');
end

%Clustering in PC space , obsolete bc has limited colors
% fh3 = figure('Units', 'Normalized', 'OuterPosition', [.00 ,0.15, 0.35,.8]);
% j = 0; figs = size(Names,2);
% for i = 1:figs; j = j+1;  subplot(figs,1,j);
%     %gscatter(Z.score(:,1)',Z.score(:,2)',Names(:,i)); for 2 d plots 
% gscatter3(Z.score(:,1)',Z.score(:,2)',Z.score(:,3)',Names(:,i)'); %if you
% %want 3D
% xlabel(['PC1 %Var=',num2str(Z.explained(1))]); ylabel(['PC2 %Var=',num2str(Z.explained(2))]); zlabel(['PC3 %Var=',num2str(Z.explained(3))]);
% end
% suptitle('Individual Cells in PC space')
% sgtitle('Individual Cells in PC space')
end


%PC vs paramter correlations, but don't show these. Need PCA matrix for
%these plots 
% for k = 1:3
% figure,
% for i = 1:size(PCA,2);
%     subplot(2,6,i)
%     scatter(Z.score(:,k),PCA(:,i),'.');
%     title(PARAMS(i))
%         xlabel(['PC' num2str(k)])
%         ylabel('standardized data')
%     end
% end
% suptitle('Individual cell correlations to PC by parameter')








