%outpus indexes of treatments so you don't have to type in well numbers
%when making graphs 


%Example: Datasheet = 'L:\Databases\ImagingExperimentSheets\Ram experiments\2019-07-19-10A184A1EREG-EGF.xlsx'
%[F,L] = datasheetIndexAR(Datasheet)
%Version 1
%need to update to include pre treatments and also second treatments..but
%who even uses this? 

%field names will replace '.' with '_' so instead of 0.1 ng/ml it will say
%0_1 ng/ml. 

%idea: add another field labeling each well in order 

%%BUGS THAT I NEED TO FIX. 1. controlNAN is emtpy 2. F and L must be exactly
%the same so add bywell to F. also still not being exactly the same 
%address second treatment problem 

function [Finalstruct, Labels] =  datasheetIndexAR_Tx1(Datasheet)

[pmd, idx] = iman_readdatasheet(Datasheet);
Celltypes = fieldnames(idx.Cell); %names the cell types
Treatment1 = fieldnames(idx.Tx1);
%Treatment2 = fieldnames(idx.Tx2);
for j = 1:length(Treatment1);
    %for k = 1:length(Treatment2)
    DosesT1.(Treatment1{j}) = unique(idx.Tx1.(Treatment1{j}).dose);%make a structure with field names and each field name has indexes of that
   % DosesT2.(Treatment2{k}) = unique(idx.Tx2.(Treatment2{k}).dose);
    %field name
   % Doses.unit.(Treatment1{j}) = idx.Tx1.(Treatment1{j}).dunit;
   % Doses.unit.(Treatment2{k}) = idx.Tx2.(Treatment2{k}).dunit;
    %Doses.(Treatment2{k}) = unique(idx.Tx2.(Treatment2{k}).dose);%for treatment2
    %end
end %%not sure why I'm doing this because this already exists as idx or pmd, maybe just easier to visualize 
%might be easier for running loops?
%% 
for i = 1:length(Celltypes);
    for j =  1:length(Treatment1);
        %for jj = 1:length(Treatment2);
            for l = 1:length(DosesT1.(Treatment1{j}));
                %for k = 1:length(DosesT2.(Treatment2{jj}));
                    Finalstruct.(Celltypes{i}) = idx.Cell.(Celltypes{i}).xy; %cell type idx
                    Finalstruct.(Treatment1{j}) = idx.Tx1.(Treatment1{j}).xy; %Tx1 idx
                    %Finalstruct.(Treatment2{jj}) = idx.Tx1.(Treatment1{jj}).xy; %Tx1 idx
                    Finalstruct.(strcat(Celltypes{i},'_',Treatment1{j})) = intersect(...
                        idx.Cell.(Celltypes{i}).xy,idx.Tx1.(Treatment1{j}).xy); %%treatments by cell type
                    A = idx.Tx1.(Treatment1{j}).dose == (DosesT1.(Treatment1{j})(l)); %logical for each dose
                    B = idx.Tx1.(Treatment1{j}).xy(A); %use logical to get xys for this dose, may include this
                    
                    if ~isempty(B); %only do next lines if B is not empty
                        a = strcat(Celltypes{i},'_',Treatment1{j},num2str((DosesT1.(Treatment1{j})(l))));
                        dosename = regexprep(a,{'\.'},{'_'});%field names can't have 0.05ng/ml so  need to convert '.' to '_'
                        Finalstruct.(dosename) =...
                            intersect(idx.Cell.(Celltypes{i}).xy,B); %xys for each dose by cell type
                    end
               %$end
            end
        %end
    end
end


clear A
clear B
clear dosename 

% create another structure for that takes the xy's from F and puts in the
% labels for each xy
fn = fieldnames(Finalstruct);
for p = 1:numel(fn);
    if( isnumeric(Finalstruct.(fn{p})) ); 
        for o = 1:length(Finalstruct.(fn{p}));
        a = (Finalstruct.(fn{p}));
        b = a(o);
       [row,col] = find(cellfun(@(x)isequal(x,b),pmd.xy));
         c{p,o} = [pmd.Tx1{row, col, 1},  pmd.Tx1{row, col, 2},pmd.Tx1{row,col,3}];
         Labels.(fn{p}) = c(p,:);
       end
    end
end
%get rid of empty arrays in Labels. 
fn = fieldnames(Labels);
for p = 1:numel(fn);
       Labels.(fn{p}) =  Labels.(fn{p})(~cellfun('isempty',Labels.(fn{p})));
end
%%make one more field name that has full treatment for each wells
A = cell2mat(pmd.xy(:));%convert xy cell array into matric 
B = sort(A);
C= B(~isnan(B))';%now I have a list of all xys 
%% 

for i = C;
     [row,col] = find(cellfun(@(x)isequal(x,i),pmd.xy));
     Labels.All{i} = [pmd.Cell{row, col, 1},pmd.Gene{row, col, 1} pmd.Tx1{row, col, 1},  pmd.Tx1{row, col, 2},pmd.Tx1{row,col,3}];
    Finalstruct.All = C;
end
