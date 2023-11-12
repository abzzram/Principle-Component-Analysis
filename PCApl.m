function [Z, PCA, PARAMS, Names ] = PCApl(d,channel,Datasheet, varargin1,varargin2)
%Principle component analysis pipeline
%function takes output of ct_data proc (d{i}.data.ekar, where i is a well)
%d should be wells by 1
% First runs ct_pulse analysis, then shows pulse detection to allow for
% validation
%then will paramterize data, run PCA, and show plots. If set to break up
%data into windows, will skip pca plots. and continue to plsr??. DONT USE
%WINDOWS FUNCTION YET, IT IS PROBLEMATIC

%If analysing multiple different experiments at once: Run PCApl for each experiment with plot
%on off, then vertcat the outsputs PCA and Names. Then run pca function on
%the concatonated matrix and then use PCAplots and PCAkclust.

%INPUTS:
%d= datastructure output of ct_dataproc. (d{i}.data.ekar, where i is a well)
%d should be wells by 1
%channel = a string of which channel to run the analysis on. ex: "ekar"
%idx = idx output from datasheetreader
%vargarin = non-defualt paramters as name value pairs.
%Datasheet = path to data sheet
%OUTPUTS:
%Z = Cell array of structures Z{g}.(coeff, score, latent, tsquared,
%explained). Where each cell is the data divided into evenly spaced
%intervals. When p.windows =0, Z will only have one cell. When p.windows
%=1, Z will have 2 intervals because the final loop will process data
%divided into two windows.
%PCA = the standarized matrix that the PCA was done. cell by paramters
%PARAMS = the names of each paramter in each column of PCA
%Names = cell by 3 matrix of labels for each cell. Cell type, drug, dose

%PARAMETERS:
%varargin1 (paramters for main function):
%'chkpulse' - Logical, 1 if you want to validate ct_pulseanalysis
%output. 'numchks' -(default = 25) numner of traces to display
%'samplingtime' - %time between samples (minutes)
%'starttime' - at which time point should the analysis start. Default
%'paramonly' - Paramterize data only (def = 0), set to 1 if you want to
%             paramterize data only
%'linear only' - Set to 1 if you only want to include linear
%terms in the model. Default is 0, this includes mean^2, max^2, mean deriv^2
%'ploton' - display P]plots. (Def =1)
%varargin2 (paramters for ct_pulse analysis):


%set default paramters
p.chkpulse = 1; %1 if you want to validate pulse detection
p.numchks = 25; %number of traces to plot during pulse validaton,
p.samplingtime = 6; %time between samples (minutes)
p.starttime = 1; %at which time point shuold the analysis start
p.windows = 0; %separate data into windows?
p.paramonly = 0;
p.linearonly = 0;
p.ploton = 1;
p.onehr = 60/p.samplingtime; %how many tps is 1 hrs/how long should windows be increasing by
p.minterval = 15; %set to 15 time points bc this minimum allowed in paramterize data
p.plotwin = 1;%is this used?
p.kmeans = 1;
q.min_rise = 0.06; %assuming erk traces, these are good defults to set
q.min_vy =0.06;
q.narm = 10;

%input option pair parsing
p.wells = 1:size(d,1);%which wells to analyze
nin = length(varargin1);
nin2 = length(varargin2);
if rem(nin,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
if rem(nin2,2) ~= 0; warning(['Additional inputs must be provided as ',...
        'option, value pairs']);  end
%Splits pairs to a structure and exhange defaults
for s = 1:2:nin;   p.(lower(varargin1{s})) = varargin1{s+1};   end
for s = 1:2:nin2;   q.(lower(varargin2{s})) = varargin2{s+1};   end

%Assemble name/value sequence for passing parameters forward
ppass = [fieldnames(p),struct2cell(p)]';
ppass2 = [fieldnames(q),struct2cell(q)]';
p.endtime = size(d{(p.wells(1))}.data.(channel),2); % time which the analysis should end

if p.windows == 1; p.windowlength = p.onehr:(p.onehr):p.endtime/2;
else p.windowlength = [p.endtime];  end

%to create labels
[pmd, idx] = iman_readdatasheet(Datasheet);


Z = cell(length(p.windowlength),1);
for g = 1:length(p.windowlength);
    %define window length
    if p.windows ==1;
        twin1 = p.starttime:(p.windowlength(g)):(p.endtime-(p.windowlength));%define first point of window
        twin1 = twin1(1:(length(twin1)-1)); %so that twin stays within length of expeirment
        twin2 =  twin1+p.windowlength(g)-1;
    else twin1 =p.starttime; twin2=p.endtime; end
    if length(twin1(1):twin2(1)) > p.minterval; %only run if window is larger than 15tp
        %conduct pulse anlaysis and parameterize data
        ctpw.(channel) =cell(max(p.wells),length(twin1)); %initialize
        paramatw = cell(max(p.wells),length(twin1));pw= cell(max(p.wells),length(twin1));%initiatlize
        for i = p.wells;
            for j = 1:length(twin1);
                ctpw.(channel){i,j} = ct_pulseanalysis(d{i}.data.(channel)(:,twin1(j):twin2(j)),ppass2{:},'MAXW',p.windowlength(g));
                %check pulse analysis
                if p.chkpulse ==1 && i == p.wells(1) && j ==1 %only check for the first window length
                    figure,
                    t = 0 ;
                    for  k = 1:(p.numchks);
                        t = t+1;
                        subplot(ceil(sqrt(p.numchks)),ceil(sqrt(p.numchks)),t);
                        hold on;
                        plot(d{i}.data.(channel)(k,twin1(j):twin2(j))');
                        scatter((ctpw.(channel){i}(k,:).pkpos),ctpw.(channel){i}(k,:).peak);
                    end
                    fprintf('Continue if peaks are detected. Shorter windows may not have as good pulse detection \n') %not sure how to fix it
                    shg; pause
                end
                %paramterize data
                [paramatw{i,j}, pw{i,j}] = ParameterizeData(d{i}.data.(channel)(:,twin1(j):twin2(j)),ctpw.(channel){i,j},'channel', channel);
            end
        end
        if p.paramonly == 0;
            if p.linearonly == 1; numpar = 8; else numpar = 11;%assemble window lengths
                Win1 = 1:numpar:(length(twin1)*numpar); Win2 = 8:numpar:(length(twin1)*numpar);
            end
            Sizes = cellfun(@numel,ctpw.(channel)); PCA = zeros(sum(Sizes(:,1)),Win2(end));
            PARAMS =cell(1,Win2(end)); Names = cell(sum(Sizes(:,1)),3); %initialize
            for j = 1:length(twin1)
                h = 0;
                for i = p.wells;
                    %create labels
                    [row,col] = find(cellfun(@(x)isequal(x,i),pmd.xy));
                    for k = 1:length(ctpw.(channel){i,j});
                        h =h+1;
                        %put each window of activity into PCA matric (one giant matrix
                        %with all windows)
                        
                        PCA(h,Win1(j):Win2(j)) = paramatw{i,j}(k,:) ;
                        PARAMS(Win1(j):Win2(j)) =  strcat(num2str(j),pw{(p.wells(1))}.params');
                        if isempty(Names{end}); %this is super sketchy, will it work for other datasheets?
                            Names(h,1) = strcat(pmd.Cell(row, col, 1),' ',pmd.Gene(row, col, 1)); %assign name
                            Names(h,2) = strrep(pmd.Tx1(row, col,1),' ','');
                            Names(h,3) = strrep(pmd.Tx1(row, col,2),' ','');
                        end
                        if p.linearonly == 0 %put in non linear parameters
                            PCA(h,Win2(j)+1) = paramatw{i,j}(k,1).^2; %mean squared
                            PCA(h,Win2(j)+2) = paramatw{i,j}(k,2).^2; %max squared
                            PCA(h,Win2(j)+3) = paramatw{i,j}(k,3).^2; %mean deriv squared
                            PARAMS(Win2(j)+1) =  {strcat(num2str(j),'mean^2')};
                            PARAMS(Win2(j)+2) =   {strcat(num2str(j),'max^2')};
                            PARAMS(Win2(j)+3) =    {strcat(num2str(j),'mean driv^2')};
                        end
                    end
                end
            end
            PCA = zscore(PCA); %standardize
            [Z{g}.coeff, Z{g}.score, Z{g}.latent, Z{g}.tsquared, Z{g}.explained] = pca(PCA);
        end
        
    else fprintf('skipping window length( %2d time points) because window is too short\n',p.windowlength(g));
    end
end
%plots
if p.ploton == 1
    [fh1, fh2,fh3] = PCAplot(Z{g},Names,PARAMS);
    if p.kmeans == 1
        [fh4] = PCAkclust(Z{g}, Names);%kmeans clustering
    end
end

end




