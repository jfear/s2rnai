clear
clc

%% file name of the expression data (requird format)
filename = 'DESeq2_normalized_read_counts_matrix_filteredLow_NetREX_Format.txt';
ImportData = importdata(filename);

%% expression data and gene symbols
Import_GeneSymbol = ImportData.textdata;
ExprssionRaw = ImportData.data;

%% normlize Expression
NumGene = size(ExprssionRaw,1);
for i = 1 : NumGene
    MinRow = min(ExprssionRaw(i,:));
    MaxRow = max(ExprssionRaw(i,:));
    ExpNormalized(i,:) = ((ExprssionRaw(i,:) - MinRow) / (MaxRow-MinRow))*2 -1;
    if(sum(isnan(ExpNormalized(i,:)))>0)
        w=1;
    end
end


%% find overlap with MKnet
load MKnet.mat
[OverlapGene, IdMKGene, IdImportGene] = intersect(GeneSymbol, Import_GeneSymbol, 'stable');
[OverlapTF, IdMKTF, IdImportTF] = intersect(TFSymbol, Import_GeneSymbol, 'stable');
MKnet = MKnet';
MKnet_Overlap = MKnet(IdMKGene, :);
MKnet_Overlap = MKnet_Overlap(:, IdMKTF);
Expression = ExpNormalized(IdImportGene,:);
GeneSymbol_Output = GeneSymbol(IdMKGene);
TFSymbol_Output = TFSymbol(IdMKTF);

%% Input parameters
[Input.NumGene Input.NumTF] = size(MKnet_Overlap);
Input.NumExp = size(Expression,2);
Input.GEMatrix = Expression;
Input.GeneGraph = (abs(corrcoef(Expression'))>0.95).*abs(corrcoef(Expression'));%sparse((GraphGE>0.5).*GraphGE);
Input.TFGraph = zeros(Input.NumTF);%sparse((GraphTF>0.3).*GraphTF);
Input.S0 = MKnet_Overlap;%S_prior_M';%double(GEData>=1) + double(GEData<=-1);
Input.Exist = (MKnet_Overlap~=0); %Input.S0>0;
Input.A0 = rand(Input.NumTF, Input.NumExp);%exp_M_norm_TF;
Input.eta = 0;
Input.lambda = 0;
Input.mu = 1; % 1 is good for most time
Input.kappa = 0.01; % none zere if have a gene-gene network as input
Input.xi = 1; % range [0,1]
Input.IterNum = 300;
Input.C = 0.5; % bound for S
Input.M = 0.5; % bound for A
Input.KeepEdge = 100000; % number of edges in the prior remains
Input.AddEdge = 200000; % number of edges to add in
% load SoldAold_xi_1_M0.5_C0.5.mat
Input.Sold = 0;
Input.Aold = 0;


%%run NetREX
% [Temp] = NetREX_EdgeControl(Input);

%% screen parameters
% TotalEdge = [250000 300000 350000];
% KeepEdge = 50000:10000:100000;
% for ii = 1 : length(KeepEdge)
%      for jj = 1 : length(TotalEdge)
%             clc
%             disp(['KeepEdge: ' num2str(KeepEdge(ii)) ' AddEdge: ' num2str(TotalEdge(jj)-KeepEdge(ii))])
%             Input.KeepEdge = KeepEdge(ii);
%             Input.AddEdge = TotalEdge(jj)-KeepEdge(ii);
%             [Temp] = NetREX_EdgeControl(Input);%DNCA_l0_xi(Input);
%             Temp.Sex = 'F';
%             Temp.KeepEdge = Input.KeepEdge;
%             Temp.AddEdge = Input.AddEdge;
%             Temp.kappa = Input.kappa;
%             Temp.xi = Input.xi;
%             eval(['save ' 'BrainNet_xi_1_mu_1_Keep(' num2str(Input.KeepEdge) ')_Add(' num2str(Input.AddEdge) ')_SA.mat Temp' ])
%      end
% end

