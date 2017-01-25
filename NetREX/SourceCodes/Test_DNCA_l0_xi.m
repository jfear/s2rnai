clear
clc

load Male_Deletion_Data_Norm[-1,+1].mat%Fly_Female_Data.mat%DongYeon_result.mat%DongYeon_male.mat%
load Genie3_male_best_net.mat
load SoldAold_xi_4_M0.2_C0.2.mat%intial_Male_SA.mat

[Input.NumGene Input.NumTF] = size(S_prior_M');
Input.NumExp = min(size(exp_M_norm));
Input.GEMatrix = exp_M_norm;
Input.GeneGraph = (abs(corrcoef(exp_M_norm'))>0.95).*abs(corrcoef(exp_M_norm'));%sparse((GraphGE>0.5).*GraphGE);
Input.TFGraph = zeros(Input.NumTF);%sparse((GraphTF>0.3).*GraphTF);
Input.S0 = (S_prior_M').*(Net_Best_GENIE3);%S_prior_M';%double(GEData>=1) + double(GEData<=-1);
Input.Exist = (S_prior_M~=0)'.*(Net_Best_GENIE3~=0); %Input.S0>0;
Input.A0 = rand(Input.NumTF, Input.NumExp);%exp_M_norm_TF;
Input.lambda = 0.3;
Input.eta = 0.5;
Input.mu = 1;
Input.kappa = 0.1;
Input.xi = 4; % avoid ||S||_F^2 being to large
Input.IterNum = 300;
Input.C = 0.2; % bound for S
Input.M = 0.2; % bound for A
Input.KeepEdge = 100000;
Input.AddEdge = 200000;
Input.Sold = 0%Male_Sold;%zeros(Input.NumGene, Input.NumTF);
Input.Aold = 0%Male_Aold;%zeros(Input.NumTF, Input.NumExp);

% CorNet = corrcoef([exp_M_norm' exp_M_norm_TF']);
% SCor = CorNet(1:length(exp_M_norm),end-(length(exp_M_norm_TF)-1):end);
% Temp.S = double(abs(SCor)>0.6);
% 
% Exist = sum(sum((Temp.S~=0).*((S_prior_M~=0))'));
% Newadd = sum(sum((Temp.S~=0).*((1-(S_prior_M~=0))')));
% 
% 
% % PPI_Fold = PPI_Enrichment(Temp.S);
% 
% disp(['Existing: ' num2str(Exist) ' NewAdd: ' num2str(Newadd)])
% 
% Temp.A = ones(Input.NumTF, Input.NumExp);
% % save Correlation_Network_06.mat Temp
% Filename = ['Male_n11_CorNet_06.txt'];
% fileh = fopen(Filename, 'w');%fopen('PairedGene_ReNCANet_0d20d8_j0d55.txt', 'w');
% Temp = Temp;
% DisReNCA = pdist(double(Temp.S~=0 ), 'jaccard');
% NumGene = length(Genesymbol);
% count = 1;
% idd = [];
% for i = 1 : NumGene-1
%     for j = i+1 : NumGene
%         if 1-DisReNCA(count) > 0.5
%             fprintf(fileh, '%d %d %f\n', j, i, 1-DisReNCA(count));
%             idd = [idd i j];
%         end
%         count = count + 1;
%     end
% end
% fclose(fileh);

% fid  = fopen('fly_male_FBID_list.txt', 'w');
% for i = 1 : Input.NumGene
%     fprintf(fid, '%s\n', GeneFBID{i});
% end

% % %%Prior network
% fileh = fopen('Test_Male_K20k_A380k.txt', 'w');%fopen('PairedGene_ReNCANet_0d20d8_j0d55.txt', 'w');
% DisReNCA = pdist(double(Temp.S~=0 ), 'jaccard');
% NumGene = length(Genesymbol);
% count = 1;
% idd = [];
% for i = 1 : NumGene-1
%     for j = i+1 : NumGene
%         if 1-DisReNCA(count) > 0.5
%             fprintf(fileh, '%d %d %f\n', j, i, 1-DisReNCA(count));
%             idd = [idd i j];
%         end
%         count = count + 1;
%     end
% end
% fclose(fileh);
%             
%%coexpressed edge
% fid = fopen('fly_male_coexpressed_genepair.txt', 'w');
% [m,n] = size(Input.GeneGraph);
% for i = 1 : m
%     for j = i+1 : n
%         if(Input.GeneGraph(i,j) ~= 0)
%             fprintf(fid, '%d %d 1\n', i,j);
%         end
%     end
% end
% fclose(fid);


alpha = [0.02 0.025 0.03 0.035 0.04]%[0.05:0.05:1]%[0.01 0.05 0.1 0.15 0.2]%[0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8] % eta-lambda
beta = [0.2 0.21 0.22 0.23 0.24]%[1:10]%[1 20 40 60 80 100]%[2 3 4 5 6 7 8 9 10] % (eta+lambda) / (eta-lambda)
kappa = [0.05];%[0.05 0.1 ]%[0.05 0.1]

TotalEdge = [200000 250000 300000 350000 400000 450000 500000];
KeepEdge = 24146;%10000:10000:100000;

for ii = 1 : length(KeepEdge)
     for jj = 1 : length(TotalEdge)
        for kk = 1 : length(kappa)
            clc
%             lambda = (beta(jj)-alpha(ii))/2;%alpha(ii)*(beta(jj)-1)/2;
%             eta = (alpha(ii)+beta(jj))/2;%alpha(ii)*(1+beta(jj))/2;
            disp(['KeepEdge: ' num2str(KeepEdge(ii)) ' AddEdge: ' num2str(TotalEdge(jj)-KeepEdge(ii)) ' Kappa: ' num2str(kappa(kk))])
            Input.KeepEdge = KeepEdge(ii);
            Input.AddEdge = TotalEdge(jj)-KeepEdge(ii);
            Input.kappa = kappa(kk);
            [Temp] = NetREX_EdgeControl(Input);%DNCA_l0_xi(Input);
            Temp.Sex = 'M';
            Temp.KeepEdge = Input.KeepEdge;
            Temp.AddEdge = Input.AddEdge;
            Temp.kappa = Input.kappa;
            Temp.xi = Input.xi;
            eval(['save ' 'DrosDel_Male_OverlapGENIE3_NetREXEC_xi_2_mu_1_Keep(' num2str(Input.KeepEdge) ')_Add(' num2str(Input.AddEdge) ')_Kappa(' num2str(kappa(kk)) ')_SA.mat Temp' ])
        end
     end
end

% %% output gene with increasing edges
% DegreeS0 = sum(Input.S0'~=0);
% DegreeS = sum(Temp.S'~=0);
% Filename = ['GenewithAddingEdges.txt'];
% fileh = fopen(Filename, 'w');
% for i = 1 : Input.NumGene
%     if(DegreeS(i)-DegreeS0(i)>500)
%         fprintf(fileh, '%s %d\n', Genesymbol{i}, DegreeS(i)-DegreeS0(i));
%     end
% end


