function [output] = NetREX_EdgeControl(Input)

% optimization: 
%              min_{S,A}: 1/2||E-SA||_F^2 + \kappa tr(S^TL_GS) + \mu ||A||_F^2
%                    s.t. ||A||_{\infty} \leq M
%                         ||S||_{\infty} \leq C
%                         ||S_0oS||_0 \leq # kept edges
%                         ||\bar{S_0}oS||_0 \leq # adding edges
%Input:
%Input.NumGene:  number of genes
%Input.NumTF:    number of TF
%Input.NumExp:   number of experiments
%Input.GEMatrix: the Gene expression data 
%Input.GeneGraph: weighted adjancey matrix of Gene network
%Input.TFGraph:   weighted adjancey matrix of TF network
%Input.S0:        initial point for GENE-TF network
%Input.KeepEdge; Input.AddEdge; Input.kappa; Input.mu
%
%Output:
%Output.S:      GENE-TF network
%Output.A:      TF activity matrix
%
%Algorithm details can be find in "Research Diary: July07"

%process data
LapT = diag(sum(Input.TFGraph)) - Input.TFGraph;
% LapTmax = eigs(LapT, 1);
LapG = diag(sum(Input.GeneGraph)) - Input.GeneGraph;
% LapGmax = eigs(LapG, 1);



% Aold = Input.A0;%2*rand(Input.NumTF, Input.NumExp) - 1;
% Aold = (abs(Aold)<=Input.M).*Aold + (abs(Aold)>Input.M)*Input.M.*sign(Aold);
% Sold = Input.S0;
% Sold = (abs(Sold)<=Input.C).*Sold + (abs(Sold)>Input.C)*Input.C.*sign(Sold);
% [Aold] = ClosedForm4A(Input, Sold);

if(sum(Input.Sold(:)) == 0 || sum(Input.Aold(:)) == 0)
    [Sold, Aold] = NCA_1(Input, Input.S0, Input.A0, LapG, 5);
    F_Sold = Sold;
    F_Aold = Aold;
    save SoldAold_xi_1_M0.5_C0.5.mat F_Sold F_Aold
else
    Sold = Input.Sold;
    Aold = Input.Aold;
end


ObjF = Objfunction(Input, Sold, Aold, LapG);
disp(['IterNum 0: OBJ ', num2str(ObjF(1)) ])
disp(sprintf('IterNum 0: OBJ: %10.5f Fitting: %10.5f Existing: %d Adding: %d ||S||: %f', ObjF(1), norm(Input.GEMatrix-Sold*Aold, 'fro'), sum(sum((Input.Exist).*Sold~=0)), sum(sum((1-Input.Exist).*Sold~=0)), norm(Sold,'fro')))
for k = 1 : Input.IterNum
    
%%    A first S next    
%      d(k) = norm(Sold*Sold', 'fro');
%      [Anew] = PALM_A(Input, Sold, Aold, d(k));
%      
%      c(k) = norm(Anew*Anew', 'fro') + Input.kappa*norm(LapG, 'fro');
%      [Snew] = PALM_S(Input, Sold, Anew, LapG, c(k));
    
%%   S first A next
%    c(k) = norm(Aold*Aold', 'fro') + 2*Input.kappa*norm(LapG, 'fro');
%    [Snew] = PALM_S(Input, Sold, Aold, LapG, c(k));
   c(k) = norm(Aold*Aold', 'fro') + 2*norm(Input.kappa*LapG, 'fro');
   [Snew] = PALM_S_EdgeControl(Input, Sold, Aold, LapG, c(k));%PALM_S_ElasticNet(Input, Sold, Aold, LapG, c(k));

   d(k) = norm(Snew*Snew', 'fro');
   [Anew] = PALM_A(Input, Snew, Aold, d(k));

%%   Accelerate...
%     [SnewT] = RegressionWithFixSupport(Input, Sold, Anew);
%     [Snew, Anew] = NCA_1(Input, Snew, Anew, 2);
%     [Anew] = ClosedForm4A(Input, Snew);
    
    
    
    ObjF = [ObjF Objfunction(Input, Snew, Anew, LapG)];
    Aold = Anew;
    Sold = Snew;
    
    if(abs(ObjF(end)-ObjF(end-1)) < 10 || isnan(ObjF(end)))
        break;
    end
    
%     disp(['IterNum ', num2str(k), ': OBJ \t', num2str(ObjF(k+1)) ])
    disp(sprintf('IterNum %d: OBJ: %10.5f Fitting: %10.5f Existing: %d Adding: %d ||S||: %f', k, ObjF(k+1), norm(Input.GEMatrix-Sold*Aold, 'fro'), sum(sum((Input.Exist).*Sold~=0)), sum(sum((1-Input.Exist).*Sold~=0)), norm(Sold,'fro')))
end


output.A = Aold;
output.S = Sold;
output.kappa = Input.kappa;
output.eta = Input.eta;
output.lambda = Input.lambda;
output.mu = Input.mu;
output.C = Input.C;
output.M = Input.M;

end

function [val] = Objfunction(Input, S, A, LapG)
    Remaining = sum(sum((Input.Exist.*S~=0)));
    Adding = sum(sum((1-Input.Exist).*S~=0));
    val = 0.5*norm(Input.GEMatrix - S*A, 'fro')^2 + (Input.eta-Input.lambda)*Remaining + (Input.eta+Input.lambda)*Adding + Input.kappa*trace(S'*LapG*S) + Input.mu*norm(A, 'fro')^2 + Input.xi*norm(S, 'fro')^2;
end

function [Anew] = PALM_A(Input, Sold, Aold, dk)
    [m,n] = size(Aold);
    Vk = Aold - (1/dk)*(Sold'*Sold*Aold-Sold'*Input.GEMatrix);
    Anewt = (1/(1+((2*Input.mu)/dk)))*Vk;
    Anew = (abs(Anewt)<=Input.M).*Anewt + (abs(Anewt)>Input.M)*Input.M.*sign(Anewt);
end

function [Anew] = ClosedForm4A(Input, Sold)
    Anewt = inv(Sold'*Sold+Input.mu*eye(Input.NumTF))*Sold'*Input.GEMatrix;
    Anew = (abs(Anewt)<=Input.M).*Anewt + (abs(Anewt)>Input.M)*Input.M.*sign(Anewt);
end


function [Snew] = PALM_S(Input, Sold, Aold, LapG, ck)
    [m,n] = size(Sold);
    Uk = Sold - (1/ck) * (Sold*Aold*Aold' + 2*Input.kappa*LapG*Sold - Input.GEMatrix*Aold');
    
    % hard thresholding
    tmp1 = sqrt(2*(Input.eta - Input.lambda)/ck);
    tmp2 = sqrt(2*(Input.eta + Input.lambda)/ck);
    SExist = Uk.*Input.Exist;
    SAdding = Uk.*(1-Input.Exist);
    
    SExist_HT = (abs(SExist)>=tmp1).*SExist;
    SAdding_HT = (abs(SAdding)>=tmp2).*SAdding;
    SnewT = SExist_HT + SAdding_HT;
    Snew = (abs(SnewT)>=Input.C)*Input.C.*sign(SnewT) + (abs(SnewT)<Input.C).*SnewT;
    
%     P = zeros(m,n);
%     P = tmp1*Input.Exist + tmp2*(1-Input.Exist);
%     SnewT = arrayfun(@(x,p) HardThreshold(x,p), Uk, P);
%     Snew = (abs(SnewT)>=Input.C)*Input.C.*sign(SnewT) + (abs(SnewT)<Input.C).*SnewT;
%     for i = 1 : m
%         for j = 1 : n
%             if(Input.Exist(i,j) ~= 0)
%                 tmp = sqrt(2*(Input.eta - Input.lambda)/ck);
%                 Snew(i,j) = sign(Uk(i,j))*min([abs(HardThreshold(Uk(i,j), tmp)), Input.C]);
%             else
%                 tmp = sqrt(2*(Input.eta + Input.lambda)/ck);
%                 Snew(i,j) = sign(Uk(i,j))*min([abs(HardThreshold(Uk(i,j), tmp)), Input.C]);
%             end
%         end
%     end
end

function [Snew] = PALM_S_ElasticNet(Input, Sold, Aold, LapG, ck)
    [m,n] = size(Sold);
    Uk = Sold - (1/ck) * (Sold*Aold*Aold' + 2*Input.kappa*LapG*Sold  - Input.GEMatrix*Aold');
    
    b = (2*Input.xi)/ck;
    
    Uk = Uk / (1+b);
    
    % hard thresholding
    tmp1 = sqrt(2*(Input.eta - Input.lambda)/ck)/sqrt(b+1); %c / sqrt(1+b)
    tmp2 = sqrt(2*(Input.eta + Input.lambda)/ck)/sqrt(b+1);
    SExist = Uk.*Input.Exist;
    SAdding = Uk.*(1-Input.Exist);
    
    SExist_HT = (abs(SExist)>=tmp1).*SExist;
    SAdding_HT = (abs(SAdding)>=tmp2).*SAdding;
    SnewT = SExist_HT + SAdding_HT;
    Snew = (abs(SnewT)>=Input.C)*Input.C.*sign(SnewT) + (abs(SnewT)<Input.C).*SnewT;
    

end

function [Snew] = PALM_S_EdgeControl(Input, Sold, Aold, LapG, ck)
    [m,n] = size(Sold);
    Uk = Sold - (1/ck) * (Sold*Aold*Aold' + 2*Input.kappa*LapG*Sold + 2*Input.xi*Sold - Input.GEMatrix*Aold');
    UkP = (abs(Uk)>=Input.C)*Input.C.*sign(Uk) + (abs(Uk)<Input.C).*Uk;
    
    % pick up # of required edges
    SExist = UkP.*Input.Exist;
    SExist_ABS = abs(SExist);
    SAdding = UkP.*(1-Input.Exist);
    SAdding_ABS = abs(SAdding);
    
    % sort and pick up
    SExist_ABS_Sort = sort(SExist_ABS(:), 'descend');
    SExist_Pick = (SExist_ABS>=SExist_ABS_Sort(Input.KeepEdge)).*SExist;
    SAdding_ABS_Sort = sort(SAdding_ABS(:), 'descend');
    SAdding_Pick = (SAdding_ABS>=SAdding_ABS_Sort(Input.AddEdge)).*SAdding;
    Snew = SExist_Pick + SAdding_Pick;
    
end


function [Snew] = RegressionWithFixSupport(Input, Sold, Aold)
    [m,n] = size(Sold);
    Snew = zeros(m,n);
    for i = 1 : m
        Index0 = find(Sold(i,:));
        Aind0 = Aold(Index0,:);
        Stemp = Input.GEMatrix(i,:)*Aind0'*inv(Aind0*Aind0'+Input.xi*eye(length(Index0)));%Input.GEMatrix(i,:)*pinv(Aind0);%
        if(sum(Stemp == 0) > 0)
            disp(['what'])
        end
        StempT = (abs(Stemp)>=Input.C)*Input.C.*sign(Stemp) + (abs(Stemp)<Input.C).*Stemp;
        Snew(i,Index0) = StempT;
    end
    
end


function Tx = HardThreshold(x, p)
    if(abs(x) > p)
        Tx = x;
    elseif(abs(x) == p)
        t = rand(1);
        Tx = (t>=0.5)*p + (t<0.5)*0;
    else
        Tx = 0;
    end
end

function [S, A] = NCA_1(Input, Sold, Aold, LapG, IterNum)

Spport = Sold;
Ir = Input.xi*eye(Input.NumGene);
for i = 1 : IterNum
    [Anew] = ClosedForm4A(Input, Sold);
    [Snew] = bartelsStewart(2*Input.kappa*LapG+2*Ir, [], [], Anew*Anew', Input.GEMatrix*Anew');%RegressionWithFixSupport(Input, Spport, Anew);
    
    disp(sprintf('Itr %d: %f', i, norm(Input.GEMatrix - Snew*Anew, 'fro')))
    
    Aold = Anew;
    Aold = (abs(Aold)<=Input.M).*Aold + (abs(Aold)>Input.M)*Input.M.*sign(Aold);
    Sold = Snew;
    Sold = (abs(Sold)<=Input.C).*Sold + (abs(Sold)>Input.C)*Input.C.*sign(Sold);
end

S = Sold;
A = Aold;

end