% % % % % % % % % % % % % % % % % 
methodname ='HyHNID'    
netname    ='PPI'  
DisSetName ='OMIM' 
fdate_start = datestr(now,'yyyy.mmm.dd-HH.MM.SS');  
fdataset        = ['data',filesep,'DataSet.PPI.OMIM.DSimHPO.mat']
load(fdataset) 
Matrix_gene_dis00      =  AdjGfD ; 
n_disgenes_eachdisease = sum(Matrix_gene_dis00,1)'; 

% % % % % % % % % % %    
rng(1111)
[n_gene_all, n_disease] = size(Matrix_gene_dis00);             
nCVTimes  = 5,   n_fold    = 5 ,   CVtype = [num2str(n_fold),'FCV'], MinSizeDisGeSet = n_fold,  

dis_IDset = find(n_disgenes_eachdisease>=MinSizeDisGeSet); 
% dis_IDset = dis_IDset(1:2)   %%%for test only %%%%%%%%%%%    
n_dis =  length( dis_IDset )  

n_disease_in_Table   = length( dis_IDset   ); 
nCV_list             = zeros( n_disease_in_Table, 1 );  	
matAUROC_nCVTimes    = cell(nCVTimes,n_disease_in_Table);
matAUPRC_nCVTimes    = cell(nCVTimes,n_disease_in_Table);
matRec5_nCVTimes     = cell(nCVTimes,n_disease_in_Table);
matPrec5_nCVTimes    = cell(nCVTimes,n_disease_in_Table);
methodset_nCVTimes   = cell(nCVTimes,n_disease_in_Table);
for  i_cv = 1:nCVTimes    
    disp(['i_cv-',num2str(i_cv) ]) 
    %
    idx_res     = 0;     
    for ii_dis = 1:n_disease_in_Table
        Matrix_gene_dis_copy = Matrix_gene_dis00 ;   
        ID_dis               = dis_IDset(ii_dis);  
        disp(['i_cv-',num2str(i_cv),', ii_dis-',num2str(ii_dis),', ID_dis-',num2str(ID_dis)]) 
        ac_gene_dis00 = Matrix_gene_dis_copy(:,ID_dis ); 
        idx_pos       = find( ac_gene_dis00 );  n_pos = length( idx_pos); 
        idx_neg       = find( ~ac_gene_dis00 ); n_neg = length( idx_neg); 
        n_fold_real   = min(n_fold, n_pos) ;   
        ind_fold_pos  = crossvalind('Kfold', n_pos, n_fold_real ) ; 
        ind_fold_neg  = crossvalind('Kfold', n_neg, n_fold_real ) ; 
        %
        matAUROC    = [];
        matAUPRC    = [];
        matRec5     = []; 
        matPrec5    = [];         
        for i_fold = 1:n_fold_real 
            tic  
            % idx_pos_train = idx_pos(ind_fold_pos~=i_fold);
            idx_pos_test    = idx_pos(ind_fold_pos==i_fold);  n_pos_test =length(idx_pos_test); 
            %   
            % idx_neg_test_WG        = idx_neg ;  n_neg_test_all = length(  idx_neg_test_WG  ) ;               
            % idx_test_pos_neg_WG    = [idx_neg_test_WG; idx_pos_test ] ;  
            % 
			idx_neg_test_X100_RC     = idx_neg(   randperm(n_neg,  min(n_pos_test*99,n_neg)  ) ); n_neg_test_x100 = length(  idx_neg_test_X100_RC  ) ; % 每个测试基因选择100个对应控制基因，总共n_pos_test*100
            idx_test_pos_neg_X100_RC = [idx_neg_test_X100_RC; idx_pos_test ] ;  
            %

            idx_test_pos_neg = idx_test_pos_neg_X100_RC; 
            AdjGfD_t                       = Matrix_gene_dis_copy; 
            AdjGfD_t(idx_pos_test,ID_dis ) = 0 ;   
			P0_G  = AdjGfD_t(:, ID_dis );     % 
			P0_D  = zeros(n_disease,1); P0_D(ID_dis)= 1; 

            % 
            TableScores = table;  
            TableScores.rand=rand(length(ac_gene_dis00),1); 
            h_para =[]; 
            [~, TableScores_t ] =A_HyHNID([],AdjGfG,AdjGfD_t,AdjDfD, P0_G,P0_D , [], [], h_para,  [], [], [], [])  ;       
            TableScores = [TableScores,TableScores_t];
            %
            methodset = TableScores.Properties.VariableNames ;  
            %
            test_real = ac_gene_dis00(idx_test_pos_neg);  
            [AUROCset, AUPRCset , Rec5set  , Prec5set ,n_method, methodnames] = getPerf(test_real,TableScores(idx_test_pos_neg,:) ) ; 
            % 
            matAUROC(i_fold,:) = AUROCset;
            matAUPRC(i_fold,:) = AUPRCset;
            matRec5(i_fold,:)  = Rec5set; 
            matPrec5(i_fold,:) = Prec5set; 
            methodset_nCVTimes{i_cv}= methodset;  
            toc
        end 
        %
        matAUROC_nCVTimes{i_cv,ii_dis} =  matAUROC;
        matAUPRC_nCVTimes{i_cv,ii_dis} =  matAUPRC;
        matRec5_nCVTimes{i_cv,ii_dis}  =  matRec5; 
        matPrec5_nCVTimes{i_cv,ii_dis} =  matPrec5;  
        % toc
    end 
    toc 
    disp('  ')
% 
end 
%  
matAUROC_nCVTimes = cat(1,matAUROC_nCVTimes{:});    
matAUPRC_nCVTimes = cat(1,matAUPRC_nCVTimes{:});    
matRec5_nCVTimes  = cat(1,matRec5_nCVTimes{:}); 
matPrec5_nCVTimes = cat(1,matPrec5_nCVTimes{:}); 
% 
matRESmean = [mean(matAUROC_nCVTimes,1);mean(matAUPRC_nCVTimes,1);mean(matRec5_nCVTimes,1);mean(matPrec5_nCVTimes,1) ]   
tbRESmean  = array2table(matRESmean, 'VariableNames', methodset_nCVTimes{1}, 'RowNames',{'AUROC','AUPRC','Rec5','Prec5'}) 
% 
%  save 
dir_results = 'results';  
if ~exist(dir_results,'dir'); mkdir(dir_results);end 
fdate_cmplt = datestr(now,'yyyy.mmm.dd-HH.MM.SS');
parastr    = sprintf('CVtype=%s_CVtime=%d_MSDGS%d', CVtype  ,  nCVTimes, MinSizeDisGeSet );   
outfile    = [dir_results,filesep,'ResPerf_',methodname,'_',netname,'_',DisSetName,'_',parastr,'_',fdate_start,'--',fdate_cmplt,'.mat'] 
save([outfile],  'tbRESmean',   'fdate_cmplt'  , '-v7.3' )   ;   

     
     
    
