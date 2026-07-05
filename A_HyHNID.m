function [ L_hyper_set, TableScores, TableScores_D ] = ...
	A_HyHNID(L_hyper_set, AdjGfG,AdjGfD,AdjDfD, P0_G,P0_D , eta_iniProRatio_Dis, W_interlayer, h_para,  MDL, MDL_NormalizationType_Individual, gFS, dFS)      
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
    % % Input % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    % AdjGfG: matrix that records associatins from (f) genes (G) to Genes (G)  
    % AdjGfD: : a matrix that record associatins from Diseases (D) to Genes (G) GfD
    % AdjDfD: matrix that records  associatins from Diseases (D) to Disease (G)  
    % P0_G: initial in Gene network
    % P0_D: initial in Disease network    
    % By Ju Xiang 
    % Email: xiang.ju@foxmail.com, xiangju@csu.edu.cn  
    %  
    % 2026-6-28
    % 
    % % % % % % % % % % % % % % % % % % % % % % % % % %   
    % diffusion constant in each network is 1， and  
    if  isempty(  h_para   ) 
        h_para.alpha_g =0.125; 
        h_para.beta_g  =0.125; 
        h_para.alpha_d =0.5; 
        h_para.beta_d  =1; 
        h_para.kg = 100; 
        h_para.kd = 20; 
    end 
    alpha_g =  h_para.alpha_g ;
    beta_g  =  h_para.beta_g ;
    alpha_d =  h_para.alpha_d ;
    beta_d  =  h_para.beta_d ;
    kg = h_para.kg ;
    kd = h_para.kd ;
    %
    if ~exist('MDL','var') || isempty( MDL ); MDL = 'M'; end 
    switch MDL
        case 'M'; NormalizationType='LaplacianNormalizationMeanDegree';  
        case 'C'; NormalizationType='ProbabilityNormalizationColumn';   
        case 'R'; NormalizationType='ProbabilityNormalizationRow';      
        case 'L'; NormalizationType='LaplacianNormalization';   
        case 'None'; NormalizationType= 'None';  
        otherwise; error('There is no definition.');
    end  
	
	if isempty( eta_iniProRatio_Dis ); eta_iniProRatio_Dis = 1 ; end 
    %if isvector( W_interlayer ); eta_iniProRatio_Dis=W_interlayer(2); W_interlayer = W_interlayer(1) ;  
	%else; eta_iniProRatio_Dis      = 0.5 ;      end 
	% 
    if isempty( W_interlayer ); W_interlayer = 0.5;  end 
	IntraDiffusionConstant = 1-W_interlayer;  
    %       
    % IntraDiffusionConstant = [1] ;       
    % IntraDiffusionConstant = [ ] ;       
    % if isempty( IntraDiffusionConstant ); IntraDiffusionConstant = 1 ; end 

    %  
    if isempty( MDL_NormalizationType_Individual ); MDL_NormalizationType_Individual = 'None'; end 	
    % MDL_NormalizationType_Individual ='R'   %% 考虑 各层网络响应的均值
    switch MDL_NormalizationType_Individual
        case 'M'; NormalizationType_Individual='LaplacianNormalizationMeanDegree';  
        case 'C'; NormalizationType_Individual='ProbabilityNormalizationColumn';   
        case 'R'; NormalizationType_Individual='ProbabilityNormalizationRow';      
        case 'L'; NormalizationType_Individual='LaplacianNormalization';   
        case 'None'; NormalizationType_Individual= 'None';  
        otherwise; error('There is no definition.');
    end  	
    %
    beta = 1;  
    %%%eta_iniProRatio_Dis      = 1;  
    
    % scheme 2          
    % %     IndividualNormalizationType = NormalizationType ;    
     A_gene_gene = AdjGfG;
     A_dis_dis   = AdjDfD; 
	if isempty( gFS ); gFS='gBMAmax'; end 
	switch  gFS
		case {'gBMA', 'gBMAmax', 'gBMAm'}
			[ sim   ] = getFuncSim_IN(AdjGfD,AdjDfD, 'BMA'   )  ;  
			A_gene_gene = sparse( sim  );  sim =[]; 
			A_gene_gene = max(A_gene_gene,AdjGfG );  
		case 'None' ; A_gene_gene = AdjGfG; 
		otherwise ; error('There is no definition.')
			
	end 
	%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
	if isempty( dFS ); dFS='dBMAmax'; end 
	switch  dFS
		case {'dBMA', 'dBMAmax', 'dBMAm'}
			[ sim   ] = getFuncSim_IN(AdjGfD',AdjGfG, 'BMA'   )  ;  
			A_dis_dis = sparse( sim  );  sim =[]; 
			A_dis_dis = max(A_dis_dis,AdjDfD );  
		case 'None' ; A_dis_dis = AdjDfD; 
		otherwise ; error('There is no definition.')
			
	end 
 
    % % construct hyperedges % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
    [Ng,Nd]= size(AdjGfD);
    if isempty(L_hyper_set)
        % gene hyperedges from gene-gene associations
        %         Hg_g = getAdjKnnColumns( A_gene_gene,  kg , 0, 0 )'; 
        Hg_g = AdjGfG; 
        Hg_g(   sub2ind( size(Hg_g), 1:Ng,1:Ng )      ) = 1; Hg_g = sparse(Hg_g);  % 使 中心节点保持在超边中 
        d_node_Hg_g = sum(Hg_g,2)+eps;
        d_edge_Hg_g = sum(Hg_g,1)'+eps; 
        L_Hg_g = speye(Ng) - diag(1./d_node_Hg_g)*Hg_g*diag(1./d_edge_Hg_g)*Hg_g'; L_Hg_g = sparse(L_Hg_g);
        % disease hyperedges from disease-disease associations 
        %         Hd_d = getAdjKnnColumns( A_dis_dis,  kd , 0, 0 )'; 
        Hd_d = AdjDfD;
        Hd_d(   sub2ind( size(Hd_d), 1:Nd,1:Nd )      ) = 1; Hd_d = sparse(Hd_d);  % 使 中心节点保持在超边中
        d_node_Hd_d = sum(Hd_d,2)+eps;
        d_edge_Hd_d = sum(Hd_d,1)'+eps;
        L_Hd_d = speye(Nd) - diag(1./d_node_Hd_d)*Hd_d*diag(1./d_edge_Hd_d)*Hd_d'; L_Hd_d = sparse(L_Hd_d);
        %
        L_hyper_set.L_Hg_g = L_Hg_g;
        L_hyper_set.L_Hd_d = L_Hd_d;
        L_hyper_set.Hg_g = Hg_g;
        L_hyper_set.Hd_d = Hd_d;

%         L_hyper  =  [ alpha_g*L_Hg_g ,         sparse(Ng,Nd)              ;  ...
%                         sparse(Nd,Ng),         alpha_d*L_Hd_d  ]    ;   
%         L_hyper_set.L_hyper = L_hyper;

        % gene hyperedges from gene-disease associations
        Hg_d = sparse(AdjGfD);  
        d_node_Hg_d = sum(Hg_d,2)+eps; 
        d_edge_Hg_d = sum(Hg_d,1)'+eps;  %size(d_edge_Hg_d) 
        L_Hg_d = speye(Ng) - diag(1./d_node_Hg_d)*Hg_d*diag(1./d_edge_Hg_d)*Hg_d'; L_Hg_d = sparse(L_Hg_d);
        % disease hyperedges from disease-gene associations 
        Hd_g = sparse(AdjGfD'); 
        d_node_Hd_g = sum(Hd_g,2)+eps; 
        d_edge_Hd_g = sum(Hd_g,1)'+eps; 
        L_Hd_g = speye(Nd) - diag(1./d_node_Hd_g)*Hd_g*diag(1./d_edge_Hd_g)*Hd_g'; L_Hd_g = sparse(L_Hd_g);
        %
        L_hyper_set.L_Hg_d = L_Hg_d;
        L_hyper_set.L_Hd_g = L_Hd_g; 
        L_hyper_set.Hg_d = Hg_d;
        L_hyper_set.Hd_g = Hd_g; 
%         %
%         L_hyper  =  [ alpha_g*L_Hg_g + beta_g*L_Hg_d ,         sparse(Ng,Nd)              ;  ...
%                           sparse(Nd,Ng),                       alpha_d*L_Hd_d +  beta_d*L_Hd_g  ]    ; 
%         %
%         L_hyper_set.L_hyper = L_hyper;

    else
        L_Hg_g = L_hyper_set.L_Hg_g;
        L_Hd_d = L_hyper_set.L_Hd_d; 
        L_Hg_d = L_hyper_set.L_Hg_d;
        L_Hd_g = L_hyper_set.L_Hd_g; 

        Hg_g = L_hyper_set.Hg_g;
        Hd_d = L_hyper_set.Hd_d;   
        Hg_d = L_hyper_set.Hg_d;
        Hd_g = L_hyper_set.Hd_g;         
    end    
    %
    L_hyper  =  [ alpha_g*L_Hg_g + beta_g*L_Hg_d ,         sparse(Ng,Nd)              ;  ...
                      sparse(Nd,Ng),                       alpha_d*L_Hd_d +  beta_d*L_Hd_g  ]    ; 
  
	 
	%  NormalizationType  D-G矩阵做标准化 效果可以（优先采用）  
    M_GfD =  ( AdjGfD );  %size(M_GfD) 
    M_DfG =  ( AdjGfD'); % size(M_DfG) 
%     M_GfG = sparse( A_gene_gene);
%     M_DfD = sparse( A_dis_dis);

    % 临时测试 每一类网络单独标准化  效果可以 
	% NormalizationType_Individual = 'ColNormalization' ; 
    if ~strcmp(NormalizationType_Individual,'None') 
        % WAdj = getNormalizedMatrix_IN(Adj, NormalizationType, SetIsolatedNodeSelfLoop  ,IsMeanDegreeForExistNodeSet )
        M_GfG = getNormalizedMatrix_IN( M_GfG , NormalizationType_Individual, true  );  %测试 gene  网络 标准化处理 
        M_DfD = getNormalizedMatrix_IN( M_DfD , NormalizationType_Individual, true  );   %测试 dis  网络 标准化处理 
        
        M_GfD =  getNormalizedMatrix_IN(M_GfD   , NormalizationType_Individual, false , false);  %   from disease space to gene space 
        M_DfG =  getNormalizedMatrix_IN(M_DfG   , NormalizationType_Individual, false , false);  %   from gene space to disease space
     end 
   % 
    
    Mall   = [   IntraDiffusionConstant.*A_gene_gene,          W_interlayer.*M_GfD              ;  ...
                          W_interlayer.*M_DfG,                 IntraDiffusionConstant.*A_dis_dis   ]    ; 
		    
    Mall = getNormalizedMatrix_IN( Mall , NormalizationType, true  );  % nnz(Mall)
    Mall = diag(sum(Mall,2))-Mall +L_hyper ;             %  nnz(Mall)
    Mall = sparse( Mall + beta*speye( size(Mall) )  );  
 
    % 
    P0_all = sparse( [ (1-eta_iniProRatio_Dis)*P0_G; eta_iniProRatio_Dis*P0_D]  ) ;   
    % % runing  % % %        
    Tf       = 5;  
    %     [t,Y_t_node] = ode23(@(t,X) (  - Mall*X)  ,[0 Tf], P0_all);
    [t,Y_t_node] = ode23(@(t,X) ( getdFdt(t,X, Mall )  )  ,[0 Tf], P0_all);
    % Xt_seq = Y_t_node'; 
   
%     Xmean = trapz( t,Y_t_node )';      
    Xmax  = max( Y_t_node , [], 1 )';      
    %  for genes    
    [N_genenode, N_disnode] = size( AdjGfD );
%     Xmean_G = Xmean(  1:N_genenode , :   );     
    Xmax_G  = Xmax(   1:N_genenode , :   );     
    % 
    TableScores = table;          
%     TableScores.(['Xmn'] )    = Xmean_G  ;      
%     TableScores.(['Xmx'] )    = Xmax_G   ;   
    TableScores.(['HyHNID'] )    = Xmax_G   ;   
   
    %
    if  nargout>2
        %%% for diseases
        % Xmean_D = Xmean(   N_genenode +1:end,                  : );     
        Xmax_D  = Xmax(    N_genenode +1:end,                  : );     
        TableScores_D = table;      
        TableScores_D.(['Xmn_d'] )    = Xmean_D  ;      
        TableScores_D.(['Xmx_d'] )    = Xmax_D   ;    
    else
        TableScores_D = [] ;                
    end 
      
end
 
% % % % % % % % % % %
function  dFt = getdFdt(t,X, Mall )

    dFt =  - Mall*X ; 


end 
 
% % % % % % % % % % % % % % % % % % % 
function WAdj = getNormalizedMatrix_IN(Adj, NormalizationType, SetIsolatedNodeSelfLoop  ,IsMeanDegreeForExistNodeSet )
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
% Adj  adjecent matrix
% % NormalizationType: 
% % 'probability normalization'  
% % 'laplacian normalization' 
% SetIsolatedNodeSelfLoop    set isolated node
% >= Matlab 2016
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     
%     if ~issparse(Adj)
%         Adj = sparse( Adj );
%     end   
    if ~exist('IsMeanDegreeForExistNodeSet','var') || isempty(IsMeanDegreeForExistNodeSet)
        IsMeanDegreeForExistNodeSet = false;
    end
 

    if ischar(NormalizationType)
    %         NormalizationType =  (NormalizationType);
        switch  lower( NormalizationType )
            case lower( { 'column','col',  ...
                    'ProbabilityNormalizationColumn','ProbabilityNormalizationCol',...
                    'ProbabilityColumnNormalization','ProbabilityColNormalization',...
                    'NormalizationColumn','NormalizationCol' , ...
                    'ColumnNormalization','ColNormalization'   })
                NormalizationName = 'ProbabilityNormalization' ;  %  'Random Walk'  
                dim =1;
            case lower({ 'row' ,'ProbabilityNormalizationRow' ,'NormalizationRow' ,'ProbabilityRowNormalization' ,'RowNormalization'   })
                NormalizationName = 'ProbabilityNormalization' ;  %  'Random Walk'  
                dim =2;
            case lower('LaplacianNormalization')
                NormalizationName = NormalizationType; 
            case lower('LaplacianNormalizationMeanDegree')
                NormalizationName = NormalizationType; 
            case lower('ColNorm2')
                NormalizationName = NormalizationType; 
            case lower('RowNorm2')
                NormalizationName = NormalizationType; 
            case lower({'none', 'None', 'NONE'})
                % NormalizationName = 'None'; 
                WAdj = Adj; 
                return; 
            otherwise
                error(['There is no type of normalization: ',char( string(NormalizationType) )] );
        end
        
    elseif isnumeric(  NormalizationType   ) 
        NormalizationName =  ( 'ProbabilityNormalization' ) ;  %  'Random Walk'  
        dim = NormalizationType; 
        
    elseif isempty( NormalizationType )
        WAdj = Adj; 
        return;  
        
    else; error('There is no defintion of NormalizationType')
    end 
    % NormalizationName = lower( NormalizationName );
    % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %     %
    matlabveryear = textscan(version('-release'),'%d') ; matlabveryear=matlabveryear{1};
    if matlabveryear>=2016 
        switch lower( NormalizationName )
            case lower( 'ProbabilityNormalization' )
                degrees = sum(Adj,dim);
                if any( degrees~=1)
                    WAdj = Adj./ ( degrees+eps  );           
                    % % WAdj = Adj./ repmat( degrees +eps,[size(Adj,1),1]); 
                else
                    WAdj = Adj; 
                end
                % 
                if SetIsolatedNodeSelfLoop  && size(Adj,1)==size(Adj,2) 
                    ii = find( ~degrees ); 
                    idx = sub2ind( size(Adj), ii,ii ); 
                    WAdj(idx) = 1;  % set to be 1 for isolated nodes, 
                end

            case lower( 'LaplacianNormalization')
                deg_rowvec = ( sum(Adj,1) ).^0.5;  
                deg_colvec = ( sum(Adj,2) ).^0.5;   
                WAdj = (Adj./(deg_colvec+eps))./(deg_rowvec+eps) ;    
                % 
                if SetIsolatedNodeSelfLoop && size(Adj,1)==size(Adj,2)
                    ii = find( ~sum(Adj,2) ) ; 
                    % size(  WAdj )
                    % size(  Adj )
                    WAdj( sub2ind( size(Adj), ii,ii ) ) = 1;  % set to be 1 for isolated nodes, 
                end

            case lower( 'LaplacianNormalizationMeanDegree')
    %             n_node = length( Adj ); 
                k_col = sum( Adj, 2 ); 
                k_row = sum( Adj, 1 );             
                if IsMeanDegreeForExistNodeSet 
                    km1 = sum(k_col)./ (nnz(k_col)+eps);  
                    km2 = sum(k_row)./ (nnz(k_row)+eps);                  
                else
                    km1 = sum(k_col)./ length(k_col);  
                    km2 = sum(k_row)./ length(k_row);                  
                end
                WAdj = Adj./( (km1.^0.5)*(km2.^0.5)  +eps) ;    
                % 
                if SetIsolatedNodeSelfLoop  && size(Adj,1)==size(Adj,2)
                    ii = find( ~sum(Adj,2) ); 
                    WAdj( sub2ind( size(Adj), ii,ii ) ) = 1;  % set to be 1 for isolated nodes, 
                end

            case lower( {'ColNorm2'} )   
                WAdj = Adj./ ( sqrt(sum( Adj.^2 ,1 )) +eps ); 

            case lower( {'RowNorm2'} )    
                WAdj = Adj./ ( sqrt(sum( Adj.^2 ,2 )) +eps ); 

            case lower( {'None','none'} )
                WAdj = Adj;   % 不做任何处理  
            otherwise
                error(['NormalizationName is wrong: ',char(string(NormalizationName) )   ]);
        end

        % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %    
    else
    %     error('errorerrorerrorerrorerrorerrorerrorerror');
        switch lower( NormalizationName )
            case lower( 'ProbabilityNormalization' )
                degrees = sum(Adj,dim);
                if any( degrees~=1)
                    % WAdj = Adj./ ( degrees+eps  ); 
                    WAdj = getMatrixOperation(Adj, degrees, './') ;
                     % % WAdj = Adj./ repmat( degrees +eps,[size(Adj,1),1]); 
                else
                    WAdj = Adj; 
                end
                % 
                if SetIsolatedNodeSelfLoop  && size(Adj,1)==size(Adj,2) 
                    ii = find( ~degrees ); 
                    idx = sub2ind( size(Adj), ii,ii ); 
                    WAdj(idx) = 1;  % set to be 1 for isolated nodes, 
                end

            case lower( 'LaplacianNormalization')
                deg_rowvec = ( sum(Adj,1) ).^0.5;  
                deg_colvec = ( sum(Adj,2) ).^0.5;   
                % %             WAdj = (Adj./(deg_colvec+eps))./(deg_rowvec+eps) ; 
                WAdj = getMatrixOperation(Adj,  deg_colvec, './') ;
                WAdj = getMatrixOperation(WAdj, deg_rowvec, './') ;
                % 
                if SetIsolatedNodeSelfLoop && size(Adj,1)==size(Adj,2)
                    ii = find( ~sum(Adj,2) ) ; 
                    % size(  WAdj )
                    % size(  Adj )
                    WAdj( sub2ind( size(Adj), ii,ii ) ) = 1;  % set to be 1 for isolated nodes, 
                end

            case lower( 'LaplacianNormalizationMeanDegree')
                n_node = length( Adj );
                km = sum( Adj(:) )./ n_node;  
                WAdj = Adj./( (km.^0.5)*(km.^0.5)  +eps) ;    
                % 
                if SetIsolatedNodeSelfLoop  && size(Adj,1)==size(Adj,2)
                    ii = find( ~sum(Adj,2) ); 
                    WAdj( sub2ind( size(Adj), ii,ii ) ) = 1;  % set to be 1 for isolated nodes, 
                end

            case lower( {'ColNorm2'} )   
                % % WAdj = Adj./ ( sqrt(sum( Adj.^2 ,1 )) +eps ); 
                WAdj = getMatrixOperation(Adj,  sqrt(sum( Adj.^2 ,1 )), './') ;

            case lower( {'RowNorm2'} )    
                % % WAdj = Adj./ ( sqrt(sum( Adj.^2 ,2 )) +eps ); 
                WAdj = getMatrixOperation(Adj,  sqrt(sum( Adj.^2 ,2 )), './') ;

            case lower( {'None','none'} )
                WAdj = Adj;   % 不做任何处理  
            otherwise
                error(['NormalizationName is wrong: ',char(string(NormalizationName) )   ]);
        end 

    end
 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [ kd,km ] = getGIPSim(Adm_interaction, gamma0_d, gamma0_m,  AvoidIsolatedNodes   )
function [ FuncSim_final , FuncSim2_final ] = getFuncSim_IN(AdjGfD,DisSim, FuncSimMethod   ) 
% Xiang  2019-11-16 
% Ref: DOSim.  
% if ~exist('gamma0_d','var') || isempty(gamma0_d)
%    gamma0_d = 1 ;  
% end
% if ~exist('gamma0_m','var') || isempty(gamma0_m)
%    gamma0_m = 1 ;  
% end
% if ~exist('AvoidIsolatedNodes','var') || isempty(AvoidIsolatedNodes)
%    AvoidIsolatedNodes = 1 ;  
% end
% % % % % % % % % % % % % % % % % % % 
[ng_all, nd_all] = size( AdjGfD  ); 
% if AvoidIsolatedNodes
nodes_g = sum(AdjGfD,2)~=0 ; 
nodes_d = sum(AdjGfD,1)~=0 ; 
ExistIsolatedNodes = ng_all~=nnz(nodes_g) || nd_all~=nnz( nodes_d )   ; 
if ExistIsolatedNodes
    Agd= ( AdjGfD( nodes_g , nodes_d ) );
    DS  = DisSim(nodes_d,nodes_d); 
else 
    Agd= ( AdjGfD );
    DS  = DisSim ; 
end 
% % % % % % % % % % % % % % % % % % % % % % % % % % 
Agd      = logical(Agd) ; 
[ng, nd] = size( Agd  );  
DS(sub2ind(size(DS),1:nd,1:nd ) )=  1 ; 
Nd_col   = sum(Agd,2);  % Number of diseases of each gene. 
FuncSim2_final =0; 
% 
if strcmpi(FuncSimMethod, 'Mean')  
    M_gSet_SumDS = zeros(ng,nd) ; 
    for  ig = 1:ng
        M_gSet_SumDS(ig,:) = sum( DS( Agd(ig,:), : ), 1 )  ;   
    end
    FuncSim = M_gSet_SumDS*Agd'./(Nd_col.*Nd_col' +eps) ;    

else   
    M_gSet_maxDS = zeros(ng,nd) ; 
    for  ig = 1:ng 
        M_gSet_maxDS(ig,:) = max( DS( Agd(ig,:), : ), [], 1 )  ;   
    end 
    %
    switch  lower( FuncSimMethod )
        case lower( 'BMA'  )
            Sim = M_gSet_maxDS*Agd'; 
            FuncSim = (Sim+Sim')./(Nd_col +Nd_col' + eps ) ;   
            %
%             FuncSim2 = zeros( ng, ng   ); 
%             for  ig = 1:ng
%                 for jg = 1:ng
%                     ds = DS(Agd(ig,:),Agd(jg,:)) ;
%                     FuncSim2(ig,jg) = ( sum( max(ds,[],2) )+sum( max(ds,[],1) ) )./( Nd_col(ig)+Nd_col(jg) + eps ) ;  
%                 end
%             end
            %   
        case lower( 'Max' )
            FuncSim = zeros( ng, ng   );
            for ig=1:ng
                FuncSim(:,ig) = sum(M_gSet_maxDS(:,Agd(ig,:) ) ,2);              
            end

        case lower(  'funSimMax' )
            Sim = M_gSet_maxDS*(Agd./Nd_col )'; 
            FuncSim = max(Sim, Sim') ;  

        case lower(  'funSimAvg' )
            Sim = M_gSet_maxDS*(Agd./Nd_col )'; 
            FuncSim = (Sim+Sim')./2 ;  

        otherwise
            error('No definition'); 

    end
    
end

% % % % % % % % % % % % % % % % % % % % 
if ExistIsolatedNodes
    FuncSim_final = eye( ng_all  ); 
    FuncSim_final(  nodes_g , nodes_g ) =  FuncSim ;      
    % % 
%     FuncSim2_final = eye( ng_all  ); 
%     FuncSim2_final(  nodes_g , nodes_g ) =  FuncSim2 ;      
% %     ss= (Sim+Sim');
% %     SSf =  zeros( ng_all  ); 
% %     SSf(  nodes_g , nodes_g ) =  ss ;      
else
    FuncSim_final  =  FuncSim ; 
%     FuncSim2_final =  FuncSim2 ;      
    
end

return;  
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [ kd,km ] = getGIPSim_IN(Adm_interaction, gamma0_d, gamma0_m,  AvoidIsolatedNodes, RemoveNonoverlapPairs   )
% Xiang  2019-11-16 
% Ref: van Laarhoven, Twan, Sander B. Nabuurs, and Elena Marchiori.
% "Gaussian interaction profile kernels for predicting drug–target interaction." 
% Bioinformatics 27, no. 21 (2011): 3036-3043.
%interaction: relation matrix between disease and miRNA, row:disease    column:miRNA  
% if ~exist('gamma0_d','var') || isempty(gamma0_d)
%    gamma0_d = 1 ;  
% end
% if ~exist('gamma0_m','var') || isempty(gamma0_m)
%    gamma0_m = 1 ;  
% end
% if ~exist('AvoidIsolatedNodes','var') || isempty(AvoidIsolatedNodes)
%    AvoidIsolatedNodes = 1 ;  
% end
% % % % % % % % % % % % % % % % % % % 
if isempty( gamma0_d  ) && isempty( gamma0_m  )
    error( 'both gamma0_d and gamma0_m are empty. No output.'   );
end
[nd_all, nm_all] = size( Adm_interaction  ); 
if AvoidIsolatedNodes
    nodes_d = sum(Adm_interaction,2)~=0 ; 
    nodes_m = sum(Adm_interaction,1)~=0 ; 
    Adm=double( Adm_interaction( nodes_d , nodes_m ) );
else
    Adm=double( Adm_interaction );
end
%
if ~exist( 'RemoveNonoverlapPairs','var'  ) || isempty( RemoveNonoverlapPairs  )
    RemoveNonoverlapPairs = true ;
end

% % % % % % % % % % % % % % % % % % % 
[nd, nm] = size( Adm  );  
SumOfSquares = sum(  Adm(:).^2 ); 
% % calculate gamad for Gaussian kernel calculation
if ~isempty( gamma0_d  )
    gamma_d = gamma0_d./( ( SumOfSquares ) / nd );
    %calculate Gaussian kernel for the similarity between disease: kd
    D   = Adm*Adm';
    dd  = diag( D  ); 
    kd2 = exp(-gamma_d*( dd+(dd'-2*D)   )  );  
    if RemoveNonoverlapPairs
        kd2(~D) = 0 ;        
    end 
    if AvoidIsolatedNodes  
        kd = zeros( nd_all  );  
        kd(  nodes_d , nodes_d ) =  kd2 ;       
    else
        kd = kd2;     
    end
    %
    D = [];
    kd2 = []; 
end

% % calculate gamam for Gaussian kernel calculation
if ~isempty( gamma0_m  )
    gamma_m = gamma0_m./( ( SumOfSquares ) / nm );
    %calculate Gaussian kernel for the similarity between miRNA: km
    E=Adm'*Adm;
    mm  = diag( E  ); 
    km2 = exp(-gamma_m*( mm+(mm'-2*E)   )  );  
    %
    if RemoveNonoverlapPairs
        km2(~E) = 0 ;        
    end    
    % 
    if AvoidIsolatedNodes   
        km = zeros( nm_all  );  
        km(  nodes_m , nodes_m ) =  km2 ;      
    else 
        km = km2;     
    end
    E = [];
    km2 = [];     
end

end





