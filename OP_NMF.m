function [W H] = opnmf(X, K, w0, initMeth, max_iter, tol, ...
    iter0, save_step, outputdir)

[Dinit,N] = size(X);

% Basic argument check
if ~isscalar(K) || ~isnumeric(K) || K<1 || K>min(Dinit,N) || K~=round(K)
    error('opnmf:badK','K should be positive integer no larger than the number of rows or columns in X');
end
if ~ismatrix(X) || ~isnumeric(X) || ~isreal(X) || any(any(X<0)) || any(any(~isfinite(X)))
    error('opnmf:badX','opnmf:X must be a matrix of non-negative values.')
end

% get rid of the background in order to boost the computational efficiency
% assumption: background corresponds to positions with stackwise mean value
% equal to zero
mean_im = mean(X,2);
data_matrix_nz = X((mean_im>0),:) ; 
X = data_matrix_nz ; clear data_matrix_nz ;

% variables
check_step = 100;
D = size(X,1);

% initialize w0
if ~exist('w0','var') || isempty(w0)
     %%disp('Initializing w0: ');
    if ~exist('initMeth', 'var') || isempty(initMeth)
        initMeth = 1;
    end
    switch initMeth
        case 0
             %%disp('random initialization ...') ;
            W = rand(D,K);
        case 1 % NNDSVD
             %%disp('NNDSVD initialization ...') ;
            [W,~] = NNDSVD(X,K,0) ;
        case 2 % NNDSVDa
             %%disp('NNDSVDa initialization ...') ;
            [W,~] = NNDSVD(X,K,1) ;
        case 3 % NNDSVDar
             %%disp('NNDSVDar initialization ...') ;
            [W,~] = NNDSVD(X,K,2) ;
        case 4 % NNDSVD using randomized SVD
             %%disp('NNDSVD initialization using random SVD calculation for efficiency ...') ;
            [W,~] = NNDSVD(X,K,3) ;
        otherwise
             %%disp('NNDSVD initialization ...') ;
            [W,~] = NNDSVD(X,K,0) ;
    end
     %%disp('done') ;
else
    W = w0;
    clear w0 ;
end

% check variables
if ~exist('max_iter', 'var') || isempty(max_iter)
    max_iter = 50000;
end
if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-5;
end
if (~exist('iter0','var') || isempty(iter0))
    iter0 = 1;
end
if (~exist('save_step','var') || isempty(save_step) )
    save_step = floor(max_iter/10) ;
end
if (exist('outputdir','var') && ~isempty(outputdir))    
    % then create the directory
    if(~strcmp(outputdir(end),'/'))
        outputdir=[outputdir '/'];
    end
    if(~exist(outputdir,'dir'))
        success = mkdir(outputdir);
        if(~success)
            error('opnmf:BadDir',['Output directory ' outputdir ' can not be created']);
        end
    end
end

% start optimization
XX = X * X';
for iter=iter0:max_iter
    W_old = W;
    if mod(iter,check_step)==0
         %%fprintf('iter=% 5d ', iter);
    end
    
    % multiplicative update rule
    W = W .* (XX*W) ./ (W*(W'*XX*W));
    W = W ./ norm(W);
    
    diffW = norm(W_old-W, 'fro') / norm(W_old, 'fro');
    if diffW<tol
         %%fprintf('converged after %d steps.\n', iter);
        break;
    end
    
    if mod(iter,check_step)==0
         %%fprintf('diff=%.10f, ', diffW);
         %%fprintf('obj=%.10f', norm(X-W*(W'*X), 'fro'));
         %%fprintf('\n');
    end
    
    % save intermediate results
    if ( exist('outputdir','var') && ~isempty(outputdir) )
        if ( mod(iter,save_step) == 0 )
             %%fprintf('Saving intermediate results ...');
            save([outputdir 'IntermResultsExtractBases.mat'],'iter','D','K','W','-v7.3') ;
             %%fprintf('done\n') ;
        end
    end
end

% Reordering the output - putting them in standard form, loosely following
% what is done in the nnmf matlab function
H = W'*X ;

hlen = sqrt(sum(H.^2,2));
if any(hlen==0)
    warning(message('opnmf:LowRank', K - sum( hlen==0 ), K));
    hlen(hlen==0) = 1;
end
clear H
Wh = bsxfun(@times,W,hlen');

% Then order by W
[~,idx] = sort(sum(Wh.^2,1),'descend'); clear Wh
W = W(:,idx);
H = W'*X ;

% put results to original dimension
WW = zeros(Dinit,K) ;
WW(mean_im>0,:) = W ; clear W;
W = WW ; clear WW
   end

function [val, idx] = yael_kmin (v, k)

[val, idx] = sort (v);

val = val (1:k, :);
idx = idx (1:k, :);
end

function [W,H] = NNDSVD(A,k,flag)
%
%----------------------check the input matrix------------------------------
if numel(find(A<0)) > 0
    error('The input matrix contains negative elements !')
end
%--------------------------------------------------------------------------

%size of the input matrix
[m,n] = size(A);

%the matrices of the factorization
W = zeros(m,k);
H = zeros(k,n);

% 1st SVD --> partial SVD rank-k to the input matrix A.
if ( flag == 3)
    % use random svd for efficient computation
    l = max(3*k,20) ;
    [U,S,V]=randpca(A,k,true,8,l);
else
    % use standard matlab svn implementation
    [U,S,V] = svds(A,k);
end

%choose the first singular triplet to be nonnegative
W(:,1)     =  sqrt(S(1,1)) * abs(U(:,1) );         
H(1,:)     =  sqrt(S(1,1)) * abs(V(:,1)'); 


% 2nd SVD for the other factors (see table 1 in our paper)
for i=2:k
    uu = U(:,i); vv = V(:,i);
    uup = pos(uu); uun = neg(uu) ;
    vvp = pos(vv); vvn = neg(vv);
    n_uup = norm(uup);
    n_vvp = norm(vvp) ;
    n_uun = norm(uun) ;
    n_vvn = norm(vvn) ;
    termp = n_uup*n_vvp; termn = n_uun*n_vvn;
    if (termp >= termn)
        W(:,i) = sqrt(S(i,i)*termp)*uup/n_uup; 
        H(i,:) = sqrt(S(i,i)*termp)*vvp'/n_vvp;
    else
        W(:,i) = sqrt(S(i,i)*termn)*uun/n_uun; 
        H(i,:) = sqrt(S(i,i)*termn)*vvn'/n_vvn;
    end
end
%------------------------------------------------------------

%actually these numbers are zeros
W(find(W<0.0000000001))=0.1;
H(find(H<0.0000000001))=0.1;

if(exist('flag','var'))
    % NNDSVDa: fill in the zero elements with the average
    if flag==1
        ind1      =  find(W==0) ;
        ind2      =  find(H==0) ;
        average   =  mean(A(:)) ;
        W( ind1 ) =  average    ;
        H( ind2 ) =  average    ;
        
        % NNDSVDar: fill in the zero elements with random values in the space [0:average/100]
    elseif flag==2
        ind1      =  find(W==0) ;
        ind2      =  find(H==0) ;
        n1        =  numel(ind1);
        n2        =  numel(ind2);
        
        average   =  mean(A(:))       ;
        W( ind1 ) =  (average*rand(n1,1)./100)  ;
        H( ind2 ) =  (average*rand(n2,1)./100)  ;
    end
end
end
%--------------------------------------------------------------------------
%end of the nndsvd function



%This function sets to zero the negative elements of a matrix
%--------------------------------------------------------------------------
function [Ap] = pos(A)
Ap = (A>=0).*A;
end
%--------------------------------------------------------------------------

%This functions sets to zero the positive elements of a matrix and takes
%the absolute value of the negative elements
%--------------------------------------------------------------------------
function [Am] = neg(A);
Am = (A<0).*(-A);
end

            
   
function ri = rand_index(p1, p2, varargin)

    % Parse the input and throw errors
    adj = 0;
    if nargin == 0
    end
    if nargin > 3
        error('Too many input arguments');
    end
    if nargin == 3
        if strcmp(varargin{1}, 'adjusted')
            adj = 1;
        else
            error('%s is an unrecognized argument.', varargin{1});
        end
    end
    if length(p1)~=length(p2)
        error('Both partitions must contain the same number of points.');
    end
    
	% Preliminary computations and cleansing of the partitions
    N = length(p1);
    [~, ~, p1] = unique(p1);
    N1 = max(p1);
    [~, ~, p2] = unique(p2);
    N2 = max(p2);
    
    % Create the matching matrix
    for i=1:1:N1
        for j=1:1:N2
            G1 = find(p1==i);
            G2 = find(p2==j);
            n(i,j) = length(intersect(G1,G2));
        end
    end
    
    % If required, calculate the basic rand index
    if adj==0
        ss = sum(sum(n.^2));
        ss1 = sum(sum(n,1).^2);
        ss2 =sum(sum(n,2).^2);
        ri = (nchoosek2(N,2) + ss - 0.5*ss1 - 0.5*ss2)/nchoosek2(N,2);
    end
    
    
    % Otherwise, calculate the adjusted rand index
    if adj==1
        ssm = 0;
        sm1 = 0;
        sm2 = 0;
        for i=1:1:N1
            for j=1:1:N2
                ssm = ssm + nchoosek2(n(i,j),2);
            end
        end
        temp = sum(n,2);
        for i=1:1:N1
            sm1 = sm1 + nchoosek2(temp(i),2);
        end
        temp = sum(n,1);
        for i=1:1:N2
            sm2 = sm2 + nchoosek2(temp(i),2);
        end
        NN = ssm - sm1*sm2/nchoosek2(N,2);
        DD = (sm1 + sm2)/2 - sm1*sm2/nchoosek2(N,2);
        ri = NN/DD;
    end 
    

    % Special definition of n choose k
    function c = nchoosek2(a,b)
        if a>1
            c = nchoosek(a,b);
        else
            c = 0;
        end
    end
end

function [val, idx] = yael_kmax (v, k)

[val, idx] = sort (v, 'descend');

val = val (1:k, :);
idx = idx (1:k, :);
end