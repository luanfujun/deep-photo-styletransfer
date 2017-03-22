function [val m1 m2 mi]=bipartite_matching(varargin)
% BIPARTITE_MATCHING Solve a maximum weight bipartite matching problem
%
% [val m1 m2]=bipartite_matching(A) for a rectangular matrix A 
% [val m1 m2 mi]=bipartite_matching(x,ei,ej,n,m) for a matrix stored
% in triplet format.  This call also returns a matching indicator mi so
% that val = x'*mi.
%
% The maximum weight bipartite matching problem tries to pick out elements
% from A such that each row and column get only a single non-zero but the
% sum of all the chosen elements is as large as possible.
%
% This function is slightly atypical for a graph library, because it will
% be primarily used on rectangular inputs.  However, these rectangular
% inputs model bipartite graphs and we take advantage of that stucture in
% this code.  The underlying graph adjency matrix is 
%   G = spaugment(A,0); 
% where A is the rectangular input to the bipartite_matching function.
%
% Matlab already has the dmperm function that computes a maximum
% cardinality matching between the rows and the columns.  This function
% gives us the maximum weight matching instead.  For unweighted graphs, the
% two functions are equivalent.
%
% Note: If ei and ej contain duplicate edges, the results of this function
% are incorrect.
%
% See also DMPERM
%
% Example:
%   A = rand(10,8); % bipartite matching between random data
%   [val mi mj] = bipartite_matching(A);
%   val

% David F. Gleich and Ying Wang
% Copyright, Stanford University, 2008-2009
% Computational Approaches to Digital Stewardship

% 2008-04-24: Initial coding (copy from Ying Wang matching_sparse_mex.cpp)
% 2008-11-15: Added triplet input/output
% 2009-04-30: Modified for gaimc library
% 2009-05-15: Fixed error with empty inputs and triple added example.

[rp ci ai tripi n m] = bipartite_matching_setup(varargin{:});

if isempty(tripi)
    error(nargoutchk(0,3,nargout,'struct'));
else    
    error(nargoutchk(0,4,nargout,'struct'));
end


if ~isempty(tripi) && nargout>3
    [val m1 m2 mi] = bipartite_matching_primal_dual(rp, ci, ai, tripi, n, m);
else
    [val m1 m2] = bipartite_matching_primal_dual(rp, ci, ai, tripi, n, m);
end

function [rp ci ai tripi n m]= bipartite_matching_setup(A,ei,ej,n,m)
% convert the input

if nargin == 1
    if isstruct(A)
        [nzi nzj nzv]=csr_to_sparse(A.rp,A.ci,A.ai);
    else
        [nzi nzj nzv]=find(A); 
    end
    [n m]=size(A);
    triplet = 0;
elseif nargin >= 3 && nargin <= 5    
    nzi = ei;
    nzj = ej;
    nzv = A;
    if ~exist('n','var') || isempty(n), n = max(nzi); end
    if ~exist('m','var') || isempty(m), m = max(nzj); end
    triplet = 1;
else    
    error(nargchk(3,5,nargin,'struct'));
end
nedges = length(nzi);

rp = ones(n+1,1); % csr matrix with extra edges
ci = zeros(nedges+n,1);
ai = zeros(nedges+n,1);
if triplet, tripi = zeros(nedges+n,1); % triplet index
else tripi = [];
end

%
% 1. build csr representation with a set of extra edges from vertex i to
% vertex m+i
%
rp(1)=0;
for i=1:nedges
    rp(nzi(i)+1)=rp(nzi(i)+1)+1;
end
rp=cumsum(rp); 
for i=1:nedges
    if triplet, tripi(rp(nzi(i))+1)=i; end % triplet index
    ai(rp(nzi(i))+1)=nzv(i);
    ci(rp(nzi(i))+1)=nzj(i);
    rp(nzi(i))=rp(nzi(i))+1;
end
for i=1:n % add the extra edges
    if triplet, tripi(rp(i)+1)=-1; end % triplet index
    ai(rp(i)+1)=0;
    ci(rp(i)+1)=m+i;
    rp(i)=rp(i)+1;
end
% restore the row pointer array
for i=n:-1:1
    rp(i+1)=rp(i);
end
rp(1)=0;
rp=rp+1;

%
% 1a. check for duplicates in the data
%
colind = false(m+n,1);
for i=1:n
    for rpi=rp(i):rp(i+1)-1
        if colind(ci(rpi)), error('bipartite_matching:duplicateEdge',...
            'duplicate edge detected (%i,%i)',i,ci(rpi)); 
        end
        colind(ci(rpi))=1;
    end
    for rpi=rp(i):rp(i+1)-1, colind(ci(rpi))=0; end % reset indicator
end


function [val m1 m2 mi]=bipartite_matching_primal_dual(...
                            rp, ci, ai, tripi, n, m)
% BIPARTITE_MATCHING_PRIMAL_DUAL                         

alpha=zeros(n,1); % variables used for the primal-dual algorithm
beta=zeros(n+m,1);
queue=zeros(n,1);
t=zeros(n+m,1);
match1=zeros(n,1);
match2=zeros(n+m,1);
tmod = zeros(n+m,1);
ntmod=0;


% 
% initialize the primal and dual variables
%
for i=1:n
    for rpi=rp(i):rp(i+1)-1
        if ai(rpi) > alpha(i), alpha(i)=ai(rpi); end
    end
end
% dual variables (beta) are initialized to 0 already
% match1 and match2 are both 0, which indicates no matches
i=1;
while i<=n
    % repeat the problem for n stages
    
    % clear t(j)
    for j=1:ntmod, t(tmod(j))=0; end
    ntmod=0;
    

    % add i to the stack
    head=1; tail=1;
    queue(head)=i; % add i to the head of the queue
    while head <= tail && match1(i)==0
        k=queue(head);
        for rpi=rp(k):rp(k+1)-1
            j = ci(rpi);
            if ai(rpi) < alpha(k)+beta(j) - 1e-8, continue; end % skip if tight
            if t(j)==0,
                tail=tail+1; queue(tail)=match2(j);
                t(j)=k;
                ntmod=ntmod+1; tmod(ntmod)=j;
                if match2(j)<1,
                    while j>0, 
                        match2(j)=t(j);
                        k=t(j);
                        temp=match1(k);
                        match1(k)=j;
                        j=temp;
                    end
                    break; % we found an alternating path
                end
            end
        end
        head=head+1;
    end
    
    if match1(i) < 1, % still not matched, so update primal, dual and repeat
        theta=inf;
        for j=1:head-1
            t1=queue(j);
            for rpi=rp(t1):rp(t1+1)-1
                t2=ci(rpi);
                if t(t2) == 0 && alpha(t1) + beta(t2) - ai(rpi) < theta,
                    theta = alpha(t1) + beta(t2) - ai(rpi);
                end
            end
        end
        
        for j=1:head-1, alpha(queue(j)) = alpha(queue(j)) - theta; end
        
        for j=1:ntmod, beta(tmod(j)) = beta(tmod(j)) + theta; end
            
        continue;
    end
        
    i=i+1; % increment i
end

val=0;
for i=1:n
    for rpi=rp(i):rp(i+1)-1
        if ci(rpi)==match1(i), val=val+ai(rpi); end
    end
end
noute = 0; % count number of output edges
for i=1:n
    if match1(i)<=m, noute=noute+1; end
end
m1=zeros(noute,1); m2=m1; % copy over the 0 array
noute=1;
for i=1:n
    if match1(i)<=m, m1(noute)=i; m2(noute)=match1(i);noute=noute+1; end
end

if nargout>3
    mi= false(length(tripi)-n,1);
    for i=1:n
        for rpi=rp(i):rp(i+1)-1
            if match1(i)<=m && ci(rpi)==match1(i), mi(tripi(rpi))=1; end
        end
    end
end



