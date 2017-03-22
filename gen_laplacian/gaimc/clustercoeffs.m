function cc=clustercoeffs(A,weighted,normalized)
% CLUSTERCOEFFS Compute undirected clustering coefficients for a graph
%
% ccfs=clustercoeffs(A) compute normalized, weighted clustering
% coefficients from a graph represented by a symmetric adjacency matrix A.
%
% ccfs=clustering(A,weighted,normalized) takes optional parameters to
% control normalization and weighted computation.  If normalized=0 or
% false, then the computation is not normalized by d*(d-1) in the
% unweighted case.  If weighted=0 or false, then the weights of the matrix
% A are ignored.  Either parameter will assume it's default value if you 
% specify an empty matrix.
%
% See also DIRCLUSTERCOEFFS
%
% Example:
%   load_gaimc_graph('clique-10');
%   cc = clustercoeffs(A) % they are all equal! as we expect in a clique

% David F. Gleich
% Copyright, Stanford University, 2008-2009

% History
% 2009-05-15: First history comment, originally written in 2008

if ~exist('normalized','var') || isempty(normalized), normalized=true; end
if ~exist('weighted','var') || isempty(weighted), weighted=true; end
donorm=1; usew=1;
if ~normalized, donorm=0; end
if ~weighted, usew=0; end

if isstruct(A)
    rp=A.rp; ci=A.ci; %ofs=A.offset;
    if usew, ai=A.ai; end
else
    if ~isequal(A,A'), error('gaimc:clustercoeffs',...
            'only undirected (symmetric) inputs allowed: see dirclustercoeffs');
    end
    if usew, [rp ci ai]=sparse_to_csr(A); 
    else [rp ci]=sparse_to_csr(A); 
    end
    if any(ai)<0, error('gaimc:clustercoeffs',...
            ['only positive edge weights allowed\n' ...
             'try clustercoeffs(A,0) for an unweighted comptuation']); 
    end
end
n=length(rp)-1;

cc=zeros(n,1); ind=false(n,1); cache=zeros(n,1); 
ew=1; ew2=1;
for v=1:n
    for rpi=rp(v):rp(v+1)-1
        w=ci(rpi); if usew, ew=ai(rpi); end
        if v~=w, ind(w)=1; cache(w)=ew^(1/3); end
    end
    curcc=0; d=rp(v+1)-rp(v);
    % run two steps of bfs to try and find triangles. 
    for rpi=rp(v):rp(v+1)-1
        w=ci(rpi); if v==w, d=d-1; continue; end % discount self-loop
        for rpi2=rp(w):rp(w+1)-1
            x=ci(rpi2); if x==w, continue; end
            if ind(x)
                if usew, ew=ai(rpi); ew2=ai(rpi2); end
                curcc=curcc+ew^(1/3)*ew2^(1/3)*cache(x);
            end
        end
    end
    if donorm&&d>1, cc(v)=curcc/(d*(d-1)); 
    elseif d>1, cc(v)=curcc;
    end
    for rpi=rp(v):rp(v+1)-1, w=ci(rpi); ind(w)=0; end % reset indicator
end


