function [sci sizes] = scomponents(A)
% SCOMPONENTS Compute the strongly connected components of a graph
%
% ci=scomponents(A) returns an index for the component number of every 
% vertex in the graph A.  The total number of components is max(ci).
% If the input is undirected, then this algorithm outputs just the 
% connected components.  Otherwise, it output the strongly connected
% components.
%
% The implementation is from Tarjan's 1972 paper: Depth-first search and 
% linear graph algorithms. In SIAM's Journal of Computing, 1972, 1, 
% pp.146-160.
%
% See also DMPERM
%
% Example:
%   load_gaimc_graph('cores_example'); % the graph A has three components
%   ci = scomponents(A)
%   ncomp = max(ci)               % should be 3
%   R = sparse(1:size(A,1),ci,1,size(A,1),ncomp); % create a restriction matrix
%   CG = R'*A*R;                  % create the graph with each component 
%                                 % collapsed into a single node.

% David F. Gleich
% Copyright, Stanford University, 2008-2009

% History
% 2008-04-11: Initial coding

if isstruct(A), rp=A.rp; ci=A.ci; %ofs=A.offset;
else [rp ci]=sparse_to_csr(A); 
end

n=length(rp)-1; sci=zeros(n,1); cn=1;
root=zeros(n,1); dt=zeros(n,1); t=0;
cs=zeros(n,1); css=0; % component stack
rs=zeros(2*n,1); rss=0; % recursion stack holds two nums (v,ri)
% start dfs at 1
for sv=1:n
    v=sv; if root(v)>0, continue; end
    rss=rss+1; rs(2*rss-1)=v; rs(2*rss)=rp(v); % add v to the stack
    root(v)=v; sci(v)=-1; dt(v)=t; t=t+1;
    css=css+1; cs(css)=v; % add w to component stack
    while rss>0
        v=rs(2*rss-1); ri=rs(2*rss); rss=rss-1; % pop v from the stack
        while ri<rp(v+1)
            w=ci(ri); ri=ri+1;
            if root(w)==0
                root(w)=w; sci(w)=-1;  dt(w)=t; t=t+1;
                css=css+1; cs(css)=w; % add w to component stack
                rss=rss+1; rs(2*rss-1)=v; rs(2*rss)=ri; % add v to the stack
                v=w; ri=rp(w); continue; 
            end
        end
        for ri=rp(v):rp(v+1)-1
            w=ci(ri); 
            if sci(w)==-1
                if dt(root(v))>dt(root(w)), root(v)=root(w); end
            end
        end
        if root(v)==v
            while css>0
                w=cs(css); css=css-1; sci(w)=cn;
                if w==v, break; end
            end
            cn=cn+1;
        end
    end
end

if nargout>1
    sizes=accumarray(sci,1,[max(sci) 1]);
end
