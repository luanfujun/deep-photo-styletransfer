function [d dt ft pred] = dfs(A,u,full,target)
% DFS Compute depth first search distances, times, and tree for a graph
%
% [d dt ft pred] = dfs(A,u) returns the distance (d), the discover (dt) and
% finish time (ft) for each vertex in the graph in a depth first search 
% starting from vertex u.
%   d = dt(i) = ft(i) = -1 if vertex i is not reachable from u
% pred is the predecessor array.  pred(i) = 0 if vertex (i)  
% is in a component not reachable from u and i != u.
%
% [...] = dfs(A,u,1) continues the dfs for all components of the graph, not
% just the connected component with u
% [...] = dfs(A,u,[],v) stops the dfs when it hits the vertex v
%
% Note 1: When target is specified, the finish time array only records the
% finish time for all vertices that actually finished.  The array will then
% be undefined on a significant portion of the graph, but that doesn't
% indicate the vertices are unreachable; they just haven't been reached
% yet.
%
% Example:
%   load_gaimc_graph('dfs_example.mat') % use the dfs example from Boost
%   d = dfs(A,1)
%
% See also BFS

% David F. Gleich
% Copyright, Stanford University, 2008-2009

% History
% 2008-04-10: Initial coding

if ~exist('full','var') || isempty(full), full=0; end
if ~exist('target','var') || isempty(full), target=0; end

if isstruct(A), rp=A.rp; ci=A.ci; 
else [rp ci]=sparse_to_csr(A); 
end

n=length(rp)-1; 
d=-1*ones(n,1); dt=-1*ones(n,1); ft=-1*ones(n,1); pred=zeros(1,n);
rs=zeros(2*n,1); rss=0; % recursion stack holds two nums (v,ri)

% start dfs at u
t=0; targethit=0;
for i=1:n
    if i==1, v=u;
    else v=mod(u+i-1,n)+1; if d(v)>0, continue; end, end
    d(v)=0; dt(v)=t; t=t+1; ri=rp(v);
    rss=rss+1; rs(2*rss-1)=v; rs(2*rss)=ri; % add v to the stack
    while rss>0
        v=rs(2*rss-1); ri=rs(2*rss); rss=rss-1; % pop v from the stack
        if v==target || targethit, 
            ri=rp(v+1); targethit=1; % end the algorithm if v is the target
        end 
        while ri<rp(v+1)
            w=ci(ri); ri=ri+1;
            if d(w)<0
                d(w)=d(v)+1; pred(w)=v;
                rss=rss+1; rs(2*rss-1)=v; rs(2*rss)=ri; % add v to the stack
                v=w; ri=rp(w);
                dt(v)=t; t=t+1; continue; % discover a new vertex!
            end
        end
        ft(v)=t; t=t+1; % finish with v
    end
    if ~full, break; end
end

    
