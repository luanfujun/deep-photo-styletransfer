function [d dt pred] = bfs(A,u,target)
% BFS Compute breadth first search distances, times, and tree for a graph
%
% [d dt pred] = bfs(A,u) returns the distance (d) and the discover time
% (dt) for each vertex in the graph in a breadth first search 
% starting from vertex u.
%   d = dt(i) = -1 if vertex i is not reachable from u
% pred is the predecessor array.  pred(i) = 0 if vertex (i)  
% is in a component not reachable from u and i != u.
%
% [...] = bfs(A,u,v) stops the bfs when it hits the vertex v
%
% Example:
%   load_gaimc_graph('bfs_example.mat') % use the dfs example from Boost
%   d = bfs(A,1)
%
% See also DFS

% David F. Gleich
% Copyright, Stanford University, 2008-20098

% History
% 2008-04-13: Initial coding

if ~exist('target','var') || isempty(full), target=0; end

if isstruct(A), rp=A.rp; ci=A.ci; 
else [rp ci]=sparse_to_csr(A); 
end

n=length(rp)-1; 
d=-1*ones(n,1); dt=-1*ones(n,1); pred=zeros(1,n);
sq=zeros(n,1); sqt=0; sqh=0; % search queue and search queue tail/head

% start bfs at u
sqt=sqt+1; sq(sqt)=u; 
t=0;  
d(u)=0; dt(u)=t; t=t+1; pred(u)=u;
while sqt-sqh>0
    sqh=sqh+1; v=sq(sqh); % pop v off the head of the queue
    for ri=rp(v):rp(v+1)-1
        w=ci(ri);
        if d(w)<0
            sqt=sqt+1; sq(sqt)=w; 
            d(w)=d(v)+1; dt(w)=t; t=t+1; pred(w)=v; 
            if w==target, return; end
        end
    end
end
