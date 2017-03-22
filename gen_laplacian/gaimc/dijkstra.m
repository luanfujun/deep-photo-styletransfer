function [d pred]=dijkstra(A,u)
% DIJKSTRA Compute shortest paths using Dijkstra's algorithm
%
% d=dijkstra(A,u) computes the shortest path from vertex u to all nodes 
% reachable from vertex u using Dijkstra's algorithm for the problem.  
% The graph is given by the weighted sparse matrix A, where A(i,j) is 
% the distance between vertex i and j.  In the output vector d, 
% the entry d(v) is the minimum distance between vertex u and vertex v.
% A vertex w unreachable from u has d(w)=Inf.
%
% [d pred]=dijkstra(A,u) also returns the predecessor tree to generate 
% the actual shorest paths.  In the predecessor tree pred(v) is the 
% vertex preceeding v in the shortest path and pred(u)=0.  Any
% unreachable vertex has pred(w)=0 as well.  
%
% If your network is unweighted, then use bfs instead.
%
% See also BFS
%
% Example:
%   % Find the minimum travel time between Los Angeles (LAX) and
%   % Rochester Minnesota (RST).
%   load_gaimc_graph('airports')
%   A = -A; % fix funny encoding of airport data
%   lax=247; rst=355;
%   [d pred] = dijkstra(A,lax);
%   fprintf('Minimum time: %g\n',d(rst));
%   % Print the path
%   fprintf('Path:\n');
%   path =[]; u = rst; while (u ~= lax) path=[u path]; u=pred(u); end
%   fprintf('%s',labels{lax}); 
%   for i=path; fprintf(' --> %s', labels{i}); end, fprintf('\n');

% David F. Gleich
% Copyright, Stanford University, 2008-200909

% History
% 2008-04-09: Initial coding
% 2009-05-15: Documentation

if isstruct(A), 
    rp=A.rp; ci=A.ci; ai=A.ai; 
    check=0;
else
    [rp ci ai]=sparse_to_csr(A); check=1;
end
if check && any(ai)<0, error('gaimc:dijkstra', ...
        'dijkstra''s algorithm cannot handle negative edge weights.'); end

n=length(rp)-1; 
d=Inf*ones(n,1); T=zeros(n,1); L=zeros(n,1);
pred=zeros(1,length(rp)-1);

n=1; T(n)=u; L(u)=n; % oops, n is now the size of the heap

% enter the main dijkstra loop
d(u) = 0;
while n>0
    v=T(1); ntop=T(n); T(1)=ntop; L(ntop)=1; n=n-1; % pop the head off the heap
    k=1; kt=ntop;                   % move element T(1) down the heap
    while 1,
        i=2*k; 
        if i>n, break; end          % end of heap
        if i==n, it=T(i);           % only one child, so skip
        else                        % pick the smallest child
            lc=T(i); rc=T(i+1); it=lc;
            if d(rc)<d(lc), i=i+1; it=rc; end % right child is smaller
        end
        if d(kt)<d(it), break;     % at correct place, so end
        else T(k)=it; L(it)=k; T(i)=kt; L(kt)=i; k=i; % swap
        end
    end                             % end heap down
    
    % for each vertex adjacent to v, relax it
    for ei=rp(v):rp(v+1)-1            % ei is the edge index
        w=ci(ei); ew=ai(ei);          % w is the target, ew is the edge weight
        % relax edge (v,w,ew)
        if d(w)>d(v)+ew
            d(w)=d(v)+ew; pred(w)=v;
            % check if w is in the heap
            k=L(w); onlyup=0; 
            if k==0
                % element not in heap, only move the element up the heap
                n=n+1; T(n)=w; L(w)=n; k=n; kt=w; onlyup=1;
            else kt=T(k);
            end
            % update the heap, move the element down in the heap
            while 1 && ~onlyup,
                i=2*k; 
                if i>n, break; end          % end of heap
                if i==n, it=T(i);           % only one child, so skip
                else                        % pick the smallest child
                    lc=T(i); rc=T(i+1); it=lc;
                    if d(rc)<d(lc), i=i+1; it=rc; end % right child is smaller
                end
                if d(kt)<d(it), break;      % at correct place, so end
                else T(k)=it; L(it)=k; T(i)=kt; L(kt)=i; k=i; % swap
                end
            end
            % move the element up the heap
            j=k; tj=T(j);
            while j>1,                       % j==1 => element at top of heap
                j2=floor(j/2); tj2=T(j2);    % parent element
                if d(tj2)<d(tj), break;      % parent is smaller, so done
                else                         % parent is larger, so swap
                    T(j2)=tj; L(tj)=j2; T(j)=tj2; L(tj2)=j; j=j2;
                end
            end  
        end
    end
end
