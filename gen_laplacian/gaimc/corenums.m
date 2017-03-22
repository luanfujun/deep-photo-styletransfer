function [d rt]=corenums(A)
% CORENUMS Compute the core number for each vertex in the graph.
%
% [cn rt]=corenums(A) returns the core numbers for each vertex of the graph
% A along with the removal order of the vertex.  The core number is the 
% largest integer c such that vertex v exists in a graph where all 
% vertices have degree >= c.  The vector rt returns the removal time 
% for each vertex.  That is, vertex vi was removed at step rt[vi].
%
% This method works on directed graphs but gives the in-degree core number.
% To get the out-degree core numbers, call corenums(A').
%
% The linear algorithm comes from:
% Vladimir Batagelj and Matjaz Zaversnik, "An O(m) Algorithm for Cores 
% Decomposition of Networks."  Sept. 1 2002.
%
% Example:
%   load_gaimc_graph('cores_example'); % the graph A has three components
%   corenums(A)
%

% See also WCORENUMS

% David F. Gleich
% Copyright, Stanford University, 2008-2009

% History
% 2008-04-21: Initial Coding

if isstruct(A), rp=A.rp; ci=A.ci; %ofs=A.offset;
else [rp ci]=sparse_to_csr(A); 
end
n=length(rp)-1;
nz=length(ci);

% the algorithm removes vertices by computing bucket sort on the degrees of
% all vertices and removing the smallest.

% compute in-degrees and maximum indegree
d=zeros(n,1); maxd=0; rt=zeros(n,1);
for k=1:nz, newd=d(ci(k))+1; d(ci(k))=newd; if newd>maxd; maxd=newd; end, end

% compute the bucket sort
dp=zeros(maxd+2,1); vs=zeros(n,1); vi=zeros(n,1); % degree position, vertices
for i=1:n, dp(d(i)+2)=dp(d(i)+2)+1; end % plus 2 because degrees start at 0
dp=cumsum(dp); dp=dp+1;
for i=1:n, vs(dp(d(i)+1))=i; vi(i)=dp(d(i)+1); dp(d(i)+1)=dp(d(i)+1)+1; end
for i=maxd:-1:1, dp(i+1)=dp(i); end

% start the algorithm
t=1;
for i=1:n
    v = vs(i); dv = d(v); rt(v)=t; t=t+1;
    for rpi=rp(v):rp(v+1)-1
        w=ci(rpi); dw=d(w);
        if dw<=dv, % we already removed w
        else % need to remove edge (v,w), which decreases d(w)
            % swap w with the vertex at the head of its degree
            pw=vi(w); % get the position of w
            px=dp(dw+1); % get the pos of the vertex at the head of dw list
            x=vs(px);
            % swap w, x
            vs(pw)=x; vs(px)=w; vi(w)=px; vi(x)=pw;
            % decrement the degree of w and increment the start of dw
            d(w)=dw-1;
            dp(dw+1)=px+1;
        end
    end
end


    
