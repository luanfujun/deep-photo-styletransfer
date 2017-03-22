function [nzi,nzj,nzv] = csr_to_sparse(rp,ci,ai,ncols)
% CSR_TO_SPARSE Convert from compressed row arrays to a sparse matrix
%
% A = csr_to_sparse(rp,ci,ai) returns the sparse matrix represented by the
% compressed sparse row representation rp, ci, and ai.  The number of
% columns of the output sparse matrix is max(max(ci),nrows).  See the call
% below.
%
% A = csr_to_sparse(rp,ci,ai,ncol) While we can infer the number of rows 
% in the matrix from this expression, you may want a
% different number of 
%
% [nzi,nzj,nzv] = csr_to_sparse(...) returns the arrays that feed the
% sparse call in matlab.  You can use this to avoid the sparse call and
% customize the behavior.
%
% This command "inverts" the behavior of sparse_to_csr.
% Repeated entries in the matrix are summed, just like sparse does.  
%
% See also SPARSE SPARSE_TO_CSR
% 
% Example:
%   A=sparse(6,6); A(1,1)=5; A(1,5)=2; A(2,3)=-1; A(4,1)=1; A(5,6)=1; 
%   [rp ci ai]=sparse_to_csr(A); 
%   A2 = csr_to_sparse(rp,ci,ai)

% David F. Gleich
% Copyright, Stanford University, 2008-2009

%  History
%  2009-05-01: Initial version
%  2009-05-16: Documentation and example

nrows = length(rp)-1;
nzi = zeros(length(ci),1);
for i=1:nrows
    for j=rp(i):rp(i+1)-1
        nzi(j) = i;
    end
end

if nargout<2,
    if nargin>3,
        nzi = sparse(nzi,ci,ai,nrows,ncols);
    else
        % we make the matrix square unless there are more columns
        ncols = max(max(ci),nrows);
        if isempty(ncols), ncols=0; end
        nzi = sparse(nzi,ci,ai,nrows,ncols);
    end
else
    nzj = ci;
    nzv = ai;
end 
