function As = convert_sparse(A)
% CONVERT_SPARSE Convert a sparse matrix to the native gaimc representation
%
% As = convert_sparse(A) returns a struct with the three arrays defining
% the compressed sparse row structure of A.
%
% Example:
%   load('graphs/all_shortest_paths_example')
%   As = convert_sparse(A)
%
% See also SPARSE_TO_CSR SPARSE

% David F. Gleich
% Copyright, Stanford University, 2008-2009

% History
% 2009-04-29: Initial coding

[rp ci ai] = sparse_to_csr(A);
As.rp = rp;
As.ci = ci;
As.ai = ai;
