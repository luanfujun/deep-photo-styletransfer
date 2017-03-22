function varargout=load_gaimc_graph(graphname)
% LOAD_GAIMC_GRAPH Loads a graph from the gaimc library
%
% load_gaimc_graph is a helper function to load a graph provided with the
% library regardless of the current working directory.  
%
% If it's called without any output arguments, it functions just like a
% load command executed on the .mat file with the graph.  If it's called
% with an output arguemnt, it functions just like a load command with
% output arguments.  It's somewhat complicated to explain because this is
% just a convinence function to make the examples work for any path, and
% not just from the gaimc root directory.
%
% Example:
%   % equivalent to load('graphs/airports.mat') run from the gaimc directory
%   load_gaimc_graph('airports') 
%   % equivalent to P=load('graphs/kt-7-2.mat') run from the gaimc directory
%   P=load_gaimc_graph('kt-7-2.mat') 
%   % so you don't have to put the path in for examples!

% David F. Gleich
% Copyright, Stanford University, 2008-2009

% History
% 2009-04-27: Initial coding


path=fileparts(mfilename('fullpath'));
if nargout==0
    evalin('caller',['load(''' fullfile(path,'graphs',graphname) ''');']);
else
    P = load(fullfile(path,'graphs',graphname));
    varargout{1} = P;
end
    
