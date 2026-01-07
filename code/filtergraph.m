function [g_simp, members, nodesize, D_simp] = filtergraph(g,d,varargin)
%FILTERGRAPH filter a graph g to produce a simpler one where nodes within
%"distance" d will be collapsed to a single node. Here we use "distance" in
%a loose sense, since "g" can be a directed graph, thus the "distance" or
%path length is not necessarily symmetric. 
%
%  [g_simp, members, nodesize, D_simp] = filtergraph(g,d,'reciprocal',1)
%
% input:
%   g: a graph or digraph (MATLAB object), to be simplified. 
%      Returned using the tknndigraph.m file
%   d: a threshold under which nodes of g will be collapsed to a new node
%   
% parameters:
%   reciprocal: whether we require both the path length from x to y and
%   from y to x to be under d in the filtering step. Default, true.
%
% output:
%   g_simp: the simplified graph (analogous to a Mapper shape graph)
%   members: members of the old graph (old nodes) in each new node
%   nodesize: # members in each new node
%   D_simp: simplified "distance" matrix. This gives the shortest distance
%   between any two group of "members" in the original graph.
%
% --------------------------- what is this ------------------------------
% Complimenting the philosophy behind a Mapper shape graph, this function
% simplifies a larger graph (directed or undirected) to a more compact one.
%
% This simplification takes three steps:
% 1. construct a filtered graph, whose connected components will be mapped
% to the nodes of the new graph. Specifically, this graph only retains
% connectivity between old node x and y, iff the shortest path length (pl)
% between x and y is less than a threshold "d" (we can either require
% pl(x,y) and pl(y,x) to be both less than d or either).
% 2. construct a new, simplified graph whose (new) nodes are
% associated with the connected components of the filtered graph from step 1
% 3. There is a link from new node x to new node y, iff there is at least
% one link from f^-1(x) to f^-1(y) in the old graph "g", and the average
% weight of links from f^-1(x) to f^-1(y) in the old graph "g" is the
% weight of the link in the simplified graph. 
% -----------------------------------------------------------------------
%{

Edited version of the one created by MZ on 9/12/2019

%}


p=inputParser();
p.addParameter('reciprocal',1);% if require pl(x,y)<0 and pl(y,x)<0 (pl = path length)
p.parse(varargin{:});
par = p.Results;

A = weightedAdj(g);% adjacency matrix
D = distances(g);% geodesic distance matrix

% -- connectivity within a distance threshold
if par.reciprocal
    A_ = (D < d) & (D' < d); % add edge iff mutual nbrs
else
    A_ = (D < d) | (D' < d); % add edge iff at least one is nbr of other
end
A_ = zerodiag(A_); % remove self-loop

% -- create graph out of nodes within said threshold
% create undirected graph corresponding to A_
if isfield(g.Nodes,'Name')
    g_ = graph(A_,g.Nodes.Name); % preserve node names if available
else 
    g_ = graph(A_);
end

% -- find connected components (define the new nodes)
% nodes in (undirected) g_ are connected if there exists a path b/w them
idx_newnodes = conncomp(g_);
if isfield(g.Nodes,'Name')
    [members, nodesize] = index2cell(idx_newnodes,g_.Nodes.Name);
else
    [members, nodesize] = index2cell(idx_newnodes);
end

% simplifyDistance() and simplifyAdj() are defined later in the script

% -- define distance between new nodes
D_simp = simplifyDistance(D,idx_newnodes); % shortest path between new nodes. 

% -- construct simplified graph
g_simp = digraph(simplifyAdj(A,idx_newnodes),'OmitSelfLoops');
end


% simplifyAdj() : Compresses the original adj matrix A into a simplified one that captures average connectivity between groups.

function A_simp = simplifyAdj(A, idx_newnodes)
% SIMPLIFYADJ simplifiy adjacency matrix A, such that A_simp(i,j) reflects
% the average connectivity between blocks of A.
    uidx = unique(idx_newnodes); % unique indices
    N_newnodes = length(uidx);
    
    A_simp = zeros(N_newnodes, N_newnodes); % initialize simplified distance matrix
    for n = 1:N_newnodes
        for m = 1:N_newnodes
            A_simp(n,m) = mean(mean(A(idx_newnodes==n, idx_newnodes==m)));% average connectivity between groups
        end
    end
end

% simplifyDistance() : compresses the original distance matrix D into a simplified one that captures the shortest topological connections between component groups. This fn captures proximity i.e. how "close" different regions are in the graph's geometry.


function D_simp = simplifyDistance(D, idx_newnodes)
% SIMPLIFYDISTANCE given a distance matrix D, and a node-assignment vector
% idx_newnodes. 
    uidx = unique(idx_newnodes);
    N_newnodes = length(uidx);
    
    D_simp = Inf(N_newnodes, N_newnodes); % initialize simplified distance matrix
    for n = 1:N_newnodes
        for m = 1:N_newnodes
            % D_simp(n,m) = min {dist(u,v) : u \in block n and v \in block m}
            D_simp(n,m) = min(min(D(idx_newnodes==n, idx_newnodes==m)));
        end
    end
end

function [members, nodesize] = index2cell(idx_newnodes,oldnodenames)
% INDEX2CELL convert a vector of labels of new nodes for each old nodes
% (idx_newnodes) to a cell array where each cell contains the names of the
% members of each new node. "nodesize" gives the size of the new nodes.
    if nargin<2 || isempty(oldenodenames)
        oldnodenames = (1:length(idx_newnodes))';
    end
    
    uidx = unique(idx_newnodes);% unique indices of new nodes
    Nidx = length(uidx);% number of new nodes
    
    members = cell(Nidx,1);% 1 cell = old nodes included in a new node
    nodesize = zeros(Nidx,1);% size of new node
    
    for n = 1:Nidx
        members{n} = oldnodenames(idx_newnodes == uidx(n));
        nodesize(n) = length(members{n});
    end
end

