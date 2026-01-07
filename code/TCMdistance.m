function [TCM_pl] = TCMdistance(g,nodet,weighted)
%TCMDISTANCE calculate a time connectivity matrix where the edge is the
%minimal path length between two sample points.
%   [TCM_pl] = TCMdistance(g,nodet,weighted)
%   input:
%       g: matlab graph, particularly simplified graph returned by filtergraph()
%       nodet: members of each node as indices of time points, i.e. "members" cell array returned 
%              by filtergraph()
%       weighted: whether to use the weight of the graph, default = false,
%       use as unweighted graph.
%   output:
%       TCM_pl: temporal connectivity using shortest path length on the
%       shape graph.
%{
created by MZ, 04-27-2019
modifications:
(8-23-2019) debug weight assignment -- no longer require input graph to be
weighted. also accommodate situation where t0~=0.
%}


%% Summary:
% Rows and columns represent time points (extracted from nodet)
% TCM_pl(t1, t2) = connectivity between time point t1 and time point t2. Hence the name "temporal" connectivity
% However, connectivity measures in distmat are based on the spatial graph.
% What is the spatial proximity between two time points, as measured through the shape graph?

% weighted parameter is optional. If not provided, set to false
if nargin<3 || isempty(weighted)
    weighted = false;
end

% if no weights provided, then use uniform weights
if ~weighted
    g.Edges.Weight=ones(height(g.Edges),1);
end

% -- extract time
t = unique(cell2mat(nodet)); % num of unique time points in all node clusters of "members"
Nt = length(t);
t_0 = min(t);

% -- construct TCM from graph distance
TCM_pl = NaN(Nt,Nt);  % initialize connectivity matrix

% store shortest path length between all pairs of nodes in distmat
% if no path exists, it returns Inf
distmat = distances(g);

% -- assign diagonal elements to 0 
% i.e. distance from any point to itself is 0
for i=1:g.numnodes
        TCM_pl(nodet{i}-t_0+1,nodet{i}-t_0+1) = 0;
end

% -- assign off diagonal elements
for i=1:g.numnodes
    for j=i+1:g.numnodes  % iterate through all pairs of nodes (i,j)
        TCM_pl(nodet{i}-t_0+1,nodet{j}-t_0+1) = nanmin(TCM_pl(nodet{i}-t_0+1,nodet{j}-t_0+1),distmat(i,j));
        TCM_pl(nodet{j}-t_0+1,nodet{i}-t_0+1) = nanmin(TCM_pl(nodet{j}-t_0+1,nodet{i}-t_0+1),distmat(j,i));
    end
end


end

