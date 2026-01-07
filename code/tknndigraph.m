function g = tknndigraph(XorD,k,tidx,varargin)
%TKNNDIGRAPH construct a directed graph based on k-nearest neighbors which
%much include its temporal neighbors. Here we do so by simply set the
%distance between consecutive time points to zero before running knn.
%   g = tknngraph(XorD,k,tidx)
% input:
%   XorD: a N-by-d matrix (X) of the coordinates of N points in d-sim
%       space, or a N-by-N distance matrix (D). If the input is X, D is
%       computed from X based on Euclidean distances. 
%   k: # nearest neighbors
%   tidx: a vector of N integers, two points x, y are considered temporal
%   neighbors iff tidx[x]+1 = tidx[y] or tidx[x]-1 = tidx[y].
% output:
%   g: matlab graph object (unweighted, undirected).
%{
created by MZ, 8-16-2019
modifcations:
(8-20-2019) add option to not enforce reciprocity.
(2-5-2020) add option to determine whether temporal link can be a spatial
link. parameter: timeExcludeSpace

%}
p = inputParser;
p.addParameter('reciprocal',true)% spatially reciprocal
p.addParameter('timeExcludeSpace',true)% whether temporal links are allowed to be spatial links
p.parse(varargin{:})
par = p.Results;

% -- check input and obtain distance matrix D
[nr,nc]=size(XorD);
if nr~=nc || any(any(XorD~=XorD'))
    D = pdist2(XorD,XorD); % pairwise L2 dist
else
    D = XorD;
end
Nn = length(D); % number of nodes

D(logical(eye(Nn))) = Inf; % exclude self-loops i.e. ensures a point can't be its own neighbor

% -- find indices for temporal links D_{i(t),i(t+1)}
% t_wafter(i) = true if i has a temporal neighbor i.e. next time point exists
t_wafter = circshift(tidx,-1,1) - 1 == tidx; 
% t_after_idx(i, i+1) = 1 if there is edge from time t to time t+1
t_after_idx = circshift(diag(t_wafter),1,2);
% If timeExcludeSpace=true then set distance b/w consecutive time points to 0.
if par.timeExcludeSpace
    D(t_after_idx) = 0;
end

% -- compute adjacency matrix
% Example:                         
% D =   [Inf  5   2   8]         Ic =      [3  2  4  1]  ← for point 1: nearest is 3, then 2
%       [5  Inf  3   1]      →             [4  3  1  2]  ← for point 2: nearest is 4, then 3
%       [2   3  Inf  6]                    [1  2  4  3]  ← for point 3: nearest is 1, then 2
%       [8   1   6  Inf]                   [2  3  1  4]  ← for point 4: nearest is 2, then 3
%

A = zeros(Nn,Nn); % initialize adj matrix
% sort rows in ascending order; only store indices Ic of sorted values
[~,Ic]=sort(D,2); 

% A is computed. Add edge to k nearest neighbors (spatial)
%      1   2   3   4
% 1 [  0   1   1   0 ]  ← point 1 connects to 2, 3
% 2 [  0   0   1   1 ]  ← point 2 connects to 3, 4
% 3 [  1   1   0   0 ]  ← point 3 connects to 1, 2
% 4 [  0   1   1   0 ]  ← point 4 connects to 2, 3
%
I = sub2ind([Nn Nn], repmat((1:Nn)',1,k), Ic(:,1:k));
A(I(:))=1;

% -- exclude or retain temporal links as spatial links 
if par.timeExcludeSpace
    A_space = A.* (~t_after_idx); % remove temporal links
else
    A_space = A;
end

% -- enforce symmetry of spatial links
if par.reciprocal
    A_space = A_space & A_space'; % edge exists if they are mutual nbrs
else
    A_space = A_space | A_space'; % at least one is a nbr of the other
end

% -- (re-)incoporate temporal links
A = t_after_idx | A_space; 

% -- convert to graph
g = digraph(A);
end

