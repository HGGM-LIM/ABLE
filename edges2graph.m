function varargout = edges2graph(varargin);
%
% Syntax :
% [Graph] = edges2graph(edgeList, mapVector);
%
% This scripts creates a sparse Graph from an edge list. The Nodes will be
% the points inside the edge list and the connections between nodes is
% given by the edge list matrix
%
% Input Parameters:
%       edgeList                : Edge list (Nx2 matrix). 
%       bool                    : Boolean variable to create a logical
%                                 graph. All edges weights will be 1.
%
% Output Parameters:
%      Graph                    : Output Graph.
%      nodeList                 : List of the Nodes (vectorized, sorted 
%                                 and unified version of the edgeList matrix)
%      ordEdges                 : EdgeList according to the obtained Graph.
%                                 It represents the connections different from
%                                 0 in the obtained Graph. It connects the 
%                                 position of the node with the original 
%                                 edgeList matrix.
%                                 edgeList = nodeList(ordEdges);
% 
%                                 
%
% See also:  
%______________________________
% Authors: Yasser Aleman Gomez 
% Radiology and Psychiatric Departments, CHUV
% February 13th 2019
% Version $1.0

%% =================== Checking Input Parameters ======================= %%
if nargin == 0
    errormsg('Please enter a correct Surface struct');
    return;
end
edgeList = varargin{1};
%% =================== End of Checking Input Parameters ================ %%
%% ========================== Main Program ============================= %%

nodeList = unique(edgeList(:));
Nv = length(nodeList); % Graph dimension
Graph = sparse(Nv,Nv); % Creating an empty sparse matrix dimension

% Reordering the nodes list
[~,ordEdges] = ismember(edgeList,nodeList);


ind = sub2ind(size(Graph),[ordEdges(:,1);ordEdges(:,2)],[ordEdges(:,2);ordEdges(:,1)]); % Edges indexes

% Creating Edges weights
if nargin == 2
    mapVector = varargin{2};
    % Logical graph
    if  length(mapVector)== 1
        Graph(ind) = 1;
    elseif length(mapVector)== size(edgeList,1)
        
        % Mean value between each point and its neighbors
        
        Graph(ind) = [mapVector;mapVector] % Connection weight
    end
    
else
    Graph(ind) = 1;
    if nargin == 3
        
        % Graph(ind) = Distance; % Grafo de Distancia

    end
end
varargout{1} = Graph;
varargout{2} = nodeList;
varargout{3} = ordEdges;

%% ==========================End of Main Program ============================= %%
return