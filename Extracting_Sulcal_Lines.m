function varargout = Extracting_Sulcal_Lines(varargin)
%
% Syntax :
%    [skelEdges, Graph] = Extracting_Sulcal_Topology(Surf, opts);
%
% This scripts extracts sulcal skeleton from the distance transform map computed for
% a surface.
%
% Input Parameters:
%        Surf                           : Surface variable (file, struct or cellarray).
%                                         number of points in its corresponding surface).
%        opts (optional input)          : Options:
%                                         opts.curvthr - Curvature threshold for defining sulcal regions
%                                         opts.dephtmap - Depth Map
%                                         opts.endpthr - Number of Endpoints threshold.
%                                          -   opts.endpthr = 0; % Means main Sulcal Line
%                                          -   opts.endpthr = 1; % All Sulcal Lines
%                                          -   opts.endpthr > 2; % opts.endpthr - 2 Secondary lines + Main Path
%
% Output Parameters:
%        skelEdges                      : Skeleton edges.
%        Graph                          : Skeleton Graph.
%
% See also: Remove_Branches Sulcal_Lines_Extraction Computing_Sulcal_Basins Surface_Checking
% Recursively_Remove_Branches Reordering_Gyral_Crowns
%__________________________________________________
% Authors: Yasser Aleman Gomez and Javier Santoja
% LIM, HUGGM
% November 13th 2014
% Version $1.0

%% ============================= Checking Inputs ======================= %%
if nargin < 2
    error('Two Inputs are mandatory');
end
Surf = varargin{1};
Pial = varargin{2};
curvmap = varargin{3};

% Surface Checking
Surf = Surface_Checking(Surf);

if ischar(curvmap) % Verifying if the curvature is a file
    if exist(curvmap,'file') % Verifying if the curvature exists
        try
            [curvmap] = read_cfiles(curvmap); % Reading surface
        catch
            error('Unrecognized curvature format');
        end
    end
end

% Curvature checking
try
    if (size(curvmap,1) ~= size(Surf.SurfData.vertices,1))
        error('Different sizes between curvature map and surface');
    end
catch
    error('Different sizes between curvature map and surface');
end

% Checking Options
if nargin < 4
    opts.dephtmap = 0;
    opts.endpthr  = 6;
    opts.filteriter = 100;
elseif nargin == 4
    opts = varargin{4};
    if isstruct(opts)
        if ~isfield(opts,'endpthr')
            opts.endpthr = 6;          % Number of Endpoints threshold.
            %                                      -   opts.endpthr = 0; % Means main Sulcal Line
            %                                      -   opts.endpthr = 1; % All Sulcal Lines
            %                                      -   opts.endpthr > 2; % opts.endpthr - 2 Secondary lines + Main Path
        end
        if ~isfield(opts,'depthmap')
            opts.depthmap = 0;         % Depth Map
        end
        if ~isfield(opts,'filteriter')
            opts.filteriter = 100;
        else
            if length(opts.depthmap)~= size(Surf.SurfData.vertices,1)
                warning('Different sizes between depth map and surface. Depth map will not be used');
                opts.depthmap = 0;
            end
        end
    end
else
    error('Unrecognized options');
end

if nargin > 4
    error('Too Many Input Parameters');
end
if nargout > 2
    error('Too Many Output Parameters');
end
%% ========================= End of Checking Inputs ==================== %%

%% =================== Computing Neighbor points ======================= %%
% load('/media/COSAS/scripts/Gyral_Crowns_and_Sulcal_Lines_Extraction/matlab.mat');

% % % % load('/media/Neuro/ibas/TestSurfOBJ/Surf.mat')
% % % % Surf = subSurf;
% % % % curvmap = Surf.Is;
Surf.Is = curvmap;
opts.sMethod = 4; opts.nSmooth = 12; smoothTesting = Geodesic_Map_Smoothing(Surf, curvmap);
Surf.Is = smoothTesting;
curvmap = smoothTesting;
% Surf = varargin{1};
% distGraph = surface2graph(Surf);
% D = graphallshortestpaths(distGraph);
distTransfmap = Surf.Is+eps;

%% filter surface
Pial = filterSurface(Pial, opts.filteriter);
%% end of filter surface

[C,~] = connected_components([Surf.SurfData.faces;repmat(size(Surf.SurfData.vertices,1),1,3)]);
non_repeated = find(histc(C(:), 1:max(C)) == 1);
C(ismember(C,non_repeated)) = 0;

components = unique(C);
num_comp = size(components,2);
Pial.Is = C';

% find the biggest component
bInd = 1;
for i = 1:num_comp
    if length(find(C == components(i))) > length(find(C == components(bInd))) && components(i) ~= 0
        bInd = i;
    end
end

[subPial,~,indCell] = Extract_Sub_Surface(Pial, components(bInd));
[U] = Skeleton_Smooth(subPial.SurfData.vertices,subPial.SurfData.faces);
[ind_endpoints] = Endpoint_Extraction(U,subPial,curvmap(indCell{1}));
endpoints = indCell{1}(ind_endpoints);

if ~exist('endpoints','var') || isempty(endpoints)
    varargout{1} = [];
    varargout{2} = [];
    return
end

vertLevel = Compute_Level_from_Boundary(Surf);
candidates = find(vertLevel == 1);
loc = dsearchn(Pial.SurfData.vertices(candidates,:),Pial.SurfData.vertices(endpoints,:));
end_points = candidates(loc);
end_points = unique(end_points);


%% test for obtaining sulcal lines
fixedPoints = end_points;
tempExcluded = [];

bEdges = boundary_faces(Surf.SurfData.faces);
bPoints = unique(bEdges);

[~,CF] = connected_components(boundary_faces(Surf.SurfData.faces));
components = zeros(length(Surf.SurfData.vertices),1);
components(bEdges(:,1)) = CF;
components(bEdges(:,2)) = CF;

queue = bPoints;

indNeighEP = find(sum(ismember(Surf.SurfData.faces,end_points),2));
indNeighEP = unique(Surf.SurfData.faces(indNeighEP,:));
fixedPoints = [ fixedPoints; indNeighEP ];

SurfOrig = Surf;

queue(ismember(queue,fixedPoints)) = [];
[Trip] = Vert_Neibp(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),size(Surf.SurfData.faces,1));

% % count = 0;
% % 
% % g = Plot_Surf(Surf);
% % set(g(1),'EdgeColor',[0.8 0.8 0.8]);
% % set(gcf, 'color', [1 1 1])
% % a = gca;
% % view(a,210,-10);
% % hold on;
% % plot3(Surf.SurfData.vertices(ind_endpoints,1),Surf.SurfData.vertices(ind_endpoints,2),Surf.SurfData.vertices(ind_endpoints,3),'.b','MarkerSize',40);
% % 
% % colorbar('off');
% % set(gca, 'XLim', [-29.9010 -12.4657])
% % set(gca, 'YLim', [67.7164 79.1301])
% % 
% % writerObj = VideoWriter('/home/afernandez/vboxwindows/presentaciones/20190412/maxcurvpath/thinning_video.mp4');
% % writerObj.FrameRate = 16;
% % 
% % % open the video writer
% % open(writerObj);
% % 
% % temp = getframe(a.Parent);
% % frame = im2frame(temp.cdata);
% % 
% % % write the frames to the video
% % writeVideo(writerObj, frame);




while ~isempty(queue)
    [~,indQueue] = min(curvmap(queue));
    currentPoint = queue(indQueue);
    queue(indQueue) = [];
    
    currentComp = components(currentPoint);
    
    neighbors = unique((Trip(currentPoint,3:Trip(currentPoint,2)+2)));
    neighbors(neighbors == currentPoint) = [];
    
    neighComps = nonzeros(components(neighbors));
    numBoundNeigh = length(neighComps);
    
    if numBoundNeigh > 2
        if length(unique(neighComps)) > 1
            fixedPoints = [fixedPoints; currentPoint];
        else
            tempExcluded = [tempExcluded; currentPoint];
        end
    else
        faces2remove = find(sum(ismember(Surf.SurfData.faces,currentPoint),2));
        
        if length(faces2remove) == 1 % check if removing the point breaks the connectivity
            currentFace = Surf.SurfData.faces(faces2remove,:);
            currentEdges = [ currentFace(1) currentFace(2); ...
                currentFace(1) currentFace(3);
                currentFace(2) currentFace(3) ];
            
            currentEdges = sort(currentEdges,2);
            boundEdges = boundary_faces(Surf.SurfData.faces);
            boundEdges = sort(boundEdges,2);
            boundEdges = unique(boundEdges,'rows');
            
            if sum(ismember(currentEdges,boundEdges,'rows')) == 3 % if all the edges are in the boundary, the face can't be removed
                fixedPoints = [fixedPoints; currentPoint; neighbors'];
                continue;
            end
            
        end
        
        SurfAux = Surf;
        SurfAux.SurfData.faces(faces2remove,:) = [];
        
        % check if we are removing a fixedPoint
        if sum(ismember(fixedPoints,unique(SurfAux.SurfData.faces(:)))) < length(fixedPoints)
            fixedPoints = [fixedPoints; currentPoint];
            continue;
        end
        
        Surf = SurfAux;
        
% %         g(1).Faces = Surf.SurfData.faces;
% %         g(1).Vertices = Surf.SurfData.vertices;
% %         
% %         temp = getframe(a.Parent);
% % %         imwrite(temp.cdata,['/home/afernandez/vboxwindows/presentaciones/20190412/maxcurvpath/thinning_' num2str(count) '.tiff']);
% %         frame = im2frame(temp.cdata);
% % 
% %         % write the frames to the video
% %         writeVideo(writerObj, frame);
% %         
% %         count = count + 1;
        
        [Trip] = Vert_Neibp(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),size(Surf.SurfData.faces,1));
        
        components(neighbors) = currentComp;
        bPoints(bPoints == currentPoint) = [];
        bPoints = unique([bPoints; neighbors']);
        
        tempExcluded(ismember(tempExcluded,neighbors)) = [];
        
        queue = bPoints;
        queue(ismember(queue,fixedPoints)) = [];
        queue(ismember(queue,tempExcluded)) = [];
    end
    
end

% % % close the writer object
% % close(writerObj);

w12 = (distTransfmap(Surf.SurfData.faces(:,1)) + distTransfmap(Surf.SurfData.faces(:,2)))/2; % Edge weight between vertex 1 and vertex 2
w23 = (distTransfmap(Surf.SurfData.faces(:,2)) + distTransfmap(Surf.SurfData.faces(:,3)))/2; % Edge weight between vertex 2 and vertex 3
w13 = (distTransfmap(Surf.SurfData.faces(:,1)) + distTransfmap(Surf.SurfData.faces(:,3)))/2; % Edge weight between vertex 1 and vertex 3
[a,location] = min([w12 w23 w13],[],2); % Finding the minimum edge weight for each face
[b,pos] = sort([w12 w23 w13]');
pos = pos(2:3,:)';  % Removing edges with the minimum weight

Face2Branches = zeros(size(pos,1),2); % Faces 2 Branches
FacPos = [1 2;2 3; 1 3];
Branches = [0 0];
Branch2Remove = [0 0];
Fac = size(pos,2); % Number of Branches per triangle
Fac2Rem = 3 - Fac; % Number of Branches from the same triangle that will be removed
for k = 1:size(pos,1) % Removing edges with the minimum weight
    
    % Branches
    TempFacpos = FacPos(pos(k,:),:);
    TempBranch = reshape(Surf.SurfData.faces(k,TempFacpos'),[2 Fac])';
    Branches = [Branches;TempBranch]; % Branches
    
    % Branches to be removed
    Pos2Rem = FacPos(find(ismember(FacPos,TempFacpos,'rows') == 0),:); % Branches that will be removed
    TempBranch = reshape(Surf.SurfData.faces(k,Pos2Rem'),[2 Fac2Rem])';
    Branch2Remove = [Branch2Remove;TempBranch];
end
ind = find(ismember(sort(Branches')',sort(Branch2Remove')','rows') == 1); % Removing
Branches(ind,:) = [];

[~,CF] = connected_components(boundary_faces(Surf.SurfData.faces));
numComp = length(unique(CF));
% ReBuilding Graph
listBranches = [[Branches(:,1);Branches(:,2)],[Branches(:,2);Branches(:,1)]];
[C,iak,~] = unique(listBranches,'rows','last');
distances = ((distTransfmap([Branches(:,1);Branches(:,2)]) +  distTransfmap([Branches(:,2);Branches(:,1)]))/2).^-1;
distances = distances(iak);
Graph = sparse(C(:,1),C(:,2),distances,length(Surf.SurfData.vertices),length(Surf.SurfData.vertices));
% Graph with no cycles and minimum total weight
[Tree, ~] = graphminspantree(Graph);
[X,Y] = find(Tree);


if numComp == 1
    Branches = [X Y];
else
    
    Branches_St = [X Y];
    listBranches = [[Branches_St(:,1);Branches_St(:,2)],[Branches_St(:,2);Branches_St(:,1)]];
    [C,iak,~] = unique(listBranches,'rows','last');
    Graph = sparse(C(:,1),C(:,2),ones(length(C),1),length(Surf.SurfData.vertices),length(Surf.SurfData.vertices));

    ind2explore = find(ismember(sort(Branches')',sort(Branches_St')','rows')==0);
    lengthVect = zeros(length(ind2explore),1);
    for st = 1:length(ind2explore)
        tempVar = Branches(ind2explore(st),:);
        stpoint = tempVar(1);
        endpoint = tempVar(2);
        [len,~] = graphshortestpath(Graph,stpoint,endpoint);
        lengthVect(st) = len;
    end
    
    ind2keep = find(lengthVect > 50);
    segments = Branches(ind2explore(ind2keep),:);
    segments = unique(sort(segments')','rows');
    
    occurences = accumarray(segments(:),segments(:)*0+1);
    indRepeated = find(occurences > 1);
    branch2Add = zeros(length(indRepeated),2);
    for i = 1:length(indRepeated)
        [X ~] = find(segments == indRepeated(i));
        curvSegments = sum(distTransfmap(segments(X,:)),2);
        [~,indBranch] = max(curvSegments);
        branch2Add(i,:) = segments(X(indBranch),:);
    end
    
    Branches = [Branches_St; branch2Add];
end

% ReBuilding Graph
listBranches = [[Branches(:,1);Branches(:,2)],[Branches(:,2);Branches(:,1)]];
[C,iak,~] = unique(listBranches,'rows','last');
distances = ((distTransfmap([Branches(:,1);Branches(:,2)]) +  distTransfmap([Branches(:,2);Branches(:,1)]))/2).^-1;
distances = distances(iak);
Graph = sparse(C(:,1),C(:,2),distances,length(Surf.SurfData.vertices),length(Surf.SurfData.vertices));

cont = sum(logical(Graph),2);

% Detecting Graph Endpoints
inde = find(cont == 1);

% Detecting Junction points
indb = find(cont > 2);

% Removing Branches that shares endpoints
indem = ismember(Branches,inde);
Branches(sum(indem,2) == 2,:) = [];

% ReBuilding Graph
listBranches = [[Branches(:,1);Branches(:,2)],[Branches(:,2);Branches(:,1)]];
[C,iak,~] = unique(listBranches,'rows','last');
distances = ((distTransfmap([Branches(:,1);Branches(:,2)]) +  distTransfmap([Branches(:,2);Branches(:,1)]))/2).^-1;
distances = distances(iak);
Graph = sparse(C(:,1),C(:,2),distances,length(Surf.SurfData.vertices),length(Surf.SurfData.vertices));

% removing Branches that not contain boundary endpoints 
Graphnew = sparse(length(Surf.SurfData.vertices),length(Surf.SurfData.vertices));
while sum(sum(Graphnew - Graph)) ~= 0%|~isempty(Branches)
    Graphnew = Graph;
    
    cont = sum(logical(Graph),2);
    
    % Detecting Graph Endpoints
    inde = find(cont == 1);
    
    % Detecting Junction points
    indb = find(cont > 2);
    
    ind2del = inde(~ismember(inde,end_points)); % detected endpoints that are not real endpoints
    branch2del = find(sum(ismember(Branches,ind2del),2)); % Edges containing endpoints that do not belong to the boundary
    
    Branches(branch2del,:) = []; % Deleting edges
    Branches = sort(Branches')';
    Branches = unique(Branches,'rows');
    Graph = sparse(length(Surf.SurfData.vertices),length(Surf.SurfData.vertices));
    if ~isempty( Branches)
        listBranches = [[Branches(:,1);Branches(:,2)],[Branches(:,2);Branches(:,1)]];
        [C,iak,~] = unique(listBranches,'rows','last');
        distances = ((distTransfmap([Branches(:,1);Branches(:,2)]) +  distTransfmap([Branches(:,2);Branches(:,1)]))/2).^-1;
        distances = distances(iak);
        Graph = sparse(C(:,1),C(:,2),distances,length(Surf.SurfData.vertices),length(Surf.SurfData.vertices));
    end
end

out_branches = Branches;


clearvars -except out_branches end_points;


% Outputs
varargout{1} = out_branches;
varargout{2} = end_points;
return;
