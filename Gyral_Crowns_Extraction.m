function varargout = Gyral_Crowns_Extraction(varargin);
%
% Syntax :
%     [WatershParcell, gCrowns] = Gyral_Crowns_Extraction(Surf, curvmap, labSulbasin, opts);
%
% This function computes the cortical surface gyral crowns using a curvature map. It is
% based on a curvature watershed restricted algorithm.
%
% Input Parameters:
%        Surf                           : Surface variable (file, struct or cellarray).
%        curvmap                        : Curvature map (file or Nx1 vector, where N is the 
%                                         number of points in its corresponding surface).
%        labSulbasin (optional input)   : Sulcal Basins (N x 1).
%        opts (optional input)          : Options:
%                                         opts.curvth - Curvature threshold for defining sulcal regions
%                                         opts.npointsth - Number of points threshold
%                                         opts.mwallind - Medial Wall Indexes
%
% Output Parameters:
%         WatershParcell                : Watershed parcellation (N x 1).
%         gCrowns                       : Gyral crowns matlab matrix (Nx4 matrix)
%                                         The first two columns contain crown edge
%                                         enpoints indexes, the third column is a
%                                         boolean variable (1 if the edge belongs to 
%                                         the principal crown line and 0 for additional 
%                                         crown edges) and the fourth column is the 
%                                         crown label according to the watershed parcellation.
%
%
% See also: Computing_Sulcal_Basins Surface_Checking Extracting_Sulcal_Line
% Recursively_Remove_Branches Reordering_Gyral_Crowns Merge_Watershed_Regions
%__________________________________________________
% Authors: Yasser Aleman Gomez
% LIM, HUGGM
% November 13th 2014
% Version $1.0


%% ============================= Checking Inputs ======================= %%
if nargin < 2
    error('Two Inputs are mandatory');
    return
end
Surf = varargin{1};
curvmap = varargin{2};

% Surface Checking
Surf = Surface_Checking(Surf); 

%Checking Curvature Map
if ischar(curvmap) % Verifying if the curvature is a file
    if exist(curvmap,'file') % Verifying if the curvature exists
        try
            [curvmap] = read_cfiles(curvmap); % Reading surface
        catch
            error('Unrecognized curvature format');
            return;
        end
    end
end
try
    if (size(curvmap,1) ~= size(Surf.SurfData.vertices,1))
        error('Different sizes between curvature map and surface');
        return;
    end
catch
    error('Different sizes between curvature map and surface');
    return;
end

if nargin < 3
    opts.curvth  = 0;       % Curvature threshold for defining sulcal regions
    opts.npointsth = 1;   % Number of points threshold
    labSulbasin = Computing_Sulcal_Basins(Surf, curvmap, opts);
elseif nargin >= 3
    labSulbasin = varargin{3};
    if ~sum(labSulbasin-floor(labSulbasin)) ==0
        error('Sulcal basins must be a postive and integer Npoints x 1 vector');
    end
end
try
    if (size(labSulbasin,1) ~= size(Surf.SurfData.vertices,1))
        error('Different sizes between sulcal basins map and surface');
        return;
    end
catch
    error('Different sizes between sulcal basins map and surface');
    return;
end

if nargin < 4
    opts.curvth  = 0;       % Curvature threshold for defining sulcal regions
    opts.npointsth = 10;   % Number of points threshold
elseif nargin == 4
    opts = varargin{4};
    if isstruct(opts)
        if ~isfield(opts,'curvth')
            opts.curvth = 0;        % Curvature threshold
        end
        if ~isfield(opts,'npointsth')
            opts.npointsth = 10;        % Number of points threshold
        end
    else
        error('Unrecognized options');
        return;
    end
end

if nargin > 4
    error('To Many Input Parameters');
    return;
end
if nargout > 2
    error('To Many Output Parameters');
    return;
end
%% ========================= End of Checking Inputs ==================== %%

disp(' ')
disp('Computing Gyral Crowns...');
tic;


%% =================== Computing Neighbor points ======================= %%
%[Surft, LabBranches] = Branch_Labelling(Surf, Slines);
[Trip] = Vert_Neibp(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),size(Surf.SurfData.faces,1));
Temp = sum(Trip);
Trip(:,Temp==0) = [];
temp = Trip(:,3:end);
indz = find(temp == 0);
temp(indz) = 1;

% ================== Creating the Labels Matrix ===================== %
% This matrix contains information about the
% labels inside each point neighborhood

Labels_Mat = labSulbasin(temp);
Labels_Mat(indz) = 0;
% =============== End of Creating the Labels Matrix ================= %

% =================== Neighbor points matrix  ======================= %
A = Indexing_Neights_V2(Trip);
T = A(:,2:end);
ind = find(T);
T(ind) = T(ind) - 2*size(Trip,1);
A(:,2:end) = T;
% ================== End of Neighbor points matrix  ================= %

indn = find(labSulbasin == 0); % Detecting Unlabeled points

if isfield(opts,'mwallind')
    indn(find(ismember(indn, opts.mwallind))) = [];
end



% % % % % % FigID = figure('numbertitle','off','name','New Figure','Color',[0 0 0], 'Position',[20 20 1500 900],'Visible','on','InvertHardcopy','off');
% % % % % % 
% % % % % % Surft = Surf;
% % % % % % Surft.Is = labSulbasin;
% % % % % % subplot(1,2,1,'Color',[0 0 0]);axis off;
% % % % % % Colors = Surf_Color(Surft);
% % % % % % Surft.SurfData.FaceVertexCData = Colors;
% % % % % % Surft.SurfData.FaceColor = 'interp';
% % % % % % Plot_Surf(Surft,'FigID',gcf);
% % % % % % view([270 0]);
% % % % % % 
% % % % % % subplot(1,2,2,'Color',[0 0 0]);
% % % % % % Plot_Surf(Surft,'FigID',gcf);axis off;
% % % % % % view([90 0]);
% % % % % % 
% % % % % % 
% % % % % % 
% % % % % % num_frames_per_second = 24;
% % % % % % aviobj = VideoWriter( '/home/yaleman/Watershed_Grow.avi'); 
% % % % % % open(aviobj);
% % % % % % frame = getframe ( gcf );
% % % % % % writeVideo(aviobj,frame);


cont = 0; %  Iterations Counter Initialization
templabSulbasin = labSulbasin; % Temporal Sulcal Basins


options.W = abs(curvmap);
[~,~,Q] = perform_fast_marching_mesh(Surf.SurfData.vertices, Surf.SurfData.faces, find(labSulbasin), options);
templabSulbasin = labSulbasin(Q);


% % % % % % % % % % ==== Recursive region growing (Watershed curvature restricted algorithm) === %
% % % % % % % % % while ~isempty(indn)
% % % % % % % % %     cont = cont + 1; % Iterations Counter
% % % % % % % % %    % disp(cont);
% % % % % % % % %     Neigh = find(sum(Labels_Mat,2) ~= 0); % Finding Labeled Neighbors
% % % % % % % % %     Neigh(ismember(Neigh,find(templabSulbasin~=0))) = [];
% % % % % % % % %     
% % % % % % % % %     [distmin,indgrow] = max(curvmap(Neigh)); % Growing point point Index
% % % % % % % % %     vert_to_grow = Neigh(indgrow); % Vertex to grow
% % % % % % % % %     
% % % % % % % % %     % Labeled with the label of the majority of its labeled neighbors
% % % % % % % % %     
% % % % % % % % %     testVar = unique(nonzeros(Labels_Mat(vert_to_grow,:)));
% % % % % % % % %     if length(testVar)~=1
% % % % % % % % %         a = 1;
% % % % % % % % %     end
% % % % % % % % %     templabSulbasin(vert_to_grow) = mode(nonzeros(Labels_Mat(vert_to_grow,:)));
% % % % % % % % %     
% % % % % % % % %     % Updating the labels matrix
% % % % % % % % %     L = A(vert_to_grow,1);
% % % % % % % % %     Labels_Mat(A(vert_to_grow,2:L+1)) =  templabSulbasin(vert_to_grow);
% % % % % % % % %     
% % % % % % % % %     % Detecting if someone is missing for labelling
% % % % % % % % %     indn = find(templabSulbasin == 0);
% % % % % % % % %     
% % % % % % % % %     if isfield(opts,'mwallind')
% % % % % % % % %         indn(find(ismember(indn, opts.mwallind))) = [];
% % % % % % % % %     end
% % % % % % % % %     
% % % % % % % % %     
% % % % % % % % %     
% % % % % % % % % % % % % % % %     Surft = Surf;
% % % % % % % % % % % % % % % %     Surft.Is = templabSulbasin;
% % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % %     subplot(1,2,1);
% % % % % % % % % % % % % % % %     cla;
% % % % % % % % % % % % % % % %     Colors = Surf_Color(Surft);
% % % % % % % % % % % % % % % %     Surft.SurfData.FaceVertexCData = Colors;
% % % % % % % % % % % % % % % %     Surft.SurfData.FaceColor = 'interp';
% % % % % % % % % % % % % % % %     Plot_Surf(Surft,'FigID',gcf);
% % % % % % % % % % % % % % % %     view([270 0]);
% % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % %     subplot(1,2,2);
% % % % % % % % % % % % % % % %     cla;
% % % % % % % % % % % % % % % %     Plot_Surf(Surft,'FigID',gcf);
% % % % % % % % % % % % % % % %     view([90 0]);
% % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % %     
% % % % % % % % % % % % % % % %     frame = getframe ( gcf );
% % % % % % % % % % % % % % % %     writeVideo(aviobj,frame);
% % % % % % % % %     
% % % % % % % % %     
% % % % % % % % %     
% % % % % % % % %     
% % % % % % % % % end
% % % % % % % % % 

% % % % % close(aviobj);

% ==== End of recursive region growing (Watershed curvature restricted algorithm) ==== %

% ============= Assigning 4000 to the medial wall values ================ %
tempSize = accumarray(templabSulbasin, templabSulbasin*0+1);
sizes2del = find( (tempSize< 10)&(tempSize~=0)); % Removing regions with less than 3 points (face)
if ~isempty(sizes2del)
    templabSulbasin(find(ismember(templabSulbasin, sizes2del))) = 0;
    Surf.Is = templabSulbasin;
    Surf2 = Surf_Corr(Surf);
    WatershParcell = Surf2.Is;
else
    WatershParcell = templabSulbasin;
end
if isfield(opts,'mwallind')
    WatershParcell(opts.mwallind) = 4000;
end
% ========== End of ssigning 4000 to the medial wall values ============= %

% ====== Detecting possible wholes in the Watershed parcellation ======== %
sts = nonzeros(unique(WatershParcell));
while ~isempty(sts)
    ind = find(WatershParcell == sts(1));
    Neigh = unique(nonzeros(Trip(ind,3:end)));
    tempParc = WatershParcell(Neigh);
    ind2rem = find(tempParc == sts(1));
    tempParc(ind2rem) = []; % Removing the labels from the same region
    cc = accumarray(tempParc,tempParc*0+1); % Detecting how many labels are present in the neighborhood
    neighLabels = find(cc~=0);
    if length(neighLabels) == 1
        WatershParcell(ind) = neighLabels;
    end
    sts(1) = [];
end
Surf.Is = WatershParcell;
Surf2 = Surf_Corr(Surf);
WatershParcell = Surf2.Is;
% == End of Detecting possible wholes in the Watershed parcellation ===== %

% ==================== Reordering Watershed labels ====================== %
[~,~, newWatershParcell] = unique(WatershParcell);
indtemp = find(WatershParcell == 4000);
if ~isempty(indtemp)
    indmedialWall = find(newWatershParcell == newWatershParcell(indtemp(1)));
    newWatershParcell(indmedialWall) = 4000;
end
WatershParcell = newWatershParcell;
% =============== End of Reordering Watershed labels ==================== %

% WatershParcell = Merge_Watershed_Regions(Surf, WatershParcell, curvmap);

% ==================== Extracting Crowns points ========================= %
temp1 = WatershParcell(temp);
temp1(indz) =  max(temp1(:))+1;
NewMat = temp1-repmat(min(temp1')',[1 size(temp1,2)]);
NewMat(indz) = 0;
a = logical(sum(logical(NewMat)')');
indc = find(a);
gCrowns = zeros(length(WatershParcell),1);
gCrowns(indc) = 1;
gCrowns = WatershParcell.*gCrowns;
% ==================== End of Extracting Crown points =================== %

% ==================== Creating Crowns Edges ============================ %
indcrown = find(gCrowns);
indcrownfaces = find(sum(ismember(Surf.SurfData.faces,indcrown),2) == 3);
crownfaces = Surf.SurfData.faces(indcrownfaces,:);
w12 = (gCrowns(crownfaces(:,1)) - gCrowns(crownfaces(:,2))); % Edge weight between vertex 1 and vertex 2
w23 = (gCrowns(crownfaces(:,2)) - gCrowns(crownfaces(:,3))); % Edge weight between vertex 2 and vertex 3
w13 = (gCrowns(crownfaces(:,1)) - gCrowns(crownfaces(:,3))); % Edge weight between vertex 1 and vertex 3
crownfaces = crownfaces(find(sum((logical(w12) +logical(w23) + logical(w13)),2)~=0),:);
crownEdges = [crownfaces(:,1) crownfaces(:,2); crownfaces(:,2) crownfaces(:,3);crownfaces(:,1) crownfaces(:,3)] ; % Edges from the intersection faces
% ================== End of Creating Crowns Edges ======================= %

% ==================== Removing some Crowns Edges ======================= %
edgeLabcrowns =  gCrowns(crownEdges);
ind2keep= find(edgeLabcrowns(:,1) - edgeLabcrowns(:,2) == 0);
indcrossEdges= find(edgeLabcrowns(:,1) - edgeLabcrowns(:,2) ~= 0);
crossEdges = crownEdges(indcrossEdges,:);
edgeLabcrowns = edgeLabcrowns(ind2keep,:);
crownEdges = [crownEdges(ind2keep,:) edgeLabcrowns(:,1)];
[~,ord] = sort(edgeLabcrowns(:,1));
gCrowns = crownEdges(ord,:);
gCrowns = [gCrowns(:,1:2) ones(size(gCrowns,1),1) WatershParcell(gCrowns(:,1))];
% ================= End of Removing some Crowns Edges  ================== %


% ==================== Removing Isolated Edges ========================== %
tempVar = unique(sort(gCrowns(:,1:2)')','rows');
cc = accumarray(tempVar(:),tempVar(:)*0+1);
rempos = find(cc == 1);
while ~isempty(rempos)
    if rempos
        [Xrem, Yrem] = find(sum(ismember(tempVar,rempos),2));
        tempVar(Xrem,:) = [];
    end
    cc = accumarray(tempVar(:),tempVar(:)*0+1);
    rempos = find(cc == 1);
end
ind2del = find(ismember(sort(gCrowns(:,1:2)')',tempVar,'rows') == 0);
gCrowns(ind2del,:) = [];
% ==================== End of Removing Isolated Edges =================== %


% =============== Reordering Gyral Crowns Edges ========================= %
gCrowns = Reordering_Gyral_Crowns(Surf, gCrowns);
% =============== End of Reordering Gyral Crowns Edges ================== %

% =============== Adding New Gyral Crowns Edges ========================= %
try
    [gCrowns] = Adding_Internal_Gyral_Crowns(Surf, curvmap, WatershParcell, gCrowns);
end

% =============== End of Adding New Gyral Crowns Edges ================== %
disp(['End of Computing Gyral Crowns ....' ]);
toc;

if isfield(opts,'mwallind')
    WatershParcell(opts.mwallind) = 4000;
end

% Outputs
varargout{1} = WatershParcell;
varargout{2} = gCrowns;
return;