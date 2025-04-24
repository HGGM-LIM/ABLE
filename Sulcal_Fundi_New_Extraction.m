function [fundiEdges,gyralCrowns,endPoints,sulcalBasins,depthMap] = Sulcal_Fundi_New_Extraction(pialSurf,pialCurv,whiteSurf,whiteCurv,annot)

DEPTH_THR = 2;
FILTER_ITER = 50;

if isstring(whiteSurf) || ischar(whiteSurf)
    Surf = Read_Surface(whiteSurf);
else
    Surf = whiteSurf;
end

if isstring(pialSurf) || ischar(pialSurf)
    Pial = Read_Surface(pialSurf);
else
    Pial = pialSurf;
end

if isstring(annot) || ischar(annot)
    annot = read_cfiles(annot); 
end

if isstring(whiteCurv) || ischar(whiteCurv)
    curv = read_cfiles(whiteCurv);
else
    curv = whiteCurv;
end

if isstring(pialCurv) || ischar(pialCurv)
    curv_pial = read_cfiles(pialCurv);
else
    curv_pial = pialCurv;
end


%curv = discrete_mean_curvature(Pial.SurfData.vertices,Pial.SurfData.faces)*-1;
%opts.nSmooth = 3;
%curv = Geodesic_Map_Smoothing(Pial,curv,opts);

disp('Depth Map Estimation');
[V,F] = signed_distance_isosurface(Pial.SurfData.vertices,Pial.SurfData.faces,'Level',10,'GridSize',400);
[V,F] = signed_distance_isosurface(V,F,'Level',-10,'GridSize',400);
[~,I,~,~] = signed_distance(V,Pial.SurfData.vertices,Pial.SurfData.faces);
indd = Pial.SurfData.faces(I,:);
indd = unique(indd(:));
depthMap = Compute_Geodesic_Distance_Transform(Pial,indd);
depthMap = Geodesic_Map_Smoothing(Surf,depthMap);

% depthMap = travelDepth(Pial.SurfData.vertices,Pial.SurfData.faces);

% load('100408-depthmap.mat'); 


disp('Brain Surface Segmentation');
[labels] = Surface_Segmentation(Pial,curv,depthMap);
opts.fill = 0;
%labels = Dilate_Surface_Label(Pial,labels,opts);

disp('Gyral Crowns Estimation');
Surf.Is = labels;
Surf.Is(Surf.Is == 0) = curv(Surf.Is == 0) > 0; % remove regions with positive curvature and low depth
Surf.Is(annot == 0) = 1;

faces2remove = find(sum(ismember(Surf.SurfData.faces,find(Surf.Is == 0)),2) > 1);
points2remove = Surf.SurfData.faces(faces2remove,:);
points2remove = unique(points2remove);
Surf.Is = Surf.Is.*0 +1;
Surf.Is(points2remove) = 0;

[subSurf,~,indCell] = Extract_Sub_Surface(Surf,0);
% filter gyral crowns
subSurfCopy = filterSurface(subSurf, FILTER_ITER);

[U] = Skeleton_Smooth(subSurfCopy.SurfData.vertices,subSurfCopy.SurfData.faces);
[ind_endpoints] = Gyral_Endpoint_Extraction(U,subSurfCopy);
[gCrowns] = Estimate_Gyral_Crowns(subSurf,curv(indCell{1}),ind_endpoints);
reLines = indCell{1}(gCrowns);
Surf.Is = Surf.Is.*0;
Surf.Is(annot == 0) = 1;
Surf.Is(reLines(:)) = 1;

Surfaux = Surf;
indZeros = find(Surf.Is == 1);
faces2remove = find(sum(ismember(Surf.SurfData.faces,indZeros),2));
Surfaux.SurfData.faces(faces2remove,:) = [];
% [subSurf,~,indCell] = Extract_Sub_Surface(Surf,0);
[C,~] = connected_components([Surfaux.SurfData.faces;repmat(size(Surfaux.SurfData.vertices,1),1,3)]);
sizeCompo = accumarray(C',C'*0+1);
sizeCompo(sizeCompo < 50) = 0;
lab2remove = find(sizeCompo == 0);
newParcell = C';
newParcell(ismember(newParcell,lab2remove)) = 0;


% load('matlab.mat');

% curv = Surf.Is;
% Surf.Is = newParcell == 0;
% [subSurf,~,indCell] = Extract_Sub_Surface(Surf,0);
% [C,~] = connected_components([subSurf.SurfData.faces;repmat(size(subSurf.SurfData.vertices,1),1,3)]);
% sizeCompo = accumarray(C',C'*0+1);
% newParcell = Surf.Is.*0;
% newParcell(indCell{1}) = C';
% structs = unique(C);

Surf.Is = newParcell;

Surfaux = Surf;
indZeros = find(newParcell ~= 0);
faces2remove = find(sum(ismember(Surf.SurfData.faces,indZeros),2));
Surfaux.SurfData.faces(faces2remove,:) = [];
gyralCrowns = boundary_faces(Surfaux.SurfData.faces);

structs = unique(newParcell);
structs(structs == 0) = [];
Nsulc = length(structs);

opts.nSmooth = 0;        % Number of smoothing iterations
opts.sMethod = 4 ;        % Smoothing method
opts.fill = 0;

points = [0];
sLines = [0 0];
[opts2{1:Nsulc}] = deal(opts);

clust = parcluster('local');
clust.NumWorkers = 4;
parpool(clust.NumWorkers);


% parfor i = 1:Nsulc % For each sulcal space
for i = 1:Nsulc % For each sulcal space
    disp(['Processing Sulcal Space ' num2str(i) ' from ' num2str(Nsulc)]);
    
    Surftt = Pial;
    Surftt.Is = newParcell;
    [Pialt,~,indCell] = Extract_Sub_Surface(Surftt, structs(i));
    Surftt = Surf;
    Surftt.Is = newParcell;
    [Surft,~,indCell] = Extract_Sub_Surface(Surftt, structs(i));
    ind = indCell{1};
    
    % skip small sulci
    if length(Surft.SurfData.vertices) < 50
        continue;
    end

    curvmapt = curv_pial(ind);         % Temporal curvature map
    depthmapt = depthMap(ind);

% try removing depth faces
    Pialaux = Pialt;
    Pialaux.Is = depthmapt;
    ind2remove = find(Pialaux.Is <= DEPTH_THR);
    faces2remove = find(sum(ismember(Pialaux.SurfData.faces,ind2remove),2));
    Pialaux.SurfData.faces(faces2remove,:) = [];
    
    if isempty(Pialaux.SurfData.faces)
        continue;
    end
    
    [~,CF] = connected_components(Pialaux.SurfData.faces);
    numComp = length(unique(CF));
    
    
    dilLab = 0;
    depth_thr_iterative = DEPTH_THR;
    
    while numComp > 1
        Pialaux = Pialt;
        depth_thr_iterative = depth_thr_iterative - 0.1;
        Pialaux.Is = depthmapt > depth_thr_iterative;
        ind2remove = find(Pialaux.Is ~= 1);
        faces2remove = find(sum(ismember(Pialaux.SurfData.faces,ind2remove),2));
        Pialaux.SurfData.faces(faces2remove,:) = [];
        [C,~] = connected_components([Pialaux.SurfData.faces;repmat(size(Pialaux.SurfData.vertices,1),1,3)]);
        sizeCompo = accumarray(C',C'*0+1);
        ind2rem = find(sizeCompo == 1);
        ind2rem = ismember(C,ind2rem);
        C(ind2rem) = 0;
        opts2{i}.reminter = 0;
        dilLab = Dilate_Surface_Label(Pialt,C',opts2{i});
        Pialaux = Pialt;
        Pialaux.Is = dilLab > 0;
        ind2remove = find(Pialaux.Is == 0);
        faces2remove = find(sum(ismember(Pialaux.SurfData.faces,ind2remove),2));
        Pialaux.SurfData.faces(faces2remove,:) = [];
        [~,CF] = connected_components(Pialaux.SurfData.faces);
        numComp = length(unique(CF));
    end
    
    if ~isempty(Pialaux.SurfData.faces) && (numComp == 1) % && (length(Pialaux.SurfData.faces)/length(Pialt.SurfData.faces) > 0.75)
        if length(dilLab) == length(depthmapt)
            Pialt.Is = dilLab > 0;
            Surft.Is = dilLab > 0;
        else
            tempMap = depthmapt > DEPTH_THR; 
         
            Pialt.Is = tempMap;
            Surft.Is = tempMap;
        end
    else
        Pialt.Is = depthmapt > DEPTH_THR;
        Surft.Is = depthmapt > DEPTH_THR;
    end
    
    if isempty(find(Pialt.Is,1))
        disp('skip');
        continue;
    end

    
    [Pialt2,~,indCell2] = Extract_Sub_Surface(Pialt,1);
    [Surft2,~,indCell2] = Extract_Sub_Surface(Surft,1);
    
    curvmapt = curvmapt(indCell2{1});
    
    opts2{i}.sulcNumber = i;
    opts2{i}.filteriter = FILTER_ITER;
    
    % Extracting Sulcal Line
    curvmapt = curvmapt + abs(min(curvmapt)) + eps;
    try
        if size(Surft2.SurfData.vertices,1) < 100 % remove very small sulci
            continue;
        end
        
        disp('Extract sulcal topology');
        [skelEdges,endpoints] = Extracting_Sulcal_Lines(Surft2,Pialt2, curvmapt, opts2{i});
        if ~isempty(endpoints)
            points = [points;ind(indCell2{1}(endpoints))] ;
        end
        if ~isempty(skelEdges)
            newLines = ind(indCell2{1}(skelEdges(:,1:2)));
            if (size(newLines, 2) == 2)
                sLines = [sLines; newLines]; % Saving Branches. Branches are saved
            else
                sLines = [sLines; newLines'];
            end
        end
    catch e
        warning(['Something is wrong with sulc ' num2str(i)]);
        warning(e.message);
    end
    
    % clearing variables
    Surftt = [];
    Pialt = [];
    Surft = [];
    ind = [];
    C = [];
    sizeCompo = [];
    indCell = [];
    curvmapt = [];
end
delete(gcp('nocreate'));
points(1,:) = [];
sLines(1,:) = [];

endPoints = points;
fundiEdges = sLines;
sulcalBasins = newParcell;

Surf = Surf_Corr(Surf);

newParcell = Surf.Is;



end
