function [sulcalWidthMap,normDepthMap,segmentationMap,nanMap,pairPoints] = Estimate_Sulcal_Width(Surf,indMedialWall,depthMap,sulcalLines,extensionLines,sulcalLabels)

Surf.Is = sulcalLabels;

labels = unique(sulcalLabels);

sulcalWidthMap = NaN(length(Surf.Is),1);
normDepthMap = Surf.Is*0;
segmentationMap = Surf.Is*0;

for j = 1:length(labels)
    
    label = labels(j);
    
    if label == 0
        continue;
    end
    
    disp(['Processing: ' num2str(j) ' of ' num2str(length(labels))]);
    
    [subSurf,~,indCell] = Extract_Sub_Surface(Surf,label);
    
    indSubSurf = indCell{1};
    subSurf.Is = subSurf.Is*0;
    
    indLine = find(ismember(indSubSurf,unique(sulcalLines(:))));
    subSurf.Is(indLine) = 1;
    indExt = find(ismember(indSubSurf,unique(extensionLines(:))));
    subSurf.Is(indExt) = 1;
    
    indLine = [indLine; indExt];
    
    if isempty(indLine) || length(indLine) < 3 || length(subSurf.Is) < 20
        continue;
    end
    
    bPoints = boundary_faces(subSurf.SurfData.faces);
    bPoints = unique(bPoints(:));
    
    depthBoundLine = Compute_Geodesic_Distance_Transform(subSurf,bPoints);
    depthLineBound = Compute_Geodesic_Distance_Transform(subSurf,indLine);
    
    normDepth = 1 - depthLineBound./(depthLineBound + depthBoundLine + eps);
    
    normDepthMap(indSubSurf) = normDepth;
    
    [sulcRegions] = Recur_Corr(subSurf,0,zeros(size(subSurf.Is)),1);
    subSurf.Is = sulcRegions;
    [subSurf] = Surf_Corr(subSurf);
    sulcRegions = subSurf.Is;
    
    segmentationMap(indSubSurf) = max(segmentationMap) + sulcRegions;
    
end
        
edge12 = normDepthMap([Surf.SurfData.faces(:,1) Surf.SurfData.faces(:,2)]); % Edge weight between vertex 1 and vertex 2
edge23 = normDepthMap([Surf.SurfData.faces(:,2) Surf.SurfData.faces(:,3)]); % Edge weight between vertex 2 and vertex 3
edge13 = normDepthMap([Surf.SurfData.faces(:,1) Surf.SurfData.faces(:,3)]); % Edge weight between vertex 1 and vertex 3
tempedge = [edge12;edge23;edge13];

[deldDist] = mean(abs(tempedge(:,1) - tempedge(:,2))); % Delta distance for point comparisons
Surf = Compute_Surface_Normals(Surf);


% Split the process: for the points where depth > 0, the width is estimated
% using the euclidean distance straightforward inside their own sulcus. In 
% regions with depth = 0, we will check the closest points in every sulci
% having in mind that the line doesn't cross the surface.

indDepthZero = find(depthMap == 0);
indDepthNonZero = find(depthMap ~= 0);

disp('Estimating sulcal width');
reverseStr = '';

pairs=[];

for i = 1:length(Surf.Is)
    
    % Display the progress
    percentDone = 100 * i / length(Surf.Is);
    msg = sprintf('Percent done: %3.1f', percentDone); 
    fprintf([reverseStr, msg]);
    reverseStr = repmat(sprintf('\b'), 1, length(msg));
    
%     idx = indDepthNonZero(i);
    idx = i;
    
    if segmentationMap(idx) == 0
        pairs=[pairs; [idx NaN NaN]];
        continue;
    end
    
    if ismember(idx,indLine)
        pairs=[pairs; [idx idx 0]];
        continue;
    end
    
    currDist = normDepthMap(idx); %find distance of my point
    points = find(((normDepthMap > currDist-deldDist) & (normDepthMap < currDist+deldDist)) ...
        & segmentationMap ~= segmentationMap(idx) ...
        & segmentationMap ~= 0 ...
        & sulcalLabels == sulcalLabels(idx));
    
    if isempty(points)
        pairs=[pairs; [idx NaN NaN]];
        continue;
    end
    
    %Calculate the euclidean distance for all ppoints
    distance_euclidean = sqrt((Surf.SurfData.vertices(idx,1) - Surf.SurfData.vertices(points(:),1)).^2 + ...
        (Surf.SurfData.vertices(idx,2) - Surf.SurfData.vertices(points(:),2)).^2 + ...
        (Surf.SurfData.vertices(idx,3) - Surf.SurfData.vertices(points(:),3)).^2);
    
    [widths,loc_point]=min(distance_euclidean);
    if widths < 35 % if distance is more than 3.5 cm discard it
        pairs=[pairs; [idx points(loc_point) widths]];
    else
        pairs=[pairs; [idx NaN NaN]];
    end
end
disp(' ');

% % for i = 1:length(indDepthZero)
% %     
% %         
% %     % Display the progress
% %     percentDone = 100 * (i + length(indDepthNonZero)) / length(Surf.Is);
% %     msg = sprintf('Percent done: %3.1f', percentDone); 
% %     fprintf([reverseStr, msg]);
% %     reverseStr = repmat(sprintf('\b'), 1, length(msg));
% %     
% %     idx = indDepthZero(i);
% %     
% %     if segmentationMap(idx) == 0
% %         pairs=[pairs; [idx NaN NaN]];
% %         continue;
% %     end
% %     
% %     if ismember(idx,indLine)
% %         pairs=[pairs; [idx idx 0]];
% %         continue;
% %     end
% %     
% %     locatedPoint = 0;
% %     
% %     currDist = normDepthMap(idx); %find distance of my point
% %     points = find(((normDepthMap > currDist-deldDist) & (normDepthMap < currDist+deldDist)) ...
% %         & segmentationMap ~= segmentationMap(idx) ...
% %         & segmentationMap ~= 0);
% %     
% %     if ~isempty(points)
% %         
% %         %Calculate the euclidean distance for all points
% %         distance_euclidean = sqrt((Surf.SurfData.vertices(idx,1) - Surf.SurfData.vertices(points(:),1)).^2 + ...
% %             (Surf.SurfData.vertices(idx,2) - Surf.SurfData.vertices(points(:),2)).^2 + ...
% %             (Surf.SurfData.vertices(idx,3) - Surf.SurfData.vertices(points(:),3)).^2);
% %         
% %         widths = 0;
% %         while ~isempty(points) && widths < 20
% %             [widths,loc_point]=min(distance_euclidean);
% %             
% %             inside = Intercept_Surf_with_LinesSegment(Surf,Surf.SurfData.vertices(idx,:),Surf.SurfData.vertices(points(loc_point),:));
% %             
% %             if inside == 0
% %                 locatedPoint = 1;
% %                 break
% %             else
% %                 clear inside;
% %                 points(loc_point) = [];
% %                 distance_euclidean(loc_point) = [];
% %             end
% %         end
% %         
% %     end
% %     
% %     if locatedPoint
% %         pairs=[pairs; [idx points(loc_point) widths]];
% %         continue
% %     else
% %         pairs=[pairs; [idx NaN NaN]];
% %     end
% %     
% % %     [widths,loc_point]=min(distance_euclidean);
% % %     pairs=[pairs; [i points(loc_point) widths]];
% % end

sulcalWidthMap = pairs(:,3);
pairPoints = pairs(:,1:2);
nanMap = isnan(sulcalWidthMap);

disp('Estimating values for unlabeled regions');
[Trip] = Vert_Neibp(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),size(Surf.SurfData.faces,1));
Temp = sum(Trip);
Trip(:,Temp==0) = [];
temp = Trip(:,3:end);
indNan = find(nanMap);
indFill = indNan(~ismember(indNan,indMedialWall));
sulcalWidthMap(indNan) = 0;
neigh = temp(indFill,:);
C = num2cell(neigh,2);
p = cell(length(C),1);

for i = 1:length(C)
    p{i} = unique(nonzeros(C{i}));
end

oldMeanVals = randn(length(indFill),1);
meanVals = zeros(length(indFill),1);

while (sum( (meanVals - oldMeanVals) < 0.001 ) ~= length(meanVals))
    oldMeanVals = meanVals;
    
    for i = 1:length(meanVals)
        meanVals(i) = mean(sulcalWidthMap(p{i}));
    end
%     tempVals = cellfun(@(p) sulcalWidthMap(p),p,'UniformOutput',false);
%     meanVals = cellfun(@mean, tempVals);
    sulcalWidthMap(indFill) = meanVals;
end

% % % p = cellfun(@nonzeros,C,'UniformOutput',false);
% % % p = cellfun(@unique,p,'UniformOutput',false);
% % % 
% % % oldMeanVals = randn(length(indFill),1);
% % % meanVals = zeros(length(indFill),1);
% % % 
% % % while (sum( (meanVals - oldMeanVals) < eps ) ~= length(meanVals))
% % %     oldMeanVals = meanVals;
% % %       
% % %     tempVals = cellfun(@(p) sulcalWidthMap(p),p,'UniformOutput',false);
% % %     meanVals = cellfun(@mean, tempVals);
% % %     sulcalWidthMap(indFill) = meanVals;
% % % end


end