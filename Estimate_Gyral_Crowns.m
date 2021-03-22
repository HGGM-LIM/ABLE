function [gCrowns] = Estimate_Gyral_Crowns(varargin)
%
% Syntax :
%       [gCrowns] = Estimate_Gyral_Crowns(GyralSurf,curvMap,fixedPoints);
%
% This function computes the gyral crowns by removing the faces with
% maximum curvature while maintaining connectivity and the points in
% fixedPoints
%
% Input Parameters:
%        GyralSurf              : Matlab surface variable with the fundi 
%                                   part removed. 
%        curvMap                : Nvert x 1 curvature map
%        fixedPoints            : Indices of the vertices that we want to
%                                   conserve
%        mWallInd               : Indices of the vertices in the medial
%                                   wall
%
% Output Parameters:
%        gCrowns                : Edges of the gyral crowns
%
%__________________________________________________
%

if nargin < 2
    error('Wrong number of arguments');
end

GyralSurf = varargin{1};
curvMap = varargin{2};

if nargin > 2
    fixedPoints = varargin{3};
else
    fixedPoints = [];
end

if nargin > 3
    mWallInd = varargin{4};   
end

indZeros = find(curvMap == 0);
faces2remove = find(sum(ismember(GyralSurf.SurfData.faces,indZeros),2));
GyralSurf.SurfData.faces(faces2remove,:) = [];

if exist('mWallInd','var')
    faces2remove = find(sum(ismember(GyralSurf.SurfData.faces,mWallInd),2));
    GyralSurf.SurfData.faces(faces2remove,:) = [];
end 

Surf = GyralSurf;

bEdges = boundary_faces(Surf.SurfData.faces);
bPoints = unique(bEdges);

[~,CF] = connected_components(boundary_faces(Surf.SurfData.faces));
components = zeros(length(Surf.SurfData.vertices),1);
components(bEdges(:,1)) = CF;
components(bEdges(:,2)) = CF;

queue = bPoints;

tempExcluded = [];

queue(ismember(queue,fixedPoints)) = [];
[Trip] = Vert_Neibp(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),size(Surf.SurfData.faces,1));

while ~isempty(queue)
    [~,indQueue] = max(curvMap(queue));
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
        Surf.SurfData.faces(faces2remove,:) = [];
        
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



gCrowns = boundary_faces(Surf.SurfData.faces);

end