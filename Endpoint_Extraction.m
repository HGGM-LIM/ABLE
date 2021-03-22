function [ind_endpoints] = Endpoint_Extraction(W,Surf,curvmap)

% % % for very small sulci, get the PCA alone
% % if length(Surf.SurfData.vertices) < 150 
% %     [~,score] = pca(W);
% %     
% %     [~,indmax] = max(score(:,1));
% %     [~,indmin] = min(score(:,1));
% %     ind_endpoints = [indmax indmin];
% %     return
% % end

vis = false;
NEIGHB_DISTANCE = 4.5;

Surft.SurfData.vertices = W;
Surft.SurfData.faces = Surf.SurfData.faces;

A = surface2graph(Surf);
dist_mat = graphallshortestpaths(A);


possible_endpoints = zeros(size(W,1),1);
principal_dirs = zeros(size(W,1),3);

num_neighbors = zeros(size(W,1),1);

point_count = zeros(size(W,1),1);
for i = 1:length(W)
    
    dist_row = dist_mat(i,:);

    
    points = W(dist_row < NEIGHB_DISTANCE,:);
    I = find(dist_row < NEIGHB_DISTANCE);
    
    num_points = size(points,1);   
    
    if num_points < 2
        continue;
    end
        
    num_neighbors(i) = size(points,1);
    
    [coef,score,~,~,explained] = pca(points);
    
    point_count(I(1:num_points)) = point_count(I(1:num_points)) + explained(1);
    
    principal_dirs(i,:) = coef(:,1)';
    
    [~,indmax] = max(score(:,1));
    [~,indmin] = min(score(:,1));
    
%     if index_I
        indices = [ I(indmin) I(indmax)];
%     else
%         indices = [ find(cumsum(neighbor_mat(i,:)) == indmin,1) find(cumsum(neighbor_mat(i,:)) == indmax,1)];
%     end
    
%     possible_endpoints(i,:) = indices;

    possible_endpoints(I(indmin)) = possible_endpoints(I(indmin)) + explained(1);
    possible_endpoints(I(indmax)) = possible_endpoints(I(indmax)) + explained(1);

%     if ismember(I(indmin),candidates)
%         possible_endpoints(I(indmin)) = possible_endpoints(I(indmin)) + 1;
%     elseif ismember(I(indmax),candidates)
%         possible_endpoints(I(indmax)) = possible_endpoints(I(indmax)) + 1;
%     end
end

% endpoint_probability = histc(possible_endpoints(:), 1:size(W,1))./point_count;
endpoint_probability = possible_endpoints./point_count;

ind_endpoints = endpoint_probability >= 0.9;

recandidates = find((endpoint_probability > 0.1) & (endpoint_probability < 0.9));
queue = recandidates;

while ~isempty(queue)
    [~,current] = max(endpoint_probability(queue));
    current = queue(current);
    queue(queue == current) = [];
    
    dist_indices = dist_mat(current,recandidates) < NEIGHB_DISTANCE*0.15;
    near_recandidates = recandidates(dist_indices);
    
    sum_prob = sum(endpoint_probability(near_recandidates));
    
    if sum_prob > 0.9
        ind_endpoints(current) = 1;
        queue(ismember(queue,near_recandidates)) = [];       
    end
end

% endpoint_probability(find(point_count)) = possible_endpoints(find(point_count))./point_count(find(point_count));
% ind_endpoints = find(endpoint_probability > 0.9);


if vis
    try
%         plot(endpoint_probability)
        g = Plot_Surf(Surf,'transpVal',0.0)
        hold on
        plot3(W(:,1),W(:,2),W(:,3),'.b','MarkerSize',25);
        plot3(W(ind_endpoints,1),W(ind_endpoints,2),W(ind_endpoints,3),'.r','MarkerSize',40);
        set(g,'EdgeColor',[0.8 0.8 0.8]);
    catch
    end
    
% %     a = gca;
% %     temp = getframe(a.Parent);
% %     imwrite(temp.cdata,['sulc' num2str(sulcNumber) '_1.tiff']);
% %     view(a,0,90)
% %     temp = getframe(a.Parent);
% %     imwrite(temp.cdata,['sulc' num2str(sulcNumber) '_2.tiff']);
% %     view(a,90,0)
% %     temp = getframe(a.Parent);
% %     imwrite(temp.cdata,['sulc' num2str(sulcNumber) '_3.tiff']);
    
% %     close;
end

% if the algorithm only detects 1 endpoint, try to perform a PCA for the
% whole sulcus
if length(find(ind_endpoints)) < 2 
    [~,score] = pca(W);
    
    [~,indmax] = max(score(:,1));
    [~,indmin] = min(score(:,1));
    ind_endpoints = [indmax indmin];
    return
end

end