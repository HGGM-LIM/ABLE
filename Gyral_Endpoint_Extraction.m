function [ind_endpoints] = Gyral_Endpoint_Extraction(U,Surf)

vis = false;
NEIGHB_DISTANCE = 20;

Surft.SurfData.vertices = U;
Surft.SurfData.faces = Surf.SurfData.faces;

Graph = surface2graph(Surf);
% dist_mat = graphallshortestpaths(A);


[Trip] = Vert_Neibp(double(Surft.SurfData.faces),size(Surft.SurfData.vertices,1),size(Surft.SurfData.faces,1));
Temp = sum(Trip);
Trip(:,Temp==0) = [];
temp = Trip(:,3:end);

cont = 0;
for i = 1:length(U)
    
    Neigh = unique(nonzeros(temp(i,:)));
    
    
    A = U(Neigh,:) - repmat(U(i,:),[length(Neigh) 1]);
    
    Anorm = normm(A);
    A = A./Anorm;
    
    
    B = sum(A(2:end,:).*repmat(A(1,:),[length(Neigh)-1 1]),2);
    
    C = find(acos(B) < pi/6); % less than 30
    
    %     ind = find(C > 0);
    if length(C) == length(B)
        cont = cont + 1;
        candidates(cont) = i;
    end
    
end

possible_endpoints = zeros(size(U,1),1);
principal_dirs = zeros(size(U,1),3);

num_neighbors = zeros(size(U,1),1);

point_count = zeros(size(U,1),1);

for i = candidates
    
    [dist_row,~,~] = graphshortestpath(Graph,i,'Directed',false);

    
    points = U(dist_row < NEIGHB_DISTANCE,:);
    I = find(dist_row < NEIGHB_DISTANCE);
        
    num_points = size(points,1);
    
    if num_points < 2
        continue;
    end
    
    point_count(I(1:num_points)) = point_count(I(1:num_points)) + 1;
        
    num_neighbors(i) = size(points,1);
    
    [coef,score] = pca(points);
    
    principal_dirs(i,:) = coef(:,1)';
    
    [~,indmax] = max(score(:,1));
    [~,indmin] = min(score(:,1));
    
        
%     if index_I
%         indices = [ I(indmin) I(indmax)];
%     else
%         indices = [ find(cumsum(neighbor_mat(i,:)) == indmin,1) find(cumsum(neighbor_mat(i,:)) == indmax,1)];
%     end
    
    if ismember(I(indmin),candidates)
        possible_endpoints(I(indmin)) = possible_endpoints(I(indmin)) + 1;
    elseif ismember(I(indmax),candidates)
        possible_endpoints(I(indmax)) = possible_endpoints(I(indmax)) + 1;
    end
    
end


ind = possible_endpoints./point_count;
ind(isnan(ind)) = 0;
ind_endpoints = find(ind);

if vis
    try
%         plot(endpoint_probability)
        g = Plot_Surf(Surf,'transpVal',0.0)
        hold on
        plot3(U(:,1),U(:,2),U(:,3),'.b','MarkerSize',25);
        plot3(U(ind_endpoints,1),U(ind_endpoints,2),U(ind_endpoints,3),'.r','MarkerSize',40);
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