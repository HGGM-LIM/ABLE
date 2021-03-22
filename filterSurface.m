function filteredSurf = filterSurface(Surf, iterations)

[Trip] = Vert_Neibp(double(Surf.SurfData.faces),size(Surf.SurfData.vertices,1),size(Surf.SurfData.faces,1));
filteredSurf = Surf;
indMat = Trip(:,[1 3:end]);
cellIndMat = num2cell(indMat',1);
cellNeighb = cellfun(@(x) unique(x(x~=0)), cellIndMat, 'UniformOutput', false);
masscenter = filteredSurf.SurfData.vertices;

% if viz
%     figure('units','normalized','outerposition',[0 0 1 1])
%     tsurf(Pial.SurfData.faces,Pial.SurfData.vertices,'FaceAlpha',0.1,'EdgeColor','k','FaceColor',[0.5 0.5 0.5],'EdgeAlpha',0.5);
%     hold on;
%     t = tsurf( ...
%         bsxfun(@plus,size(Pial.SurfData.faces,1)*(0:2),(1:size(Pial.SurfData.faces,1))'), ...
%         Pial.SurfData.vertices(Pial.SurfData.faces,:),'FaceColor',[0.7 0.8 1.0],'EdgeColor','k','EdgeAlpha',0.5, ...
%         'SpecularStrength',0.1,'FaceLighting','phong');
%     tr = tsurf([1 1 1],Pial.SurfData.vertices,'EdgeColor','r','LineWidth',2);
%     hold off;
%     axis equal;
%     view(-62,10);
%     camproj('persp');
%     light('Position',[-100.0,1.0,-1.0],'Style','infinite');
%     set(gca,'Visible','off');
%     set(gcf,'Color','w');
%     drawnow;
% end

for j = 1:iterations
    for i = 1:size(filteredSurf.SurfData.vertices,1)
        masscenter(i,:) = mean(masscenter(cellNeighb{i},:));
    end
    
%     if viz
%         set(t,'Vertices',masscenter(Pial.SurfData.faces,:));
%         drawnow;
%     end
end

filteredSurf.SurfData.vertices = masscenter;

end