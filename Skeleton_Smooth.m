function [U] = Skeleton_Smooth(V,F)
  % SKELETON_EXTRACTION Compute the "skeleton" of a surface mesh following "Skeleton
  % Extraction by Mesh Contraction" [Au et al. 2008]. The final combinatorially
  % simplification is a bit different than the QSlim-based approach in [Au et
  % al. 2008] without the thinning
  % 
  % [U] = Skeleton_Smooth(V,F)
  % 
  % Inputs:
  %   V  #V by 3 list of mesh vertex positions
  %   F  #F by 3 list of mesh indices
  % Outputs:
  %   U  #V by 3 list of skeleton vertex positions, before thinning
  %

  delta = 1e-3*max(max(V)-min(V))^2;
  U = V;
  
  viz = false;
  if viz
      figure('units','normalized','outerposition',[0 0 1 1])
      tsurf(F,V,'FaceAlpha',0.1,'EdgeColor','k','FaceColor',[0.5 0.5 0.5],'EdgeAlpha',0.5);
      hold on;
      t = tsurf( ...
          bsxfun(@plus,size(F,1)*(0:2),(1:size(F,1))'), ...
          U(F,:),'FaceColor',[0.7 0.8 1.0],'EdgeColor','k','EdgeAlpha',0.5, ...
          'SpecularStrength',0.1,'FaceLighting','phong');
      tr = tsurf([1 1 1],U,'EdgeColor','r','LineWidth',2);
      hold off;
      axis equal;
      view(-62,10);
      camproj('persp');
      light('Position',[-100.0,1.0,-1.0],'Style','infinite');
      set(gca,'Visible','off');
      set(gcf,'Color','w');
      drawnow;
%       a = gca;
%       temp = getframe(a.Parent);
%       imwrite(temp.cdata,['/media/DATOS/afernandez/matlab-utils/MyToolboxes/Sulcal_Width_Computation/figures/skeletonSmooth_0.tiff']);
  end

  b = [];

  iter = 1;
  h = avgedge(V,F);

  while true
    L = cotmatrix(U,F);
    M = massmatrix(U,F,'barycentric');
    b = union(find(any(abs(L)>1e+5,2)),b);
    U_prev = U;
    U = min_quad_with_fixed(-(M-delta*L),2*M*U,b,U(b,:));
    %U = U/sqrt(sum(doublearea(U,F))*0.5);
    %U = bsxfun(@minus,U,centroid(U,F));
    %D = max(internalangles(U,F),[],2)>0.9*pi & ...
    %  doublearea(U,F)<1e-5;
    D = all(ismember(F,b),2);
    if viz
        set(t,'Vertices',U(F,:));
        set(tr,'Faces',F(D,:),'Vertices',U);
        xlabel(sprintf('%d',numel(b)));
        drawnow;
        
%         a = gca;
%         temp = getframe(a.Parent);
%         imwrite(temp.cdata,['/home/afernandez/vboxwindows/presentaciones/20190412/skeleton/skeletonSmooth_' num2str(iter) '.tiff']);
    end
    if max(max(abs(U-U_prev))) < 1e-2*h
      break;
    end
    iter =iter + 1;
  end
  
end
