function [labels] = Surface_Segmentation(Surf,curvmap,depthMap)

labels = (curvmap > 0) & (depthMap > 1);

% NUM_CLUSTERS = 2;
% ZSCORE_CORRECTION = 5;
% FCM_THRESHOLD = 0.5;
% 
% logDepth = log(depthMap+1);
% indzeros = logDepth == 0;
% logDepth(indzeros) = [];
% zDepthmap = zscore(logDepth);
% auxCurv = curvmap(~indzeros);
% auxDepth = depthMap(~indzeros);
% zCurvmap = zscore(auxCurv);
% 
% % calculate outlier correction for curvature maps
% curvmap_deviation = max(zCurvmap) - min(zCurvmap);
% count = 1;
% while curvmap_deviation > 1.05*ZSCORE_CORRECTION*2
%     disp(['Outlier correction for curv map, iteration ' num2str(count)]);
%     zCurvmap(zCurvmap > ZSCORE_CORRECTION) = ZSCORE_CORRECTION;
%     zCurvmap(zCurvmap < -ZSCORE_CORRECTION) = -ZSCORE_CORRECTION;
%     
%     zCurvmap = zscore(zCurvmap);
%     
%     curvmap_deviation = max(zCurvmap) - min(zCurvmap);
%     count = count + 1;
% end
% 
% 
% ind = ~indzeros;
% Surfaux = Surf;
% indfaces = ismember(Surf.SurfData.faces(:,1),find(ind)) & ismember(Surf.SurfData.faces(:,2),find(ind)) & ismember(Surf.SurfData.faces(:,3),find(ind));
% Surfaux.SurfData.faces = Surf.SurfData.faces(indfaces,:);
% [C,~] = connected_components([Surfaux.SurfData.faces;repmat(size(Surfaux.SurfData.vertices,1),1,3)]);
% sizeCompo = accumarray(C',C'*0+1);
% comp2del = find(sizeCompo <= 50);
% C(find(ismember(C,comp2del))) = 0;
% 
% listComp = unique(C);
% listComp(listComp == 0) = [];
% 
% map = zeros(1,length(indzeros));
% Surf.Is = C';
% Surfaux.Is = C';
% 
% for i = listComp
%     ind = C(~indzeros) == i;
%     
%     zDMaux = zscore(zDepthmap(ind));
%     zCMaux = zscore(zCurvmap(ind));
%     
%     [fcm_center,U,~] = fcm_fix_seed([zDMaux zCMaux],NUM_CLUSTERS);
%     [~,I] = max(fcm_center(:,2));
%     
% %     [U,k_centroid] = kmeans([zDMaux zCMaux],NUM_CLUSTERS);
% %     [~,I] = max(k_centroid(:,2));
% 
% %     [subSurf,~,reind] = Extract_Sub_Surface(Surf,i);
% %  
% %     sulcSeg = lacalentada(subSurf,curvmap(reind{1}),depthMap(reind{1}));
% %     auxMap(~ind) = 0;
% %     auxMap(ind) = sulcSeg;
% 
%     auxMap(~ind) = 0;
%     auxMap(ind) = U(I,:) > FCM_THRESHOLD;
%     
%     map(~indzeros) = map(~indzeros) + auxMap;
% end
% 
% 
% % [~,U,~] = fcm_fix_seed([zDepthmap zCurvmap],NUM_CLUSTERS);
% % a = kmeans([zDepthmap zCurvmap],NUM_CLUSTERS);
% 
% % sts = unique(a);
% % means = [];
% % for i = 1:length(sts)
% %     ind = find(a == sts(i));
% %     means(i) = mean(auxCurv(ind));
% % end
% 
% % [~,I] = max(means);
% % map(~indzeros) = (a == sts(I));
% 
% % means = [];
% % for i = 1:NUM_CLUSTERS
% %     ind = find(U(i,:) > FCM_THRESHOLD);
% %     means(i) = mean(auxCurv(ind));
% % end
% % 
% % [~,I] = max(means);
% % map(~indzeros) = (U(I,:) > FCM_THRESHOLD);
% 
% % % numComp = [];
% % % numNonZeros = [];
% % % 
% % % for i = 0:0.005:1
% % %     mapAux(indzeros) = 0;
% % %     mapAux(~indzeros) = (U(I,:) > i);
% % %     mapAux = mapAux';
% % %     ind = (mapAux == 1);
% % %     
% % %     Surfaux = Surf;
% % %     indfaces = ismember(Surf.SurfData.faces(:,1),find(ind)) & ismember(Surf.SurfData.faces(:,2),find(ind)) & ismember(Surf.SurfData.faces(:,3),find(ind));
% % %     Surfaux.SurfData.faces = Surf.SurfData.faces(indfaces,:);
% % %     [C,~] = connected_components([Surfaux.SurfData.faces;repmat(size(Surfaux.SurfData.vertices,1),1,3)]);
% % %     
% % %     sizeCompo = accumarray(C',C'*0+1);
% % %     comp2del = find(sizeCompo <= 50);
% % %     C(find(ismember(C,comp2del))) = 0;
% % %     
% % %     numNonZeros = [numNonZeros; length(nonzeros(C))];    
% % %     numComp = [numComp; length(unique(C))];
% % % 
% % %     clear mapAux;
% % % end
% % % 
% % % 
% % % [~,dynThreshold] = max(numComp.*numNonZeros/max(numNonZeros)*max(numComp));
% % % dynThreshold = dynThreshold * 0.005;
% % % disp(['Setting dynamic threshold at ' num2str(dynThreshold)]);
% % % map(~indzeros) = (U(I,:) > dynThreshold);
% 
% map = map';
% 
% map = Erode_Surface_Label(Surf,map);
% ind = (map > 0);
% 
% % relabel
% Surfaux = Surf;
% indfaces = ismember(Surf.SurfData.faces(:,1),find(ind)) & ismember(Surf.SurfData.faces(:,2),find(ind)) & ismember(Surf.SurfData.faces(:,3),find(ind));
% Surfaux.SurfData.faces = Surf.SurfData.faces(indfaces,:);
% [C,~] = connected_components([Surfaux.SurfData.faces;repmat(size(Surfaux.SurfData.vertices,1),1,3)]);
% 
% non_repeated = find(histc(C(:), 1:max(C)) == 1);
% C(ismember(C,non_repeated)) = 0;
% 
% tempVar = nonzeros(C(:));
% sizeCompo = accumarray(tempVar,tempVar*0+1);
% comp2del = find(sizeCompo <= 10);
% C(find(ismember(C,comp2del))) = 0;
% 
% labels = C';

end
