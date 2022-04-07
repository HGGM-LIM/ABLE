function Launch_Morph_Processing(pialSurfLHFile,pialCurvLHFile,whiteSurfLHFile,whiteCurvLHFile,annotLHFile, ...
    pialSurfRHFile,pialCurvRHFile,whiteSurfRHFile,whiteCurvRHFile,annotRHFile, ...
    outputFile)



disp('Estimating Left Hemisphere');

[sulcalLines,gyralCrowns,endPoints,sulcalBasins,depthMap] = Sulcal_Fundi_New_Extraction(pialSurfLHFile,pialCurvLHFile,whiteSurfLHFile,whiteCurvLHFile,annotLHFile);
morph.left.sulci = sulcalLines;
morph.left.gyri = gyralCrowns;
morph.left.endpoints = endPoints;
morph.left.sulcal_basins = sulcalBasins;
morph.left.depthmap = depthMap;


disp('Estimating Right Hemisphere');

[sulcalLines,gyralCrowns,endPoints,sulcalBasins,depthMap] = Sulcal_Fundi_New_Extraction(pialSurfRHFile,pialCurvRHFile,whiteSurfRHFile,whiteCurvRHFile,annotRHFile);
morph.right.sulci = sulcalLines;
morph.right.gyri = gyralCrowns;
morph.right.endpoints = endPoints;
morph.right.sulcal_basins = sulcalBasins;
morph.right.depthmap = depthMap;



save(outputFile,'morph');

end