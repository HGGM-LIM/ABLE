# ABLE: Automated Brain Lines Extraction Based on Laplacian Skeletonization

Work in Progress

This MATLAB tool takes brain MR images processed with FreeSurfer (?h.white, ?h.pial, ?h.curv, ?h.curv.pial, and ?h.aparc.annot) and extracts gyral crown and sulcal fundi lines in the given mesh as pairs of connected vertices in the surface.

## Usage
The process can be invoked using

    Launch_Morph_Processing(pialSurfLHFile, pialCurvLHFile, whiteSurfLHFile, whiteCurvLHFile, annotLHFile, pialSurfRHFile, pialCurvRHFile, whiteSurfRHFile, whiteCurvRHFile, annotRHFile, outputFile)
 
The output is a MATLAB structure containing the lines, depth map, sulcal basins used for the segmentation, and endpoints obtained for both right and left hemispheres of the given subject.

## Installation
First clone the repository with

    git clone --recursive https://github.com/albertofpena/ABLE.git
    
Then, gptoolbox must be compiled, for that

    cd ABLE/dependencies/gptoolbox/mex
    mkdir build
    cd build
    cmake ..
    make
    
Further instructions can be found in https://github.com/alecjacobson/gptoolbox/blob/master/mex/README.md

Finally, just add the dependencies folder to the MATLAB path.

## Docker (Experimental)
An early version of a docker container has been created. The container is accessible executing:
    docker push albertofpena/able:traveldepth
    
This version estimates the depth map using Mindboggle's Travel Depth. The documentation for this part is still lacking and not all the parameters have been exposed. Further improvements can be expected in the following weeks.
