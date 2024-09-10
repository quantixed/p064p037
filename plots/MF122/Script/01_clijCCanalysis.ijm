/*
 * The aim of this macro is to use take each binary stack (in this case LC3 images segemented with labkit)
 * Find the objects and get the statistics. Note that the output is in isotropic pixels so needs conversion
 * by xy_pixel_size ^ 3
 * 
 * Uses the 'connected components analysis (box)' function from the clij2 plugin suite to label objects from the masks.
 * Requires the user to select the 'masks' folder generated from the LC3 labkit segmentation script which needs to be run before this one.
 * 
 * Written by Mary Fesenko and Stephen Royle
 */

// clean up first
run("Close All");
run("Clear Results");

// Assumes that all images to be analysed are in a subfolder called objectMaps in mainFolder
#@ File(style="directory", label="Select main folder") mainFolder

mainFolder = mainFolder + "/";

NEWresultsFolder = mainFolder + "NEW_clijCC_results/";
NEW_LC3labelsFolder = mainFolder + "NEW_LC3labelsFolder/";

File.makeDirectory(NEWresultsFolder);
File.makeDirectory(NEW_LC3labelsFolder);


run("CLIJ2 Macro Extensions", "cl_device=");

setBatchMode(true);

fileList = getFileList(mainFolder);


for (i = 0; i < fileList.length; i++) {
	currentFile = fileList[i];
	if (!endsWith(currentFile, "_mask.tiff")) continue;
	
	maskPath = mainFolder + currentFile;
	NEWresultsPath = NEWresultsFolder + replace(currentFile, "_mask.tiff", "_results.csv");

	print("Analysing: " + i + " - " + currentFile);
	open(maskPath);
	
	namefile=getTitle();
	namefilewoext=File.getNameWithoutExtension(namefile);

	Ext.CLIJ2_clear();
	stack = getTitle();
	Ext.CLIJ2_push(stack);
	
	// read out image size
	Ext.CLIJ2_getDimensions(stack, width, height, depth);
	print("Image size:", width, height, depth );
	
	// read voxel size from original image using ImageJ functions
	selectWindow(stack);
	getVoxelSize(voxel_width, voxel_height, voxel_depth, unit);
	print("Voxel size:", voxel_width, voxel_height, voxel_depth);
	
	// assume voxel_depth is bigger than voxel_width and voxel_height, and that voxel_width and voxel_height are equal
	factor = 1 / voxel_width;
	voxel_width = voxel_width * factor;
	voxel_height = voxel_height * factor;
	voxel_depth = voxel_depth * factor;
	
	// create another stack of different size
	Ext.CLIJ2_create3D(isotropic_stack, width * voxel_width, height * voxel_height , depth * voxel_depth, 16);
	
	// Scale the image so that it becomes [isotropic](https://en.wikipedia.org/wiki/Isotropy)
	Ext.CLIJ2_scale3D(stack, isotropic_stack, voxel_width, voxel_height, voxel_depth, false);
	

	// connected components labeling
	connected_stack = "connected_components_labeling";
	Ext.CLIJ2_connectedComponentsLabelingBox(isotropic_stack, connected_stack);
	run("glasbey_on_dark");
	
	//Make the object labels from the 'connected components analysis' into a mask and save this file for QC
	Ext.CLIJ2_mask(connected_stack, connected_stack, label_masks);
	// visualize labeled image
	Ext.CLIJ2_pull(label_masks);
	run("glasbey_on_dark");
	resetMinAndMax;
	saveAs("Tiff",NEW_LC3labelsFolder+namefilewoext+"_LC3labels.tiff");
		
	//Record various measurments from the labelled objects and save the results table
	Ext.CLIJ2_statisticsOfLabelledPixels(connected_stack, connected_stack);
	saveAs("Results", NEWresultsPath);
	
	// clean up
	run("Close All");
	run("Clear Results");
}
setBatchMode(false);
print("Done!");