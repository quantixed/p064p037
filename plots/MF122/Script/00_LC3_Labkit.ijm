/*
 * This script is designed to segment autophagosomes in 3D images using a pre-trained Labkit model.
 * The script assumes that the images are in a folder and that the classifier is saved in a separate file.
 * Written by Mary Fesenko 
 */

// close everything
close("*");
print("\\Clear");
run("Clear Results");

setBatchMode("show");

// Assumes that all images to be analysed are in the same folder
#@ File(style="directory", label="Select main folder") mainFolder
#@ String(style="text field", label="File extension") fileExtension
#@ File(style="file", label="Select classifier file") classifierPath


//Creates folders to store masks (for segmentation QC), "object maps" (for 3D object counter outputs QC) and results tables from 3D object counter
mainFolder = mainFolder + "/";
maskFolder = mainFolder + "masks/";

File.makeDirectory(maskFolder);

// Loop over the list and process the images in the targetFolder
fileList = getFileList(mainFolder);
for (i = 0; i < fileList.length; i++) {
	currentFile = fileList[i];
	if (endsWith(currentFile, fileExtension)) {
		maskPath = maskFolder + replace(currentFile, fileExtension, "_mask.tiff");
		resultsPath = resultsFolder + replace(currentFile, fileExtension, "_results.csv");
		objectMapPath = objectMapFolder + replace(currentFile, fileExtension, "_objectMap.tiff");
		
		print("Analysing: " + mainFolder+currentFile);
		
		// Open the image using the Bioformats importer
		run("Bio-Formats", "open=["+mainFolder+currentFile+"] color_mode=Default quiet rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");	
		
		//Normalization before Labkit segmentation
		run("Enhance Contrast...", "saturated=0.35 normalize process_all use");
		run("Apply LUT", "stack");
        run("Grays");
		
		// Runs labkit on this image, using the classifier selected at the start (classifier needs to be trained and saved separately)
		run("Segment Image With Labkit", "input=dogImage segmenter_file=[" +classifierPath + "] use_gpu=false");
		
		
		// Turn the labkit segmentation result into a binary mask 
		// Close shapes and fill holes to reduce "breaks" in the autophagosome ROIs
		// Save the mask for QC
		setThreshold(1, 255);
		run("Convert to Mask", "background=Dark calculate black");
		run("Close-", "stack");
		run("Fill Holes", "stack");
		save(maskPath);
		
		
		//close everything and start again on the next image in mainFolder
		close("*");
		close("*.csv");
		run("Clear Results");
	}
}

print("Segmentation is complete!");