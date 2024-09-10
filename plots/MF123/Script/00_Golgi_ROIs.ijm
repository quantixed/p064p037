/*
 * This script is designed to segment the Golgi apparatus in images of ATG9-labelled cells.
 * The script assumes that the images are in a folder and that the classifier is saved in a separate file.
 * Written by Mary Fesenko 
 */


// close everything
close("*");
print("\\Clear");
run("Clear Results");

setBatchMode("hide");

#@ File(style="directory", label="Select main folder (registered images)") mainFolder
#@ File    (label = "Labkit classifier", style = "file") LabkitTraining

mainFolder = mainFolder + "/";
resultsFolder = mainFolder + "results/";
maskFolder = mainFolder + "masks/";

File.makeDirectory(resultsFolder);
File.makeDirectory(maskFolder);

// get a list with all the files inside mainFolder
folderlist = getFileList(mainFolder);

// get image list in the targetFolder
list = getFileList(mainFolder);

file_extension = ".tif";

// loop over the list and process the images in the targetFolder
for (i = 0; i < list.length; i++) {
	currentFile = list[i];
	if (endsWith(currentFile, file_extension)){

		resultsPath = resultsFolder + replace(currentFile, file_extension, "_results.csv");
		
		roiManager("reset");
		print("Currently processing: " + list[i] + "\n");

		open(mainFolder+currentFile);
		
		namefile=getTitle();
		namefilewoext=File.getNameWithoutExtension(namefile);
		
		// Max intensity projection
		run("Z Project...", "projection=[Max Intensity]");
		
		//Channel splitting and renaming
		run("Split Channels");

		selectWindow("C1-MAX_"+currentFile); 		
		rename("ATG9_to_measure");
	
		selectWindow("C2-MAX_"+currentFile); 
		rename("golgi_for_roi");
	
		//Median filter to reduce backgrund graininess
		run("Median...", "radius=5");
	
		//Normalization before Labkit segmentation
		run("Enhance Contrast", "saturated=0.35");
        run("Apply LUT", "stack");
        
        //Labkit segmentation
		run("Segment Image With Labkit", "input=golgi_for_roi segmenter_file=["+ LabkitTraining + "] use_gpu=false");
		
		//Convert Labkit segmentation result into a mask
		setThreshold(1, 255);
		run("Convert to Mask");
		
		//Save Golgi masks for segmentation QC
		saveAs("Tiff",maskFolder+namefilewoext+"_mask.tiff");
		
		//checks some objects have been detected after segmentation
		//if yes, uses the Golgi mask to create an ROI, moves the Golgi ROI onto the ATG9 channel, 
		//records mean signal intensity inside Golgi ROI and saves the result
		max_value_check = getValue("Max");
		
		if (max_value_check>0) {		
			run ("Create Selection");
			roiManager("Add"); 
			run("Set Measurements...", "area mean redirect=None decimal=3");
			selectWindow("ATG9_to_measure");
			roiManager("select", 0);
			run("Measure");
			saveAs("Results", resultsPath);
		} else {
			print(currentFile + "There's nothing here. Just a dark empty void.");
		}
		
		//close everything and start again on the next image in mainFolder
		close("*");
		run("Clear Results");
		roiManager("reset");
	}
}

print("I look into the finance box just to check my status \nI look into the microscope, I see the Golgi Apparatus \nGolgi, oh, woe is me, you can't even see the sea \nGolgi, Golgi, Goh-oh-olgi, Golgi segmentation is complete \n - inspired by the song 'Golgi Apparatus' by Phish"); 