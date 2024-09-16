/*
 * This script is a version of the MitoMaskExpansionRoutine for live cells
 * It requires that the mitochondrial masks (for movies) have been generated
 * and that NanoJ processed images are available
 * It will measure bg, mito, mito_adj and cell for 2 channels for all timepoints
 */

// close everything
close("*");
print("\\Clear");
run("Clear Results");

setBatchMode("hide");

#@ File(style="directory", label="Select main folder (registered images)") mainFolder

// Create custom results table
Table.create("DATA");

// contains registered movies
mainFolder = mainFolder + "/";
// contains masks folder, outputs saved here incl timepoints from another script
outputFolder = mainFolder + "results/";
// contains masks that have an ROI marking the cell
maskFolder = outputFolder + "masks/";

// get a list with all the files inside mainFolder
list = getFileList(mainFolder);

// loop over the list
for (i = 0; i < list.length; i++) {
	if (endsWith(list[i], ".tiff")){
		roiManager("reset");
		print("Currently processing: " + list[i] + "\n");

		currentFile = list[i];

		// open the image
		run("Bio-Formats", "open="+mainFolder+currentFile+" color_mode=Default quiet rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
		namefile = getTitle();
		namefilewoext = File.getNameWithoutExtension(namefile);
		if(!File.exists(maskFolder + namefilewoext + "_mask.tiff")){
			close("*");
			continue;
		}
	
		run("Split Channels");

		selectWindow("C1-"+namefile); 		
		rename("mito_"+namefilewoext);
		mito="mito_"+namefilewoext;
	
		selectWindow("C3-"+namefile); 
		rename("TPD54_"+namefilewoext);
		TPD54="TPD54_"+namefilewoext;
	
		selectWindow("C2-"+namefile); 
		rename("ATG9_"+namefilewoext);
		ATG9="ATG9_"+namefilewoext;
	
		// Dan's background subtraction code
		imageCalculator("Add create", ATG9,mito);
		firstImageCreated = getTitle();
		imageCalculator("Add create", "Result of "+ATG9,TPD54);
		secondImageCreated = getTitle();
		close(firstImageCreated);
		selectWindow(secondImageCreated);
		if(nSlices() > 1) {
			run("Z Project...", "projection=[Max Intensity]");
			thirdImageCreated = getTitle();
			close(secondImageCreated);
			secondImageCreated = thirdImageCreated;
		}
		
		x = 0;
		y = 0;
		while ( y < getWidth()) {
			makeRectangle(x, y, 50, 50);
			roiManager("Add");
			y=y+50;									
		}
		y=0;
		x=50;
		while ( x < getWidth()) {
			makeRectangle(x, y, 50, 50);
			roiManager("Add");
			x=x+50;
			
		}
		y=50;
		x=x-50;
		while ( y < getWidth()) {
			makeRectangle(x, y, 50, 50);
			roiManager("Add");
			y=y+50;
			
		}
		y=y-50;
		x=50;
		while ( x < (getWidth()-50)) {
			makeRectangle(x, y, 50, 50);
			roiManager("Add");
			x=x+50;
			
		}
//		ns = nSlices();
		run("Set Measurements...", "mean redirect=None decimal=0");
		roiManager("multi-measure measure_all");
		
		for(n=0;n<nResults;n++){
			Int=getResult("Mean", n);
				if(n==0){
				Min=Int;
				MinBox = 0;
				}
				if(n>0){
					if (Int<Min) {
						Min = Int;
						MinBox = n;				// gives the row number
//						MinBox = n%ns;          // gives the row number if it's a stack
					}
				}
		}
		
		run("Clear Results");
		close(secondImageCreated);
		
		// add Background ROI to ROI manager
		selectWindow(TPD54); 
		roiManager("select", MinBox);
		run("Clear Results");
		roiManager("deselect");
		roiManager("delete");
		roiManager("Add"); // 0 = bg
		
		// END of BG subtraction spec
					
		// Now import mito mask
		open(maskFolder + namefilewoext + "_mask.tiff");
		rename("mask_" + mito);
		mask_mito = "mask_" + mito;
		// take the ROI that should be in file - no check here!
		roiManager("Add"); // 1 = cell(s)
		
		// check that we have the correct options set
		run("Options...", "iterations=1 count=1 black do=Nothing");
		// this line is not needed because the images are already correct 0,255
		// Usage of ROI in this step is essential, if it is not present the dilation will be a disaster
//		run("Convert to Mask", "method=Default background=Default black");
		// this file is the mito mask
		
		// we will make a 1x and an 8x dilated version for subtraction
		selectWindow(mask_mito);
		// deselect ROI otherwise the duplicated movie will be the wrong size
		run("Select None");
		run("Duplicate...", "title=mask1e duplicate");
		// select cell ROI again
		roiManager("Select", 1);
		run("Dilate", "stack");
		// now deselect again
		run("Select None");
		run("Duplicate...", "title=mask8e duplicate");
		// and reselect!
		roiManager("Select", 1);
		for (k = 0; k < 7; k++) {
			run("Dilate", "stack");
		}
		// now make mask of mito adjacent regions
		imageCalculator("Subtract create stack", "mask8e","mask1e");
		selectWindow("Result of mask8e");
		rename("adj_" + mito);
		adj_mito = "adj_" + mito;
		
		// Create mito and adj ROI
		selectWindow(mask_mito);
		nS = nSlices();
		for (kk = 1; kk <= nS; kk++) {
			setSlice(kk);
			run("Create Selection");
			Roi.setPosition(kk);
			roiManager("Add"); // 2 = mito
		}
		
		selectWindow(adj_mito);
		for (kk = 1; kk <= nS; kk++) {
			setSlice(kk);
			run("Create Selection");
			Roi.setPosition(kk);
			roiManager("Add"); // 3 = mito_adj
		}
	
		// we have
		// ii channels (2)
		// jj ROIs (2 + (2 * frames))
		// kk frames (variable)
		nROI = roiManager("count");
		
		for (ii = 0; ii < 2; ii++) {
			if (ii == 0) {
				selectWindow(TPD54);
			} else {
				selectWindow(ATG9); 
			}
			for (jj = 0; jj < 2; jj++) {
				
				for (kk = 1; kk <= nS; kk++) {
					roiManager("Select", jj);
					setSlice(kk);

					index = Table.size;
					// Add results to table
					Table.set("Mean",index,getValue("Mean"));
					Table.set("Area",index,getValue("Area"));
//						Table.set("ROI",index,jj);
					if (jj < 2) {
						// if roi is bg or cell, add that integer
						Table.set("ROI",index,jj);
					} else {
						//
						Table.set("ROI",index,2 + floor((jj - 2) / nS));
					}
					if (ii == 0) {
						Table.set("Channel",index,"TPD54");
					} else {
						Table.set("Channel",index,"ATG9");
					}
					Table.set("Frame",index,kk);
					Table.set("Filename", index, namefilewoext);
				}
			}
		}
		
		for (ii = 0; ii < 2; ii++) {
			if (ii == 0) {
				selectWindow(TPD54);
			} else {
				selectWindow(ATG9); 
			}
			for (jj = 2; jj < nROI; jj++) {
				roiManager("Select", jj);

				index = Table.size;
				// Add results to table
				Table.set("Mean",index,getValue("Mean"));
				Table.set("Area",index,getValue("Area"));
//						Table.set("ROI",index,jj);
				if (jj < 2) {
					// if roi is bg or cell, add that integer
					Table.set("ROI",index,jj);
				} else {
					//
					Table.set("ROI",index,2 + floor((jj - 2) / nS));
				}
				if (ii == 0) {
					Table.set("Channel",index,"TPD54");
				} else {
					Table.set("Channel",index,"ATG9");
				}
				Table.set("Frame",index,getSliceNumber());
				Table.set("Filename", index, namefilewoext);
			}
		}
		
		// Save results for this image in the results folder
		resultsfilename = outputFolder + namefilewoext + "_results.csv";
		Table.save(resultsfilename);
		
		close("*");		
	}
}

// Save table
outputtablepath = outputFolder + "allResults.csv";

selectWindow("DATA");
Table.save(outputtablepath);

setBatchMode("exit and display");