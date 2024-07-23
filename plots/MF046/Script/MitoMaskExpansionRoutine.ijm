// close everything
close("*");
print("\\Clear");
run("Clear Results");

setBatchMode("hide");

#@ File(style="directory", label="Select main folder (registered images)") mainFolder

// Create custom results table
Table.create("DATA");

mainFolder = mainFolder + "/";
outputFolder = mainFolder + "results/";
maskFolder = outputFolder + "masks/";

// get a list with all the files inside mainFolder
folderlist = getFileList(mainFolder);

// loop over the list
for (j = 0; j < folderlist.length; j++) {
	// there is no check here that we are finding folders
	targetFolder = folderlist[j];
	// get image list in the targetFolder
	list = getFileList(mainFolder + targetFolder);
	print("Processing directory: " + targetFolder + "\n");

	// loop over the list
	for (i = 0; i < list.length; i++) {
		if (endsWith(list[i], ".tiff")){
			roiManager("reset");
			print("Currently processing: " + list[i] + "\n");
			
			currentFile = list[i];
			// open the image
			run("Bio-Formats", "open="+mainFolder+targetFolder+currentFile+" color_mode=Default quiet rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
			namefile = getTitle();
			namefilewoext = File.getNameWithoutExtension(namefile);
		
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
			imageCalculator("Add create", "Result of "+ATG9,TPD54);
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
							MinBox = n;				
						}
					}
			}
		
			run("Clear Results");
			
			// TPD54 Background
			selectWindow(TPD54); 
			roiManager("select", MinBox);
			roiManager("measure");
			TPD54Background = getResult("Mean", 0);
		
			// ATG9 Background						
			run("Clear Results");
			selectWindow(ATG9); 
			roiManager("select", MinBox);
			roiManager("measure");
			ATG9Background=getResult("Mean", 0);
			
			run("Clear Results");
			roiManager("deselect");
			roiManager("delete");
		
			// END of BG subtraction calc
						
			// Now import mito mask
			open(maskFolder + targetFolder + namefilewoext + "_mask.tiff");
			rename("mask_" + mito);
			mask_mito = "mask_" + mito;
			// check that we have the correct options set
			run("Options...", "iterations=1 count=1 black do=Nothing");
			run("Convert to Mask");
			// this file is the mito mask
			
			// we will make a 1x and an 8x dilated version for subtraction
//			selectWindow(mask_mito);
			run("Duplicate...", "title=mask1e");
			run("Dilate");
			run("Duplicate...", "title=mask8e");
			for (k = 0; k < 7; k++) {
				run("Dilate");
			}
			// now make mask of mito adjacent regions
			imageCalculator("Subtract create", "mask8e","mask1e");
			selectWindow("Result of mask8e");
			rename("adj_" + mito);
			adj_mito = "adj_" + mito;
			
			// Create mito and adj ROI
			selectWindow(mask_mito);
			run("Create Selection");
			roiManager("Add");
			selectWindow(adj_mito);
			run("Create Selection");
			roiManager("Add");
		
			//measure TPD54 signal in each selection
			selectWindow(TPD54); 
			for (k = 0; k < 2; k++) {
				roiManager("Select", k);
				
				index = Table.size;
				// Add results to table
				Table.set("Mean",index,getValue("Mean"));
				Table.set("BG", index, TPD54Background);
				Table.set("Mean-BG",index,getValue("Mean")-TPD54Background);
				Table.set("Area",index,getValue("Area"));
				Table.set("Slice",index,k); // 0 with be mito, 1 will be adj
				Table.set("Channel",index,"TPD54");
				Table.set("Filename", index, namefilewoext);
				Table.set("Condition", index, targetFolder);
			}
			
			//measure ATG9 signal in each selection
			selectWindow(ATG9); 
			for (k = 0; k < 2; k++) {
				roiManager("Select", k);
				
				index = Table.size;
				// Add results to table
				Table.set("Mean",index,getValue("Mean"));
				Table.set("BG", index, ATG9Background);
				Table.set("Mean-BG",index,getValue("Mean")-ATG9Background);
				Table.set("Area",index,getValue("Area"));
				Table.set("Slice",index,k); // 0 with be mito, 1 will be adj
				Table.set("Channel",index,"ATG9");
				Table.set("Filename", index, namefilewoext);
				Table.set("Condition", index, targetFolder);
			}
	
			close("*");		
		}
	}			
}

// Save table
outputtablepath = outputFolder + "allResults.csv";

selectWindow("DATA");
Table.save(outputtablepath);

//print("Yo, it's the Bad Boy Chiller Crew \nHold tight, S Dog on this one \nRKSwitch on production \nWe're not no divvies 'round here mush, geet up \nHold tight, all the boys inside \nYo \n[BBCC EDIT] \n");
setBatchMode("exit and display");