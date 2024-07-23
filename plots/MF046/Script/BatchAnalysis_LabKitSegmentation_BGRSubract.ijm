// close everything
close("*");
print("\\Clear");
run("Clear Results");

setBatchMode("hide");

#@ File    (label = "Labkit classifier", style = "file") LabkitTraining
#@ File(style="directory", label="Select main folder") mainFolder
#@ String(label="Folder prefix") folderPrefix

// Create custom results table
Table.create("DATA");

mainFolder = mainFolder + "/";
outputFolder = mainFolder+"results/";
maskoutputfolder = outputFolder+"masks/";

File.makeDirectory(outputFolder);
File.makeDirectory(maskoutputfolder);

// get a list with all the files inside mainFolder
folderlist = getFileList(mainFolder);

// loop over the list
for (j = 0; j < folderlist.length; j++) {
	// check if folder has same prefix as folderPrefix
		if (startsWith(folderlist[j], folderPrefix)){
				targetFolder = folderlist[j];
				
				// get image list
				list = getFileList(mainFolder+targetFolder);
				
				print("Processing directory: " + targetFolder + "\n");
				
				// mask output folder target
				targetmaskoutputfolder = maskoutputfolder+targetFolder;
				File.makeDirectory(targetmaskoutputfolder);
				
				// loop over the list
				for (i = 0; i < list.length; i++) {
						if (endsWith(list[i], ".tiff")){
								roiManager("reset");
								print("Currently processing: " + list[i] + "\n");
								
								currentFile = list[i];
								
								// open the image
								run("Bio-Formats", "open="+mainFolder+targetFolder+currentFile+" color_mode=Default quiet rois_import=[ROI manager] view=Hyperstack stack_order=XYCZT");
								
								namefile=getTitle();
								namefilewoext=File.getNameWithoutExtension(namefile);
								
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
								TPD54Background=getResult("Mean", 0);
								
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
								
				
								//Mito segmentation
								selectWindow(mito); 
								run("Duplicate...","duplicate title=labkitChannel duplicate");
								
								
								run("Segment Image With Labkit", "input=labkitChannel segmenter_file="+ LabkitTraining + " use_gpu=false");
								setMinAndMax(0, 1);
								setAutoThreshold("Huang dark");
								run("Convert to Mask");	
								
								//Save mito masks for segmentation QC
								saveAs("tiff", targetmaskoutputfolder+namefilewoext+"_mask.tiff");
								
								
								rename("mask_"+mito);
								mask_mito="mask_"+mito;
								
								
								//Create mito ROI for each frame
								run("Create Selection");
								print(nSlices + "\n ");
								for (k = 1; k <= nSlices; k++) {
									setSlice(k);
									run("Create Selection");
								    Roi.setPosition(k);
								    roiManager("Add");
								}
								
								//measure TPD54 signal in mito selection
								selectWindow(TPD54); 
								for (m = 1; m <= nSlices; m++) {
									setSlice(m);
									roiManager("Select", m-1);
									
									index = Table.size;
									// Add results to table
									Table.set("Mean",index,getValue("Mean"));
									Table.set("BG", index, TPD54Background);
									Table.set("Mean-BG",index,getValue("Mean")-TPD54Background);
									Table.set("Area",index,getValue("Area"));
									Table.set("Slice",index,m);
									Table.set("Channel",index,"TPD54");
									Table.set("Filename", index, namefilewoext);
									Table.set("Condition", index, targetFolder);
								}
								
								//measure ATG9 signal in mito selection
								selectWindow(ATG9); 
								for (n = 1; n <= nSlices; n++) {
									setSlice(n);
									roiManager("Select", n-1);
									
									index = Table.size;
									// Add results to table
									Table.set("Mean",index,getValue("Mean"));
									Table.set("BG", index, ATG9Background);
									Table.set("Mean-BG",index,getValue("Mean")-ATG9Background);
									Table.set("Area",index,getValue("Area"));
									Table.set("Slice",index,n);
									Table.set("Channel",index,"ATG9");
									Table.set("Filename", index, namefilewoext);
									Table.set("Condition", index, targetFolder);
								}
								
								close("*");								
						}
				}		
				print("==================\n");				
		}
}

// Save table
outputtablepath = outputFolder + "results.csv";

selectWindow("DATA");
saveAs("Results", outputtablepath);


print("Yo, it's the Bad Boy Chiller Crew \nHold tight, S Dog on this one \nRKSwitch on production \nWe're not no divvies 'round here mush, geet up \nHold tight, all the boys inside \nYo \n[BBCC EDIT] \n");
setBatchMode("exit and display");