/*
* Created by Johannes Heyn, Spring 2020
* Last edited: 12 January 2022
* Macro used to determine the front and rear position of objects, i.e. cells or nuclei, by semi-manual segmentation
*
* Protocol
* See protocol_KymoSegmentation.md
*
*
*Ideas for improvement:
*Save image that the macro is working on in a 'temp' folder and move it to the actual folder in the end if successful.
*Empty temp folder after macro execution in any case.
*/


macro "kymoSegmentation [a]"{
	output = getArgument();
	isCalled = false;
	if (output != "") {
		isCalled = true;
	}

	if(is("Batch Mode")) {
		setBatchMode(false);
	}
	macroStatus = segmentation();

	if (isCalled) {
		continueSegmentation = getBoolean("Do you wish to continue?");
		if (continueSegmentation == false) {
			macroStatus = "-1";
		}
		return macroStatus;
	}
}


function segmentation() {
	returnSegmentation = "-1";
	imageName = getTitle();
	newImageName = replace(imageName, "\\.tif", "");

	// Create directory where output is saved
	dirOrigin = getDir("image");
	dirSave = dirOrigin + File.separator + "segmentation";
	if (!File.isDirectory(dirSave)) {
		File.makeDirectory(dirSave);
		print("Created directory: " + dirSave);
	}

	// Duplicate image
	titlekymo = newImageName +  ".tif";
	run("Duplicate...", "title=" + titlekymo);
	selectWindow(titlekymo);
	close("\\Others");
	pathSaveImage = dirSave + File.separator + titlekymo;
	checkPathSave(pathSaveImage);
	saveAs("Tiff", pathSaveImage);
	setOption("BlackBackground", true);

	// Start segmentation
	run("Gaussian Blur...", "sigma=2");
	setAutoThreshold("Default dark");
	setForegroundColor(255, 255, 255);
	run("Convert to Mask");
	run("Fill Holes");

	// Add pre-selection to ROI Manager
	roiManager("reset");
	roiManager("Show All with labels");
	run("Analyze Particles...", "size=1000-Infinity add");
	close("*");

	open(pathSaveImage);
	roiManager("show all");

	run("Brightness/Contrast...");
	run("Enhance Contrast", "saturated=0.35");
	roiManager("Show None");
	roiManager("Show All");
	roiManager("select", 0);
	manualSegmentation();

	// Save ROI
	roiname = newImageName + ".roi";
	pathSaveROI = dirSave + File.separator + roiname;
	checkPathSave(pathSaveROI);
	saveAs("selection", pathSaveROI);
	print("Saved " + pathSaveROI);

	// Save binarized image
	pathSaveImageBin = dirSave + File.separator + newImageName + "_binarized.tif";
	saveBinIm(pathSaveImageBin, pathSaveROI);

	close("ROI Manager");
	close("*");

	// Control results
	open(pathSaveImage);
	selectWindow(titlekymo);
	open(pathSaveROI);
	keepImage = getBoolean("Do you want to keep the image?");
	if (!keepImage) {
		startAgain = getBoolean("Do you want to start the segmentation of this image again?");
		close("*");
		isDeleted = File.delete(pathSaveImage);
		if (!isDeleted) {
			print("Error: The kymograph couldn't be deleted.");
		}
		isDeleted = File.delete(pathSaveROI);
		if (!isDeleted) {
			print("Error: The roi-file couldn't be deleted.");
		}
		isDeleted = File.delete(pathSaveImageBin);
		if (!isDeleted) {
			print("Error: The binarized image couldn't be deleted.");
		}
		else {
			print("Deleted successfully.");
		}
		if (startAgain) {
			open(dirOrigin + File.separator + imageName);
			segmentation();
		}
	}
	returnSegmentation = "0";
	return returnSegmentation;
}



function sqr(x) {
	return x*x;
}

function append(arr, value) {
	arr2 = newArray(arr.length+1);
	for (i=0; i<arr.length; i++)
	arr2[i] = arr[i];
	arr2[arr.length] = value;
	return arr2;
}


function list(name, a) {
	print(name);
	for (i=0; i<a.length; i++)
	print("   "+a[i]);
}

function checkPathSave(pathSave) {
	fPathSave = pathSave;
	fPath = File.getDirectory(fPathSave);
	if (!File.isDirectory(fPath)) {
		exit(fPath + " is not a directory");
	}
	if (File.exists(fPathSave)){
		fFileName = File.getName(fPathSave);
		msg = "A file with the name " + fFileName + " already exists." +
		" Do you wish to continue anyway?";
		ttl = "Waiting...";
		Dialog.create(ttl);
		Dialog.addMessage(msg);
		Dialog.show();
	}
	return fPathSave;
}

function saveBinIm(pathSaveImageBin, pathROI) {
	newImName = File.getNameWithoutExtension(pathSaveImageBin);
	dirSave = File.getParent(pathSaveImageBin);
	imageBinarizedName = newImName + "_binarized.tif";

	getDimensions(width, height, channels, slices, frames);
	newImage(imageBinarizedName, "8-bit black", width, height, slices);
	open(pathROI);

	setForegroundColor(255, 255, 255);
	run("Fill", "slice");

	checkPathSave(pathSaveImageBin);
	saveAs("Tiff", pathSaveImageBin);
	print("Saved " + pathSaveImageBin);
}

function manualSegmentation() {
	// User updates ROI using the brush tool.
	// Tool size has to be changed by hand; to do that double click on symbol
	setTool("brush");
	waitForUser("Manual Segmentation", "Draw ROI using brush tool. Then hit OK. " +
	"\n - " +
	"\nUsing the ROI manager you can update ROIs " +
	"\nby selecting them first in the ROI manager " +
	"\nthen extending them using the brush tool before clicking on update. " +
	"\nDelete ROIs that don't belong to the kymograph. ");
	if (roiManager("count") < 1) {
		print("No ROI detected");
		manualSegmentation();
	} else {
		roiManager("update");
		if (roiManager("count") == 1) {
			roiManager("select", 0);
		} else {
		if (roiManager("count") > 1) {
			roiManager("deselect");
			roiManager("combine");
			roiManager("select", 0);
	}}}
}
