/*
  Created by Johannes Heyn, 04/01/2022

  Macro to accelerate the selection of rectangular areas in a raw kymograph that contain a region of interest (ROI).
  Pre-processing step before the segmentation using kymoSegmentation.ijm
  1. Open an image that is named similar to "Kymograph16-0_BF.tif" (only the last sub-string "_BF.tif" is important).
  1. Run the macro.
  1. It will save all sub-kymos that the user selected in a folder "/selection" which is inside the directory of the kymograph.

  job.ijm-compatible
*/




if(is("Batch Mode")) {
	setBatchMode(false);
}


imageName = getTitle();
dirSource = getDirectory("image");

// If image displays BF, also open nuc. If images displays nuc, skip.
if (matches(imageName, ".*_nuc.*")) { // better practice: exclude nuc-images alltogether in the job list
	print("Image displays nucleus. Will skip.");
	exit
	}
else{
	if (matches(imageName, ".*_BF.*")) {
		imageBF = imageName;
		titleRawBF = getTitle();
		imageNuc = replace(imageName, "_BF", "_nuc");
		pathImNuc = dirSource + File.separator + imageNuc;
		if (File.exists(pathImNuc)) {
				open(pathImNuc);
				titleRawNuc = getTitle();
			} else {
				exit("The file " + pathImNuc + " does not exist.");
			}
		}
	else {
		exit("Could not identify the image as nucleus or BF image.");
		}
	}

dirSave = dirSource + File.separator + "selection";
if (!File.isDirectory(dirSave)) {
	File.makeDirectory(dirSave);
	print("Created directory: " + dirSave);
	}

run("Merge Channels...", "c1=" + titleRawNuc + " c4=" + titleRawBF + " create");

selectWindow("Composite");
IdComposite = getImageID();
run("Brightness/Contrast...");
run("Enhance Contrast", "saturated=0.35");

counter = 1;
selectArea(IdComposite, counter);


function checkPath(fPathSave) {
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

function selectArea(IdComposite, n) {

	selectImage(IdComposite);

	setTool("rectangle");
	waitForUser("Manual Segmentation", "Draw ROI using rectangle. Then hit OK." +
	"\nIf the image doesn't contain any ROI, close it before clicking OK.");
	if(!isActive(IdComposite)){return(-1)};
	// what happens if user doesn't draw anything and still clicks OK?
	// apperantely, the last ROI is applied
	// is there a way to check whether the user has actually done anything?

	run("Duplicate...", "duplicate");
	run("Split Channels");

	selectWindow("C1-Composite-1");
	IdNuc = getImageID();
	run("Grays");

	pathSaveNuc = dirSave + File.separator + replace(imageNuc, "\\.tif", "-" + toString(n) + "\\.tif");
	checkPath(pathSaveNuc);
	saveAs("Tiff", pathSaveNuc);
	close();

	selectWindow("C2-Composite-1");
	IdBF = getImageID();
	run("Grays");

	pathSaveBF = dirSave + File.separator + replace(imageBF, "\\.tif", "-" + toString(n) + "\\.tif");
	checkPath(pathSaveBF);
	saveAs("Tiff", pathSaveBF);
	close();

	continueSelection = getBoolean("Do you wish to select another area within the raw kymo?");

	if (continueSelection) {
		selectArea(IdComposite, n+1);
	} else {
		return 0;
	}
}
