/*
 * Created by Johannes Heyn, Spring 2020
 * Last edited: 11. March 2022
 * Macro used to create kymographs
 *
 * Requires plugin "KymoResliceWide", https://imagej.net/plugins/kymoreslicewide
 *
 * Protocol
 * 1. open stack
 * 2. draw line along first (= closest to the upper edge of image) cell trajectory
 * 			 or import line positions from linesXY.csv file, see MATLAB script 'lines2csv'
 * 3. press "0" to activate myKymo
 *
 * Notes: Goal is to create a script that runs through all files and imports corresponding line files automatically
 */


macro "myKymo [0]" {
	/*
	 * saves Kymograph and Roi of selection in video
	 * has built-in test mode to check the lines positions'
	 */

	// sets parameters correctly
	run("Line Width...", "line=5");
	setTool("line");

	// coordinates of first line and the y-shift.

	import_coordinates = getBoolean("Do you want to import the coordinates from a csv file?");
	
	del_d = 136;
	binning = 2048 / getHeight();
	del_d = del_d / binning;

	if (import_coordinates) {
	// pay attention to binning factor.
		path_lines = File.openDialog("Choose the corresponding csv file:");//automatise by extracting ID and opening file automatically?
		x1 = readLines(path_lines, 0);
		x2 = readLines(path_lines, 1);
		y1 = readLines(path_lines, 2);
		y2 = readLines(path_lines, 3);
	} else { // Gets user's line and uses its coordinates
		getSelectionCoordinates(xSelec, ySelec);
		x1_0 = xSelec[0];
		x2_0 = xSelec[1];
		y1_0 = ySelec[0];
		y2_0 = ySelec[1];
		del_y = del_d * sqrt(1 + sqr((y1_0 - y2_0) / (x1_0 - x2_0)));

		// Create list of coordinates
		y1 = newArray(1);
		y2 = newArray(1);
		y1[0] = y1_0;
		y2[0] = y2_0;

		n = 0;
		while ((y1[n] + del_y < getHeight()) && (y2[n] + del_y < getHeight())) {
			y1 = append(y1, round(y1_0 + n*del_y));
			y2 = append(y2, round(y2_0 + n*del_y));
			n = n + 1;
		}

		x1 = newArray(1);
		x2 = newArray(1);
		x1[0] = x1_0;
		x2[0] = x2_0;
		for (i = 0; i < y1.length-1; i++) {
			x1 = append(x1, x1_0);
			x2 = append(x2, x2_0);
		}
	}


	imageName = getTitle();
	isTest = getBoolean("Is this a test to check the line's positions?");
	if (!isTest) {
		pos = getNumber("This is the n-th position of the experiment. n:", 1);
		isNucleus = getBoolean("Does the kymograph display the nucleus?");
		dir1 = getDirectory("Choose a Directory");
		if (isNucleus) {
			isNucleusStr = "_nuc";
		} else {
			isNucleusStr = "_BF";
		}

	}

	for (i = 0; i < y1.length; i++) {

		selectWindow(imageName);
		makeLine(x1[i], y1[i], x2[i], y2[i]);
		if (!isTest) {
			if (isNucleus) {
				roiname="roiOfKymo" +d2s(pos,0) + "-" + d2s(i,0) + ".roi";
				saveAs("selection", dir1 + roiname);
			}

			run("KymoResliceWide ", "intensity=Maximum rotate ignore");

			kymoname = "Kymograph" + d2s(pos,0) + "-" + d2s(i,0) + isNucleusStr + ".tif";
			savepath1 = dir1 + kymoname;
			saveAs("Tiff", savepath1);
		} else {
			wait(1000);
			print("line: " + i, ", y1: " + y1[i], ", y2: " + y2[i]);
		}
	}
}





// function definitions

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

function readLines(path, colNr) {
	// returns column  in a csv file
	// returns array type num
	// identified by colNr (0-indexing)

	x1 = newArray(0);
	requires("1.42j");
	lineseparator = "\n";
	cellseparator = ",\t";

	// copies the whole RT to an array of lines
	lines=split(File.openAsString(path), lineseparator);

	// recreates the columns headers
	labels=split(lines[0], cellseparator);
	if (labels[0]==" ")
	  k=1; // it is an ImageJ Results table, skip first column
	else
	  k=0; // it is not a Results table, load all columns
	for (j=k; j<labels.length; j++)
	  setResult(labels[j],0,0);

	for (i=0; i<lines.length; i++) {
	 // print("Line: " + lines[i]);
	  items=split(lines[i], cellseparator);
	  if (colNr >items.length) {
	    	exit("colNr out of bound. The file contains less columns than colNr.");
	    }
	  for (j=k; j<items.length; j++)
	    if (j==colNr) {
	    	item = parseFloat(items[j]);
	    	if (isNaN(item)) {
	    		exit("Parsing Error. Item could not be converted to type float. \nItem: " + items[j]);
	    	}
	    	x1 = append(x1, item);
	    }

	}
	return x1;
}
