/*
 * Created by Johannes Heyn, Spring 2020
 * Last edited: 15. September 2021
 * Testing: 07/12/2021
 * Status: passed
 * 
 * Macro used to create kymographs. Compatible with job.ijm
 * 
 * Requires plugin "KymoResliceWide", https://imagej.net/plugins/kymoreslicewide
 * 
 * saves Kymograph and Roi of selection
 * has built-in test mode to check the lines positions'
 *
 * Protocol
 * 1. detect pattern via MATLAB script 'line_detection', save lines folder inside the directory of the stacks
 * 2. Invoke this macro via job.ijm
 *
 *
 * Notes: 
 * 	- Add intensity measurement for each line
 * 	- add test mode to make sure the coordinates really fit
 *
 */



import_coordinates = 1;
isTest = 0;

imageName = getTitle();
imageNameL = imageName.toLowerCase;
identIdx1 = imageNameL.indexOf("xy");
identIdx2 = imageNameL.indexOf("_", identIdx1);
if (identIdx2<=0) {
	identIdx2 = imageNameL.indexOf(".tif", identIdx1);
}
ident = imageName.substring(identIdx1 + 2, identIdx2);

args = getArgument();
myArgList = args2list(args);

imp_coordArg = parseArgs(myArgList, "imp_coord");
if (imp_coordArg != -1) {
	import_coordinates = imp_coordArg;
}

isTestArg = parseArgs(myArgList, "test");
if (isTestArg != -1) {
	isTest = isTestArg;
}

if (imageNameL.matches(".*bf.*")) {
	isNucleus = 0;
} else if (imageNameL.matches(".*_texasred.*") || imageNameL.matches(".*_mcherry.*")) {
	isNucleus = 1;
} else {
	exit("Can't infer whether stack displays nucleus or BF.");
}



	// sets parameters correctly
	run("Line Width...", "line=5");
	setTool("line");



	del_d = 68;

	if (import_coordinates) {
	
		path_lines = getDirectory("image") + "lines" + File.separator;
		if (File.isDirectory(path_lines) != 1) {
			exit("'" + path_lines + "' does not exist");
		}

	listRaw = getFileList(path_lines);
	listRaw = Array.sort(listRaw);
	lines_file = "";
	for (i = 0; i < listRaw.length; i++) {
		if(endsWith(listRaw[i], ".csv") &&  matches(listRaw[i], ".*" + ident + ".*")){
			lines_file = listRaw[i];
			break;
		}
	}
	if (lines_file == "") {
		exit("No lines file found for position " + ident);
	}

		path_lines += lines_file;
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



	if (!isTest) {
		// save all kymos of one stack in same folder named XY+ident
		imageDir = getDirectory("image");
		dir1 = imageDir + "XY" + ident;
		if(File.isDirectory(dir1)!= 1) {
			File.makeDirectory(dir1);
		}
		dir1 += File.separator;

		if (isNucleus) {
			isNucleusStr = "_nuc";
		} else {
			isNucleusStr = "_BF";
		}
	}

	for (i = 0; i < y1.length-1; i++) {
		selectWindow(imageName);
		makeLine(x1[i], y1[i], x2[i], y2[i]);
		if (!isTest) {
			if (isNucleus) {
				roiname="roiOfKymo" +d2s(ident,0) + "-" + d2s(i,0) + ".roi";
				saveAs("selection", dir1 + roiname);
			}

			run("KymoResliceWide ", "intensity=Maximum rotate ignore");

			kymoname = "Kymograph" + d2s(ident,0) + "-" + d2s(i,0) + isNucleusStr + ".tif";
			savepath1 = dir1 + kymoname;
			saveAs("Tiff", savepath1);
		} else {
			wait(1000);
			print("line: " + i, "y1: " + y1[i], "y2: " + y2[i]);
		}
	}
	makeLine(x1[i], y1[i], x2[i], y2[i]);


return 0;




function args2list(args) {
    if (lengthOf(args)==0) {
        argList2 = 0;
        } 
    else {
    	argList2 = split(args, "");
    }
    // dirty work around
    argList = newArray("0","0");
    argList = Array.concat(argList, argList2);
    return argList;
}

function parseArgs(argList, keyword) {
	// returns argument corresponding to 'keyword' in 'argList' of the form: ["keyword=argument", ...]
	// return -1 if 'keyword' not in 'argList'
    for (i=0; i<argList.length; i++) {
        arg = argList[i];
        if (arg.matches(keyword + "=.*")) {
        	print("It's a match!");
        	print(myKeyword);
        	print(arg);
        	arg = arg.replace(keyword + "=", "");
        	arg = arg.trim();
        	return arg;
        }
    }
    return -1;
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

function readLines(path, colNr) {
	// returns column  in a csv file at 'path' as array type num
	// identified by 'colNr' (0-indexing)

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
