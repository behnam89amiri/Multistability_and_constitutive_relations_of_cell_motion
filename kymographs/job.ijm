/*
  Created by Johannes Heyn, Spring 2020
  Macro to run other macro for all files ending with specified suffix in a specified folder.
  To run for all subfolders, please uncomment the appropriate lines below.
  To stop running script, press Esc or space bar really hard!
  If that doesn't work, press Ctrl + C.

	Check whether file is actually empty
  
*/

requires("1.52t");

// fetches user input
#@ File (label = "Input directory", style = "directory") input
#@ File (label = "Output directory", style = "directory") output
#@ String (label = "File suffix", value = ".tif") suffix
#@ String (label = "Identifier", value = "c2") identifier
#@ String (label = "Excluder", value = "XY01") exclude
#@ String (label = "Macro", value = "E:/05_Skripte/ImageJ/Slicer.ijm") macroFile
#@ Integer (label = "Starting from file number", value = 1) startPoint
#@ boolean (label = "Use Batch Mode") batchMode
#@ boolean (label = "Use BioFormats") BioFormats
#@ boolean (label = "Stop after displaying file list") fileListOnly


exclude = String.trim(exclude);
exclude = replace(exclude, " ", "|");
exclude = ".*(" + exclude + ").*";
startPoint -= 1;

if (batchMode) {
	setBatchMode(true);
}


print("");
print(getTimeString());
print("Starting job.ijm macro...");
//run("Monitor Memory..."); //leads to a stop
//run("Monitor Events..."); //opens Event monitors irrespective of already open ones
processFolder(input);
print("Finished macro.");



if(is("Batch Mode")) {
	setBatchMode(false);
}

// function to scan folders/subfolders to find files with correct suffix & identifier
function processFolder(input) {
	print("Executing processFolder...");
	listRaw = getFileList(input);
	listRaw = Array.sort(listRaw);
	list = newArray(0);

	//returns list of identified files
	print("List of files:");
	for (i = startPoint; i < listRaw.length; i++) {
		if(matches(listRaw[i], ".*" + exclude + ".*")) {
			continue;
		}
		if(endsWith(listRaw[i], suffix) &  matches(listRaw[i], ".*" + identifier + ".*")){
			list = append(list, listRaw[i]);
			print((i+1) + ": " + listRaw[i]);
		}
	}
	
	 // if you only want to see what files jobs will work on, stop here:
	 if(fileListOnly) {
	 	return 0;	
	 }
	 

	startTime = getTime();
	for (i = 0; i < list.length; i++) {
		close("*");
		print("");
		print("Working on file " + (i+1) + " / " + list.length);
	    interruptMacro = isKeyDown("space");
    	if (interruptMacro == true) {
        	print("Interrupted by user.");
        	setKeyDown("none");
       	 	break;
    	}
    // subfolders, uncomment if necessary
		// if(File.isDirectory(input + File.separator + list[i]))
		//     processFolder(input + File.separator + list[i]);
	
		processFile(input, output, list[i]);

		fraction = (i+1) / (list.length);
		myProgress(fraction);
		timePeriod = getTime() - startTime;
		timePeriod = timePeriod / (i+1); //average over previous cycles
		timePeriod = timePeriod / 60000; //convert to minutes
		timeLeft = timePeriod * (list.length - i - 1);
		if (timeLeft < 60){
			print("Time left: " + round(timeLeft*10)/10 + "min");
		} else {
			timeLeft /= 60;
			print("Time left: " + round(timeLeft*10)/10 + "h");
		}
  	}
}

function processFile(input, output, file) {
	// open file and evaluate macro on it
	print("Processing: " + input + File.separator + file);

	if(BioFormats) {
		run("Bio-Formats Importer", "open=" + input + File.separator + file + " color_mode=Default view=Hyperstack stack_order=XYCZT");
	} else {
		open(input + File.separator + file);
	}

  returnStatusMacro = runMacro(macroFile, output);
  print("Macro returned: " + returnStatusMacro);
  
  close("*");
  if (returnStatusMacro == "-1") {
  	print("job.ijm aborted by user");
  	exit("job.ijm aborted by user"); 
  }

}

// might be redundand since it's already defined in macros/Library.txt
function append(arr, value) {
   arr2 = newArray(arr.length+1);
   for (i=0; i<arr.length; i++)
      arr2[i] = arr[i];
   arr2[arr.length] = value;
   return arr2;
}

function myProgress(fraction) {
	//prints progress bar to log
	//parameters: 0 <= fraction <= 1

	if (fraction < 0) {
		exit("error: myProgress -> fraction < 0");
	}
	if (fraction > 1) {
		exit("error: myProgress -> fraction > 1");
	}

	showProgress(fraction);

	s="[";
	for (j=0; j<50; j++) {
		if (j < fraction * 50) {
			s = s + "/";
		} else {
			s = s + " ";
		}
	}
	s = s + "]";
	print(s);
}


function getTimeString() {
	 MonthNames = newArray("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec");
     DayNames = newArray("Sun", "Mon","Tue","Wed","Thu","Fri","Sat");
     getDateAndTime(year, month, dayOfWeek, dayOfMonth, hour, minute, second, msec);
     TimeString ="Date: "+DayNames[dayOfWeek]+" ";
     if (dayOfMonth<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+dayOfMonth+"-"+MonthNames[month]+"-"+year+"\nTime: ";
     if (hour<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+hour+":";
     if (minute<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+minute+":";
     if (second<10) {TimeString = TimeString+"0";}
     TimeString = TimeString+second;
     return TimeString;
}
