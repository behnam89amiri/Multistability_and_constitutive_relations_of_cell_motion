# Protocol
How to segment kymographs using the ImageJ macro 'kymoSegmentation.ijm'.

## Requirements
* ImageJ
* kymoSegmentation.ijm script installed
* Kymograph as a tif-file (for the creation of kymographs see kymoCreation.ijm)

## Step by step guide for a single kymograph
1. Open the two corresponding channels of the kymograph, i.e. nucleus- and BF-channel
1. Select the area containing the region of interest (ROI) on both the nucleus and BF image. Kymographs must only intersect with the left and right border of the selection
1. ctrl + shift + 'd' duplicates the frame of interest
1. Save frame using ctrl + 's'
1. Select the nucleus file and click shift + 'e'
1. Duplicate (ctrl + shift + 'd') and save the selection (ctrl + 's')
1. Run the macro 'kymoSegmentation' either via the ImageJ editor or, if already installed, via keyboard shortcut [a]
1. Adjust brightness / contrast to best visualise the edges
1. Select ROI with selection brush tool
1. Zoom in and out onto edges using crtl + mouse wheel
1. Click OK. The macro saves ROI in the .roi-file format to folder 'segmentation' in the directory of the kymograph
1. To convert the segmentation from .roi to a mask, use 

## Working on a batch of kymographs
If working on a batch of kymographs, you can use job.ijm to apply macros in a repeated fashion on all kymos within a folder.
1. Start job.ijm to sequentally open all kymographs that are being pre-processed, i.e. where the area containing the ROI needs to be selected
1. Start job.ijm to segment all pre-processed kymographs one by one. 

## To dos
[ ] Labkit
[ ] Open both channels at the same time and select the same ROI
[ ] Apply filters to visualise the image the best
