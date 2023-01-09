# Nucleus Tracking
## Summary
main contributor: Christoph Schreiber, christoph_schreiber@gmx.de   
second contact person: Johannes Heyn, j.heyn@physik.uni-muenchen.de

This project builds on work that was published in PNAS: [On the adhesion-velocity relation and length adaptation of motile cells on stepped fibronectin lanes](https://www.pnas.org/doi/10.1073/pnas.2009959118).

Using tiff-files as a stack (2D x time) these scripts allow to:
* extract the intensity of the fluorescently labelled pattern (if the pattern consists
 of only straight lines)
* determine the position of fluorescently labelled particles (e.g. nuclei)
* track those particles over time

## Procedure
The image analysis is done using line_detection.m and lines_track.m.

### Line detection

The line detection algorithm is based on the Hough transformation. To find the lines the distance and width of the lines is used as an input.

###    Cell tracking on lines

Position of cell nuclei that are located on the microlanes detected before are tracked. Filter rules are applied to get single cell tracks in the end.
