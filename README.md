# bapps
Batch PET Preprocessing in SPM8

Manuel Schütze - Feb 2015

- Description:
This batch scripts runs different SPM8 preprocessing routines on pre-aligned nifti PET images in a folder.

- How to use:
1. First trasform your dicom images into nifti images (eg using MRIcron).
2. Rename images to patient IDs.
3. Using SPM8's display function, align images so that the AC is at the origin and the PC is about the same axis. Put all images in a folder inside your project folder (something like Nifti_aligned). 
4. Navigate to the project folder in Matlab (one above the folder you created).
5. Modify options in this script and run it by typing bapps.

Output folder with different steps of the preprocessing are created automatically according to defined names. Make sure to use check reg in SPM8 to check the results.

You can only run some of the steps/routines, as long as the previous step has been run.
