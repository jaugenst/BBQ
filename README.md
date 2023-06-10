# BBQ

## BBQ = Bacterial Burden Quantification

Here are presented 2 alternative quantitative analysis workflows to quantify the bacterial burden inside macrophages from fixed samples stained by immnunofluorescence. The quantification is based on the total bacterial fluorescence per infected cell, or by bacterial volume calculation per infected cells. Each workflow is presented in 2 separated jupyter notebooks. The workflows are optimized for multi-channel z-stacks images, pre-configured for images aquired on a Zeiss microscope generating "czi" files. The code and workflow can be adapted to your own type of image and file format. 

## Installation
It is highly recommended to install all the softwares into a virtual environment so it doesn't interfere with other potential python installations for other softwares. Virtual environment are easily created and managed through Anacaonda. It is also recommended to start by the installation of PyimageJ as it is more dependent on the installation of OpenJDK and Maven and the proper configuration of JAVA_HOME and MAVEN_HOME variables. Cellpose can then be installed inside the created environment (see instruction's URL below).
Requirements and installation instructions :
> PyImageJ : https://github.com/imagej/pyimagej
> Cellpose : https://github.com/mouseland/cellpose
> ImageJ / Fiji : https://imagej.net/software/fiji/downloads
(Fiji is preferable as it containis more built-in functions). A local installation allow more flexibility to use custom built macros of imagej / fiji that could be created and tested outside the pyimagej environment.

The use of jupyter notebooks is recommended : 
> Jupyter notebooks and Jupyter lab installation instructions : https://jupyter.org/

The jupyter interface allows the separation of code that is intuitive and allow the running of block of code separately. This allows to follow and directly visualize the progression of the image treatement step by step when trying to configure the workflow.

## Workflows description

### Bacterial burden quantification by integrated density of fluorescence per cells
This workflow is derived from a method presented __[here](https://doi.org/10.1083/jcb.201603040)__ and __[here](https://doi.org/10.1371/journal.ppat.1010020)__. The improvement proposed here lies in the nuclei detection using Cellpose, followed by a fluorescence quantification by integrating the total bacterial fluorescence per infected cells. The results can also provide a rapid and precise quantification of the infectivity on your sample as the number of detected nuclei provide the number of cells per fields, and the total number of bacterial fluorescence values is the number of infected cells.

The workflow is as follow :
1. Nuclei channel extraction. From the multi-channel z-stack image, a smooth z-projection "max" image is generated and saved with the extension "_blue.tif". 
2. The nuclei channel is segmented using Cellpose. This software contains a pretrained model called "nuclei" adapted for the segmentation of nuclei stained with DAPI (or you favorite nucleus marker). In case of unsatisfactory segmentation due to different nucleus shape or size, the pretrained model can be easily used to train a custom model to improve the results (see Cellpose documentation and tutorials).
3. Bacterial channel extraction. From the multi-channel z-stack image, a z-projection "sum" image is generated and saved. 
4. Voronoi segmentation and cell ROIs generation. The ROIs from cellpose output are generated and used to perform voronoi segmentation. A new set of ROIs is then generated based on the result and the RoiSet is saved.
5. The "sum" image of the bacterial channel is opened and the fluorescence signal is thresholded to remove the background. 
6. The voronoi ROIs are applied on the "sum" image of the bacterial channel. The ROIs are iteratively selected and the RawIntDen is measured using the "Measure" function of ImageJ. The results are pooled for each field of view of each conditions and saved as a "csv" file. 

