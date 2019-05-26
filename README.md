# READ-ME

The code in this repository can be used to reproduce the “Repressilator-based inference of growth at a single-cell level” as detailed in Riglar et al, “Bacterial variability in the mammalian gut captured by a single-cell synthetic oscillator” (currently under review and released as a preprint doi: https://www.biorxiv.org/content/10.1101/472720v1). 

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.

### Prerequisites

MATLAB v2017a or later (may work on previous versions of MATLAB, but has not been tested by the authors).

```
Todo: provide a list of toolboxes that are required to run the scripts:
https://in.mathworks.com/help/matlab/ref/matlab.codetools.requiredfilesandproducts.html
```

### Installing

1. Add Riglar_RINGS_run.m and Riglar_RINGS_fit.m to the same folder.
2. Add the folder to the MATLAB path

### Example Usage

**Riglar_RINGS_run.m**

This function will prompt you to select a parent folder containing folders named ‘YFP’ and ‘CFP’. Both YFP and CFP should contain a number of parallel subfolders, corresponding to different timepoints. 

The code loops over all subfolders, loads the paired YFP and CFP images, and calls the fitting function: Riglar_RINGS_fit.m

Image data is in the form of single fluorescent jpg images of individual bacterial colonies, centered upon the colony. For this study, colonies were identified and exported as 200x200 pixel images from macroscope images of agar plates using FIJI (ImageJ).

Notes:

1. This function is hard-coded to expect a certain directory structure. Please abide by the following rules: **(i)** The plate names cannot contain ‘.’ or ‘-‘.  **(ii)** The code expects that for every YFP image there is a corresponding CFP image, and the file names differ only by replacing ‘YFP’ with ‘CFP’

2. The method outputs a csv file with all fitted theta_0 values for CFP, organized by folders.  Theta_0_YFP can be calculated as follows: theta_0_YFP = theta_0_CFP + colorPhaseShift.

3. The csv file is saved within the folder chosen at the prompt at the completion of the script. 

4. A simple visualization of the datapoints [0:2pi], y-axis, arranged by timepoint, x-axis is also returned.

Default Parameters:

1. _Default parameters for E. coli MC4100 (eg LPT320)_  
```Riglar_RINGS_run()```  
or  
```Riglar_RINGS_run(‘rMax’, 1.3, ‘expectedSlope’, 0.39, ‘slopeTol’, 0.3, ‘colorPhaseShift’, 1.5)```

2. _Default parameters for E. coli MG1655 (eg. PAS715)_  
```Riglar_RINGS_run(‘rMax’, 1.3, ‘expectedSlope’, 0.34, ‘slopeTol’, 0.4, ‘colorPhaseShift’, 0.9)```

3. _Default parameters for S. Typhimurium (eg. PAS716)_  
```Riglar_RINGS_run(‘rMax’, 1.1, ‘expectedSlope’, 0.43, ‘slopeTol’, 0.4, ‘colorPhaseShift’, 1.0)```

**Riglar_RINGS_fit.m**

This code is the workhorse that does the model fitting, and evaluates whether or not the fit parameters are outside an acceptable range. 

There are 4 self-identified failure modes.
1. Center of colony too far from center of image.
2. Slope is not within tolerance.
3. Uncertainty in fitted slope is not within tolerance.
4. CFP and YFP amplitude have opposite sign.  Fit is unrealiable.

The primary failure mode seems to be due to strong arcs in the raw data that the model fits, but which are discontinuous. The model then spans different bright arcs with different rings, often leading to an offset in x-y, and the wrong slope. 

Set plotFlag = 1 to visualize the fit output. This creates and opens an image for each datapoint, so should be used with caution on large datasets. 

The standard usage is to call this function through Riglar_RINGS_run.m, but you can also call it directly by passing it an image, and additional parameters as name-value pairs. Arguments can be given to Riglar_RINGS_run and will be passed through to Riglar_RINGS_fit. For a list of arguments, see argument parsing in Riglar_RINGS_fit.m.

Default Parameters:

1. ```Riglar_RINGS_fit(im, ‘rMax’, 1.3, ‘expectedSlope’, 0.39, ‘slopeTol’, 0.3, ‘colorPhaseShift’, 1.5)```  
The first parameter (im) is the 2 color image, with shape: height x width x 2.

### Visualization

Example successful fit visualization for (top) YFP and (bottom) CFP.  (Left) normalized image, (middle) model fit and (right) residual.

![alt text](https://github.com/HMS-IDAC/repressilator-colony-rings/blob/master/successful_fit.png "Logo Title Text 1")


Example unsuccessful fit visualization for (top) YFP and (bottom) CFP.  (Left) normalized image, (middle) model fit and (right) residual.

![alt text](https://github.com/HMS-IDAC/repressilator-colony-rings/blob/master/unsuccessful_fit.png "Logo Title Text 1")

### Test Cases

A small example dataset in provided in the folder [sample-data](https://github.com/HMS-IDAC/repressilator-colony-rings/tree/master/sample-data), comprising the first 20 images of each dataset of E. coli MC4100 repressilator bacteria (LPT320), from 3 points (2h, 4h, 6h) of the time course displayed in Figure 1 F-H. 

The script can be run through the MATLAB command window, using the LPT320 default values, as follows:

```Riglar_RINGS_run```

At the first prompt, the data folder (sample-data) is selected.

At the final prompt, the desired save directory for the csv output (theta_0.csv) is selected.

Notes:

1. Alternative parameters can be passed to Riglar_RINGS_fit using the construction  
```Riglar_RINGS_run(‘rMax’, 1.3, ‘expectedSlope’, 0.39, ‘slopeTol’, 0.3, ‘colorPhaseShift’, 1.5)```

2. Fits can be inspected by setting ‘plotFlag’ to 1  
```Riglar_RINGS_run(‘plotFlag’, 1)```

3. Included in the LPT320_t02h timepoint is an example of a failed fit (File 2h_LPT320_2_h_1-636265086740000720-YFP.jpg5). This datapoint will return a NaN value in the csv output file. 

Expected outputs:

1. On screen, visualization of datapoints fitted. 
2. On screen (if plotFlag = 1), figures showing normalized data, model fit and residuals. 1 per datapoint fit. 
3. Saved at prompt, a csv file with all fit theta_0 values arranged in columns by timepoint

Expected runtime on standard desktop:

< 10 seconds (with ‘plotFlag’ = 0)  
< 1 min (with ‘plotFlag’ = 1)

## Authors

```
David T Riglar, David L Richmond, Laurent Potvin-Trottier, Andrew A Verdegaal, Alexander D Naydich, Somenath Bakshi, Emanuele Leoncini, Johan Paulsson, Pamela A Silver
```

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE) file for details
