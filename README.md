# cellshapefft


Development ended here, merged into cjwlab/DSLMDataAnalysis project as directory ``cellshapefft/``

-------

Forked from [mdurade/Coarse-grained-anisotropy-and-size-using-FFT](https://github.com/mdurande/coarse-grained-anisotropy-and-size-using-FFT) , original code 
as developed and described in the article: 
    [Fast determination of cell anisotropy and size in epithelial tissue images using Fourier Transform](https://doi.org/10.1103/physreve.99.062401)

Abstract of the article:

Mechanical strain and stress play a major role in biological processes such as wound healing or
morphogenesis. To assess this role quantitatively, fixed or live images of tissues are acquired at a
cellular precision in large fields of views. To exploit these data, large number of cells have to be
analyzed to extract cell shape anisotropy and cell size. Most frequently, this is performed through
detailed individual cell contour determination, using so-called segmentation computer programs,
complemented if necessary by manual detection and error corrections. However, a coarse grained and
faster technique can be recommended in at least three situations. First, when detailed information
on individual cell contours is not required, for instance in studies which require only coarse-grained
average information on cell anisotropy. Second, as an exploratory step to determine whether full
segmentation can be potentially useful. Third, when segmentation is too difficult, for instance due to
poor image quality or a too large cell number. We developed a user-friendly, Fourier Transform based
image analysis pipeline. It is fast (typically 104 cells per minute with a current laptop computer) and
suitable for time, space or ensemble averages. We validate it on one set of artificial images and on
two sets of fully segmented images, from Drosophila pupa and chicken embryo; the pipeline results
are robust. Perspectives include in vitro tissues, non-biological cellular patterns such as foams, and
xyz stacks.

This fork contains further developments by [@cjwlab](https://github.com/cjwlab):

Original code by Durande (2017)
    
Further modification by Guillermo Serrano Najera (2019)
* Converted into a class for convenience
* Plotting functions modified to be faster with very large images
* Parallelization def_analysis and plotting
* Array multiplication (.*) for masking (much faster than original)

Modifications Jason Klebes (2023)
* Developments specific to our DSLM experiments: option to supress a particular stretched Gaussian noise background, which biases the anisotropy calculation.
* integrating our expReader, should work as installed from github with sister projects (utilities).
* renamed confusing variables
* color-coding outputs by "quality" (avg Fourier space intensity)
* Rearrange for parfor loops

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development.

### Prerequisites

A directory of images of cells over time.
Images cen be .png , .tiff .tif , .jpg .jpeg, .jp2.
Because of expReader, there should be at least five images with numbered names of the same length and pattern, 
for example "img_0001.png", "img_0002.png", ... .

### Installing
* Install matlab 
* Download package ``utilities`` containing class ``expReader`` and make sure it's on the MATLABPATH, or otherwise have an ``expReader`` file on path.
* Make sure ffmpeg is installed and on path (comes up with system command ``ffmpeg``), if not install a copy of ffmpeg.  If locally installed and not on path, the full ffmpeg_path can be input manually.
* Download this repository.  This repository contains script ``call_cellshape_fft.m``, class ``cellshapefft``, and function collections ``spectrum_analysis.m``, ``deformation_ellipse.m``, ``deformation_matrix.m``, ``visualization_ellipse.m``, ``visualization_strain.m``.  Function ``periodicDecomposition.m`` is contained in the utilities package, another copy is shipped here.

### How to use it in practice
To run the code:
* The file to execute is ``call_cellshape_fft.m``.  Inside this file, we create the ``param`` struct holding all user input paramteres and an object ``cellshapefft``. For fields set to empty, default values, which can be seen in ``cellshapefft.m``, will be used. We then initialize a ``cellshapefft`` object and run the analysis by calling method "full_analysis".
        
```
    Example:
    
    param = struct();
    
    param.pathin = '\\lfs.lifesci.dundee.ac.uk\lfs\cjw\Jason\cable detection\exp0186\focused_one';
    param.pathout = ['\\lfs.lifesci.dundee.ac.uk\lfs\cjw\Jason\cable detection\exp0186\results_actin_', ...
    char(datetime('now', Format='dd-MM-yyyy''_''HH-mm-ss'))];

    param.ffmpeg_path = [];
    param.timePoints = [1:48];
    param.time_avg = 1;
    param.method = 'matrix'; %choose between ellipse-fitting and inertia matrix
    param.workers = 16;
    param.chunk_size = 96; %ideally multiple of number of workers
    param.tileSize = 100; %size of tiles.
    param.overlap= 0.5;
    param.cut = 5; %masking a circle in the middle of Fourier spectrum
    param.propor = 0.02; 
    param.sigma = 0.8;
    param.strel = 4;
    
    %scale, color of lines to draw in visualization:
    param.scale = [];
    param.col='green'; 

    param.writeSpectra = []; %to choose whether spectrum data is saved
    param.stretchedNoise = []; %correction particular to interpolated DSLM data
    
    obj = cellshapefft(param);
    obj.full_analysis();

```

## Outputs
In the given output director:
* TIF stacks of images 
* A copy of the code
* data files such as ``cellshapefft_results_SangS_<t1>_<t2>.mat``, containing the objects ``M``, ``S`` and ``angS`` ("matrix" route) or ``a``,``b``, and ``phi`` (ellipse route) in arrays of shape [N , nTiles, chunk_size] for the ``chunk_size`` timepoints in a chunk.  Where N = [2,2] for inertia matrix, 2 for fields (S, angS) and 3 for fields (a,b,phi)
    * ``M`` inertia matrix.
    * ``S``, ``angS`` cell strain deviator amplitude ( 1/2*(log(L1/L2)) in paper) , and angle.
    * ``a``, ``b``, ``phi`` : (real space) ellipse major axis, minor axis, and angle.
