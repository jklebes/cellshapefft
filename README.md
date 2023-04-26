# cellshapefft

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

Modifications Jason Klebes (2022)
* Developments specific to our DSLM experiments: attempts to suppress scan stripe noise, which biases the anisotropy calculation.
* splitting off some functions (perdecomp) for re-use.
* Should work as installed from github with sister projects (+utilities).
* notebooks, new GUI app for calibrating

## Getting Started

These instructions will get you a copy of the project up and running on your local machine for development.

### Prerequisites

You will need a version of matlab older than 2015. (Because MATLAB's GUI designer GUIDE was ended in 2015.  With newer MATLAB the app still runs as usual, but the GUI is not edittable in the same way - jk)
Images should be either .png or .tiff .tif (stacks are okay) 

### Installing
* Install matlab 
* Download package ``+utilities`` containing class ``expReader`` and make sure it's on the MATLABPATH, or otherwise have an ``expReader`` file on path.
* Make sure ffmpeg is installed and on path (comes up with terminal command ``ffmpeg``), if not install a copy of ffmpeg.  If locally installed and not on path, the full ffmpeg_path can be input manually.
* Download this repository.  This repository contains script ``call_cellshape_fft.m``, class ``cellshapefft``, and function collections ``spectrum_analysis.m``, ``deformation_ellipse.m``, ``deformation_matrix.m``, ``visualization_ellipse.m``, ``visualization_strain.m``.  Function ``perdecomp.m`` is contained in the +utilities package, another copy is shipped here.

### How to use it in practice
To run the code:
* The file to execute is ``call_cellshape_fft.m``.  Inside this file, we create the ``param`` struct holding all user input paramteres and an object ``cellshapefft``. We then run the analysis by calling method "full_analysis" or "full_analysis_chunks".  Or the user can otherwise build a ``param`` struct, create a  ``cellshapefft``  object, and call the ``cellshapefft.full_analysis`` method.
        
```
    Example:
    
    param = struct();
    param.pathin = 'C:PATH\TO\IMAGES';
    param.contour = 'C:PATH\TO\MASK';
    param.pathout = 'C:PATH\TO\RESULTS';
    param.siz = [];
    param.time_points = [];
    param.chunk_size = [];
    param.tleng = [];
    param.timestep = [];
    param.tileSize = [];
    param.fres = [];
    param.cut = [];
    param.propor = [];
    param.sigma = [];
    param.scale = [];
    param.strel = [];
    param.register = [];
    param.regsize = [];
    param.workers = [];

    obj = cellshapefft(param);
    obj.full_analysis_chunks;

```

## Outputs

* a ``results_TIMESTAMP`` director is created in the output path
* TIF stacks of images ``Deformation_map...`` and ``Deformation_map_onim...``
* A snapshot of the GUI if used
* ``Param.mat``, saving input parameters
* ``Results.mat``, saving results

``Results.mat`` contains a structure ``Results``, values ``nTiles`` (number of subimages), ``tileCoords`` (positions of subimages), ``numX`` and ``numY`` (dimensions of subimages array), ``ci`` (time average window?).  Main results are in ``Results.im_regav`` structure in the form of arrays ``M``, ``S``, ``angS``, ``a``, ``b``, ``phi``.
* ``a``, ``b``, ``phi`` : (real space) ellipse major axis, minor axis, and angle.
