SNPInt-GPU
==========

SNPInt-GPU is a software providing several methods for statistical epistasis testing. SNPInt-GPU supports GPU acceleration, but can also be used without GPU hardware. The software implements logistic regression (as in PLINK epistasis testing), BOOST, log-linear regression, mutual information (MI) and information gain (IG) for pairwise testing as well as mutual information and information gain for third-order tests. Optionally, r^2 scores for testing for linkage disequilibrium (LD) can be calculated on-the-fly.

Installation instructions:
--------------------------

SNPInt-GPU requires a Linux-based operating system (e.g.\ Ubuntu) with an up-to-date C++ compiler and CUDA libraries. Furthermore, Zlib, OpenMP, BOOST and TBB libraries as well as the Cmake are required.
We tested our software on an Ubuntu 19.04 system with GCC-8, CUDA 10.1, BOOST 1.65.1 and TBB 4.4. Please ensure that you have installed the required software before continuing with the installation. Please also note, that the CUDA libraries are still required, even if you don't want to use GPU acceleration.

The step by step compilation procedure is as follows. 

1.  Either clone or download the software sources from GitHub. We recommend cloning the repository using Git:

        $ git clone https://github.com/ikmb/snpint-gpu.git  

2.  Switch into the source directory and generate the Makefile in a subfolder named "build" (will be created if necessary):

        $ cd snpint-gpu
        $ cmake -DCMAKE_BUILD_TYPE=Release -B"build"
  
3.  As a last step, switch into the build folder and type "make" or "make -j" on a multi-core system to finally build the binary. It will be located in that build folder and can directly be executed.

        $ cd build
        $ make -j

