# Dynamic Endocast By Bone -- Details and Instructions

This document describes the workflow and functions in detail, and provides instructions for its use.

## dynamicEndocastByBone function details

Inputs: 
 datapath
 XYZfile
 RBTfile
 refbone
 freezeIncrement
 savepath
 savefile
 objFolder

1. Install the latest version of XMALab found [here](https://bitbucket.org/xromm/xmalab/).
2. Follow [DLC documentation](https://github.com/AlexEMG/DeepLabCut/blob/master/docs/installation.md) to install Anaconda environment. We recommend using the provided [easy install](https://github.com/AlexEMG/DeepLabCut/blob/master/conda-environments/README.md).
3. If you haven't already, clone or download this repository, and locate the XROMM_DLCTools [functions](/functions/xrommtools.py) and [Demo Jupyter Notebook](/templates/XROMM_Pipeline_Demo.ipynb)
4. **Copy** the xrommtools.py file into the DLC Anaconda environment's utils folder, which can be found where Anaconda was installed- here: ...\Anaconda3\envs\dlc-windowsGPU\Lib\site-packages\deeplabcut\utils

      Or on a Mac:

      .../opt/anaconda3/envs/dlc-macOS-CPU/lib/python3.6/site-packages/deeplabcut/utils

![](https://user-images.githubusercontent.com/53494838/74692595-9ccb9080-51ad-11ea-9906-e6b841238ad7.png)

5.  **Importantly,** Follow DLC documentation to verify successful installation before proceeding.
