# Motion-triggered average
"Motion-triggered average (MTA)" is an algorithm to characterize pseudo-simultaneous dynamic changes in arbitrary cellular deformation and molecular activities. Using MTA, one can extract a pseudo-simultaneous time series from individually observed activities of intracellular molecules. This software was used in the analysis of the following paper:
Kunida et al., Cell Reports, 2023, [https://doi.org/10.1016/j.celrep.2023.112071](https://doi.org/10.1016/j.celrep.2023.112071)
 
# DEMO
 
You can learn how to use MTA by just run:

* Step1_SingleCellMTA.m
* Step2_CellCommonMTA.m

These codes perform single-cell MTA and cell-common MTA using sample data. At the top of the codes, specify the data files and parameters as desired.

# Features

You can run MTA using only Matlab itself without using a particular toolbox.

# Requirement

* Matlab R2022b

MTA probably works with older versions of Matlab, but we have not checked MTA with older one.

# Installation

Just download and run it on Matlab.

# Usage

```
>> cd 'MTA folder'
>> Step1_SingleCellMTA
>> Step2_CellCommonMTA
```

# Note
 
Ensure that the folder you put your data in is consistent with the settings in the code. 
 
# Author

* Yuichi Sakumura, Nara Institute of Science and Technology, Japan (saku@bs.naist.jp)
* Kazuyuki Kunida, Fujita Health University, Japan (katsuyuki.kunida@fujita-hu.ac.jp)


# License

"MTA" is under [GNU GPL](https://en.wikipedia.org/wiki/GNU_General_Public_License). 

