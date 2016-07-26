# MMICA
Maximal Mutual Information Component Analysis

The details on the algorithm and application can be found in __Detecting Microbial Dysbiosis Associated with Pediatric Crohn Disease Despite the High Variability of the Gut Microbiota__ [](http://dx.doi.org/10.1016/j.celrep.2015.12.088) published on ___Cell Reports___ in 2016.


# Usage
Users need to write an IO to input the micrbiome data into the code and analysis the output results. The input data should be a `numpy` array with size `[N_OTUs, N_samples]`, and the first `Nc` samples are control ones. The code for MMICA was written as functions which can be easily understood. To do it, one can implement the function `readCell()` to convert the microbiome data in `biom`-format or other formats to a `numpy` array.

`misfunc` includes a function `showarray()`, which prints array in a better look. It is not a part of MMICA, and one can implement `showarray()` to be simply as
```python
def showarray(x):
	print x
```

# Requirement
* python >= 2.7
* scipy

# Citation
Please consider to cite
```bibtex
@article{wang_detecting_2016,
	title = {Detecting {Microbial} {Dysbiosis} {Associated} with {Pediatric} {Crohn} {Disease} {Despite} the {High} {Variability} of the {Gut} {Microbiota}},
	volume = {14},
	url = {http://www.sciencedirect.com/science/article/pii/S2211124715015442},
	doi = {10.1016/j.celrep.2015.12.088},
	number = {4},
	journal = {Cell Reports},
	author = {Wang, Feng and Kaplan, Jess L. and Gold, Benjamin D. and Bhasin, Manoj K. and Ward, Naomi L. and Kellermayer, Richard and Kirschner, Barbara S. and Heyman, Melvin B. and Dowd, Scot E. and Cox, Stephen B. and Dogan, Haluk and Steven, Blaire and Ferry, George D. and Cohen, Stanley A. and Baldassano, Robert N. and Moran, Christopher J. and Garnett, Elizabeth A. and Drake, Lauren and Otu, Hasan H. and Mirny, Leonid A. and Libermann, Towia A. and Winter, Harland S. and Korolev, Kirill S.},
	month = feb,
	year = {2016},
	pages = {945--955}
}
```
