refsig 1.0.3
---------------------

This package works as tool for identifying the referential signal from a set of unipolar iEEG data.

Package contains 3 methods for obtaining the reference, but second method works the best and should be used in practice (see the syntax below). Methods 1 works really slow and the results are not satisfying. Method Avg is the basic algorithm.


For import, use:     

**-> from refsig import ref**

Each method is initialized by:  

**-> ref.m1(unipol**, N_iterations = 20, p = 0.2)    

**-> ref.m2(unipol)** <--this one should be used!   

**-> ref.avg(unipol)**  

	where unipol is a m*k set of iEEG data, where m is the number of channels
	and k is the number of samples. 

Based on:  

HU, S., STEAD, M., WORRELL, G. A. Automatic identification and removal of scalp reference signal for intracranial EEGs based on independent component analysis. IEEE Transactions on Biomedical Engineering, 2007.
