refsig 1.0.4
---------------------

The refsig package works as tool for identifying the referential signal from a set of unipolar iEEG data.

Package contains 3 methods for obtaining the reference:

(1) The first method creates the referential signal for the given iEEG data set using Method I from [1].
It does so using a comparative method based on correlation coefficient among the ICs of unipolar and bipolar montage.
The method does not calculate the referential signal, however it merely locates one unipolar independent component,
which has the lowest correlation coefficient with any bipolar IC.

(2) Second method creates the referential signal of an iEEG set of data
using method II described in [1].
The difference from method I is that this method is solely based on
calculating the reference instead of locating it amongst unipolar ICs.
This method has proven to be the most accurate and fastest by the experimental results in [1]
as well as in our own testing. Therefore this method should be used for obtaining the reference.
    
(3) The last method is a calculation of mere average reference from the given set of iEEG data. 
This method is the fastest one, however not as accurate as the methods based on ICA. 

For import, use:     

**-> from refsig import getref**

Each method is then initialized by:  

**(1) -> getref.find(unipol**, N_iterations = 20, p = 0.2)    

**(2) -> getref.calc(unipol)** <--this one should be used!   

**(3) -> getref.avg(unipol)**  

	where unipol is a m*k set of iEEG data, where m is the number of channels
	and k is the number of samples. 

All methods return the referential signal as a linear vector 1xk.

Based on:  

[1] HU, S., STEAD, M., WORRELL, G. A. Automatic identification and removal of scalp reference signal for intracranial EEGs based on independent component analysis. IEEE Transactions on Biomedical Engineering, 2007.
