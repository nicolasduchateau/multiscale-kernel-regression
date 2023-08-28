Code for the Multiscale kernel regression, available under the license [CeCILL-B](http://www.cecill.info/licences/Licence_CeCILL-B_V1-en.html)
  
# Multiscale kernel regression - Matlab + Python

Release 1.0 = MATLAB / Author: Nicolas Duchateau (CREATIS Lyon, France) / February 2018

Release 2.0 = PYTHON / Authors: Benoit Freiche, Fei Zheng, Nicolas Duchateau (CREATIS Lyon, France) / August 2023

Links to the corresponding publications at: <br/> https://www.creatis.insa-lyon.fr/~duchateau/#publications

------------------------------------------------------------------------------------------------------------------------
**NOTICE:**

This code is made open-access. Comments and bug reports are welcome, as well as feedback on its possible improvements.

**Published reports of research using this code (or a modified version) may cite the following articles that describes the method:**

**- *exact matching* case:**

*Bermanis A, Averbuch A, Coifman RR. Multiscale data sampling and function extension. Applied and Computational Harmonic Analysis, 2013;34(1):15-29.*
https://doi.org/10.1016/j.acha.2012.03.002

**- *inexact matching* case, corresponding to the present MATLAB implementation:**

*Duchateau N, De Craene M, Sitges M, Caselles V. Adaptation of multiscale function extension to inexact matching: Application to the mapping of individuals to a learnt manifold. In: Proceedings of SEE International Conference on Geometric Science of Information (GSI). Springer LNCS, 2013;8085:578-86.*
https://doi.org/10.1007/978-3-642-40020-9_64

------------------------------------------------------------------------------------------------------------------------
**ARCHIVE CONTENT [MATLAB]**

**testSIN.m** = example to launch the code on a sinusoid at different frequencies (to test the regression at different scales)

**interpolateMulti.m** = launching the multi-scale routine (algorithm #4 from adapted Bermanis et al. 2013)

**MINEXACT_Bermanis_ACHA_2013_ALGO_03.m** = computations at a single scale (algorithm #3 adapted from Bermanis et al. 2013)

------------------------------------------------------------------------------------------------------------------------
**ARCHIVE CONTENT [PYTHON]**

**regression.py** = class MultiScaleKernelRidge + example to launch the code on a sinusoid at different frequencies (to test the regression at different scales)
