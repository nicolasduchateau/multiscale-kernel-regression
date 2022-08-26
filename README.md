# multiscale-kernel-regression
multiscale kernel regression - Matlab

**Instructions for MATLAB code for multi-scale kernel regression.**

Release 1.0 = MATLAB scripts

Author: Nicolas Duchateau (nicolas.duchateau@creatis.insa-lyon.fr)

Date: February, 2018

Website for updates: links to the corresponding publications at: <br/> https://www.creatis.insa-lyon.fr/~duchateau/#publications

------------------------------------------------------------------------------------------------------------------------
**COPYRIGHT NOTICE:**

Multiscale kernel regression
(c) 2018 Nicolas Duchateau

This code is made open-access. Comments and bug reports are welcome, as well as feedback on its possible improvements.

Published reports of research using this code (or a modified version) may cite the following article that describes the method in the **exact matching** case: <br/> *Bermanis A, Averbuch A, Coifman RR. Multiscale data sampling and function extension. Applied and Computational Harmonic Analysis, 2013;34(1):15-29.*

The present MATLAB implementation is the one detailed in: <br/> *Duchateau N, De Craene M, Sitges M, Caselles V. Adaptation of multiscale function extension to inexact matching: Application to the mapping of individuals to a learnt manifold. In: Proceedings of SEE International Conference on Geometric Science of Information (GSI). Springer LNCS, 2013;8085:578-86.*

A complementary list of related publications / applications of the algorithm can be found at: <br/> https://www.creatis.insa-lyon.fr/~duchateau/#publications

You are free to modify, extend or distribute this code, as long as this copyright notice is included whole and unchanged.

------------------------------------------------------------------------------------------------------------------------
**ARCHIVE CONTENT:**

**testSIN.m** = example to launch the code on a sinusoid at different frequencies (to test the regression at different scales)

**interpolateMulti.m** = launching the multi-scale routine (algorithm #4 from adapted Bermanis et al. 2013)

**MINEXACT_Bermanis_ACHA_2013_ALGO_03.m** = computations at a single scale (algorithm #3 adapted from Bermanis et al. 2013)
