########################################################################
Title:   The image intra-class correlation coefficient (I2C2) for replication studies
Authors: H. Shou, A. Eloyan, S. Lee, V. Zipunnikov, A.N. Crainiceanu, Mary Beth Nebel
         B. Caffo, M.A. Lindquist and C.M. Crainiceanu 
Contact: Department of Biostatistics
         Bloomberg School of Public Health 
         Johns Hopkins University
         615 N. Wolfe St., Baltimore, MD 21205
Email:   hshou@jhsph.edu

#######################################################################################
Code Files and Descriptions

"trace_only_unba.R":  The main function that calculate I2C2 based on the trace 
                      method described in the paper. The code works for balanced 
					  and unbalanced design. 

"I2C2_inference.R":   The functions that provide inference for I2C2. I2C2.mcCI builds 
                      confidence interval of the estimated I2C2 by bootstrapping subjects;
                      I2C2.mcNulldist gives the null distribution of I2C2 by permuting the 
                      imaging observations.	If the algorithm runs on a machine with multiple
                      cores or CPUs, one can set 'ncores' greater than 1 and utilize parallel
                      computing by 'multicore' package.					  
 
"demo_vn10.R":        An example code file illustrating how I2C2 can be estimated and made 
                      inference on. A beanplot to visualize the I2C2 and its null distribution
					  is provided. The example uses Kirby21 ventricular RAVENS images which
					  are included in the zip file. 
                       
"vn10.RData":         Data file containing Kirby21 ventricular RAVENS images were thresholded 
                      at sum > 15.
"mask10.RData":       Masks indicating ventriculars in "vn10.RData"				  
##############################################################################                      
                       
