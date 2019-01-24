In order to fully reproduce the paper results, R (v3.5.0), Stan (v2.17.3) are required. 

####################################################################################################
# Codes used for the article entitled "A bifactor generalized partial credit model with flexible   #
# link functions: a novel approach in survey data", by Marcelo A. da Silva,                        #
# Anne C. Huggins-Manley, José A. Mazzon, and Jorge L. Bazán                                       #                                                  
#                                                                                                  #
# Real data is also provided and can be reproduced.                                                #
# Date: 2018-11-08                                                                                 #
####################################################################################################


Folder "Application_JAS" contains the files needed to reproduce the results presented in Section 5 
of the article. Inside this folder we provide commented R codes, as well as Stan model files, which 
are called when the file R is run. This file must be opened in software R and run by the user. Each 
model considered in the manuscript (bifac-GPC model with logit link function, bifac-GPC model with 
probit link function, and bifac-GPC model with clog-log link function) has two Stan files: one to 
find the parameter estimates and the values of the model comparison criteria and the other to 
simulate replicas of the data matrix y for the purpose of calculating the Bayesian p-value using 
the PPMC method, as described in the manuscript.

To correctly reproduce the Application with mobile banking data presented in the article, the 
following steps must be performed:
1. Install packages "loo" in software R.
2. Install packages "rstan" in software R (for information about how to install RStan package, 
   access the website: http://mc-stan.org).
2. Open the file "cloglog.R" (or the logit.R or probit.R files) in software R.
3. Properly configure the file path in line 15 of code.
4. Select all the code and run it.

After running all of the R code (which may take a few hours), the following CSV files will be 
generated in the configured directory (item 3):
1. A CSV file containing the model parameter estimates.
2. A CSV file containing the values of the model comparison criteria (DIC, WAIC and LOO).
3. A CSV file containing the Bayesian p-value.

####################################################################################################