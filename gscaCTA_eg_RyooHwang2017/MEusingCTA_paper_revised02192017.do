/***** MEusingCTA_paper revised by considering correlation and standardized variables *****/

. clear

cd "D:\Documents\UVa_Tenure\Projects\ProjectA6_GSCA\MEusingCTA\Analysis"

use ECLS_CTA_SocialSkill.dta

rename X4PRNCON p4con
rename X4PRNSOC p4int
rename X4PRNSAD p4sad
rename X4PRNIMP p4imp

rename X4TCHCON t4con
rename X4TCHPER t4int
rename X4TCHEXT t4extb
rename X4TCHINT t4intb

/*** New approach with correlation ***/
egen p4con1 = std(p4con)
egen p4int1 = std(p4int)
egen p4sad1 = std(p4sad)
egen p4imp1 = std(p4imp)
egen t4con1 = std(t4con)
egen t4int1 = std(t4int)
egen t4extb1 = std(t4extb)
egen t4intb1 = std(t4intb)
/*** End of standardizing variables ***/

/* It is required to replace the model-implied covariance matrix with the one obtained from GSCA */

sem (SS_P -> p4con1 p4int1 p4sad1 p4imp1) (SS_T -> t4con1 t4int1 t4extb1 t4intb1)
estat framework, fitted
mat sigma0 = r(Sigma)

/* Model-implied covariance matrix - Model 0 */

matrix mivc0 = ///
(1.62244,	0.37713,	0.56165,	0.55989,	0.17444,	0.16919,	0.16060,	0.10059\ ///
0.37713,	1.22850,	0.34030,	0.33924,	0.10569,	0.10251,	0.09731,	0.06095\ ///
0.56165,	0.34030,	1.50680,	0.50521,	0.15740,	0.15267,	0.14492,	0.09077\ ///
0.55989,	0.33924,	0.50521,	1.50363,	0.15691,	0.15219,	0.14447,	0.09049\ ///
0.17444,	0.10569,	0.15740,	0.15691,	1.83781,	0.81261,	0.77137,	0.48314\ ///
0.16919,	0.10251,	0.15267,	0.15219,	0.81261,	1.78816,	0.74817,	0.46861\ ///
0.16060,	0.09731,	0.14492,	0.14447,	0.77137,	0.74817,	1.71021,	0.44483\ ///
0.10059,	0.06095,	0.09077,	0.09049,	0.48314,	0.46861,	0.44483,	1.27862)

/* Repalce sigma1 with mivc0 */
forvalues i = 1/8 {
	forvalues j = 1/8 {
		 matrix sigma0[`i',`j']= mivc0[`i',`j']
	}
}
      
/* Confirmatory Tetrad Analysis */
net install tetrad, from (https://github.com/sbauldry/tetrad/raw/master)

tetrad p4con1 p4int1 p4sad1 p4imp1 t4con1 t4int1 t4extb1 t4intb1, icm1(sigma0) 



/* It is required to replace the model-implied covariance matrix with the one obtained from GSCA */

qui sem (SS_P -> p4con1 p4int1 p4sad1 p4imp1) (SS_T -> t4con1 t4int1 t4extb1 t4intb1) (SS_P -> SS_T)
qui estat framework, fitted
mat sigma1 = r(Sigma)

/* Model-implied covariance matrix - Model 1 */

matrix mivc1 = ///
(1.62635,	0.38628,	0.55270,	0.57485,	0.26539,	0.25827,	0.24575,	0.15570\ ///
0.38628,	1.23823,	0.34086,	0.35452,	0.16367,	0.15928,	0.15156,	0.09602\ ///
0.55270,	0.34086,	1.48771,	0.50726,	0.23418,	0.22790,	0.21686,	0.13739\ ///
0.57485,	0.35452,	0.50726,	1.52758,	0.24357,	0.23703,	0.22555,	0.14290\ ///
0.26539,	0.16367,	0.23418,	0.24357,	1.97788,	0.81318,	0.77378,	0.49024\ ///
0.25827,	0.15928,	0.22790,	0.23703,	0.81318,	1.97193,	0.75302,	0.47709\ ///
0.24575,	0.15156,	0.21686,	0.22555,	0.77378,	0.75302,	1.96186,	0.45397\ ///
0.15570,	0.09602,	0.13739,	0.14290,	0.49024,	0.47709,	0.45397,	1.90414)

/* Repalce sigma1 with mivc0 */
forvalues i = 1/8 {
	forvalues j = 1/8 {
		 matrix sigma1[`i',`j']= mivc1[`i',`j']
	}
}
      
/* Confirmatory Tetrad Analysis */
tetrad p4con1 p4int1 p4sad1 p4imp1 t4con1 t4int1 t4extb1 t4intb1, icm1(sigma1) 


/* It is required to replace the model-implied covariance matrix with the one obtained from GSCA */

qui sem (SS_P -> p4con1 p4int1 p4sad1 p4imp1) (SS_T -> t4con1 t4int1 t4extb1 t4intb1) (SS -> SS_P@1 SS_T@1)
qui estat framework, fitted
mat sigma2 = r(Sigma)

/* Model-implied covariance matrix - Model 2 */

matrix mivc2 = ///
(1.62287,	0.38019,	0.55809,	0.56377,	0.46559,	0.45202,	0.42925,	0.26965\ ///
0.38019,	1.23207,	0.34065,	0.34412,	0.28419,	0.27591,	0.26201,	0.16459\ ///
0.55809,	0.34065,	1.50004,	0.50514,	0.41716,	0.40501,	0.38461,	0.24161\ ///
0.56377,	0.34412,	0.50514,	1.51028,	0.42141,	0.40914,	0.38852,	0.24407\ ///
0.46559,	0.28419,	0.41716,	0.42141,	1.83697,	0.81259,	0.77165,	0.48474\ ///
0.45202,	0.27591,	0.40501,	0.40914,	0.81259,	1.78891,	0.74917,	0.47062\ ///
0.42925,	0.26201,	0.38461,	0.38852,	0.77165,	0.74917,	1.71143,	0.44691\ ///
0.26965,	0.16459,	0.24161,	0.24407,	0.48474,	0.47062,	0.44691,	1.28075)

/* Repalce sigma1 with mivc0 */
forvalues i = 1/8 {
	forvalues j = 1/8 {
		 matrix sigma2[`i',`j']= mivc2[`i',`j']
	}
}
      
/* Confirmatory Tetrad Analysis */
tetrad p4con1 p4int1 p4sad1 p4imp1 t4con1 t4int1 t4extb1 t4intb1, icm1(sigma2) 


/* It is required to replace the model-implied covariance matrix with the one obtained from GSCA */

sem (SS_P1 -> p4con1 p4int1) (SS_P2 -> p4sad1 p4imp1) (SS_T1 -> t4con1 t4int1) (SS_T2 -> t4extb1 t4intb1)
estat framework, fitted
mat sigma3 = r(Sigma)

/* Model-implied covariance matrix - Model 3 */

matrix mivc3 = ///
(1.61034,	0.61275,	0.28146,	0.28183,	0.15020,	0.15020,	0.10178,	0.10142\ ///
0.61275,	1.61518,	0.28257,	0.28294,	0.15079,	0.15079,	0.10219,	0.10182\ ///
0.28146,	0.28257,	1.65031,	0.65116,	0.13739,	0.13739,	0.13214,	0.13167\ ///
0.28183,	0.28294,	0.65116,	1.65202,	0.13757,	0.13757,	0.13231,	0.13184\ ///
0.15020,	0.15079,	0.13739,	0.13757,	1.89913,	0.89913,	0.49485,	0.49309\ ///
0.15020,	0.15079,	0.13739,	0.13757,	0.89913,	1.89913,	0.49485,	0.49309\ ///
0.10178,	0.10219,	0.13214,	0.13231,	0.49485,	0.49485,	1.65282,	0.65050\ ///
0.10142,	0.10182,	0.13167,	0.13184,	0.49309,	0.49309,	0.65050,	1.64818)

/* Repalce sigma1 with mivc0 */
forvalues i = 1/8 {
	forvalues j = 1/8 {
		 matrix sigma3[`i',`j']= mivc3[`i',`j']
	}
}
      
/* Confirmatory Tetrad Analysis */
tetrad p4con1 p4int1 p4sad1 p4imp1 t4con1 t4int1 t4extb1 t4intb1, icm1(sigma3) 


/* CTA with Model 0 va. Model 3 */

tetrad p4con1 p4int1 p4sad1 p4imp1 t4con1 t4int1 t4extb1 t4intb1, icm1(sigma0) icm2(sigma3) reps(5) seed(1234) tlist(1)

/* CTA with Model 1 va. Model 3 */

tetrad p4con1 p4int1 p4sad1 p4imp1 t4con1 t4int1 t4extb1 t4intb1, icm1(sigma1) icm2(sigma3) reps(5) seed(1234) tlist(1)

/* CTA with Model 2 va. Model 3 */

tetrad p4con1 p4int1 p4sad1 p4imp1 t4con1 t4int1 t4extb1 t4intb1, icm1(sigma2) icm2(sigma3) reps(5) seed(1234) tlist(1)

/* Confirmatory Tetrad Analysis for Model 0*/
tetrad p4con1 p4int1 p4sad1 p4imp1 t4con1 t4int1 t4extb1 t4intb1, icm1(sigma0) seed(1) tlist(1)

/* Confirmatory Tetrad Analysis for Model 1*/
tetrad p4con1 p4int1 p4sad1 p4imp1 t4con1 t4int1 t4extb1 t4intb1, icm1(sigma1) seed(1) tlist(1)

/* Confirmatory Tetrad Analysis for Model 2*/
tetrad p4con1 p4int1 p4sad1 p4imp1 t4con1 t4int1 t4extb1 t4intb1, icm1(sigma2) seed(1) tlist(1)

/* Confirmatory Tetrad Analysis for Model 3*/
tetrad p4con1 p4int1 p4sad1 p4imp1 t4con1 t4int1 t4extb1 t4intb1, icm1(sigma3) seed(1) tlist(1)


/******* Additional example *******/
/* CTA with Model 1 va. Model 2 */
tetrad p4con1 p4int1 p4sad1 p4imp1 t4con1 t4int1 t4extb1 t4intb1, icm1(sigma0) icm2(sigma2) reps(5) seed(1234) tlist(1)




