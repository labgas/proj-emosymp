*SCRIPT TO ANALYZE HEAD MOVEMENT PARAMETERS OF PROJ-EMOSYMP
AUTHORS: LUKAS VAN OUDENHOVE
@DARTMOUTH, APRIL 2021

-----------
DATA IMPORT
-----------;

%let path=C:\Users\lukas\Dropbox (Dartmouth College)\fMRI_emotion_Giao\BIDS\qc;
%let outputdir=C:\Users\lukas\Dropbox (Dartmouth College)\fMRI_emotion_Giao\BIDS\qc;
proc import
datafile="&path\group_bold.tsv"
out=work.emosymp_hmp dbms=tab;
run;


*------------------
CHECK DISTRIBUTIONS
-------------------;
* entire sample included in fmri analyses;
proc univariate data=emosymp_hmp;
where subj_included = 1 AND run_included = 1;
var dvars_std fd_mean;
histogram dvars_std fd_mean / normal (mu=est sigma=est) lognormal (theta=est zeta=est sigma=est);
run;

* per group;
proc sort data=emosymp_hmp;
by patient;
run;

proc univariate data=emosymp_hmp;
by patient;
where subj_included = 1 AND run_included = 1;
var dvars_std fd_mean;
histogram dvars_std fd_mean / normal (mu=est sigma=est) lognormal (theta=est zeta=est sigma=est);
run;


*------------------
TRANSFORM VARIABLES
-------------------;
data emosymp_hmp;
set emosymp_hmp;
z=0;
run;
*adds variable z with all zeros, needed in proc transreg;

proc transreg data=emosymp_hmp maxiter=0 nozeroconstant;
where subj_included = 1 AND run_included = 1;
model BoxCox(dvars_std/parameter=0) = identity(z);
run;
*check lambda in output, in this case -1.25
parameter is constant to make all values positive if there are negative values, hence parameter = |minimum|, see below;

proc transreg data=emosymp_hmp maxiter=0 nozeroconstant;
where subj_included = 1 AND run_included = 1;
model BoxCox(fd_mean/parameter=0) = identity(z);
run;
* lambda = 0.25;

data emosymp_hmp;
set emosymp_hmp;
where subj_included = 1 AND run_included = 1;
bc_dvars_std = (dvars_std**-1.25 -1)/-1.25;
run;
*boxcox formula, -2.5 is lambda, -1 is constant, parameter = 0;

data emosymp_hmp;
set emosymp_hmp;
where subj_included = 1 AND run_included = 1;
bc_fd_mean = (fd_mean**0.25 -1)/0.25;
run;

*check normality of box-cox transformed variables;
proc univariate data=emosymp_hmp;
where subj_included = 1 AND run_included = 1;
var bc_dvars_std bc_fd_mean;
histogram bc_dvars_std bc_fd_mean / normal (mu=est sigma=est);
run;
*distributions normal;


*-----------
MIXED MODELS
------------;
proc mixed data=emosymp_hmp;
where subj_included = 1 AND run_included = 1;
class subject run patient;
model bc_dvars_std = patient / solution residual influence;
repeated run / subject=subject type=cs r rcorr;
lsmeans patient;
run;

proc mixed data=emosymp_hmp;
where subj_included = 1 AND run_included = 1;
class subject run patient;
model bc_fd_mean = patient / solution residual influence;
repeated run / subject=subject type=cs r rcorr;
lsmeans patient;
run;

proc mixed data=emosymp_hmp;
where subj_included = 1 AND run_included = 1;
class subject run patient;
model bc_dvars_std = patient | run / solution residual influence;
repeated run / subject=subject type=cs r rcorr;
lsmeans patient;
lsmeans run / diff=all;
run;

proc mixed data=emosymp_hmp;
where subj_included = 1 AND run_included = 1;
class subject run patient;
model bc_fd_mean = patient | run / solution residual influence;
repeated run / subject=subject type=cs r rcorr;
lsmeans patient;
lsmeans run / diff=all;
run;
