
%let progname = tlc31.sas;
* input: chol.dat
* output:
* xref:
* does: descriptive statistics

*******************************************************************;
title1 "&progname Cholesterol";
filename chol "~/BIOS667/data/cholesterol.dat";
*******************************************************************;

data A;
  infile chol firstobs=2; *starts on 2nd row;
  input group id y1 y2 y3 y4 y5; *variable headers;
run;

proc sort data=A;
	by group;
run;

title 'Covariance matrix, summary stats per treatment group (original vars)';
proc corr data = A  cov ; *prints summary statistics, covariance matrix ;
  var y1-y5; *correlation coefficients for these vars;
  by group; *per trmt group;
  run;
title;

*****************************;
* transpose for boxplot means;
proc sort data = A out=AD; 
	by id; 
run;

proc transpose data=AD out=D(rename=(col1=level)) name=year;
	var y1-y5;
	by id group;
run;

proc sort data=D;
	by group;
run;

* boxplot graph means;
title 'Average cholesterol levels over time by group';
proc sgplot data=D;
	vbox level / category=year group=group connect=mean legendlabel='Group';
	xaxis 
		values = ("y1" "y2" "y3" "y4" "y5")
		valuesdisplay = ("Baseline" "6" "12" "20" "24"); *re-label x-axis tick marks;
	yaxis grid
		values = (120 to 450 by 30);
	label year='Month';
	label level='Serum cholesterol (mg/dL)';
run;
title;

*****************************************;
* transpose for series graph means;
title2 "Summary: mean, sd, #obs";
proc summary data = A  nway ;
  class group; *separate by treatment group -> 2 rows; 
  var  y1-y5; * choose 4 vars; 
  output out  = B 
    mean = y1-y5
; * output data set with summary stats per timepoint
* note n balanced, equal number of samples per timepoint ;
  run;

* transpose for series graph means;
proc transpose data=B out=C(rename=(col1=cholmean)) name=Month;
	var y1-y5;
	by group;
run;

* series graph means;
title 'Average cholesterol levels over time by group';
proc sgplot data=C;
	series x=Month y=cholmean / group=group  markers markerattrs=(size=5pt);
	yaxis grid;
	label Month = 'Month';
	label cholmean = 'Serum cholesterol (mg/dL)';
	xaxis 
		values = ("y1" "y2" "y3" "y4" "y5")
		valuesdisplay = ("Baseline" "6" "12" "20" "24"); *re-label x-axis tick marks;
run;
title;


****************************:
*  3bc) Chapter 5, Section 5.4 Interaction effect;
* (T5.3-5.5) Analysis of Response Profiles of data on Blood Lead Levels";

data chol_long (drop = y1-y5); 
	set A;
baseline=y1;
month = 0; chol = y1; output;
month = 6; chol = y2; output;
month = 12; chol = y3; output;
month = 20; chol = y4; output;
month = 24; chol = y5; output;
run;

proc sort data = chol_long;
by group month;
run;

* Interaction contrast treatment effect over time;
proc mixed data = chol_long;
class month(ref='0') id group(ref='1'); 		  *set month 0 as baseline, group=1 as baseline;
model chol = month group month * group / s chisq; *fit a model with interaction mixed vars;
repeated month / type = un subject=id r; 		  *longitudinal analyis, repeated over time, unstructure cov matrix, subject=id;
run;

*************************************************;
* 3d) Chapter 5, Section 5.7;

* Create dummy variables for group and time;
data chol_ind;
  set chol_long;
  treat = (group = 2); *indicator, if group=2(high dose), then treat=1;
  t0 = (month = 0); * if month=0, then t0=1;
  t6 = (month = 6);
  t12 = (month = 12);
  t20 = (month = 20);
  t24 = (month =24);
 run;

*(T5.7) 3d) Analysis of Response Profiles assuming equal mean at Baseline";
proc mixed data = chol_ind order=data;
     class id month;
     model chol = t6 t12 t20 t24 treat*t6 treat*t12 treat*t20 treat*t24 / s chisq;
     repeated month / type=un subject=id r;
     contrast '4 DF Test of Interaction' 
          treat*t6 1, treat*t12 1, treat*t20 1, treat*t24 1 / e chisq;
run;

******************************************;
*3e) Quadratic function of longitudinal analysis; 
data quad;
	set chol_long;
	x=(group=2); 		*indicator variable, if group=2(high dose) then x=1;
	
	t=month; 			*time;
	t2=month*month; 	*time^2;
						
						*Interaction terms;
	xt=x*month; 		*indicator*time;
	xt2=x*t2;			*indicator*time^2;
run;

proc mixed data = quad order=data;
     class id month;
     model chol = t t2 xt xt2 / s chisq;
     repeated month / type=un subject=id r;
	 contrast 'Quadratic test of interaction' xt 1, xt2 1 / e chisq;
run;

* 3f/g) Changes from baseline two-way ANOVA; 
title2 "3. (T5.8) Analysis of Response Profiles of Changes from Baseline";
data diffc;
  set chol_ind;
  if (month > 0);		  		*remove month=0, keep other months;
  change = chol-baseline; 		*mean diff (y-baseline);
  mbaseline = baseline - 230 ;	*mean-centered baseline; *reference group becomes month 6;
run;
 

proc mixed data=diffc  order=data;
     class id month;
     model change = treat t12 t20 t24  treat*t12 treat*t20 treat*t24/ s chisq;
     repeated month / type=un subject=id r;
     contrast '4DF Test of Main Effect and Interaction' 
          treat 1,treat*t12 1, treat*t20 1, treat*t24 1 / chisq;
run;
title;


* 3hi) Adjusted Changes from baseline;
title2 "4. (T5.9) Analysis of Response Profiles of Adjusted Changes from Baseline";
proc mixed data = diffc order=data;
     class id month;
     model change = mbaseline treat t12 t20 t24  treat*t12 treat*t20 treat*t24/ s chisq;
     repeated month / type=un subject=id r;
     contrast '4DF Test of Main Effect and Interaction' 
          treat 1,treat*t12 1, treat*t20 1, treat*t24 1 / chisq;

run;

* 3j) LRT same model as 3hi;
* full model with interaction;
title2 "7A. Full likelihood ratio tests for BETA are valid";

proc mixed data = diffc order=data method=ML;
     class id month group;
     model change = mbaseline treat t12 t20 t24  treat*t12 treat*t20 treat*t24/ s chisq;
     repeated month / type=un subject=id r=1,100  group=group;
  	 contrast '4DF Test of Main Effect and Interaction' 
          treat 1,treat*t12 1, treat*t20 1, treat*t24 1 / chisq;
run;

*3j) LRT 
* reduced model without interaction;
title2 "7B. Full likelihood ratio tests for BETA are valid";
proc mixed data = diffc order=data method=ML;
     class id month group;
     model change = mbaseline t12 t20 t24/ s chisq;
     repeated month / type=un subject=id r=1,100  group=group;
run;

