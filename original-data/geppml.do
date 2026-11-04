*************************************
*This file contains the GE PPML code*
*************************************
**********
* HEADER *
**********
clear all
cd "C:\Users\muham\EGEI Dissertation\afcfta\code\GEPPML_data_code\GEPPML_stata"
set mem 1g
set matsize 8000
cap log close
log using geppml, text replace
set more off
pause on

******************
*I. Prepapre Data*
******************
use ge_ppml_data_orig, replace 
*******************************
*1. Create aggregate variables*
*******************************
bysort exporter: egen output=sum(trade)
bysort importer: egen expndr=sum(trade)
****************************************
*2. Chose a country for reference group*
****************************************
replace exporter="ZZZ" if exporter=="DEU"
replace importer="ZZZ" if importer=="DEU"
gen expndr_deu0=expndr if importer=="ZZZ"
egen expndr_deu=mean(expndr_deu0)
************************
*3 Create Fixed Effects*
************************
qui tab exporter, gen(exp_fe_)
qui tab importer, gen(imp_fe_)
******************************
*4. Set additional parameters*
******************************
scalar sigma=7
save ge_ppml_data, replace

**************************
*II. GE Analysis in Stata*
**************************
*****************************
*Step 1: `Baseline' Scenario*
*****************************
***************************************
*Step 1.a: Estimate `Baseline' Gravity*
***************************************
* Load data
use ge_ppml_data, clear
* Define number of countries
qui tab exporter
local NoC=r(r)
* Equation (9)
ppml trade LN_DIST CNTG BRDR exp_fe_* imp_fe_1-imp_fe_`=`NoC'-1', iter(30) noconst
* Predict trade in the baseline (bsln)
predict trade_bsln, mu
* Save estimates
scalar DIST_est=_b[LN_DIST]
scalar BRDR_est=_b[BRDR]
scalar CNTG_est=_b[CNTG]
* Create trade costs
gen t_ij_bsln=exp(_b[LN_DIST]*LN_DIST+_b[CNTG]*CNTG+_b[BRDR]*BRDR)
gen t_ij_ctrf=exp(_b[LN_DIST]*LN_DIST+_b[CNTG]*CNTG+_b[BRDR]*BRDR*0)
replace t_ij_ctrf=t_ij_bsln if exporter==importer
gen t_ij_ctrf_1=log(t_ij_ctrf)
******************************************
*Step 1.b: Construct `Baseline' GE Indexes*
******************************************
forvalues i=1(1)`=`NoC'-1'{
qui replace exp_fe_`i'=exp_fe_`i'*exp(_b[exp_fe_`i'])
qui replace imp_fe_`i'=imp_fe_`i'*exp(_b[imp_fe_`i'])
}
qui replace exp_fe_`NoC'=exp_fe_`NoC'*exp(_b[exp_fe_`NoC'])
qui replace imp_fe_`NoC'=imp_fe_`NoC'*exp(0)
egen all_exp_fes_0=rowtotal(exp_fe_1-exp_fe_`NoC')
egen all_imp_fes_0=rowtotal(imp_fe_1-imp_fe_`NoC') 
gen temp=all_exp_fes_0 if exporter==importer
bysort importer: egen all_exp_fes_0_TB=mean(temp)
drop temp
* Equation (7)
gen omr_bsln=output*expndr_deu/(all_exp_fes_0)
* Equation (8)
gen imr_bsln=expndr/(all_imp_fes_0*expndr_deu)
* Calculate real GDP
gen rGDP_bsln_temp=output/(imr_bsln^(1/(1-sigma))) if exporter==importer
bysort exporter: egen rGDP_bsln=sum(rGDP_bsln_temp)
* Calculate share of domestic expenditure 
gen acr_bsln=trade_bsln/expndr if exporter==importer
* Calcualte baseline bilateral exports
gen exp_bsln=trade_bsln if exporter!=importer
* Calculate domestic expenditure for each country
gen exp_bsln_acr=trade_bsln if exporter==importer
* Calculate total exports for each country
bysort exporter: egen tot_exp_bsln=sum(exp_bsln)

********************************
*Step 2: `Conditional' Scenario*
********************************
******************************************
*Step 2.a: Estimate `Conditional' Gravity*
******************************************
* Re-Create Fixed Effects
drop exp_fe_* imp_fe_*
qui tab exporter, gen(exp_fe_)
qui tab importer, gen(imp_fe_)
* Equation (11)
ppml trade exp_fe_* imp_fe_1-imp_fe_`=`NoC'-1', iter(30) noconst offset(t_ij_ctrf_1)
* Predict trade in `Conditional' GE counterfactual (cndl)
predict trade_cndl, mu
**********************************************
*Step 2.b: Construct `Conditional' GE Indexes*
**********************************************
forvalues i=1(1)`=`NoC'-1'{
qui replace exp_fe_`i'=exp_fe_`i'*exp(_b[exp_fe_`i'])
qui replace imp_fe_`i'=imp_fe_`i'*exp(_b[imp_fe_`i'])
}
qui replace exp_fe_`NoC'=exp_fe_`NoC'*exp(_b[exp_fe_`NoC'])
qui replace imp_fe_`NoC'=imp_fe_`NoC'*exp(0)
egen all_exp_fes_1=rowtotal(exp_fe_1-exp_fe_`NoC')
egen all_imp_fes_1=rowtotal(imp_fe_1-imp_fe_`NoC')
* Equation (7)
gen omr_cndl=output*expndr_deu/(all_exp_fes_1)
* Equation (8)
gen imr_cndl=expndr/(all_imp_fes_1*expndr_deu)
* Calcualte baseline bilateral exports
gen exp_cndl=trade_cndl if exporter!=importer
* Calculate total exports for each country
bysort exporter: egen tot_exp_cndl=sum(exp_cndl)
* Calculate change in total exports for each country
gen tot_exp_cndl_ch=(tot_exp_cndl-tot_exp_bsln)/tot_exp_bsln*100
* Calculate real GDP
gen rGDP_cndl_temp=output/(imr_cndl^(1/(1-sigma))) if exporter==importer
bysort exporter: egen rGDP_cndl=sum(rGDP_cndl_temp)

***********************************
*Step 3: `Full Endowment' Scenario*
***********************************
*********************************************
*Step 3.a: Estimate `Full Endowment' Gravity*
*********************************************
* Define variables for loop
local i=3
local diff_all_exp_fes_sd=1
local diff_all_exp_fes_max=1
gen double trade_1_pred=trade_cndl
gen double output_bsln=output
gen double expndr_bsln=expndr
gen double phi=expndr/output if exporter==importer
gen double temp=all_exp_fes_0 if exporter==importer
bysort importer: egen double all_exp_fes_0_imp=mean(temp)
drop temp* 

gen double expndr_temp_1=phi*output if exporter==importer
bysort importer: egen double expndr_1=mean(expndr_temp_1)
drop *temp*
gen double expndr_deu01_1=expndr_1 if importer=="ZZZ"
egen double expndr_deu1_1=mean(expndr_deu01_1)
egen double expndr_deu_1=mean(expndr_deu01_1)
gen double temp=all_exp_fes_1 if exporter==importer
bysort importer: egen double all_exp_fes_1_imp=mean(temp)
drop *temp*
gen double p_full_exp_0=0
gen double p_full_exp_1=(all_exp_fes_1/all_exp_fes_0)^(1/(1-sigma))
gen double p_full_imp_1=(all_exp_fes_1_imp/all_exp_fes_0_imp)^(1/(1-sigma))
gen double imr_full_1=expndr_1/(all_imp_fes_1*expndr_deu_1)
gen double imr_full_ch_1=1
gen double omr_full_1=output*expndr_deu_1/(all_exp_fes_1)
gen double omr_full_ch_1=1

* Start loop with convergence criteria that the mean and the standard devition of the difference in each of the factory-gate prices between two subsequent iterations is smaller than the pre-defined tolerance criterium of 0.001
while (`diff_all_exp_fes_sd'>0.001) | (`diff_all_exp_fes_max'>0.001) {
* Equation (14)
gen double trade_`=`i'-1'=trade_`=`i'-2'_pred*p_full_exp_`=`i'-2'*p_full_imp_`=`i'-2'/(omr_full_ch_`=`i'-2'*imr_full_ch_`=`i'-2')
drop exp_fe_* imp_fe_*
qui tab exporter, gen (exp_fe_)
qui tab importer, gen (imp_fe_)
* Equation (11)
capture glm trade_`=`i'-1' exp_fe_* imp_fe_*, family(poisson) offset(t_ij_ctrf_1) noconst irls iter(30) 
* Predict trade flows
predict trade_`=`i'-1'_pred, mu
forvalues j=1(1)`NoC'{
qui replace exp_fe_`j'=exp_fe_`j'*exp(_b[exp_fe_`j'])
qui replace imp_fe_`j'=imp_fe_`j'*exp(_b[imp_fe_`j'])
}
egen double all_exp_fes_`=`i'-1'=rowtotal(exp_fe_1-exp_fe_`NoC') 
egen double all_imp_fes_`=`i'-1'=rowtotal(imp_fe_1-imp_fe_`NoC')
* Update output
bysort exporter: egen double output_`=`i'-1'=total(trade_`=`i'-1'_pred)
* Update expenditure
bysort importer: egen double expndr_check_`=`i'-1'=total(trade_`=`i'-1'_pred)
gen double expndr_deu0_`=`i'-1'=expndr_check_`=`i'-1' if importer=="ZZZ"
egen double expndr_deu_`=`i'-1'=mean(expndr_deu0_`=`i'-1')
gen double temp=all_exp_fes_`=`i'-1' if exporter==importer
bysort importer: egen double all_exp_fes_`=`i'-1'_imp=mean(temp)
drop temp* 
* Update factory-gate prices
gen double p_full_exp_`=`i'-1'=((all_exp_fes_`=`i'-1'/all_exp_fes_`=`i'-2')/(expndr_deu_`=`i'-1'/expndr_deu_`=`i'-2'))^(1/(1-sigma))
gen double p_full_imp_`=`i'-1'=((all_exp_fes_`=`i'-1'_imp/all_exp_fes_`=`i'-2'_imp)/(expndr_deu_`=`i'-1'/expndr_deu_`=`i'-2'))^(1/(1-sigma))
* Equation (7)
gen double omr_full_`=`i'-1'=output_`=`i'-1'/all_exp_fes_`=`i'-1'
gen double omr_full_ch_`=`i'-1'=omr_full_`=`i'-1'/omr_full_`=`i'-2'
* Update expenditure
gen double expndr_temp_`=`i'-1'=phi*output_`=`i'-1' if exporter==importer
bysort importer: egen double expndr_`=`i'-1'=mean(expndr_temp_`=`i'-1')
* Equation (8)
gen double imr_full_`=`i'-1'=expndr_`=`i'-1'/(all_imp_fes_`=`i'-1'*expndr_deu_`=`i'-1')
gen double imr_full_ch_`=`i'-1'=imr_full_`=`i'-1'/imr_full_`=`i'-2'
* Convergence criteria in terms of changes in factory-gate prices
gen double diff_p_full_exp_`=`i'-1'=p_full_exp_`=`i'-2'-p_full_exp_`=`i'-3'
sum diff_p_full_exp_`=`i'-1'
local diff_all_exp_fes_sd=r(sd)
local diff_all_exp_fes_max=abs(r(max))
local i=`i'+1
}
*************************************************
*Step 3.a: Construct `Full Endowment' GE Indexes*
*************************************************
* Calculate p^c/p
drop output expndr expndr_deu0
rename expndr_deu expndr_deu_bsln
bysort exporter: egen double output=sum(trade_`=`i'-2'_pred)
gen double expndr_temp=phi*output if exporter==importer
bysort importer: egen double expndr=mean(expndr_temp)
gen double expndr_deu0=expndr if importer=="ZZZ"
egen double expndr_deu=mean(expndr_deu0) 
forvalues j=1(1)`NoC'{
qui replace imp_fe_`j'=imp_fe_`j'*exp(_b[imp_fe_`j'])
}
gen double p_full=((all_exp_fes_`=`i'-2'/all_exp_fes_0)/(expndr_deu/expndr_deu_bsln))^(1/(1-sigma))
* Calculate output
gen output_full=p_full*output_bsln
* Equation (7)
gen double omr_full=output_full*expndr_deu/(all_exp_fes_`=`i'-2')
* Equation (8)
egen double all_imp_fes=rowtotal(imp_fe_1-imp_fe_`NoC')
gen double imr_full=expndr/(all_imp_fes_`=`i'-2'*expndr_deu)
* Calculate real GDP
gen double rGDP_full_temp=p_full*output_bsln/(imr_full^(1/(1-sigma))) if exporter==importer
bysort exporter: egen double rGDP_full=sum(rGDP_full_temp)
* Calculate expenditure
gen double expndr_full_temp=phi*output_full if exporter==importer
bysort importer: egen double expndr_full=mean(expndr_full_temp)
* Calcualte baseline bilateral exports
gen trade_full=(output_full*expndr_full*t_ij_ctrf)/(imr_full*omr_full)
gen double exp_full=trade_full if exporter!=importer
* Calculate total exports for each country
bysort exporter: egen double tot_exp_full=sum(exp_full)
* Calculate change in total exports for each country
gen double tot_exp_full_ch=(tot_exp_full-tot_exp_bsln)/tot_exp_bsln*100
* Save results
save full_static_all, replace

******************************************
*Prepare All Indexes at the Country Level*
******************************************
collapse imr_full imr_bsln imr_cndl, by(importer)
rename importer country
gen imr_full_ch=(imr_full-imr_bsln)/imr_bsln*100
gen imr_cndl_ch=(imr_cndl-imr_bsln)/imr_bsln*100
save imrs_all, replace
use full_static_all, clear
collapse omr_full omr_cndl omr_bsln rGDP_full rGDP_cndl rGDP_bsln tot_exp_* p_full* output_bsln acr*, by(exporter)
gen omr_full_ch=(omr_full-omr_bsln)/omr_bsln*100
gen omr_cndl_ch=(omr_cndl-omr_bsln)/omr_bsln*100
gen rGDP_full_ch=(rGDP_full-rGDP_bsln)/rGDP_bsln*100
gen rGDP_cndl_ch=(rGDP_cndl-rGDP_bsln)/rGDP_bsln*100
rename exporter country
joinby country using imrs_all
save all_indexes_geppml, replace
log close
