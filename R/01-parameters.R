# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
if (interactive()) {

	# Input
	root <- here::here()
	outfile <- file.path(root, "output", "parameters.qs2")

	# Packages
	source(file.path(root, "R", "modify_attach.R"))
	modify_attach(qs2, include.only = "qs_save")

} else {

	# Input
	args <- commandArgs(trailingOnly = TRUE)
	outfile  <- args[1]

	# Packages
	library(qs2, include.only = "qs_save")
}

# -------------------------------------------------------------------------
# Generate some of the data frame input
# -------------------------------------------------------------------------

# Age distribution
# TODO - Check these bounds are correct!
#      - They differ to others used.
age_dist <- data.frame(
	lower_closed = c(16, 36, 46,  66),
	upper_closed = c(35, 45, 65, 100),
	freq         = c(14, 32, 32,  22) / 100
)
age_dist <- transform(
	age_dist,
	age_point = round((lower_closed + upper_closed) / 2),
	age_range = sprintf("A_%dto%d", lower_closed, upper_closed)
)
age_dist <- age_dist[c("age_range", "lower_closed", "upper_closed", "age_point", "freq")]

# TPT eff
tpt_data <- data.frame(
	regimen = c("6H"  , "3HP", "3HR" ),
	eff     = c( 0.6  ,  0.64,  0.67 ),
	ae      = c( 0.027,  0.01,  0.068)
)

parameters <- list(

	# general inputs ----------------------------------------------------------
	dt = 1,
	n_samples = 1,
	t_end = 86,
	t_hor = 86,

	# cohort inputs -----------------------------------------------------------
	cohort_size = 10000,

	# Age eligibility for TPT (0 or 1)
	# TODO - these need checking with JV to ensure precise
	#      - working on assumption we want [16-35), [35-45), [45-65) and [65-100)
	age_eligibility = c(a16_35 = 1, a16_45 = 0, a16_65 = 0, a16_100 = 0),

	# Burden by country of origin distribution (% 0 to 100) (See load_WHO script)
	# TODO - check aforementioned script
	burden = c(`0-40` = 0.4, `40-100` = 8.6, `100-150` = 7, `150+` = 84),

	# Age distribution
	age_dist = age_dist,

	# TPT eff
	tpt_data = tpt_data,

	# TB epi inputs -----------------------------------------------------------
	#r0                  = 0.2, # TB cases arising from each new case
	#p_ptbld             = 3.05, # prev by age X * RR of COPD in PTBLD (Byrne IJID)
	#rr_ptbld_mu         = 1.14 - 1, # Romanowski (2.91) vs Menzies (1.14)
	#r_tbdur             = 1 / (66 / 365), # TB duration (ETS)
	#r_etbdur            = 1 / (91 / 365), # eTB duration (ETS)
	#r_txdur             = 1 / (230 / 365), # TX duration (ETS)
	#p_tbhosp            = 0.1, #0.345, # Fraction inpatint TB (experts)
	r_reactivation      = 0.0826094004, # Start value from Menzies (this is calibrated)
	r_reactivation_slow = 0.000594, # Start value from Menzies (this is calibrated)
	r_reactivation_fn_f = 0.0826094004, # Start value from Menzies (this is calibrated)
	r_reactivation_fn_s = 0.000594, # Start value from Menzies (this is calibrated)
	r_slow              = 0.8721665962, # Start value from Menzies (this is calibrated)
	ptbld_switch        = 1,

	# Model C  -----------------------------------------------------------
	## This parameters only work for model C, currently not in use
	# Prevalence and dist of recent infection
	p_y1 = 0.26 * 0.5,
	p_y2 = 0.26 * 0.5,
	p_y3 = (1 - 0.26) * 0.5,

	# Initial state from Horton et al
	p_inf   = 0.83,
	p_min   = 0.14,
	p_subc  = 1 - (0.83 + 0.14),
	p_selfc = 0,

	yr_fast      = 5,
	step_decline = 1,

	# LTBI cascade input ------------------------------------------------------
	#test = c("T-SPOT.TB", "QuantiFERON", "Tuberculin Skin Test (5mm)", "Tuberculin Skin Test (10mm)", "Tuberculin Skin Test (15mm)"),
	test           = "Tuberculin Skin Test (5mm)", # "C-TB"
	tpt            = "3HR", #c("3HR", "3HP", "6H"),
	true_pos_tspo  = 1,
	true_neg_tspo  = 1,
	true_pos_qft   = 0.99319,
	true_neg_qft   = 0.9694,
	true_pos_tst5  = 1,
	true_neg_tst5  = 0.70188,
	true_pos_ctb   = 0.99319,#0.99434,
	true_neg_ctb   = 0.9694,#0.95706,
	
	
	
	#true_prev      = 0.2084517, # -TSPOT positivity in predict, taken as gold standard

	# LTBI Cascade components (choose 0 to 1)
	tst_return     = 1-0.15,  # Only appplies if using TST test or C-TB 
	test_2nd       = sqrt(0.45) , # Proportion attends second test appointment, costs staff time
	tpt_start      = sqrt(0.45), #0.66, # Loutet. changed to Latent tuberculosis testing and treatment programme for migrants
	tpt_completion = 0.75,#0.95, # Loutet

	# fraction not completing enjoying full benefits
	tpt_lfup_eff = 0,



	# Epi parameters distributions --------------------------------------------
	true_prev_dist	     = "PERT", # c("Beta", "PERT"),
	true_prev_beta	     = 0.20,#84517,
	true_prevsd_beta	   = 0.02,
	true_prev_pert	     = 0.20,#84517,
	true_prevmin_pert	   = 0.20,#84517*0.8,
	true_prevmax_pert	   = 0.20,#84517*1.2,
	true_prevlam_pert	   = 4,

	tpt_eff_dist	       = "PERT", # c("Beta", "PERT"),
	tpt_eff_beta	       = 0.67,
	tpt_effsd_beta	     = 0.1,
	tpt_eff_pert	       = 0.67,
	tpt_effmin_pert	     = 0.60,
	tpt_effmax_pert	     = 0.90,
	tpt_efflam_pert	     = 4,

	tpt_ae_dist	         = "PERT", # c("Beta", "PERT"),
	tpt_ae_beta	         = 1,
	tpt_aesd_beta	       = 0.2,
	tpt_ae_pert	         = 1,
	tpt_aemin_pert	     = 0.8,
	tpt_aemax_pert	     = 1.2,
	tpt_aelam_pert	     = 4,

	r_tbdur_dist	       = "Normal", # c("Beta", "PERT"),
	r_tbdur_normal	     = 78,
	r_tbdursd_normal	   = 0.5,
	r_tbdur_pert	       = 78,
	r_tbdurmin_pert	     = 78*0.8,
	r_tbdurmax_pert	     = 78*1.2,
	r_tbdurlam_pert	     = 4,

	r_etbdur_dist	     = "Normal", # c("Beta", "PERT"),
	r_etbdur_normal	     = 117,
	r_etbdursd_normal	 = 0.76,
	r_etbdur_pert	     = 117,
	r_etbdurmin_pert	 = 117*0.8,
	r_etbdurmax_pert	 = 117*1.2,
	r_etbdurlam_pert	 = 4,

	r_txdur_dist	     = "Normal", # c("Beta", "PERT"),
	r_txdur_normal	     = 230,
	r_txdursd_normal	 = 0.5,
	r_txdur_pert	     = 230,
	r_txdurmin_pert	     = 230*0.8,
	r_txdurmax_pert	     = 230*1.2,
	r_txdurlam_pert	     = 4,

	p_tbhosp_dist	     = "PERT", # c("Beta", "PERT"),
	p_tbhosp_normal	     = 0.1,
	p_tbhospsd_normal	 = 0.02,
	p_tbhosp_pert	     = 0.1,
	p_tbhospmin_pert	 = 0.1*0.8,
	p_tbhospmax_pert	 = 0.1*1.2,
	p_tbhosplam_pert	 = 4,

	p_ptbld_dist	     = "PERT", # c("Beta", "PERT"),
	p_ptbld_normal	     = 3.05,
	p_ptbldsd_normal	 = 0.2,
	p_ptbld_pert	     = 3.05,
	p_ptbldmin_pert	     = 3.05*0.8,
	p_ptbldmax_pert	     = 3.05*1.2,
	p_ptbldlam_pert	     = 4,

	rr_ptbld_mu_dist	 = "PERT", # c("Beta", "PERT"),
	rr_ptbld_mu_normal	 = 0.14,
	rr_ptbld_musd_normal = 0.02,
	rr_ptbld_mu_pert	 = 0.14,
	rr_ptbld_mumin_pert  = 0.14*0.8,
	rr_ptbld_mumax_pert  = 0.14*1.2,
	rr_ptbld_mulam_pert  = 4,


	r0_dist	            = "PERT", # c("Beta", "PERT"),
	r0_normal	        = 0.2,
	r0sd_normal         = 0.02,
	r0_pert	            = 0.2,
	r0min_pert          = 0.2*0.8,
	r0max_pert          = 0.2*1.2,
	r0lam_pert          = 4,


	# Cost distributions ------------------------------------------------------
	p_dr                  = 0.035, # discount
	adversecost_dist	  = "PERT", # c("PERT", "Gamma"),
	cost_adverse_gamma	  = 1500,
	cost_adversesd_gamma  = 100,
	cost_adverse_pert	  =	1742,
	cost_adversemin_pert  =	1742 * 0.75,
	cost_adversemax_pert  =	1742 * 1.25,
	cost_adverselam	      =	4,

	testcost_qf_dist	  =	"Gamma", # c("Gamma", "PERT"),
	cost_test_qf_gamma	  =	(25 + 31.5), # 62.7,
	cost_testsd_qf_gamma  =	10,
	cost_test_qf_pert	    =	(25 + 31.5),
	cost_testmin_qf_pert  =	(25 + 31.5) * 0.75,
	cost_testmax_qf_pert  =	(25 + 31.5) * 1.25,
	cost_testqf_lam	      =	4,

	
	testcost_tst_dist	  =	"Gamma", # c("Gamma", "PERT"),
	cost_test_tst_gamma	  =	60.7,
	cost_testsd_tst_gamma =	10,
	cost_test_tst_pert	  =	60.7,
	cost_testmin_tst_pert =	40,
	cost_testmax_tst_pert =	110,
	cost_testtst_lam	    =	4,
	
	
	testcost_ctb_dist	  =	"Gamma", # c("Gamma", "PERT"),
	cost_test_ctb_gamma	  =	20,
	cost_testsd_ctb_gamma =	3.5,
	cost_test_ctb_pert	  =	20,
	cost_testmin_ctb_pert =	5,
	cost_testmax_ctb_pert =	30,
	cost_testctb_lam	    =	4,
	
	testcost_staff_ctb_dist	=	"Gamma", # c("Gamma", "PERT"),
	cost_test_staff_ctb_gamma	  =	31.5,
	cost_testsd_staff_ctb_gamma =	5,
	cost_test_staff_ctb_pert	  =	31.5,
	cost_testmin_staff_ctb_pert =	20,
	cost_testmax_staff_ctb_pert =	40,
	cost_teststaffctb_lam	    =	4,

	
	testcost_tsp_dist	  =	"Gamma", # c("Gamma", "PERT"),
	cost_test_tsp_gamma	  =	67.61,
	cost_testsd_tsp_gamma =	10,
	cost_test_tsp_pert	  =	67.61,
	cost_testmin_tsp_pert =	40,
	cost_testmax_tsp_pert =	100,
	cost_testtsp_lam	    =	4,

	test2_negcost_tst_dist	   =	"Gamma", # c("Gamma", "PERT"),
	cost_test2_neg_tst_gamma	 =	19.2,
	cost_test2_negsd_tst_gamma =	5,
	cost_test2_neg_tst_pert	   =	19.2,
	cost_test2_negmin_tst_pert =	10,
	cost_test2_negmax_tst_pert =	30,
	cost_test2_negtst_lam	     =	4,
	
	test2_negcost_qf_dist	     =	"Gamma", # c("Gamma", "PERT"),
	cost_test2_neg_qf_gamma	   =	0,
	cost_test2_negsd_qf_gamma  =	0,
	cost_test2_neg_qf_pert	   =	0,
	cost_test2_negmin_qf_pert  =	0,
	cost_test2_negmax_qf_pert  =	0,
	cost_test2_negqf_lam	    = 	0,
	
	test2_negcost_ctb_dist	   =	"Gamma", # c("Gamma", "PERT"),
	cost_test2_neg_ctb_gamma	 =	19.2,
	cost_test2_negsd_ctb_gamma =	5,
	cost_test2_neg_ctb_pert	   =	19.2,
	cost_test2_negmin_ctb_pert =	10,
	cost_test2_negmax_ctb_pert =	30,
	cost_test2_negctb_lam	     =	4,
	
	test2_negcost_tsp_dist	  =	"Gamma", # c("Gamma", "PERT"),
	cost_test2_neg_tsp_gamma	 =	0,
	cost_test2_negsd_tsp_gamma =	0,
	cost_test2_neg_tsp_pert	   =	0,
	cost_test2_negmin_tsp_pert =	0,
	cost_test2_negmax_tsp_pert =	0,
	cost_test2_negtsp_lam	     =	0,
	
	
	test2_poscost_tst_dist	   =	"Gamma", # c("Gamma", "PERT"),
	cost_test2_pos_tst_gamma	 =	19.2,
	cost_test2_possd_tst_gamma =	5,
	cost_test2_pos_tst_pert	   =	19.2,
	cost_test2_posmin_tst_pert =	10,
	cost_test2_posmax_tst_pert =	30,
	cost_test2_postst_lam	     =	4,
	
	test2_poscost_qf_dist	     =	"Gamma", # c("Gamma", "PERT"),
	cost_test2_pos_qf_gamma	   =	37.7,
	cost_test2_possd_qf_gamma  =	5,
	cost_test2_pos_qf_pert	   =	37.7,
	cost_test2_posmin_qf_pert  =	20,
	cost_test2_posmax_qf_pert  =	50,
	cost_test2_posqf_lam	    = 	4,
	
	test2_poscost_ctb_dist	   =	"Gamma", # c("Gamma", "PERT"),
	cost_test2_pos_ctb_gamma	 =	50.9,
	cost_test2_possd_ctb_gamma =	5,
	cost_test2_pos_ctb_pert	   =	50.9,
	cost_test2_posmin_ctb_pert =	20,
	cost_test2_posmax_ctb_pert =	70,
	cost_test2_posctb_lam	     =	4,
	
	test2_poscost_tsp_dist	  =	"Gamma", # c("Gamma", "PERT"),
	cost_test2_pos_tsp_gamma	 =	0,
	cost_test2_possd_tsp_gamma =	0,
	cost_test2_pos_tsp_pert	   =	0,
	cost_test2_posmin_tsp_pert =	0,
	cost_test2_posmax_tsp_pert =	0,
	cost_test2_postsp_lam	     =	0,
	
	
	
	positivecost_dist	    =	"PERT", # c("PERT", "Gamma"),
	cost_positive_gamma	  =	78.76, #49, Using inflated Pareek cost of per LTBI diagnosed minus test costs
	cost_positivesd_gamma =	20, #10,
	cost_positive_pert	  =	78.76,
	cost_positivemin_pert = 78.76 * 0.75,
	cost_positivemax_pert = 78.76 * 1.25,
	cost_positivelam	    = 4,

	tptcost_6h_dist	      = "Gamma", # c("Gamma", "PERT"),
	cost_tpt_6h_gamma	  =	408,
	cost_tptsd_6h_gamma	  =	80,
	cost_tpt_6h_pert	  =	408,
	cost_tptmin_6h_pert	  =	80,
	cost_tptmax_6h_pert	  =	600,
	cost_tpt6h_lam	      =	4,

	tptcost_3hr_dist	  =	"Gamma", # c("Gamma", "PERT"),
	cost_tpt_3hr_gamma	  =	192.8,
	cost_tptsd_3hr_gamma  =	20,
	cost_tpt_3hr_pert	  =	192.8,
	cost_tptmin_3hr_pert  =	192.8 * 0.75,
	cost_tptmax_3hr_pert  =	192.8 * 1.25,
	cost_tpt3hr_lam	      =	4,

	tptcost_3hp_dist	  =	"Gamma", # c("Gamma", "PERT"),
	cost_tpt_3hp_gamma	  =	560,
	cost_tptsd_3hp_gamma  =	80,
	cost_tpt_3hp_pert	  =	560,
	cost_tptmin_3hp_pert  =	350,
	cost_tptmax_3hp_pert  =	700,
	cost_tpt3hp_lam	      =	4,

	tbtxcost_dist	      = "PERT", # c("PERT", "Gamma"),
	cost_tbtx_gamma	      =	2597, # Cavany and Abubabkar range1565,#2597, changed to Abubakar inflated, also in Abubakar 2018
	cost_tbtxsd_gamma	  =	150, #300,
	cost_tbtx_pert	      =	1644, #2081,
	cost_tbtxmin_pert	  =	1361,
	cost_tbtxmax_pert	  = 1928,
	cost_tbtxlam	      =	4,

	tbdx_ptbcost_dist	   =	"PERT", # c("PERT", "Gamma"),
	cost_tbdx_ptb_gamma	   =	766.95, # IDEA
	cost_tbdx_ptbsd_gamma  =	50,
	cost_tbdx_ptb_pert     =	766.95,
	cost_tbdx_ptbmin_pert  =	766.95 * 0.75,
	cost_tbdx_ptbmax_pert  =	766.95 * 1.25,
	cost_tbdx_ptblam	   =	4,

	tbdx_eptbcost_dist	    =	"PERT", # c("PERT", "Gamma"),
	cost_tbdx_eptb_gamma	=	1113.63, # IDEA
	cost_tbdx_eptbsd_gamma  =	100,
	cost_tbdx_eptb_pert     =	1113.63,
	cost_tbdx_eptbmin_pert  =	1113.63 * 0.75,
	cost_tbdx_eptbmax_pert  =	1113.63 * 1.25,
	cost_tbdx_eptblam	    =	4,

	tbhocost_dist	      =	"PERT", # c("PERT", "Gamma"),
	cost_tbho_gamma	      =	5207,
	cost_tbhosd_gamma	  =	500,
	cost_tbho_pert	      =	2677, #5207,  Pareek inpatient 2010 inflaated to 2022
	cost_tbhomin_pert	  =	2677 * 0.75,
	cost_tbhomax_pert	  = 2677 * 1.25,
	cost_tbholam	      =	4,

	postcost_dist	      =	"PERT", # c("PERT", "Gamma"), # cost of COPD
	cost_post_gamma	      =	779,
	cost_postsd_gamma	  =	30,
	cost_post_pert	      =	1275,
	cost_postmin_pert	  =	861,
	cost_postmax_pert	  = 1690,
	cost_postlam	      =	4,


	# QOL input ---------------------------------------------------------------
	ptbqol_dist	     = "Beta", # c("Beta", "PERT"),
	qol_ptb_beta	 = 0.72,
	qol_ptbsd_beta	 = 0.02,
	qol_ptb_pert	 = 0.72,
	qol_ptbmin_pert	 = 0.65,
	qol_ptbmax_pert	 = 0.79,
	qol_ptblam_pert	 = 4,

	eptbqol_dist	 = "Beta", # c("Beta", "PERT"),
	qol_eptb_beta	 = 0.8, # 0.72, using Hassan et al QoL, for 40% (0.5 QoL) Linadenitis and 60% pleuritis (1 QoL)
	qol_eptbsd_beta	 = 0.02,
	qol_eptb_pert	 = 0.80,
	qol_eptbmin_pert = 0.71,
	qol_eptbmax_pert = 0.89,
	qol_eptblam_pert = 4,

	ptxqol_dist	     = "PERT", # c("PERT", "Beta"),
	qol_ptx_beta	 = 0.88, # 0.8182, # dale et al based on Bauer et al
	qol_ptxsd_beta	 = 0.13,
	qol_ptx_pert	 = 0.8182,
	qol_ptxmin_pert	 = 0.7,
	qol_ptxmax_pert	 = 0.85,
	qol_ptxlam_pert	 = 4,

	postqol_dist	 = "Beta", # c("Beta", "PERT"),
	qol_post_beta	 = 0.947, # Quaife et al , based on COPD prevalence
	qol_postsd_beta	 = 0.01,
	qol_post_pert	 = 0.947,
	qol_postmin_pert = 0.92,
	qol_postmax_pert = 0.99,
	qol_postlam_pert = 4,

	aeqol_dist	     = "PERT", #c("PERT", "Beta"),
	qol_ae_beta	     = 0.25,
	qol_aesd_beta	 = 0.05,
	qol_ae_pert	     = 0.004794521, # Campbell
	qol_aemin_pert	 = 0.004794521 * 0.75,
	qol_aemax_pert	 = 0.004794521 * 1.25,
	qol_aelam_pert	 = 4,

	tbhoqol_dist	 = "PERT", #c("PERT", "Beta"),
	qol_tbho_gamma	 = 0.015,
	qol_tbhosd_gamma = 0.01,
	qol_tbho_pert	 = 0.015,
	qol_tbhomin_pert = 0.015 * 0.75,
	qol_tbhomax_pert = 0.015 * 1.25,
	qol_tbholam	     = 4
)


# -------------------------------------------------------------------------
# Save the parameters -----------------------------------------------------
# -------------------------------------------------------------------------
qs_save(parameters, outfile)
