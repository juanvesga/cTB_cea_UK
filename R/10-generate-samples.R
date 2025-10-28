# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
if (interactive()) {

	# Input
	root <- here::here()
	infile_p_cfr_ets     <- file.path(root, "output", "p_cfr_ets.qs2")
	infile_p_ptb_ets     <- file.path(root, "output", "p_ptb_ets.qs2")
	infile  <- file.path(root, "output", "parameters.qs2")
	outfile <- file.path(root, "output", "samples.qs2")

	# Packages
	source(file.path(root, "R", "modify_attach.R"))
	modify_attach(qs2,     include.only = c("qs_read", "qs_save"))
	modify_attach(freedom, include.only = "rpert")

} else {

	# Input
	args <- commandArgs(trailingOnly = TRUE)
	infile  <- args[1]
	outfile  <- args[2]

	# Packages
	library(qs2,     include.only = c("qs_read", "qs_save"))
	library(freedom, include.only = "rpert")
}

# -------------------------------------------------------------------------
# Load parameters ---------------------------------------------------------
# -------------------------------------------------------------------------
parameters <- qs_read(infile)

# -------------------------------------------------------------------------
# Helper functions --------------------------------------------------------
# -------------------------------------------------------------------------
# Truncated normal
rtruncnorm <- function(n, mu, sigma, low, high) {
	# find quantiles that correspond the the given low and high levels.
	p_low <- pnorm(low, mu, sigma)
	p_high <- pnorm(high, mu, sigma)

	# draw quantiles uniformly between the limits and pass these
	# to the relevant quantile function.
	if(n>1){
		outres<-qnorm(runif(n, p_low, p_high), mu, sigma)
	}else{
		outres<-1
	}
	return(outres)

}


rnorm_ptb <- function(n, mu, sigma) {

	if(n>1){
		outres<-rnorm(n, mu, sigma)
	}else{
		outres<-mu
	}

	return(outres)

}

# Sample QALY and Costs (will take mean for n_samples ==1)
.gamma_or_pert_samples <- function(n, dist, gamma_mean, gamma_sd, pert_mean, pert_min, pert_max, pert_lambda) {
	if (n == 1) {
		switch(dist, Gamma = gamma_mean, PERT = pert_mean, stop())
	} else {
		switch(
			dist,
			Gamma = rgamma(n, shape = (gamma_mean^2) / (gamma_sd^2), rate = gamma_mean / (gamma_sd^2)),
			PERT = rpert(n, x.min = pert_min, x.max = pert_max, x.mode = pert_mean, lambda = pert_lambda),
			stop()
		)
	}
}

.estBetaParams <- function(mu, var) {
	alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
	beta <- alpha * (1 / mu - 1)
	list(alpha = alpha, beta = beta)
}


.beta_or_pert_samples <- function(n, dist, beta_mean, beta_sd, pert_mean, pert_min, pert_max, pert_lambda) {
	if (n == 1) {
		switch(dist, Beta = beta_mean, PERT = pert_mean, stop())
	} else {
		switch(
			dist,
			Beta = {
				mu <- beta_mean
				var <- beta_sd^2
				alpha <- ((1 - mu) / var - 1 / mu) * mu ^ 2
				beta <- alpha * (1 / mu - 1)
				rbeta(n, shape1 = alpha, shape2 = beta)
			},
			PERT = rpert(n, x.min = pert_min, x.max = pert_max, x.mode = pert_mean, lambda = pert_lambda),
			stop()
		)
	}
}

.normal_or_pert_samples <- function(n, dist, normal_mean, normal_sd, pert_mean, pert_min, pert_max, pert_lambda) {
	if (n == 1) {
		switch(dist, Normal = normal_mean, PERT = pert_mean, stop())
	} else {
		switch(
			dist,
			Normal = rnorm(n, normal_mean, normal_sd),
			PERT = rpert(n, x.min = pert_min, x.max = pert_max, x.mode = pert_mean, lambda = pert_lambda),
			stop()
		)
	}
}

# -------------------------------------------------------------------------
# Generate samples --------------------------------------------------------
# -------------------------------------------------------------------------
samples <- with(parameters, list(

	# IGRA positivity
	sim_true_prev = .beta_or_pert_samples(
		n_samples,
		dist        = true_prev_dist,
		beta_mean   = true_prev_beta,
		beta_sd     = true_prevsd_beta,
		pert_mean   = true_prev_pert,
		pert_min    = true_prevmin_pert,
		pert_max    = true_prevmax_pert,
		pert_lambda = true_prevlam_pert
	),

	# TPT efficacy
	sim_tpt_eff = .beta_or_pert_samples(
		n_samples,
		dist        = tpt_eff_dist,
		beta_mean   = tpt_eff_beta,
		beta_sd     = tpt_effsd_beta,
		pert_mean   = tpt_eff_pert,
		pert_min    = tpt_effmin_pert,
		pert_max    = tpt_effmax_pert,
		pert_lambda = tpt_efflam_pert
	),

	# TPT adverse events
	sim_tpt_ae = .beta_or_pert_samples(
		n_samples,
		dist        = tpt_ae_dist,
		beta_mean   = tpt_ae_beta,
		beta_sd     = tpt_aesd_beta,
		pert_mean   = tpt_ae_pert,
		pert_min    = tpt_aemin_pert,
		pert_max    = tpt_aemax_pert,
		pert_lambda = tpt_aelam_pert
	),

	# PTB duration
	sim_r_tbdur = .normal_or_pert_samples(
		n_samples,
		dist        = r_tbdur_dist,
		normal_mean   = r_tbdur_normal,
		normal_sd     = r_tbdursd_normal,
		pert_mean   = r_tbdur_pert,
		pert_min    = r_tbdurmin_pert,
		pert_max    = r_tbdurmax_pert,
		pert_lambda = r_tbdurlam_pert
	),

	# EPTB duration
	sim_r_etbdur = .normal_or_pert_samples(
		n_samples,
		dist        = r_etbdur_dist,
		normal_mean   = r_etbdur_normal,
		normal_sd     = r_etbdursd_normal,
		pert_mean   = r_etbdur_pert,
		pert_min    = r_etbdurmin_pert,
		pert_max    = r_etbdurmax_pert,
		pert_lambda = r_etbdurlam_pert
	),

	# TB Tx duration
	sim_r_txdur = .normal_or_pert_samples(
		n_samples,
		dist        = r_txdur_dist,
		normal_mean   = r_txdur_normal,
		normal_sd     = r_txdursd_normal,
		pert_mean   = r_txdur_pert,
		pert_min    = r_txdurmin_pert,
		pert_max    = r_txdurmax_pert,
		pert_lambda = r_txdurlam_pert
	),


	# TB hospitalised
	sim_p_tbhosp = .normal_or_pert_samples(
		n_samples,
		dist        = p_tbhosp_dist,
		normal_mean   = p_tbhosp_normal,
		normal_sd     = p_tbhospsd_normal,
		pert_mean   = p_tbhosp_pert,
		pert_min    = p_tbhospmin_pert,
		pert_max    = p_tbhospmax_pert,
		pert_lambda = p_tbhosplam_pert
	),

	# PTBLD
	sim_p_ptbld = .normal_or_pert_samples(
		n_samples,
		dist        = p_ptbld_dist,
		normal_mean   = p_ptbld_normal,
		normal_sd     = p_ptbldsd_normal,
		pert_mean   = p_ptbld_pert,
		pert_min    = p_ptbldmin_pert,
		pert_max    = p_ptbldmax_pert,
		pert_lambda = p_ptbldlam_pert
	),

	# PTBLD mort RR
	sim_rr_ptbld_mu = .normal_or_pert_samples(
		n_samples,
		dist        = rr_ptbld_mu_dist,
		normal_mean   = rr_ptbld_mu_normal,
		normal_sd     = rr_ptbld_musd_normal,
		pert_mean   = rr_ptbld_mu_pert,
		pert_min    = rr_ptbld_mumin_pert,
		pert_max    = rr_ptbld_mumax_pert,
		pert_lambda = rr_ptbld_mulam_pert
	),


	# Secondary cases
	sim_r0 = .normal_or_pert_samples(
		n_samples,
		dist        = r0_dist,
		normal_mean   = r0_normal,
		normal_sd     = r0sd_normal,
		pert_mean   = r0_pert,
		pert_min    = r0min_pert,
		pert_max    = r0max_pert,
		pert_lambda = r0lam_pert
	),


	# Campaign costs
	sim_adverse_cost = .gamma_or_pert_samples(
		n_samples,
		dist        = adversecost_dist,
		gamma_mean  = cost_adverse_gamma,
		gamma_sd    = cost_adversesd_gamma,
		pert_mean   = cost_adverse_pert,
		pert_min    = cost_adversemin_pert,
		pert_max    = cost_adversemax_pert,
		pert_lambda = cost_adverselam
	),

	# QFT
	sim_test_cost_qfn = .gamma_or_pert_samples(
		n_samples,
		dist        = testcost_qf_dist,
		gamma_mean  = cost_test_qf_gamma,
		gamma_sd    = cost_testsd_qf_gamma,
		pert_mean   = cost_test_qf_pert,
		pert_min    = cost_testmin_qf_pert,
		pert_max    = cost_testmax_qf_pert,
		pert_lambda = cost_testqf_lam
	),

	# TSPOT
	sim_test_cost_tspo = .gamma_or_pert_samples(
		n_samples,
		dist        = testcost_tsp_dist,
		gamma_mean  = cost_test_tsp_gamma,
		gamma_sd    = cost_testsd_tsp_gamma,
		pert_mean   = cost_test_tsp_pert,
		pert_min    = cost_testmin_tsp_pert,
		pert_max    = cost_testmax_tsp_pert,
		pert_lambda = cost_testtsp_lam
	),

	# TST
	sim_test_cost_tst = .gamma_or_pert_samples(
		n_samples,
		dist        = testcost_tst_dist,
		gamma_mean  = cost_test_tst_gamma,
		gamma_sd    = cost_testsd_tst_gamma,
		pert_mean   = cost_test_tst_pert,
		pert_min    = cost_testmin_tst_pert,
		pert_max    = cost_testmax_tst_pert,
		pert_lambda = cost_testtsy_lam
	),
	
	sim_test_cost_ctb = .gamma_or_pert_samples(
	  n_samples,
	  dist        = testcost_ctb_dist,
	  gamma_mean  = cost_test_ctb_gamma,
	  gamma_sd    = cost_testsd_ctb_gamma,
	  pert_mean   = cost_test_ctb_pert,
	  pert_min    = cost_testmin_ctb_pert,
	  pert_max    = cost_testmax_ctb_pert,
	  pert_lambda = cost_testctb_lam
	),
	
	sim_test_staff_cost_ctb = .gamma_or_pert_samples(
	  n_samples,
	  dist        = testcost_staff_ctb_dist,
	  gamma_mean  = cost_test_staff_ctb_gamma,
	  gamma_sd    = cost_testsd_staff_ctb_gamma,
	  pert_mean   = cost_test_staff_ctb_pert,
	  pert_min    = cost_testmin_staff_ctb_pert,
	  pert_max    = cost_testmax_staff_ctb_pert,
	  pert_lambda = cost_teststaffctb_lam
	),
	
	# QFT
	sim_test2_neg_cost_qfn = .gamma_or_pert_samples(
	  n_samples,
	  dist        = test2_negcost_qf_dist	,
	  gamma_mean  = cost_test2_neg_qf_gamma,
	  gamma_sd    = cost_test2_negsd_qf_gamma,
	  pert_mean   = cost_test2_neg_qf_pert,
	  pert_min    = cost_test2_negmin_qf_pert,
	  pert_max    = cost_test2_negmax_qf_pert,
	  pert_lambda = cost_test2_negqf_lam
	),
	
	sim_test2_pos_cost_qfn = .gamma_or_pert_samples(
	  n_samples,
	  dist        = test2_poscost_qf_dist	,
	  gamma_mean  = cost_test2_pos_qf_gamma,
	  gamma_sd    = cost_test2_possd_qf_gamma,
	  pert_mean   = cost_test2_pos_qf_pert,
	  pert_min    = cost_test2_posmin_qf_pert,
	  pert_max    = cost_test2_posmax_qf_pert,
	  pert_lambda = cost_test2_posqf_lam
	),
	
	# TSPOT
	sim_test2_neg_cost_tspo = .gamma_or_pert_samples(
	  n_samples,
	  dist        = test2_negcost_tsp_dist,
	  gamma_mean  = cost_test2_neg_tsp_gamma,
	  gamma_sd    = cost_test2_negsd_tsp_gamma,
	  pert_mean   = cost_test2_neg_tsp_pert,
	  pert_min    = cost_test2_negmin_tsp_pert,
	  pert_max    = cost_test2_negmax_tsp_pert,
	  pert_lambda = cost_test2_negtsp_lam
	),
	
	sim_test2_pos_cost_tspo = .gamma_or_pert_samples(
	  n_samples,
	  dist        = test2_poscost_tsp_dist,
	  gamma_mean  = cost_test2_pos_tsp_gamma,
	  gamma_sd    = cost_test2_possd_tsp_gamma,
	  pert_mean   = cost_test2_pos_tsp_pert,
	  pert_min    = cost_test2_posmin_tsp_pert,
	  pert_max    = cost_test2_posmax_tsp_pert,
	  pert_lambda = cost_test2_postsp_lam
	),
	
	# TST
	sim_test2_neg_cost_tst = .gamma_or_pert_samples(
	  n_samples,
	  dist        = test2_negcost_tst_dist,
	  gamma_mean  = cost_test2_neg_tst_gamma,
	  gamma_sd    = cost_test2_negsd_tst_gamma,
	  pert_mean   = cost_test2_neg_tst_pert,
	  pert_min    = cost_test2_negmin_tst_pert,
	  pert_max    = cost_test2_negmax_tst_pert,
	  pert_lambda = cost_test2_negtst_lam
	),
	
	sim_test2_pos_cost_tst = .gamma_or_pert_samples(
	  n_samples,
	  dist        = test2_poscost_tst_dist,
	  gamma_mean  = cost_test2_pos_tst_gamma,
	  gamma_sd    = cost_test2_possd_tst_gamma,
	  pert_mean   = cost_test2_pos_tst_pert,
	  pert_min    = cost_test2_posmin_tst_pert,
	  pert_max    = cost_test2_posmax_tst_pert,
	  pert_lambda = cost_test2_postst_lam
	),
	
	# CTB
	
	sim_test2_neg_cost_ctb = .gamma_or_pert_samples(
	  n_samples,
	  dist        = test2_negcost_ctb_dist,
	  gamma_mean  = cost_test2_neg_ctb_gamma,
	  gamma_sd    = cost_test2_negsd_ctb_gamma,
	  pert_mean   = cost_test2_neg_ctb_pert,
	  pert_min    = cost_test2_negmin_ctb_pert,
	  pert_max    = cost_test2_negmax_ctb_pert,
	  pert_lambda = cost_test2_negctb_lam
	),
	
	sim_test2_pos_cost_ctb = .gamma_or_pert_samples(
	  n_samples,
	  dist        = test2_poscost_ctb_dist,
	  gamma_mean  = cost_test2_pos_ctb_gamma,
	  gamma_sd    = cost_test2_possd_ctb_gamma,
	  pert_mean   = cost_test2_pos_ctb_pert,
	  pert_min    = cost_test2_posmin_ctb_pert,
	  pert_max    = cost_test2_posmax_ctb_pert,
	  pert_lambda = cost_test2_posctb_lam
	),

	# Positive test assessment
	sim_positive_cost = .gamma_or_pert_samples(
		n_samples,
		dist        = positivecost_dist,
		gamma_mean  = cost_positive_gamma,
		gamma_sd    = cost_positivesd_gamma,
		pert_mean   = cost_positive_pert,
		pert_min    = cost_positivemin_pert,
		pert_max    = cost_positivemax_pert,
		pert_lambda = cost_positivelam
	),

	# TPT cost based on regimen
	sim_tpt_cost = switch(
		tpt,
		"3HR" = .gamma_or_pert_samples(
			n_samples,
			dist        = tptcost_3hr_dist,
			gamma_mean  = cost_tpt_3hr_gamma,
			gamma_sd    = cost_tptsd_3hr_gamma,
			pert_mean   = cost_tpt_3hr_pert,
			pert_min    = cost_tptmin_3hr_pert,
			pert_max    = cost_tptmax_3hr_pert,
			pert_lambda = cost_tpt3hr_lam
		),
		"3HP" = .gamma_or_pert_samples(
			n_samples,
			dist        = tptcost_3hp_dist,
			gamma_mean  = cost_tpt_3hp_gamma,
			gamma_sd    = cost_tptsd_3hp_gamma,
			pert_mean   = cost_tpt_3hp_pert,
			pert_min    = cost_tptmin_3hp_pert,
			pert_max    = cost_tptmax_3hp_pert,
			pert_lambda = cost_tpt3hp_lam
		),
		"6H" = .gamma_or_pert_samples(
			n_samples,
			dist        = tptcost_6h_dist,
			gamma_mean  = cost_tpt_6h_gamma,
			gamma_sd    = cost_tptsd_6h_gamma,
			pert_mean   = cost_tpt_6h_pert,
			pert_min    = cost_tptmin_6h_pert,
			pert_max    = cost_tptmax_6h_pert,
			pert_lambda = cost_tpt6h_lam
		)
	),

	# TB Dx cost
	sim_tbdx_ptb_cost = .gamma_or_pert_samples(
		n_samples,
		dist        = tbdx_ptbcost_dist,
		gamma_mean  = cost_tbdx_ptb_gamma,
		gamma_sd    = cost_tbdx_ptbsd_gamma,
		pert_mean   = cost_tbdx_ptb_pert,
		pert_min    = cost_tbdx_ptbmin_pert,
		pert_max    = cost_tbdx_ptbmax_pert,
		pert_lambda = cost_tbdx_ptblam
	),

	sim_tbdx_eptb_cost = .gamma_or_pert_samples(
		n_samples,
		dist        = tbdx_eptbcost_dist,
		gamma_mean  = cost_tbdx_eptb_gamma,
		gamma_sd    = cost_tbdx_eptbsd_gamma,
		pert_mean   = cost_tbdx_eptb_pert,
		pert_min    = cost_tbdx_eptbmin_pert,
		pert_max    = cost_tbdx_eptbmax_pert,
		pert_lambda = cost_tbdx_eptblam
	),

	# TB Tx cost
	sim_tbtx_cost = .gamma_or_pert_samples(
		n_samples,
		dist        = tbtxcost_dist,
		gamma_mean  = cost_tbtx_gamma,
		gamma_sd    = cost_tbtxsd_gamma,
		pert_mean   = cost_tbtx_pert,
		pert_min    = cost_tbtxmin_pert,
		pert_max    = cost_tbtxmax_pert,
		pert_lambda = cost_tbtxlam
	),

	# TB hospitalised cost
	sim_tbho_cost = .gamma_or_pert_samples(
		n_samples,
		dist        = tbhocost_dist,
		gamma_mean  = cost_tbho_gamma,
		gamma_sd    = cost_tbhosd_gamma,
		pert_mean   = cost_tbho_pert,
		pert_min    = cost_tbhomin_pert,
		pert_max    = cost_tbhomax_pert,
		pert_lambda = cost_tbholam
	),

	# Post-TB lung disease
	sim_post_cost = .gamma_or_pert_samples(
		n_samples,
		dist        = postcost_dist,
		gamma_mean  = cost_post_gamma,
		gamma_sd    = cost_postsd_gamma,
		pert_mean   = cost_post_pert,
		pert_min    = cost_postmin_pert,
		pert_max    = cost_postmax_pert,
		pert_lambda = cost_postlam
	),

	# TB QoL
	sim_tb_qol = .beta_or_pert_samples(
		n_samples,
		dist        = ptbqol_dist,
		beta_mean   = qol_ptb_beta,
		beta_sd     = qol_ptbsd_beta,
		pert_mean   = qol_ptb_pert,
		pert_min    = qol_ptbmin_pert,
		pert_max    = qol_ptbmax_pert,
		pert_lambda = qol_ptblam_pert
	),

	# EPTB Qol
	sim_eptb_qol = .beta_or_pert_samples(
		n_samples,
		dist        = eptbqol_dist,
		beta_mean   = qol_eptb_beta,
		beta_sd     = qol_eptbsd_beta,
		pert_mean   = qol_eptb_pert,
		pert_min    = qol_eptbmin_pert,
		pert_max    = qol_eptbmax_pert,
		pert_lambda = qol_eptblam_pert
	),

	sim_tx_qol = .beta_or_pert_samples(
		n_samples,
		dist        = ptxqol_dist,
		beta_mean   = qol_ptx_beta,
		beta_sd     = qol_ptxsd_beta,
		pert_mean   = qol_ptx_pert,
		pert_min    = qol_ptxmin_pert,
		pert_max    = qol_ptxmax_pert,
		pert_lambda = qol_ptxlam_pert
	),

	# POstTB QoL
	sim_post_qol = .beta_or_pert_samples(
		n_samples,
		dist        = postqol_dist,
		beta_mean   = qol_post_beta,
		beta_sd     = qol_postsd_beta,
		pert_mean   = qol_post_pert,
		pert_min    = qol_postmin_pert,
		pert_max    = qol_postmax_pert,
		pert_lambda = qol_postlam_pert
	),

	# AE TPT QoL
	sim_ae_qol = .beta_or_pert_samples(
		n_samples,
		dist        = aeqol_dist,
		beta_mean   = qol_ae_beta,
		beta_sd     = qol_aesd_beta,
		pert_mean   = qol_ae_pert,
		pert_min    = qol_aemin_pert,
		pert_max    = qol_aemax_pert,
		pert_lambda = qol_aelam_pert
	),

	# TB hospitalised QoL
	sim_tbho_qol = .beta_or_pert_samples(
		n_samples,
		dist        = tbhoqol_dist,
		beta_mean   = qol_tbho_beta,
		beta_sd     = qol_tbhosd_beta,
		pert_mean   = qol_tbho_pert,
		pert_min    = qol_tbhomin_pert,
		pert_max    = qol_tbhomax_pert,
		pert_lambda = qol_tbholam
	)
))

# Add CFR samples

cfr    <- qs_read(infile_p_cfr_ets)

# Sample from truncated normal (to avoid negatives)
rands <- rtruncnorm(parameters$n_samples, 1, cfr$st_err, low = 0, high = Inf)


samples$sim_mu_tbtx=cfr$ptb$y%*%t(rands)


samples$sim_mu_etbtx=cfr$eptb$y%*%t(rands)


# Add Prob PTB samples

prob_ptb  <- qs_read(infile_p_ptb_ets)

primer <- rnorm_ptb (parameters$n_samples,
					 prob_ptb$p_ptb$y[1],
					 prob_ptb$st_err)


samples$sim_p_ptb= prob_ptb$ratio%*%t(primer)



# -------------------------------------------------------------------------
# Save the samples --------------------------------------------------------
# -------------------------------------------------------------------------
qs_save(samples, outfile)
