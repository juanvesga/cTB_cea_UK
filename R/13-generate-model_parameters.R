# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
if (interactive()) {

	# Input
	root                 <- here::here()
	infile_parameters    <- file.path(root, "output", "parameters.qs2")
	infile_conmat        <- file.path(root, "output", "conmat.qs2")
	infile_dfage         <- file.path(root, "output", "dfage.qs2")
	infile_ptbld_age     <- file.path(root, "output", "ptbld_age.qs2")
	#infile_p_ptb_ets     <- file.path(root, "output", "p_ptb_ets.qs2")
	infile_AE_rates_NICE <- file.path(root, "output", "AE_rates_NICE.qs2")
	infile_qale          <- file.path(root, "output", "qale.qs2")
	infile_age_rates     <- file.path(root, "output", "age_rates.qs2")
#	infile_p_cfr_ets     <- file.path(root, "output", "p_cfr_ets.qs2")
	infile_samples       <- file.path(root, "output", "samples.qs2")
	outfile              <- file.path(root, "output", "model_parameters.qs2")

	# Packages
	source(file.path(root, "R", "modify_attach.R"))
	modify_attach(data.table, include.only = "as.data.table")
	modify_attach(qs2, include.only = c("qs_read", "qs_save"))

} else {

	# Input
	args                 <- commandArgs(trailingOnly = TRUE)
	infile_parameters    <- args[1]
	infile_conmat        <- args[2]
	infile_dfage         <- args[3]
	infile_ptbld_age     <- args[4]
	infile_AE_rates_NICE <- args[5]
	infile_qale          <- args[6]
	infile_age_rates     <- args[7]
	infile_samples       <- args[8]
	outfile              <- args[9]

	# Packages
	library(data.table, include.only = "as.data.table")
	library(qs2, include.only = c("qs_read", "qs_save"))
}


# -------------------------------------------------------------------------
# Load parameters ---------------------------------------------------------
# -------------------------------------------------------------------------
parameters          <- qs_read(infile_parameters)

COHORT_SIZE         <- parameters$cohort_size
ELIGIBILITY         <- parameters$age_eligibility
TEST                <- parameters$test
TEST_2ND            <- parameters$test_2nd
TST_RETURN          <- parameters$tst_return
TPT_START           <- parameters$tpt_start
TPT_COMPLETION      <- parameters$tpt_completion
TPT_LFUP_EFF        <- parameters$tpt_lfup_eff
#TRUE_PREV           <- parameters$true_prev
#P_PTBLD             <- parameters$p_ptbld
TRUE_POS_TSPO       <- parameters$true_pos_tspo
TRUE_NEG_TSPO       <- parameters$true_neg_tspo
TRUE_POS_QFT        <- parameters$true_pos_qft
TRUE_NEG_QFT        <- parameters$true_neg_qft
TRUE_POS_TST5       <- parameters$true_pos_tst5
TRUE_NEG_TST5       <- parameters$true_neg_tst5
TRUE_POS_CTB        <- parameters$true_pos_ctb
TRUE_NEG_CTB        <- parameters$true_neg_ctb
TPT                 <- parameters$tpt
TPT_DATA            <- parameters$tpt_data
YR_FAST             <- parameters$yr_fast
STEP_DECLINE        <- parameters$step_decline
DT                  <- parameters$dt
#R0                  <- parameters$r0
R_REACTIVATION_SLOW <- parameters$r_reactivation_slow
R_REACTIVATION      <- parameters$r_reactivation
R_REACTIVATION_FN_F <- parameters$r_reactivation_fn_f
R_REACTIVATION_FN_S <- parameters$r_reactivation_fn_s
R_SLOW              <- parameters$r_slow
#RR_PTBLD_MU         <- parameters$rr_ptbld_mu
#R_TBDUR             <- parameters$r_tbdur
#R_ETBDUR            <- parameters$r_etbdur
#R_TXDUR             <- parameters$r_txdur
#P_TBHOSP            <- parameters$p_tbhosp
P_DR                <- parameters$p_dr
T_END               <- parameters$t_end
T_HOR               <- parameters$t_hor
N_SAMPLES           <- parameters$n_samples
PTBLD_SWITCH        <- parameters$ptbld_switch

# -------------------------------------------------------------------------
# Load other data ---------------------------------------------------------
# -------------------------------------------------------------------------
con_mat   <- qs_read(infile_conmat) # contact matrix for distributing secondary TB cases
dfage     <- qs_read(infile_dfage)
ptbld_age <- qs_read(infile_ptbld_age)
#prob_ptb  <- qs_read(infile_p_ptb_ets)
prob_ae   <- qs_read(infile_AE_rates_NICE)
qale      <- qs_read(infile_qale)
df_rates  <- qs_read(infile_age_rates)
#cfr_tx    <- qs_read(infile_p_cfr_ets)
samples   <- qs_read(infile_samples)

# -------------------------------------------------------------------------
# stuff
# -------------------------------------------------------------------------

# Get Age eligibility --------------------------------------------------------
age_elig <- numeric(length = 85)
a0 <- 15
age_elig[(16 - a0):(100 - a0)] <- ELIGIBILITY[["a16_100"]]
age_elig[(16 - a0):( 65 - a0)] <- ELIGIBILITY[["a16_65"]]
age_elig[(16 - a0):( 45 - a0)] <- ELIGIBILITY[["a16_45"]]
age_elig[(16 - a0):( 35 - a0)] <- ELIGIBILITY[["a16_35"]]

# Set TBI positivity full Method ---------------------------------------------------------

# init_testpos<-tbipos$n_testpos*age_elig
# init_testneg<-tbipos$n_testneg*age_elig
#
#
# # Re-scale to get 10k of the elected age groups
# pos_f<- init_testpos/sum(init_testpos+init_testneg)
# neg_f<- init_testneg/sum(init_testpos+init_testneg)
#
# init_testpos<-pos_f*par$cohort_size
# init_testneg<-neg_f*par$cohort_size
#

# Set positivity test all 25yrold and predict -----------------------------


# Load PostTB data --------------------------------------------------------------
prev_ptbld <- ptbld_age$prev

# load ETS data on proportions developing eTB and pTB --------------------------------------------------------------
# ptb <- approx(
# 	x = prob_ptb$mid_age,
# 	y = prob_ptb$est,
# 	xout = seq(16,100),
# 	method = "linear"
# )

# Load adverse events data --------------------------------------------------------------
ae_6h <- approx(
	x = prob_ae$age,
	y = prob_ae$ae6H,
	xout = seq(16,100),
	method = "linear"
)

ae_3hr <- approx(
	x=prob_ae$age,
	y =prob_ae$ae3HR,
	xout = seq(16,100),
	method = "linear"
)

ae_3hp <- approx(
	x = prob_ae$age,
	y = prob_ae$ae3HP,
	xout = seq(16,100),
	method = "linear"
)

# Set TPT efficacy according to Regimen selection-----------------------------
#tpt_eff <- with(TPT_DATA, eff[regimen == TPT])
tpt_ae <- switch(TPT, "6H" = ae_6h, "3HP" = ae_3hp, "3HR" = ae_3hr, stop())


# Briggs QALY approach (not in use)-------------------------------------------
dQALY <- qale$dQALY


# Load CFR for TB, Tb treatment ----------------------------------------------
# dptb <- subset(cfr_tx, tb_type2 == "pulmonary")
# detb <- subset(cfr_tx, tb_type2 == "extrapulmonary")
#
# cfr_eptbtx <- approx(
# 	x = detb$mid_age,
# 	y = detb$prop,
# 	xout = seq(16,100),
# 	method = "linear"
# )
#
# cfr_ptbtx <- approx(
# 	x = dptb$mid_age,
# 	y = dptb$prop,
# 	xout = seq(16,100),
# 	method = "linear"
# )
#
# #Smooth CFR
# steps <- (unique(df_rates$cfr))
# y <- c(0,steps,steps[3])
# x<-c(
# 	1,
# 	floor( sum( df_rates$cfr == steps[1L] ) / 2 ),
# 	floor( median( which(df_rates$cfr == steps[2L]) ) ),
# 	floor( median( which(df_rates$cfr == steps[3L]) ) ),
# 	max( which( df_rates$cfr == steps[3L] ) )
# )
#
# smooth_cfr <- approx(x, y, xout = seq_along(df_rates$cfr))


# Get competing hazard with CFR

#Currently those set to 0 because we assume no TB deaths before treatment
# mutb <- smooth_cfr$y *0#par$r_tbdur*df_rates$cfr/(1-df_rates$cfr)
#
# muetb <- smooth_cfr$y *0#par$r_etbdur * cfr_eptbtx$y /(1- cfr_eptbtx$y)
#
# mutbtx <- cfr_ptbtx$y #par$r_txdur*cfr_ptbtx$y/(1-cfr_ptbtx$y)
#
# muetbtx <- cfr_eptbtx$y #par$r_txdur*cfr_eptbtx$y/(1-cfr_eptbtx$y)


# creates a step function to reduce progreesion overtiime (currently not used)------------
shape_f <- numeric(length = 85) + 1
stopifnot(length(YR_FAST) == 1L, YR_FAST <= 85, YR_FAST >= 1)
shape_f[YR_FAST:85] <- STEP_DECLINE



## Output list
model_pars <- list(
	dt                  = DT,
	cohort_size         = COHORT_SIZE,
	contact_matrix      = con_mat,
	#r0                  = R0,
	N_age               = 85, # TODO - This should come in via parameters along with some checks
	mu                  = df_rates$mu,
	mu_tb               = df_rates$mu*0,# df_rates$cfr,
	#mu_tbtx             = mutbtx,#cfr_ptbtx$y,
	#mu_etbtx            = muetbtx,#cfr_eptbtx$y,
	#p_ptb               = ptb$y,
	qol                 = df_rates$qol,
	dQALY               = dQALY,
	#init_testpos        = init_testpos,
	#init_testneg        = init_testneg,
	tpt_start           = TPT_START,
	tpt_completion      = TPT_COMPLETION,
	age_elig            = age_elig,
	#r_tpteff            = tpt_eff,
	tpt_ae               = tpt_ae$y,
	age_rate            = 1,
	r_reactivation_slow = R_REACTIVATION_SLOW,
	r_reactivation      = R_REACTIVATION,
	r_reactivation_fn_f = R_REACTIVATION_FN_F,
	r_reactivation_fn_s = R_REACTIVATION_FN_S,
	r_slow              = R_SLOW,
	prev_ptbld          = prev_ptbld,
	#rr_ptbld_mu         = RR_PTBLD_MU,
	#r_tbdur             = R_TBDUR,
	#r_etbdur            = R_ETBDUR,
	#r_txdur             = R_TXDUR,
	#p_tbhosp            = P_TBHOSP,
	disc_rate           = P_DR,
	rate_shape          = shape_f,
	n_shape_steps       = 85,
	true_pos_tspo       = TRUE_POS_TSPO,
	true_neg_tspo       = TRUE_NEG_TSPO,
	true_pos_qft        = TRUE_POS_QFT,
	true_neg_qft        = TRUE_NEG_QFT,
	true_pos_tst5        = TRUE_POS_TST5,
	true_neg_tst5        = TRUE_NEG_TST5,
	true_pos_ctb        = TRUE_POS_CTB,
	true_neg_ctb        = TRUE_NEG_CTB,
	tst_return          = TST_RETURN,
	test_2nd            = TEST_2ND,
	tpt_lfup_eff        = TPT_LFUP_EFF,
	t_end               = T_END,
	t_hor               = T_HOR,
	n_samples           = N_SAMPLES,
	dfage_value         =dfage$value,
	ptbld_switch        = PTBLD_SWITCH
)


res <- list(model_pars = model_pars, samples = samples)

# -------------------------------------------------------------------------
# Save list ---------------------------------------------------------------
# -------------------------------------------------------------------------
qs_save(res, outfile)
