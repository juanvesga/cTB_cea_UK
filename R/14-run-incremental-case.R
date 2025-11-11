# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------

root <- here::here()
infile_model_parameters <- file.path(root, "output", "model_parameters.qs2")
model_parameters    <- qs_read(infile_model_parameters)

# cascade params
uptake_val <- 0.45; increm_val <- 0.01
# unit price
price_vals <- seq(10,26,by=4)
# C-Tb sensit/specif: TRUE positive/negative rates
pars_true_pos_ctb <- model_parameters$model_pars$true_pos_ctb #
pars_true_neg_ctb <- model_parameters$model_pars$true_neg_ctb #

# differential return rate for negative result patients. By default set to 1
neg_test_diff_return <- 1
val_sim_test2_neg_cost_ctb <- model_parameters$samples$sim_test2_neg_cost_ctb*neg_test_diff_return

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  # create folder names
  # string abt unit price
  price_str <- ifelse(length(price_vals)>1,
            paste0("_pricemin",
                min(price_vals),"max",max(price_vals),"incr",unique(diff(price_vals))),
            paste0("_price",unique(price_vals)) )
  # string abt test performance
  test_perf_str <- ifelse(pars_true_pos_ctb>model_parameters$model_pars$true_pos_qft,                        "_ctb_high_sensit", "")
  print(test_perf_str)
  
  # create string if negative test return rate is different
  neg_test_ret_str <- ifelse(neg_test_diff_return==1,
              "",paste0("_negretrate",neg_test_diff_return))
  
  # subfolder name together
  subfolder_name <- with( list(uptake_val,increm_val,price_vals), 
            paste0("uptakeval",uptake_val,
              "incr",increm_val, price_str, 
              test_perf_str, neg_test_ret_str))
    # format(Sys.time(), "%Y_%m_%d_%H_%M")

### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ### ###
  
  # Input
  infile_pars_discrete2   <- file.path(root, "output", "pars_cohort_discrete2_new.qs2")
  infile_model_0          <- file.path(root, "models", "cohort_discrete0_new.R")
  infile_model_1          <- file.path(root, "models", "cohort_discrete1_new.R")
  infile_model_2          <- file.path(root, "models", "cohort_discrete2_new.R")
  infile_model_3          <- file.path(root, "models", "cohort_discrete3_new.R")
  out_folder <- file.path(root, "results", subfolder_name )
  dir.create(out_folder)
  
  # Packages
  source(file.path(root, "R", "modify_attach.R"))
  modify_attach(qs2,         include.only = c("qs_read", "qs_save"))
  modify_attach(matrixStats, include.only = c("colCumsums", "rowCumsums"))
  


# -------------------------------------------------------------------------
# Helper functions --------------------------------------------------------
# -------------------------------------------------------------------------



# Simulate Main function --------------------------------------------------

simulate_main <- function(generator, pars, samples, 
                  modeltype, price_range=price_vals, # seq(8,24,2),
                  overall_uptake=uptake_val      ) {
  
  
  t_hor <- pars$t_hor
  t_end <- pars$t_end
  igra <-  c("C-TB")
  ptbld_vals <- c(1)
  c_ctb <- price_range # seq(8,24,2)
  
  
  # Create memory to allocate values
  # TODO - check these dims
  resmat <- list()
  icer_sims<-matrix(0,nrow=length(igra)*length(ptbld_vals)*length(c_ctb),ncol=pars$n_samples)
  incQALY_sims<-matrix(0,nrow=length(igra)*length(ptbld_vals)*length(c_ctb),ncol=pars$n_samples)
  incCost_sims<-matrix(0,nrow=length(igra)*length(ptbld_vals)*length(c_ctb),ncol=pars$n_samples)
  
  # load parameters from calibration for Model B
  fitted <- qs_read(infile_pars_discrete2)
  
  pars$r_reactivation      <- with(fitted, x[X == "r_reactivation"])
  pars$r_reactivation_fn_f <- with(fitted, x[X == "r_reactivation_fn_f"])
  pars$r_reactivation_fn_s <- with(fitted, x[X == "r_reactivation_fn_s"])
  pars$r_reactivation_slow <- with(fitted, x[X == "r_reactivation_slow"])
  pars$r_slow              <- with(fitted, x[X == "r_slow"])
  
  
  
  # Test level --------------------------------------------------------------
  pars0<-pars
  samples0<-samples
  counter <- 0
  for (ig in igra){
    
    
    
    for (ic in c_ctb){
      
      
      
      
    # overall_uptake <- 0.64
    
    pars$test <- ig  # Set test
    pars0$test <- "QuantiFERON"  # Set test
    samples0$sim_test_cost <- samples$sim_test_cost_qfn
    samples0$sim_test2_neg_cost <- samples$sim_test2_neg_cost_qfn
    samples0$sim_test2_pos_cost <- samples$sim_test2_pos_cost_qfn
    pars0$true_pos <- pars$true_pos_qft
    pars0$true_neg <- pars$true_neg_qft
    pars0$p_2nd  <- sqrt(overall_uptake)
    pars0$tpt_start<- sqrt(overall_uptake)
    
    # true positives and negatives according to test-------------------
    if (ig == "T-SPOT.TB"){
      samples$sim_test_cost <- samples$sim_test_cost_tspo
      samples$sim_test2_neg_cost <- samples$sim_test2_neg_cost_tspo
      samples$sim_test2_pos_cost <- samples$sim_test2_pos_cost_tspo
      pars$true_pos <- pars$true_pos_tspo
      pars$true_neg <- pars$true_neg_tspo
      pars$p_2nd  <- pars$test_2nd
    } else if (ig == "QuantiFERON") {
      samples$sim_test_cost <- samples$sim_test_cost_qfn
      samples$sim_test2_neg_cost <- samples$sim_test2_neg_cost_qfn
      samples$sim_test2_pos_cost <- samples$sim_test2_pos_cost_qfn
      pars$true_pos <- pars$true_pos_qft
      pars$true_neg <- pars$true_neg_qft
      pars$p_2nd  <- pars$test_2nd
      
    } else if (ig == "Tuberculin Skin Test (5mm)") { 
      samples$sim_test_cost <- samples$sim_test_cost_tst
      samples$sim_test2_neg_cost <- samples$sim_test2_neg_cost_tst
      samples$sim_test2_pos_cost <- samples$sim_test2_pos_cost_tst
      pars$true_pos <- pars$true_pos_tst5
      pars$true_neg <- pars$true_neg_tst5
      pars$p_2nd  <- pars$tst_return
    } else if (ig == "C-TB") {
        
        samples$sim_test_cost      <- samples$sim_test_staff_cost_ctb + ic
      samples$sim_test2_neg_cost     <- val_sim_test2_neg_cost_ctb # samples$sim_test2_neg_cost_ctb
      samples$sim_test2_pos_cost     <- samples$sim_test2_pos_cost_ctb
      pars$true_pos <- pars_true_pos_ctb
      pars$true_neg <- pars_true_neg_ctb
      pars$p_2nd  <- pars$tst_return
    }  else {
      stop("Currently we are only considering T-SPOT and QuantiFERON.")
    }
    
    # define cascade -------------------------------------------------------------
    pars0$cascade <-
      pars0$p_2nd * # minus attrition those not attending 2nd test visit
      pars0$tpt_start * # fraction positive starting tpt
      (pars0$tpt_completion + (1 - pars0$tpt_completion) * pars0$tpt_lfup_eff)  # GFraction starting completig tpt
    
    
    pars$cascade <-
      pars$p_2nd * # minus attrition those not attending 2nd test visit
      pars$tpt_start * # fraction positive starting tpt
      (pars$tpt_completion + (1 - pars$tpt_completion) * pars$tpt_lfup_eff)  # GFraction starting completig tpt
    
    
    # PTBLD Level ---------------------------------------------------------
    
    
    for (ld in seq_along(ptbld_vals)){

      p0 <- pars0      
      p2 <- pars
      

      p0$ptbld_switch <- ptbld_vals[ld]
      p2$ptbld_switch <- ptbld_vals[ld]
      
      # run model
      runs <- simulate_cohort(generator,p0, p2, samples0,samples)
      
      full_pars<-runs$input
      
      
      # Process results and get ICER estimate
      icer <- get_icer(full_pars, runs, 0)
      
      # Fill overall results in premade matrix
      counter<-counter+1
      id<-counter
      
      mat       <- icer$out_table
      mat$igra  <- ig
      mat$ptbld <- ptbld_vals[ld]
        mat$c_ctb <- ic
      resmat[[counter]]<-mat
      
      # ICERS , QALy and costs simulations for CEAC
      icer_sims[id,]<-icer$icer[t_hor,]
      incQALY_sims[id,]<-icer$delta_qaly[t_hor,]
      incCost_sims[id,]<-icer$delta_cost[t_hor,]
      
    }
    
    
    }
    
    
    
    
  }
  res <- do.call(rbind, resmat)
  list(tab = res,
       icer_sims=icer_sims,
       incQALY_sims=incQALY_sims,
       incCost_sims=incCost_sims)
}


# Simulate cohort function ------------------------------------------------


simulate_cohort <- function(generator, 
                    input0,input, 
                    samples0,samples) {
  
  t_end <- input$t_end
  t     <- seq(1, t_end, length.out = t_end)
  n     <- length(samples$sim_tb_qol)
  
  # Create memory to store results (0 for baseline, 1 for intervention)
  arr <- array(0, c(n, length(t), 85))
  mat <- matrix(0, nrow = length(t), ncol = n)
  
  pop_age0         <- arr
  N_tot0           <- mat
  new_TBcases0     <- mat
  new_secTBcases0  <- mat
  new_eTBcases0    <- mat
  new_TBcases_neg0 <- mat
  new_TBcases_pos0 <- mat
  new_TBdeaths0    <- mat
  new_PTBLDdeaths0 <- mat
  new_deaths0      <- mat
  new_tpt0         <- mat
  qaly0            <- mat
  qaly_LE0         <- mat
  qaly_TB0         <- mat
  qaly_TBtx0       <- mat
  qaly_TBhosp0     <- mat
  qaly_postTB0     <- mat
  qaly_tptae0      <- mat
  qaly_pop0        <- mat
  qaly_all0        <- mat
  costs0           <- mat
  cost_test0       <- mat
  cost_tpt0        <- mat
  cost_tptae0      <- mat
  cost_tb0         <- mat
  cost_postTB0     <- mat
  shape0           <- mat
  
  pop_age1         <- arr
  N_tot1           <- mat
  new_TBcases1     <- mat
  new_secTBcases1  <- mat
  new_eTBcases1    <- mat
  new_TBcases_neg1 <- mat
  new_TBcases_pos1 <- mat
  new_TBdeaths1    <- mat
  new_PTBLDdeaths1 <- mat
  new_deaths1      <- mat
  new_tpt1         <- mat
  qaly1            <- mat
  qaly_LE1         <- mat
  qaly_TB1         <- mat
  qaly_TBtx1       <- mat
  qaly_TBhosp1     <- mat
  qaly_postTB1     <- mat
  qaly_tptae1      <- mat
  qaly_pop1        <- mat
  qaly_all1        <- mat
  costs1           <- mat
  cost_test1       <- mat
  cost_tpt1        <- mat
  cost_tptae1      <- mat
  cost_tb1         <- mat
  cost_postTB1     <- mat
  shape1           <- mat
  
  res <- list()
  
  # If we have samples > 1 we need this loop too load specific paramaters for each iteration
  for (jj in 1:n) {
    input0$qol_ptb_weight   <- samples0$sim_tb_qol[jj]
    input0$qol_eptb_weight  <- samples0$sim_eptb_qol[jj]
    input0$qol_ptbld_weight <- samples0$sim_post_qol[jj]
    input0$qol_tx_weight    <- samples0$sim_tx_qol[jj]
    input0$qol_tptae_loss   <- samples0$sim_ae_qol[jj]
    input0$qol_tbhosp_loss  <- samples0$sim_tbho_qol[jj]
    input0$c_tptae          <- samples0$sim_adverse_cost[jj]
    input0$c_test           <- samples0$sim_test_cost[jj]
    input0$c_test2_pos      <- samples0$sim_test2_pos_cost[jj]
    input0$c_test2_neg      <- samples0$sim_test2_neg_cost[jj]
    input0$c_posi_test      <- samples0$sim_positive_cost[jj]
    input0$c_tpt            <- samples0$sim_tpt_cost[jj]
    input0$c_tbdx_ptb       <- samples0$sim_tbdx_ptb_cost[jj]
    input0$c_tbdx_eptb      <- samples0$sim_tbdx_eptb_cost[jj]
    input0$c_tbtx           <- samples0$sim_tbtx_cost[jj]
    input0$c_tbhosp         <- samples0$sim_tbho_cost[jj]
    input0$c_post           <- samples0$sim_post_cost[jj]
    
    
    input$qol_ptb_weight   <- samples$sim_tb_qol[jj]
    input$qol_eptb_weight  <- samples$sim_eptb_qol[jj]
    input$qol_ptbld_weight <- samples$sim_post_qol[jj]
    input$qol_tx_weight    <- samples$sim_tx_qol[jj]
    input$qol_tptae_loss   <- samples$sim_ae_qol[jj]
    input$qol_tbhosp_loss  <- samples$sim_tbho_qol[jj]
    input$c_tptae          <- samples$sim_adverse_cost[jj]
    input$c_test           <- samples$sim_test_cost[jj]
    input$c_test2_pos      <- samples$sim_test2_pos_cost[jj]
    input$c_test2_neg      <- samples$sim_test2_neg_cost[jj]
    input$c_posi_test      <- samples$sim_positive_cost[jj]
    input$c_tpt            <- samples$sim_tpt_cost[jj]
    input$c_tbdx_ptb       <- samples$sim_tbdx_ptb_cost[jj]
    input$c_tbdx_eptb      <- samples$sim_tbdx_eptb_cost[jj]
    input$c_tbtx           <- samples$sim_tbtx_cost[jj]
    input$c_tbhosp         <- samples$sim_tbho_cost[jj]
    input$c_post           <- samples$sim_post_cost[jj]
    
    
    # Pass sampled Epi parameters
    input0$init_testpos   <- round(input0$cohort_size * input0$dfage_value * samples0$sim_true_prev[jj])
    input0$init_testneg   <- (input0$cohort_size * input0$dfage_value) - input0$init_testpos
    input0$r_tpteff       <- samples0$sim_tpt_eff[jj]
    input0$r_tptae        <- input0$tpt_ae * samples0$sim_tpt_ae[jj]
    input0$r_tbdur        <- 1 / (samples0$sim_r_tbdur[jj]/ 365)  # TB duration (ETS)
    input0$r_etbdur       <- 1 / (samples0$sim_r_etbdur[jj]/ 365) # eTB duration (ETS)
    input0$r_txdur        <- 1 / (samples0$sim_r_txdur[jj]/ 365)  # TX duration (ETS)
    input0$p_tbhosp       <- samples0$sim_p_tbhosp[jj]
    input0$p_ptbld        <- input0$prev_ptbld * samples0$sim_p_ptbld[jj] * input0$ptbld_switch
    input0$rr_ptbld_mu    <- samples0$sim_rr_ptbld_mu[jj]
    input0$r0             <- samples0$sim_r0[jj]
    input0$mu_tbtx        <- samples0$sim_mu_tbtx[,jj]
    input0$mu_etbtx       <- samples0$sim_mu_etbtx[,jj]
    input0$p_ptb          <- samples0$sim_p_ptb[,jj]
    
    
    input$init_testpos   <- round(input$cohort_size * input$dfage_value * samples$sim_true_prev[jj])
    input$init_testneg   <- (input$cohort_size * input$dfage_value) - input$init_testpos
    input$r_tpteff       <- samples$sim_tpt_eff[jj]
    input$r_tptae        <- input$tpt_ae * samples$sim_tpt_ae[jj]
    input$r_tbdur        <- 1 / (samples$sim_r_tbdur[jj]/ 365)  # TB duration (ETS)
    input$r_etbdur       <- 1 / (samples$sim_r_etbdur[jj]/ 365) # eTB duration (ETS)
    input$r_txdur        <- 1 / (samples$sim_r_txdur[jj]/ 365)  # TX duration (ETS)
    input$p_tbhosp       <- samples$sim_p_tbhosp[jj]
    input$p_ptbld        <- input$prev_ptbld * samples$sim_p_ptbld[jj] * input$ptbld_switch
    input$rr_ptbld_mu    <- samples$sim_rr_ptbld_mu[jj]
    input$r0             <- samples$sim_r0[jj]
    input$mu_tbtx        <- samples$sim_mu_tbtx[,jj]
    input$mu_etbtx       <- samples$sim_mu_etbtx[,jj]
    input$p_ptb          <- samples$sim_p_ptb[,jj]
    
    
    full_input<-input
    full_input0<-input0
    
    # Remove unneeded values to quieten odin
    nms <- c(
      "true_pos_tst5",
      "true_neg_tst5",
      "true_pos_ctb",
      "true_neg_ctb",
      "true_pos_tspo",
      "true_neg_tspo",
      "true_pos_qft",
      "true_neg_qft",
      "tst_return",
      "test_2nd",
      "tpt_lfup_eff",
      "t_end",
      "t_hor",
      "n_samples",
      "test",
      "init_testpos_y1",
      "init_testpos_y2",
      "init_testpos_y3",
      "init_testpos_y10",
      "cohort_size",
      "tpt_ae",
      "prev_ptbld",
      "dfage_value",
      "ptbld_switch"
    )
    
    # Create parameters object for baseline
    
    
    input1 <- input
    input0[nms] <- NULL
    input1[nms] <- NULL
    
    
    # generate baseline model
    mod0 <- cohort_generator$new(user = input0)
    
    # Run baseline model
    y0 <- mod0$run(t) #  baseline
    
    # generate intervention model
    mod1 <- cohort_generator$new(user = input1)
    
    # Run intervention model
    y1 <- mod1$run(t) # intervention
    
    # Store results
    res$y0[[jj]] <- y0
    res$y1[[jj]] <- y1
    
    for (a in 1:input$N_age) {
      col_txt <- paste0("N_tot_age[", a, "]")
      pop_age0[jj, , a] <- y0[, col_txt]
      pop_age1[jj, , a] <- y1[, col_txt]
    }
    
    N_tot0[, jj]           <- y0[, "N_tot"]
    new_TBcases0[, jj]     <- y0[, "new_TBcases"]
    new_secTBcases0[, jj]  <- y0[, "new_secondaryTB"]
    new_eTBcases0[, jj]    <- y0[, "new_eTBcases"]
    new_TBdeaths0[, jj]    <- y0[, "new_TBdeaths"]
    new_TBcases_pos0[, jj] <- y0[, "new_allTBcases_pos"]
    new_TBcases_neg0[, jj] <- y0[, "new_allTBcases_neg"]
    new_PTBLDdeaths0[, jj] <- y0[, "new_PTBLDdeaths"]
    new_deaths0[, jj]      <- y0[, "new_deaths"]
    new_tpt0[, jj]         <- y0[, "new_tpt"]
    qaly_LE0[, jj]         <- y0[, "qaly_LE_loss"]
    qaly_TB0[, jj]         <- y0[, "qaly_TB"]
    qaly_TBtx0[, jj]       <- y0[, "qaly_TBtx"]
    qaly_TBhosp0[, jj]     <- y0[, "qaly_TBhosp_loss"]
    qaly_postTB0[, jj]     <- y0[, "qaly_postTB"]
    qaly_tptae0[, jj]      <- y0[, "qaly_tptae_loss"]
    qaly_pop0[, jj]        <- y0[, "qaly_pop"]
    qaly_all0[, jj]        <- y0[, "qaly_all"]
    costs0[, jj]           <- y0[, "costs"]
    cost_test0[, jj]       <- y0[, "cost_test"]
    cost_tpt0[, jj]        <- y0[, "cost_tpt"]
    cost_tptae0[, jj]      <- y0[, "cost_tptae"]
    cost_tb0[, jj]         <- y0[, "cost_tb"]
    cost_postTB0[, jj]     <- y0[, "cost_post"]
    
    N_tot1[, jj]           <- y1[, "N_tot"]
    new_TBcases1[, jj]     <- y1[, "new_TBcases"]
    new_secTBcases1[, jj]  <- y1[, "new_secondaryTB"]
    new_eTBcases1[, jj]    <- y1[, "new_eTBcases"]
    new_TBdeaths1[, jj]    <- y1[, "new_TBdeaths"]
    new_TBcases_pos1[, jj] <- y1[, "new_allTBcases_pos"]
    new_TBcases_neg1[, jj] <- y1[, "new_allTBcases_neg"]
    new_PTBLDdeaths1[, jj] <- y1[, "new_PTBLDdeaths"]
    new_deaths1[, jj]      <- y1[, "new_deaths"]
    new_tpt1[, jj]         <- y1[, "new_tpt"]
    qaly_LE1[, jj]         <- y1[, "qaly_LE_loss"]
    qaly_TB1[, jj]         <- y1[, "qaly_TB"]
    qaly_TBtx1[, jj]       <- y1[, "qaly_TBtx"]
    qaly_TBhosp1[, jj]     <- y1[, "qaly_TBhosp_loss"]
    qaly_postTB1[, jj]     <- y1[, "qaly_postTB"]
    qaly_tptae1[, jj]      <- y1[, "qaly_tptae_loss"]
    qaly_pop1[, jj]        <- y1[, "qaly_pop"]
    qaly_all1[, jj]        <- y1[, "qaly_all"]
    costs1[, jj]           <- y1[, "costs"]
    cost_test1[, jj]       <- y1[, "cost_test"]
    cost_tpt1[, jj]        <- y1[, "cost_tpt"]
    cost_tptae1[, jj]      <- y1[, "cost_tptae"]
    cost_tb1[, jj]         <- y1[, "cost_tb"]
    cost_postTB1[, jj]     <- y1[, "cost_post"]
  }
  
  testpos <- sum(input$model_pars$init_testpos)
  testneg <- sum(input$model_pars$init_testneg)
  
  base <- list(
    pop_age = pop_age0,
    N_tot = N_tot0,
    new_TBcases = new_TBcases0,
    new_secTBcases = new_secTBcases0,
    new_eTBcases = new_eTBcases0,
    new_TBdeaths = new_TBdeaths0,
    new_TBcases_pos = (new_TBcases_pos0),
    new_TBcases_neg = (new_TBcases_neg0),
    new_PTBLDdeaths = new_PTBLDdeaths0,
    new_deaths = new_deaths0,
    new_tpt = new_tpt0,
    qaly_LE = qaly_LE0,
    qaly_TB = qaly_TB0,
    qaly_TBtx = qaly_TBtx0,
    qaly_TBhosp = qaly_TBhosp0,
    qaly_postTB = qaly_postTB0,
    qaly_tptae = qaly_tptae0,
    qaly_pop = qaly_pop0,
    qaly_all = qaly_all0,
    costs = costs0,
    cost_test = cost_test0,
    cost_tpt = cost_tpt0,
    cost_tptae = cost_tptae0,
    cost_tb = cost_tb0,
    cost_post = cost_postTB0,
    shape = shape0
  )
  
  itv <- list(
    pop_age = pop_age1,
    N_tot = N_tot1,
    new_TBcases = new_TBcases1,
    new_secTBcases = new_secTBcases1,
    new_eTBcases = new_eTBcases1,
    new_TBdeaths = new_TBdeaths1,
    new_TBcases_pos = new_TBcases_pos1 ,
    new_TBcases_neg = new_TBcases_neg1,
    new_PTBLDdeaths = new_PTBLDdeaths1,
    new_deaths = new_deaths1,
    new_tpt = new_tpt1,
    qaly_LE = qaly_LE1,
    qaly_TB = qaly_TB1,
    qaly_TBtx = qaly_TBtx1,
    qaly_TBhosp = qaly_TBhosp1,
    qaly_postTB = qaly_postTB1,
    qaly_tptae = qaly_tptae1,
    qaly_pop = qaly_pop1,
    qaly_all = qaly_all1,
    costs = costs1,
    cost_test = cost_test1,
    cost_tpt = cost_tpt1,
    cost_tptae = cost_tptae1,
    cost_tb = cost_tb1,
    cost_post = cost_postTB1,
    shape = shape1
  )
  
  out <- list(
    y0 = res$y0,
    y1 = res$y1,
    base = base,
    itv = itv,
    input=full_input,
    input0=full_input0
  )
  
  return(out)
}


# Get ICER function -------------------------------------------------------


get_icer <- function(par, sims, plot_switch=1){
  
  
  
  t_hor <-par$t_hor
  tt <- par$t_end
  n  <- par$n_samples
  
  # Discount matrix (currently not in use because discounting occurs inside the model)
  i <- par$disc_rate
  v<- 1 / (1+i)
  discount_factors <- v^(0:tt)
  discount_factors <- discount_factors[1:tt]
  
  # Get basline and intervention objects
  base <- sims$base
  itv  <- sims$itv
  
  # ICER --------------------------------------------------------------------
  qaly_base  <- colCumsums(base$qaly_all)
  qaly_itv   <- colCumsums(itv$qaly_all)
  delta_qaly <- qaly_base - qaly_itv # Note we do baseline - intervention because we are using QALY losses not qaly gains
  
  # Briggs Method=(not in use)============================
  # Base QALYs
  
  #
  # baseline_qaly<- sum(par$dQALY*(par$init_testpos+par$init_testneg))
  #
  # tmp<-  base$qaly_pop #* discount_factors
  # #
  # #
  #  cum_base_qaly_vector<- colCumsums(tmp)
  # # #
  # # #qaly_base<-baseline_qaly-cum_base_qaly_vector
  # #
  #  qaly_base<-cum_base_qaly_vector
  # # #
  # # # # ITV QALYS
  # # #
  #  tmp<-itv$qaly_pop#* discount_factors
  # # #
  #  cum_itv_qaly_vector<- colCumsums(tmp)
  # # #
  # # # qaly_itv<-baseline_qaly-cum_itv_qaly_vector
  # qaly_itv<-cum_itv_qaly_vector
  #
  # Incremental QALYs
  #delta_qaly<- abs(qaly_itv-qaly_base)
  
  
  # Costs differences  ------------------------------------------------------
  costs_base <- colCumsums(base$costs)
  costs_itv  <- colCumsums(itv$costs)
  
  # Incremental costs
  delta_cost <- costs_itv - costs_base
  
  # ICER --------------------------------------------------------------------
  icer  <- delta_cost / delta_qaly
  
  # disaggregated costs and qalys -------------------------------------------
  
  # QALYs by categories
  # TODO - Shall we split this up. There's a lot going on here.
  delta_cost_cats <- c(
    postTB = mean(t(rowCumsums(t(itv$cost_post)))[tt,] - t(rowCumsums(t(base$cost_post)))[tt,]),
    
    TB = mean(t(rowCumsums(t(itv$cost_tb)))[tt,] - t(rowCumsums(t(base$cost_tb)))[tt,]),
    
    Test = mean(t(rowCumsums(t(itv$cost_test)))[tt,] - t(rowCumsums(t(base$cost_test)))[tt,]),
    
    TPT = mean(t(rowCumsums(t(itv$cost_tpt)))[tt,] - t(rowCumsums(t(base$cost_tpt)))[tt,]),
    
    tptAE = mean(t(rowCumsums(t(itv$cost_tptae)))[tt,] - t(rowCumsums(t(base$cost_tptae)))[tt,])
  )
  
  # costs by categories
  delta_qaly_cats<-c(
    LE = mean(t(rowCumsums(t(itv$qaly_LE)))[tt,]-t(rowCumsums(t(base$qaly_LE)))[tt,]),
    
    TB = mean(t(rowCumsums(t(itv$qaly_TB)))[tt,]-t(rowCumsums(t(base$qaly_TB)))[tt,]),
    
    TBhosp = mean(t(rowCumsums(t(itv$qaly_TBhosp)))[tt,]-t(rowCumsums(t(base$qaly_TBhosp)))[tt,]),
    
    TBtx = mean(t(rowCumsums(t(itv$qaly_TBtx)))[tt,]-t(rowCumsums(t(base$qaly_TBtx)))[tt,]),
    
    postTB = mean(t(rowCumsums(t(itv$qaly_postTB)))[tt,]-t(rowCumsums(t(base$qaly_postTB)))[tt,]),
    
    tptAE = mean(t(rowCumsums(t(itv$qaly_tptae)))[tt,]-t(rowCumsums(t(base$qaly_tptae)))[tt,])
  )
  
  ptbld_qaly_loss <- abs(
    mean(t(rowCumsums(t(itv$qaly_postTB)))[tt,] - t(rowCumsums(t(base$qaly_postTB)))[tt,])
  )
  
  
  
  
  # Other epi outputs -------------------------------------------------------
  
  pop<- drop(base$N_tot)
  
  
  totalTBcases <- mean(colSums(base$new_TBcases + base$new_eTBcases))
  
  totalsecondaryTBcases <- mean(colSums(base$new_secTBcases))
  
  totalTBdeaths <- mean(colSums(base$new_TBdeaths))
  
  PTBLDdeaths <- mean(colSums(base$new_PTBLDdeaths) - colSums(itv$new_PTBLDdeaths))
  
  deaths <- base$new_PTBLDdeaths
  
  tbdeaths <- mean(colSums(base$new_TBdeaths))
  
  alldeaths <- mean(colSums(base$new_deaths))
  
  
  TBcasesaverted<-mean(
    colSums(base$new_TBcases + base$new_eTBcases) - colSums(itv$new_TBcases + itv$new_eTBcases)
  )
  
  TBdeathsaverted <- mean(colSums(base$new_TBdeaths) - colSums(itv$new_TBdeaths))
  
  meanICER <- mean(icer[t_hor,])
  
  icer_range <- quantile(icer[t_hor,], probs = c(0.025, 0.975), names = FALSE)
  
  meanQALY <- mean(delta_qaly[t_hor,])
  
  meanCost <- mean(delta_cost[t_hor,])
  
  # When using > 1 samples find fraction of ICERs below 20k or 30k
  icer_f20k <- length(which(icer[t_hor,] < 20000)) / n
  
  icer_f30k <- length(which(icer[t_hor,] < 30000)) / n
  
  c_eff_bin<- if( (icer_f20k > 0.5) &&  (icer_f30k > 0.9)) 1 else 0
  
  
  # Results table -----------------------------------------------------------
  outp <- c(
    totalTBcases,
    TBcasesaverted,
    totalTBdeaths,
    TBdeathsaverted,
    PTBLDdeaths,
    meanCost,
    delta_cost_cats,
    mean(qaly_base[t_hor,]),
    mean(qaly_itv[t_hor,]),
    mean(delta_qaly[t_hor,]),
    delta_qaly_cats,
    meanICER,
    icer_range[1],
    icer_range[2],
    icer_f20k*1e2,
    icer_f30k*1e2
  )
  
  labs <- c(
    "Total TB cases baseline",
    "Total TB cases averted",
    "Total TB deaths baseline",
    "Total TB deaths averted",
    "PostTBLD deaths",
    "Incremental cost",
    "delta cost ptbld",
    "delta cost TB",
    "delta cost test",
    "delta cost tpt",
    "delta cost tpt AE",
    "QALY loss base",
    "QALY loss itv",
    "Incremental QALY",
    "delta qaly LE",
    "delta qaly TB",
    "delta qaly TBhosp",
    "delta qaly TBtx",
    "delta qaly ptbld",
    "delta qaly tpt AE",
    "ICER",
    "ICER_low",
    "ICER up",
    "Sims under 20k(%)",
    "Sims under 30k(%)"
  )
  
  names(outp) <- labs
  
  # Final list object
  res <- list(
    delta_qaly      = delta_qaly,
    delta_cost      = delta_cost,
    icer            = icer,
    delta_qaly_cats = delta_qaly_cats,
    delta_cost_cats = delta_cost_cats,
    out_table       = as.list(outp),
    c_eff_bin       = c_eff_bin,
    meanICER        = meanICER
  )
  
} # end of get_icer function

  

# -------------------------------------------------------------------------
# Load parameters ---------------------------------------------------------
# -------------------------------------------------------------------------
pars_discrete2      <- qs_read(infile_pars_discrete2)
nsamples            <- model_parameters$model_pars$n_samples

# -------------------------------------------------------------------------
# Other inputs ------------------------------------------------------------
# -------------------------------------------------------------------------

# Models
model_files <- c(
  infile_model_2        # Model B with suppressive effect of TPT
)


# Sample QALY and Costs
s <- model_parameters$samples

# Go through scenarios ----------------------------------------------------

# cascade combinations

tpt_start=seq(uptake_val,1,by=increm_val)
test_2nd =seq(uptake_val,1,by=increm_val)

casc_scens <- data.frame(
  tpt_start= rep( tpt_start, each = length(test_2nd)) ,
  test_2nd = rep( test_2nd , length(tpt_start))
)

results_list<-list()
icer_list<-list()
qaly_list<-list()
cost_list<-list()
inb_list<-list()

# Choose model type
for(i in seq_along(model_files)) {
  
  # Generate ODIN object ----------------------------------------------------
  model_f <- model_files[i]
  cohort_generator<-odin::odin(model_f)
  
  
  modeltype <- basename(tools::file_path_sans_ext(model_f))
  
  p <- model_parameters$model_pars
  main_tab_list <- vector("list", nrow(casc_scens)  )
  icer_tab_list <- vector("list", nrow(casc_scens)  )
  qaly_tab_list <- vector("list", nrow(casc_scens) )
  cost_tab_list <- vector("list", nrow(casc_scens)  )
  inb_tab_list  <- vector("list", nrow(casc_scens)  )
  
  ii<-0
  
  for(jj in 1:nrow(casc_scens)){
    ii <- ii+1
    
    if (ii %% 5 == 0) { print(paste0(ii,"/",nrow(casc_scens))) }
    
    p$tpt_start <- casc_scens$tpt_start[jj]
    p$tst_return <- casc_scens$test_2nd[jj]
    
    # Call main function
    runs_main <- simulate_main(cohort_generator, p, s, modeltype, 
                        price_range=price_vals, overall_uptake=uptake_val )
    
    
    #INB = (ΔE × K) − ΔC
    threshold <- 20000
    
    tab<-as.data.frame(runs_main$tab)
    icer<-as.data.frame(runs_main$icer_sims)
    qaly<-as.data.frame(runs_main$incQALY_sims)
    cost<-as.data.frame(runs_main$incCost_sims)
    inb <-as.data.frame( (runs_main$incQALY_sims*threshold) - runs_main$incCost_sims)
    tab$inb<-as.data.frame( (runs_main$incQALY_sims*threshold) - runs_main$incCost_sims)
    
    
    tab$start   <-toString(casc_scens$tpt_start[jj])
    tab$test_2nd<-toString(casc_scens$test_2nd[jj])
    
    tab$model_struct="B"
    if(i==1){
      tab$tpt_effect= "0 suppressive"} else {tab$tpt_effect="1 curative"}
    
    
    icer$start<-tab$start
    qaly$start<-tab$start
    cost$start<-tab$start
    inb$start<-tab$start
    
    icer$test_2nd<-tab$test_2nd
    qaly$test_2nd<-tab$test_2nd
    cost$test_2nd<-tab$test_2nd
    inb$test_2nd<-tab$test_2nd
    
    icer$model_struct<-tab$model_struct
    qaly$model_struct<-tab$model_struct
    cost$model_struct<-tab$model_struct
    inb$model_struct<-tab$model_struct
    
    icer$tpt_effect<-tab$tpt_effect
    qaly$tpt_effect<-tab$tpt_effect
    cost$tpt_effect<-tab$tpt_effect
    inb$tpt_effect<-tab$tpt_effect
    
    icer$igra<-tab$igra
    qaly$igra<-tab$igra
    cost$igra<-tab$igra
    inb$igra <-tab$igra
    
    icer$ptbld<-tab$ptbld
    qaly$ptbld<-tab$ptbld
    cost$ptbld<-tab$ptbld
    inb$ptbld <-tab$ptbld
    
    icer$c_ctb<-tab$c_ctb
    qaly$c_ctb<-tab$c_ctb
    cost$c_ctb<-tab$c_ctb
    inb$c_ctb <-tab$c_ctb
    
    
    main_tab_list[[ii]] <- tab
    icer_tab_list[[ii]] <- icer
    qaly_tab_list[[ii]] <- qaly
    cost_tab_list[[ii]] <- cost
    inb_tab_list[[ii]] <- inb
    
    
    
  }
  
  main_tab_mat <- do.call("rbind",main_tab_list)
  main_tab_mat[] <- lapply(  main_tab_mat, unlist)
  results_list[[i]]<-main_tab_mat
  
  icer_tab_mat <- do.call("rbind",icer_tab_list)
  icer_tab_mat[] <- lapply(  icer_tab_mat, unlist)
  icer_list[[i]]<-icer_tab_mat
  
  qaly_tab_mat <- do.call("rbind",qaly_tab_list)
  qaly_tab_mat[] <- lapply(  qaly_tab_mat, unlist)
  qaly_list[[i]]<-qaly_tab_mat
  
  cost_tab_mat <- do.call("rbind",cost_tab_list)
  cost_tab_mat[] <- lapply(  cost_tab_mat, unlist)
  cost_list[[i]]<-cost_tab_mat
  
  inb_tab_mat   <- do.call("rbind",inb_tab_list)
  inb_tab_mat[] <- lapply(  inb_tab_mat, unlist)
  inb_list[[i]] <- inb_tab_mat
  
  
  main_tab_mat <- do.call(rbind, main_tab_list)
  main_tab_mat[] <- lapply(main_tab_mat, unlist)
  filename <- paste0(modeltype, "_results.csv")
  write.csv(main_tab_mat, file.path(out_folder, filename), row.names = FALSE)
  print(paste("run model", i))
}

# -------------------------------------------------------------------------
# Save main results ---------------------------------------------------------------
# -------------------------------------------------------------------------

if(nsamples==1){
  
  # Main results
  results_mat<-do.call("rbind",results_list)
  results_mat[] <- lapply(results_mat, unlist)
  results_mat<-results_mat[with(results_mat, order(model_struct,tpt_effect, igra, ptbld,  test_2nd, start,c_ctb)), ]
  
  results_summary<-cbind(
    results_mat[,"model_struct"],
    results_mat[,"tpt_effect"],
    results_mat[,"igra"],
    results_mat[,"ptbld"],
    results_mat[,"test_2nd"],
    results_mat[,"start"],
    results_mat[,"c_ctb"],
    results_mat[,"Incremental cost"],
    results_mat[,"Incremental QALY"],
    results_mat[,"ICER"],
    results_mat[,"inb"])
  
  
  if(model_parameters$model_pars$t_hor==5){
    
    file_name<-paste0("all_results_single_5yr.csv")
    write.csv( results_mat,file.path(out_folder, file_name), row.names = FALSE)
    
    file_name<-paste0("summary_results_5yr.csv")
    write.csv( results_summary,file.path(out_folder, file_name), row.names = FALSE)
    
    
    
  }else{
    
    
    file_name<-paste0("all_results_single_incremental_case.csv")
    write.csv( results_mat,file.path(out_folder, file_name), row.names = FALSE)
    
    file_name<-paste0("summary_results_incremental_case.csv")
    write.csv( results_summary,file.path(out_folder, file_name), row.names = FALSE)
  }
  
  
  
  
}else	{
  
  # Main results
  results_mat<-do.call("rbind",results_list)
  results_mat[] <- lapply(results_mat, unlist)
  results_mat<-results_mat[with(results_mat, order(model_struct,tpt_effect, igra, ptbld,  test_2nd, start)), ]
  file_name<-paste0("all_results_incremental_case.csv")
  write.csv( results_mat,file.path(out_folder, file_name), row.names = FALSE)
  
  
  # ICER all sims results
  icer_mat<-do.call("rbind",icer_list)
  icer_mat[] <- lapply(icer_mat, unlist)
  icer_mat<-icer_mat[with(icer_mat, order(model_struct,tpt_effect, igra, ptbld,  test_2nd, start)), ]
  file_name<-paste0("all_icer_incremental_case.csv")
  write.csv( icer_mat,file.path(out_folder, file_name), row.names = FALSE)
  
  # ICER all sims results
  inb_mat<-do.call("rbind",inb_list)
  inb_mat[] <- lapply(inb_mat, unlist)
  inb_mat<-icer_mat[with(inb_mat, order(model_struct,tpt_effect, igra, ptbld,  test_2nd, start)), ]
  file_name<-paste0("all_inb_incremental_case.csv")
  write.csv( icer_mat,file.path(out_folder, file_name), row.names = FALSE)
  
  
  # Cost all sims cost
  cost_mat<-do.call("rbind",cost_list)
  cost_mat[] <- lapply(cost_mat, unlist)
  cost_mat<-cost_mat[with(cost_mat, order(model_struct,tpt_effect, igra, ptbld, test_2nd, start)), ]
  file_name<-paste0("all_cost_incremental_case.csv")
  write.csv( cost_mat,file.path(out_folder, file_name), row.names = FALSE)
  
  
  # QALY all sims qaly
  qaly_mat<-do.call("rbind",qaly_list)
  qaly_mat[] <- lapply(qaly_mat, unlist)
  qaly_mat<-qaly_mat[with(qaly_mat, order(model_struct,tpt_effect, igra, ptbld, test_2nd, start)), ]
  file_name<-paste0("all_qaly_incremental_case.csv")
  write.csv( qaly_mat,file.path(out_folder, file_name), row.names = FALSE)

}

# move QS2 files into a sub-sub-folder
dir.create(file.path("output",subfolder_name), showWarnings=F)
file.rename(list.files("output",pattern="\\.qs2$", full.names=T),
            file.path("output",subfolder_name,
              basename(list.files("output", pattern="\\.qs2$", full.names=T))))
