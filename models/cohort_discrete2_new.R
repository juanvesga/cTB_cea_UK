
# time deffinition --------------------------------------------------------

dt <- user()
initial(time) <- 0
update(time) <- (step + 1) * dt

tt<-step #as.integer(t)




# Initial conditions ---------------------------------------------------
baseline<- if (cascade==0) 0 else 1

initial(Igranegf0[])     <- if (cascade==0) init_testneg[i] * frac_recent else
  init_testneg[i] * frac_recent * true_neg + init_testneg[i] * frac_recent * (1-true_neg)* (1-p_cascade)

initial(Igranegs0[])     <- if (cascade==0) init_testneg[i] * (1-frac_recent) else
  init_testneg[i] * (1-frac_recent) * true_neg + init_testneg[i] * (1-frac_recent) * (1-true_neg)* (1-p_cascade)

initial(Igranegf1[])     <- if (cascade==0) 0 else init_testneg[i] * frac_recent * (1-true_neg) * p_cascade

initial(Igranegs1[])     <- if (cascade==0) 0 else init_testneg[i] * (1-frac_recent) * (1-true_neg) * p_cascade

initial(Igraposf0[])     <- if (cascade==0) init_testpos[i] * frac_recent else
  init_testpos[i] * frac_recent * (1-true_pos) + init_testpos[i] * frac_recent * true_pos * (1-p_cascade)

initial(Igraposs0[])     <- if (cascade==0) init_testpos[i] * (1-frac_recent) else
  init_testpos[i] * (1-frac_recent) * (1-true_pos) + init_testpos[i] * (1-frac_recent) * true_pos* (1-p_cascade)

initial(Igraposf1[])     <- if (cascade==0) 0 else init_testpos[i] * frac_recent * true_pos * p_cascade

initial(Igraposs1[])     <- if (cascade==0) 0 else init_testpos[i] * (1-frac_recent) * true_pos * p_cascade
initial(pTB[])        <- 0
initial(eTB[])        <- 0
initial(pTX[])        <- 0
initial(eTX[])        <- 0
initial(R[])          <- 0
initial(LD[])         <- 0

# Core equations ----------------------------------------------------------

# IGRA negative fast progression no TPT
#Igraneg=IGRA negative ; f = fast prgression to TB; 0 = no TPT

update(Igranegf0[1]) <-
  Igranegf0[1] -
  n_ageoIgranegf0[1]-
  n_Igranegf0_Igranegs0[1] -
  n_Igranegf0_tb[1] -
  n_muIgranegf0[1]

update(Igranegf0[2:N_age]) <-
  Igranegf0[i] +
  n_ageoIgranegf0[i-1]-
  n_ageoIgranegf0[i]-
  n_Igranegf0_Igranegs0[i] -
  n_Igranegf0_tb[i] -
  n_muIgranegf0[i]


# IGRA negative slow progression no TPT
#Igraneg=IGRA negative ; s = slow prgression to TB; 0 = no TPT
update(Igranegs0[1]) <-
  Igranegs0[1] -
  n_ageoIgranegs0[1]-
  n_Igranegs0_tb[1] -
  n_muIgranegs0[1]

update(Igranegs0[2:N_age]) <-
  Igranegs0[i] +
  n_ageoIgranegs0[i-1]+
  n_Igranegf0_Igranegs0[i-1]-
  n_ageoIgranegs0[i] -
  n_Igranegs0_tb[i] -
  n_muIgranegs0[i]



# IGRA negative fast progression no TPT
#Igraneg=IGRA negative ; f = fast prgression to TB; 1 = TPT

update(Igranegf1[1]) <-
  Igranegf1[1] -
  n_ageoIgranegf1[1]-
  n_Igranegf1_Igranegs1[1] -
  n_Igranegf1_tb[1] -
  n_muIgranegf1[1]

update(Igranegf1[2:N_age]) <-
  Igranegf1[i] +
  n_ageoIgranegf1[i-1]-
  n_ageoIgranegf1[i]-
  n_Igranegf1_Igranegs1[i] -
  n_Igranegf1_tb[i] -
  n_muIgranegf1[i]


# IGRA negative slow progression no TPT
#Igraneg=IGRA negative ; s = slow prgression to TB; 1 = TPT

update(Igranegs1[1]) <-
  Igranegs1[1] -
  n_ageoIgranegs1[1]-
  n_Igranegs1_tb[1] -
  n_muIgranegs1[1]

update(Igranegs1[2:N_age]) <-
  Igranegs1[i] +
  n_ageoIgranegs1[i-1]+
  n_Igranegf1_Igranegs1[i-1]-
  n_ageoIgranegs1[i] -
  n_Igranegs1_tb[i] -
  n_muIgranegs1[i]


# IGRA positive fast progression no TPT
#Igrapos=IGRA positive ; f = fast prgression to TB; 0 = no TPT

update(Igraposf0[1]) <-
  Igraposf0[1] -
  n_ageoIgraposf0[1]-
  n_Igraposf0_Igraposs0[1] -
  n_Igraposf0_tb[1] -
  n_muIgraposf0[1]

update(Igraposf0[2:N_age]) <-
  Igraposf0[i] +
  n_ageoIgraposf0[i-1]-
  n_ageoIgraposf0[i]-
  n_Igraposf0_Igraposs0[i] -
  n_Igraposf0_tb[i] -
  n_muIgraposf0[i]


# IGRA positive slow progression no TPT
#Igrapos=IGRA negative ; s = slow prgression to TB; 0 = no TPT
update(Igraposs0[1]) <-
  Igraposs0[1] -
  n_ageoIgraposs0[1]-
  n_Igraposs0_tb[1] -
  n_muIgraposs0[1]

update(Igraposs0[2:N_age]) <-
  Igraposs0[i] +
  n_ageoIgraposs0[i-1]+
  n_Igraposf0_Igraposs0[i-1]-
  n_ageoIgraposs0[i] -
  n_Igraposs0_tb[i] -
  n_muIgraposs0[i]



# IGRA positive fast progression no TPT
#Igrapos=IGRA negative ; f = fast prgression to TB; 1 = TPT

update(Igraposf1[1]) <-
  Igraposf1[1] -
  n_ageoIgraposf1[1]-
  n_Igraposf1_Igraposs1[1] -
  n_Igraposf1_tb[1] -
  n_muIgraposf1[1]

update(Igraposf1[2:N_age]) <-
  Igraposf1[i] +
  n_ageoIgraposf1[i-1]-
  n_ageoIgraposf1[i]-
  n_Igraposf1_Igraposs1[i] -
  n_Igraposf1_tb[i] -
  n_muIgraposf1[i]


# IGRA positive slow progression no TPT
#Igrapos=IGRA negative ; s = slow prgression to TB; 1 = TPT

update(Igraposs1[1]) <-
  Igraposs1[1] -
  n_ageoIgraposs1[1]-
  n_Igraposs1_tb[1] -
  n_muIgraposs1[1]

update(Igraposs1[2:N_age]) <-
  Igraposs1[i] +
  n_ageoIgraposs1[i-1]+
  n_Igraposf1_Igraposs1[i-1]-
  n_ageoIgraposs1[i] -
  n_Igraposs1_tb[i] -
  n_muIgraposs1[i]

# Pulmonary TB
update(pTB[1]) <-
  pTB[1]  -
  n_ageopTB[1]-
  n_ptb_ptx[1] -
  n_mupTB[1] -
  n_mupTB_tb[1]

update(pTB[2:N_age]) <-
  pTB[i] +
  n_ageopTB[i-1]+
  n_ptb[i-1]  -
  n_ageopTB[i] -
  n_ptb_ptx[i] -
  n_mupTB[i] -
  n_mupTB_tb[i]


# extra-Pulmonary TB
update(eTB[1]) <-
  eTB[1]  -
  n_ageoeTB[1]-
  n_etb_etx[1] -
  n_mueTB[1]

update(eTB[2:N_age]) <-
  eTB[i] +
  n_ageoeTB[i-1]+
  n_etb[i-1] -
  n_ageoeTB[i]-
  n_etb_etx[i] -
  n_mueTB[i]


# Pulmonary TB on treatment
update(pTX[1]) <-
  pTX[1] -
  n_ageopTX[1]-
  n_ptx_ld[1] -
  n_ptx_r[1] -
  n_mupTX[1] -
  n_mupTX_tx[1]

update(pTX[2:N_age]) <-
  pTX[i] +
  n_ageopTB[i-1] +
  n_ptb_ptx[i-1] -
  n_ageopTX[i]-
  n_ptx_ld[i] -
  n_ptx_r[i] -
  n_mupTX[i] -
  n_mupTX_tx[i]

# extra-Pulmonary TB on treatment
update(eTX[1]) <-
  eTX[1] -
  n_ageoeTX[1]-
  n_etx_r[1] -
  n_mueTX[1]-
  n_mueTX_tx[1]

update(eTX[2:N_age]) <-
  eTX[i] +
  n_ageoeTX[i-1]+
  n_etb_etx[i-1] -
  n_ageoeTX[i]-
  n_etx_r[i] -
  n_mueTX[i] -
  n_mueTX_tx[1]

# post pulmonary TB lung disease
update(LD[1]) <-
  LD[1]  -
  n_ageoLD[1]-
  n_muLD[1] -
  n_muLD_ld[1]

update(LD[2:N_age]) <-
  LD[i] +
  n_ageoLD[i-1] +
  n_ptx_ld[i-1] -
  n_ageoLD[i]-
  n_muLD[i] -
  n_muLD_ld[i]

# Recovered
update(R[1]) <-
  R[1] -
  n_ageoR[1]-
  n_muR[1]

update(R[2:N_age]) <-
  R[i] +
  n_ageoR[i-1]+
  n_ptx_r[i-1] +
  n_etx_r[i-1] -
  n_ageoR[i]-
  n_muR[i]




dim(Igranegf0)<- N_age
dim(Igranegs0)<- N_age
dim(Igranegf1)<- N_age
dim(Igranegs1)<- N_age
dim(Igraposf0)<- N_age
dim(Igraposs0)<- N_age
dim(Igraposf1)<- N_age
dim(Igraposs1)<- N_age
dim(pTB)    <- N_age
dim(eTB)    <- N_age
dim(pTX)    <- N_age
dim(eTX)    <- N_age
dim(R)      <- N_age
dim(LD)     <- N_age

# Parametes ---------------------------------------------------

shape_fac <- if (as.integer(step) >= n_shape_steps)
  rate_shape[n_shape_steps]*
  (1- (as.integer(step)-n_shape_steps)) else
    rate_shape[as.integer(step) + 1]



p_mu[]             <- mu[i]#1 - exp(-mu[i] * dt)
p_mu_tb[]          <- mu_tb[i]#1 - exp(-mu_tb[i]* dt)
p_mu_tbtx[]        <- mu_tbtx[i]#1 - exp(-mu_tbtx[i]* dt)
p_mu_etbtx[]       <- mu_etbtx[i]#1 - exp(-mu_etbtx[i]* dt)
p_reac             <- 1  -exp(-r_reactivation*dt*shape_fac)
p_reac_fn_f        <- 1  -exp(-r_reactivation_fn_f*dt*shape_fac)
p_reac_fn_s        <- 1  -exp(-r_reactivation_fn_s*dt*shape_fac)
p_reac_slow        <- 1  -exp(-r_reactivation_slow*dt*shape_fac)
p_slow             <- 1  -exp(-r_slow*dt)
p_cascade          <- 1 - exp(-cascade* dt)
p_tptae[]          <- 1 - exp(-r_tptae[i]* dt)
p_tpteff           <- r_tpteff* dt
p_mu_ptbld[]       <- 1 - exp(- mu[i] * rr_ptbld_mu * dt)
p_tbdur            <- 1 - exp(-r_tbdur* dt)
p_etbdur           <- 1 - exp(-r_etbdur* dt)
p_txdur            <- 1 - exp(-r_txdur* dt)
p_aging            <- 1


# Hard coded

frac_recent<- 0.09#0.052 # From Schwalb Lancet

dim(p_mu)<-N_age
dim(p_mu_tb)<-N_age
dim(p_mu_tbtx)<-N_age
dim(p_mu_etbtx)<-N_age
dim(p_tptae)<-N_age
dim(p_mu_ptbld)<-N_age


r0<-user()
N_age <- user()
contact_matrix[,]<-user()
mu[]<-user()
mu_tb[]<-user()
mu_tbtx[]<-user()
mu_etbtx[]<-user()
init_testpos[]<-user()
init_testneg[]<-user()
age_elig[]<-user()
r_reactivation <-user()
r_reactivation_fn_f <-user()
r_reactivation_fn_s <-user()
r_reactivation_slow <-user()
rate_shape[]<-user()
n_shape_steps<-user()
r_slow <-user()
true_neg <- user()
true_pos <- user()
cascade<-user()
p_2nd  <-user()
tpt_start<-user()
tpt_completion<-user()
r_tpteff <-user()
r_tptae[] <-user()
age_rate <-user()
p_ptb[] <-user()
p_ptbld[] <-user()
rr_ptbld_mu<-user()
r_tbdur<-user()
r_etbdur<-user()
r_txdur<-user()
p_tbhosp<-user()
qol[]<-user()
dQALY[]<-user()
qol_ptb_weight<-user()
qol_eptb_weight<-user()
qol_ptbld_weight<-user()
qol_tx_weight<-user()
qol_tptae_loss<-user()
qol_tbhosp_loss<-user()
c_tptae<- user()
c_test <- user()
c_test2_pos <- user()
c_test2_neg <- user()
c_posi_test<-user()
c_tpt<-user()
c_tbdx_ptb<-user()
c_tbdx_eptb<-user()
c_tbtx<-user()
c_tbhosp<-user()
c_post<-user()
disc_rate<-user(0.035)



dim(mu)<-N_age
dim(contact_matrix) <- c(N_age, N_age)
dim(mu_tb)<-N_age
dim(mu_tbtx)<-N_age
dim(mu_etbtx)<-N_age
dim(init_testpos)<-N_age
dim(init_testneg)<-N_age
dim(age_elig)<-N_age
dim(rate_shape)<-n_shape_steps
dim(r_tptae)<-N_age
dim(p_ptb)<-N_age
dim(p_ptbld)<-N_age
dim(qol)<-N_age
dim(dQALY)<-N_age

# Transitions -------------------------------------------------------------

# Background mortality

n_muIgranegf0[]     <- Igranegf0[i] * p_mu[i]
n_muIgranegs0[]     <- Igranegs0[i] * p_mu[i]
n_muIgranegf1[]     <- Igranegf1[i] * p_mu[i]
n_muIgranegs1[]     <- Igranegs1[i] * p_mu[i]
n_muIgraposf0[]     <- Igraposf0[i] * p_mu[i]
n_muIgraposs0[]     <- Igraposs0[i] * p_mu[i]
n_muIgraposf1[]     <- Igraposf1[i] * p_mu[i]
n_muIgraposs1[]     <- Igraposs1[i] * p_mu[i]
n_mupTB[]           <- pTB[i]  * p_mu[i]
n_mueTB[]           <- eTB[i]  * p_mu[i]
n_mupTX[]           <- pTX[i]  * p_mu[i]
n_mueTX[]           <- eTX[i]  * p_mu[i]
n_muLD[]            <- LD[i]   * p_mu[i]
n_muR[]             <- R[i]    * p_mu[i]



# Transitions Igra neg fast no tpt
n_Igranegf0_tb[]   <- (Igranegf0[i] - n_muIgranegf0[i]) * p_reac_fn_f

n_Igranegf0_Igranegs0[] <- (Igranegf0[i] - n_muIgranegf0[i] - n_Igranegf0_tb[i]) * p_slow

n_ageoIgranegf0[] <- Igranegf0[i] - n_muIgranegf0[i] - n_Igranegf0_tb[i] - n_Igranegf0_Igranegs0[i]


# Transitions Igra neg slow no tpt
n_Igranegs0_tb[] <- (Igranegs0[i] - n_muIgranegs0[i]) * p_reac_fn_s

n_ageoIgranegs0[] <- Igranegs0[i] - n_muIgranegs0[i] - n_Igranegs0_tb[i]


# Transitions Igra neg fast tpt
n_Igranegf1_tb[]   <- (Igranegf1[i] - n_muIgranegf1[i]) * p_reac_fn_f * (1-p_tpteff)

n_Igranegf1_Igranegs1[] <- (Igranegf1[i] - n_muIgranegf1[i] - n_Igranegf1_tb[i]) * p_slow

n_ageoIgranegf1[] <- Igranegf1[i] - n_muIgranegf1[i] - n_Igranegf1_tb[i] - n_Igranegf1_Igranegs1[i]


# Transitions Igra neg slow  tpt
n_Igranegs1_tb[] <- (Igranegs1[i] - n_muIgranegs1[i]) * p_reac_fn_s  * (1-p_tpteff)

n_ageoIgranegs1[] <- Igranegs1[i] - n_muIgranegs1[i] - n_Igranegs1_tb[i]



# Transitions Igra pos fast no tpt
n_Igraposf0_tb[]   <- (Igraposf0[i] - n_muIgraposf0[i]) * p_reac

n_Igraposf0_Igraposs0[] <- (Igraposf0[i] - n_muIgraposf0[i] - n_Igraposf0_tb[i]) * p_slow

n_ageoIgraposf0[] <- Igraposf0[i] - n_muIgraposf0[i] - n_Igraposf0_tb[i] - n_Igraposf0_Igraposs0[i]


# Transitions Igra pos slow no tpt
n_Igraposs0_tb[] <- (Igraposs0[i] - n_muIgraposs0[i]) * p_reac_slow

n_ageoIgraposs0[] <- Igraposs0[i] - n_muIgraposs0[i] - n_Igraposs0_tb[i]


# Transitions Igra pos fast tpt
n_Igraposf1_tb[]   <- (Igraposf1[i] - n_muIgraposf1[i]) * p_reac * (1-p_tpteff)

n_Igraposf1_Igraposs1[] <- (Igraposf1[i] - n_muIgraposf1[i] - n_Igraposf1_tb[i]) * p_slow

n_ageoIgraposf1[] <- Igraposf1[i] - n_muIgraposf1[i] - n_Igraposf1_tb[i] - n_Igraposf1_Igraposs1[i]


# Transitions Igra pos slow  tpt
n_Igraposs1_tb[] <- (Igraposs1[i] - n_muIgraposs1[i]) * p_reac_slow  * (1-p_tpteff)

n_ageoIgraposs1[] <- Igraposs1[i] - n_muIgraposs1[i] - n_Igraposs1_tb[i]


# Overall number developing  primary TB (transitioned from TBI)


new_primary_tb[]<-
  n_Igranegf0_tb[i] +
  n_Igranegs0_tb[i] +
  n_Igranegf1_tb[i] +
  n_Igranegs1_tb[i] +
  n_Igraposf0_tb[i] +
  n_Igraposs0_tb[i] +
  n_Igraposf1_tb[i] +
  n_Igraposs1_tb[i]


# Overall secondary cases (arising from active TB in the model , with a R0=r0 )
# and allocated by age using a contact matrix from POLYMOD (Cij)
new_secondary_tb_mat[,]<-new_primary_tb[i]*contact_matrix[i,j]*r0

new_secondary_tb[]<-sum(new_secondary_tb_mat[,i])

# Overall Pulmonary TB
n_ptb[]      <- ( new_primary_tb[i] + new_secondary_tb[i]) *p_ptb[i]

# Pulmonary TB transitions
n_mupTB_tb[] <- ((pTB[i]-n_mupTB[i]) * p_mu_tb[i])

n_ptb_ptx[]   <- ((pTB[i]-n_mupTB[i]-n_mupTB_tb[i]) * p_tbdur)

n_ageopTB[]   <-  ((pTB[i]-n_mupTB[i]-n_mupTB_tb[i] -  n_ptb_ptx[i]) * p_aging)

# Overall extra-Pulmonary TB
n_etb[]      <-  (new_primary_tb[i] + new_secondary_tb[i] ) - n_ptb[i]

# extra-Pulmonary TB transitions
n_etb_etx[]   <- ((eTB[i]-n_mueTB[i]) * p_etbdur)

n_ageoeTB[]   <-  eTB[i]-n_mueTB[i] -  n_etb_etx[i]

# Pulmonary TB on treatment transitions
n_mupTX_tx[] <- ((pTX[i]-n_mupTX[i]) * p_mu_tbtx[i])

n_ptx_out[]  <- ((pTX[i]-n_mupTX[i]-n_mupTX_tx[i]) * p_txdur)

n_ptx_ld[]   <- (n_ptx_out[i] * p_ptbld[i])

n_ptx_r[]    <- n_ptx_out[i] - n_ptx_ld[i]

n_ageopTX[]  <- ((pTX[i]-n_mupTX[i]-n_mupTX_tx[i]-n_ptx_out[i]) * p_aging)


# extra Pulmonary TB on treatment transitions
n_mueTX_tx[] <- ((eTX[i]-n_mupTX[i]) * p_mu_etbtx[i])

n_etx_r[]    <- ((eTX[i]-n_mueTX[i]-n_mueTX_tx[i] ) * p_txdur)

n_ageoeTX[]  <- ((eTX[i]-n_mueTX[i]-n_mueTX_tx[i]-n_etx_r[i]) * p_aging)

# postPulmonary TB transitions
n_muLD_ld[]  <- ((LD[i]-n_muLD[i]) * p_mu_ptbld[i])

n_ageoLD[]   <- ((LD[i]-n_muLD[i] - n_muLD_ld[i]) *p_aging)

# Recovered transitions
n_ageoR[]   <- ((R[i]-n_muR[i]) * p_aging)





dim(n_Igranegf0_tb)<-N_age
dim(n_Igranegf0_Igranegs0)<-N_age
dim(n_Igranegs0_tb)<-N_age
dim(n_Igranegf1_tb)<-N_age
dim(n_Igranegf1_Igranegs1)<-N_age
dim(n_Igranegs1_tb)<-N_age
dim(n_Igraposf0_tb)<-N_age
dim(n_Igraposf0_Igraposs0)<-N_age
dim(n_Igraposs0_tb)<-N_age
dim(n_Igraposf1_tb)<-N_age
dim(n_Igraposf1_Igraposs1)<-N_age
dim(n_Igraposs1_tb)<-N_age
dim(new_primary_tb)<-N_age
dim(new_secondary_tb)<-N_age
dim(new_secondary_tb_mat)<-c(N_age,N_age)
dim(n_ptb)      <- N_age
dim(n_mupTB_tb) <- N_age
dim(n_ptb_ptx)   <- N_age
dim(n_etb)      <- N_age
dim(n_etb_etx)   <- N_age
dim(n_mupTX_tx) <- N_age
dim(n_mueTX_tx) <- N_age
dim(n_ptx_out)  <- N_age
dim(n_ptx_ld)   <- N_age
dim(n_ptx_r)    <- N_age
dim(n_etx_r)    <- N_age
dim(n_muLD_ld)  <- N_age

dim(n_muIgranegf0)    <-N_age#
dim(n_muIgranegs0)    <-N_age#
dim(n_muIgranegf1)    <-N_age#
dim(n_muIgranegs1)    <-N_age#
dim(n_muIgraposf0)    <-N_age#
dim(n_muIgraposs0)    <-N_age#
dim(n_muIgraposf1)    <-N_age#
dim(n_muIgraposs1)    <-N_age#
dim(n_mupTB)    <- N_age
dim(n_mueTB)    <- N_age
dim(n_mupTX)  <- N_age
dim(n_mueTX)  <- N_age
dim(n_muR)      <- N_age
dim(n_muLD)  <- N_age


dim(n_ageoIgranegf0)<-N_age
dim(n_ageoIgranegs0)<-N_age
dim(n_ageoIgranegf1)<-N_age
dim(n_ageoIgranegs1)<-N_age
dim(n_ageoIgraposf0)<-N_age
dim(n_ageoIgraposs0)<-N_age
dim(n_ageoIgraposf1)<-N_age
dim(n_ageoIgraposs1)<-N_age
dim(n_ageopTB)    <- N_age
dim(n_ageoeTB)    <- N_age
dim(n_ageopTX)  <- N_age
dim(n_ageoeTX)  <- N_age
dim(n_ageoR)      <- N_age
dim(n_ageoLD)  <- N_age





# Get output objects ------------------------------------------------------------------

n_tot[] <-
  Igranegf0[i]+
  Igranegs0[i]+
  Igranegf1[i]+
  Igranegs1[i]+
  Igraposf0[i]+
  Igraposs0[i]+
  Igraposf1[i]+
  Igraposs1[i]+
  pTB[i] +
  eTB[i] +
  pTX[i] +
  eTX[i] +
  LD[i] +
  R[i]

new_deaths_age[]<-
  n_muIgranegf0[i]+
  n_muIgranegs0[i]+
  n_muIgranegf1[i]+
  n_muIgranegs1[i]+
  n_muIgraposf0[i]+
  n_muIgraposs0[i]+
  n_muIgraposf1[i]+
  n_muIgraposs1[i]+
  n_mupTB[i]+
  n_mueTB[i]+
  n_mupTX[i]+
  n_mueTX[i]+
  n_muR[i]+
  n_muLD[i]+
  n_mupTB_tb[i]+
  n_mupTX_tx[i]+
  n_mueTX_tx[i]+
  n_muLD_ld[i]



new_allTBcases_pos_age[] <-
  n_Igraposf0_tb[i]+
  n_Igraposs0_tb[i]+
  n_Igraposf1_tb[i]+
  n_Igraposs1_tb[i]



new_allTBcases_neg_age[] <-
  n_Igranegf0_tb[i]+
  n_Igranegs0_tb[i]+
  n_Igranegf1_tb[i]+
  n_Igranegs1_tb[i]


new_allTBcases_tpt_age[] <-
  n_Igranegf1_tb[i]+
  n_Igranegs1_tb[i]+
  n_Igraposf1_tb[i]+
  n_Igraposs1_tb[i]


new_TBcases_age[] <- n_ptb[i]


new_TBtx_age[] <- n_ptb_ptx[i] + n_etb_etx[i]


new_eTBcases_age[] <- n_etb[i]


new_TBdeaths_age[]<- n_mupTB_tb[i] + n_mupTX_tx[i] + n_mueTX_tx[i]

new_PTBLDdeaths_age[]<- n_muLD_ld[i]

new_tpt_age[]<- if(tt==1) init_testpos[i] * true_pos * p_cascade * age_elig[i] +
  init_testneg[i] * (1-true_neg) * p_cascade * age_elig[i] else 0



dr<- 1/((1+disc_rate)^tt)


# QALYs -------------------------------------------------------------------

PopQALY[]<-
  dr * (    Igranegf0[i]*qol[i]+
            Igranegs0[i]*qol[i]+
            Igranegf1[i]*qol[i]+
            Igranegs1[i]*qol[i]+
            Igraposf0[i]*qol[i]+
            Igraposs0[i]*qol[i]+
            Igraposf1[i]*qol[i]+
            Igraposs1[i]*qol[i]+
            pTB[i]*(qol_ptb_weight/ref_qol)+
            eTB[i]*(qol_eptb_weight/ref_qol)*qol[i]+
            pTX[i]*(qol_tx_weight/ref_qol)*qol[i]+
            eTX[i]*(qol_tx_weight/ref_qol)*qol[i]+
            LD[i]*(qol_ptbld_weight)*qol[i]+
            R[i]*qol[i])
ref_qol<-0.91
# TB states QALYs
qaly_TB_age[]<-
  dr * (
    (pTB[i]*qol[i] - pTB[i]*(qol_ptb_weight/ref_qol) * qol[i])+# min(pTB[i]*qol_ptb_weight,pTB[i]*qol[i]) )+
      (eTB[i]*qol[i] - eTB[i]*(qol_eptb_weight/ref_qol)* qol[i]))

qaly_TBtx_age[]<-
  dr * (
    (pTX[i]*qol[i] - pTX[i]*(qol_tx_weight/ref_qol  )* qol[i])+#min(pTX[i]*qol_tx_weight,pTX[i]*qol[i]) )+
      (eTX[i]*qol[i] - eTX[i]*(qol_tx_weight/ref_qol  )* qol[i]) )#min(eTX[i]*qol_tx_weight,eTX[i]*qol[i]) ))

# PTLD states QALYs
qaly_postTB_age[]<- dr*(LD[i]*qol[i]- LD[i]*(qol_ptbld_weight)*qol[i])# min(LD[i]*qol_ptbld_weight,LD[i]*qol[i])


# Life Expectancy Qaly loss
qaly_LE_age[]<- dr * (
  n_mupTB_tb[i]+
    n_mupTX_tx[i]+
    n_mueTX_tx[i]+
    n_muLD_ld[i])*dQALY[i]



# TB hospitalisation event Qaly loss
qaly_TB_hosp_age[]<-dr * (qol_tbhosp_loss* new_TBcases_age[i] * p_tbhosp)



# TPT AE QALY loss
qaly_tptae_age[]<-  dr * (qol_tptae_loss*r_tptae[i]* new_tpt_age[i] )


QALY_ALL_age[]<-
  qaly_TB_age[i]+
  qaly_TBtx_age[i]+
  qaly_postTB_age[i]+
  qaly_TB_hosp_age[i]+
  qaly_LE_age[i]+
  qaly_tptae_age[i]



# Costs


oneoffcosts_age[]<-
  if(tt==1)
    age_elig[i]*(init_testpos[i]  * c_test  +  init_testpos[i]  * p_2nd * c_test2_pos +
       init_testneg[i] * c_test +  init_testneg[i]  * p_2nd * c_test2_neg)  else 0

oneoffcosts<-   sum(oneoffcosts_age) * dr



starttptcosts_age[]<-
  if(tt==1)
    new_tpt_age[i]  * (c_tpt+c_posi_test) * age_elig[i]  else 0

starttptcosts<-   sum(starttptcosts_age) * dr





# need to add cost at start of those not completing
cost_tpt<-starttptcosts#sum(new_tpt_age) * c_tpt

cost_tptae_age[]<- dr * new_tpt_age[i] * r_tptae[i] * c_tptae



cost_tb <- dr * (   sum(n_ptb_ptx) * (c_tbdx_ptb + c_tbtx) +
				  	sum(n_etb_etx) * (c_tbdx_eptb + c_tbtx) +
                   sum(new_TBcases_age) * p_tbhosp * c_tbhosp)

cost_ptbld<- dr * (sum(LD)*c_post)


# Output vectors
output(N_tot_age[]) <- n_tot[i]
output(N_tot) <- sum(n_tot)

output(new_TBcases)<-  sum(new_TBcases_age)
output(new_eTBcases)<- sum(new_eTBcases_age)
output(new_secondaryTB)<-sum(new_secondary_tb)
output(new_ptbld_cases)<-sum(n_ptx_ld)

output(new_allTBcases_pos)<-sum(new_allTBcases_pos_age)
output(new_allTBcases_neg)<-sum(new_allTBcases_neg_age)
output(new_allTBcases_tpt)<-sum(new_allTBcases_tpt_age)

output(new_TBtx)<- sum(new_TBtx_age)

output(new_TBdeaths)<- sum(new_TBdeaths_age)
output(new_PTBLDdeaths)<-sum(new_PTBLDdeaths_age)
output(new_deaths)<- sum(new_deaths_age)

output(new_TBdeaths_age)   <- new_TBdeaths_age
output(new_PTBLDdeaths_age)<- new_PTBLDdeaths_age
output(new_deaths_age)     <- new_deaths_age



output(new_tpt)    <- sum(new_tpt_age)
output(qaly_LE_loss)<-sum(qaly_LE_age)
output(qaly_TB)<-sum(qaly_TB_age)
output(qaly_TBtx)<-sum(qaly_TBtx_age)
output(qaly_TBhosp_loss)<-sum(qaly_TB_hosp_age)
output(qaly_postTB)<-sum(qaly_postTB_age)
output(qaly_tptae_loss)<-sum(qaly_tptae_age)
output(qaly_pop)<-sum(PopQALY) - (sum(qaly_TB_hosp_age) + sum(qaly_tptae_age))
output(qaly_all)<-sum(QALY_ALL_age)
output(costs)<- oneoffcosts+cost_tpt+sum(cost_tptae_age)+cost_tb+cost_ptbld
output(cost_test)<- oneoffcosts
output(cost_tpt)<-cost_tpt
output(cost_tptae)<-sum(cost_tptae_age)
output(cost_tb)<-cost_tb
output(cost_post)<-cost_ptbld
output(ptbld_pyr)<-sum(LD)




# Define dimensions
dim(new_allTBcases_neg_age)<-N_age
dim(new_allTBcases_pos_age)<-N_age
dim(new_allTBcases_tpt_age)<-N_age
dim(new_TBcases_age)<-N_age
dim(new_eTBcases_age)<-N_age
dim(new_PTBLDdeaths_age)<-N_age
dim(new_TBdeaths_age)<-N_age
dim(new_TBtx_age)<-N_age
dim(new_tpt_age)<-N_age
dim(new_deaths_age)<-N_age
dim(N_tot_age)<-N_age
dim(n_tot)<-N_age
dim(QALY_ALL_age)<-N_age
dim(qaly_LE_age)<-N_age
dim(qaly_TB_age)<-N_age
dim(qaly_TBtx_age)<-N_age
dim(qaly_postTB_age)<-N_age
dim(qaly_tptae_age)<-N_age
dim(qaly_TB_hosp_age)<-N_age
dim(cost_tptae_age)<-N_age
dim(oneoffcosts_age)<-N_age
dim(starttptcosts_age)<-N_age
dim(PopQALY)<-N_age


