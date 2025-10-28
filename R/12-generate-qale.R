# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
if (interactive()) {

	# Input
	root              <- here::here()
	infile_parameters <- file.path(root, "output", "parameters.qs2")
	infile_qaly_input <- file.path(root, "output", "qaly_input.qs2")
	outfile           <- file.path(root, "output", "qale.qs2")

	# Packages
	source(file.path(root, "R", "modify_attach.R"))
	modify_attach(qs2,     include.only = c("qs_read", "qs_save"))
	modify_attach(data.table, include.only = c("setDT", "setDF", "rbindlist", "setnames"))

} else {

	# Input
	args              <- commandArgs(trailingOnly = TRUE)
	infile_parameters <- args[1]
	infile_qaly_input <- args[2]
	outfile           <- args[3]

	# Packages
	library(qs2,     include.only = c("qs_read", "qs_save"))
	library(data.table, include.only = c("setDT", "setDF", "rbindlist", "setnames"))
}

# -------------------------------------------------------------------------
# Constants ---------------------------------------------------------------
# -------------------------------------------------------------------------
# TODO - I'm assuming we won't change smr and can take it as 1
# TODO - I'm assuming we won't change country and can take it as UK
#country      <- "UK"
#smr          <- 1
qcm          <- 1
time_horizon <-100 #parameters$time_horizon Seee full LE


# -------------------------------------------------------------------------
# Load parameters ---------------------------------------------------------
# -------------------------------------------------------------------------
parameters <- qs_read(infile_parameters)
qaly_input <- qs_read(infile_qaly_input)

R <- parameters$p_dr


# -------------------------------------------------------------------------
# stuff
# -------------------------------------------------------------------------
q.male    <- qaly_input$q.male
q.female  <- qaly_input$q.female
qol       <- qaly_input$qol
age_bands <- qaly_input$age_bands

l_x_est <- function(dat){
	y <- dat[c("Age", "UK")]
	names(y) <- c("x","q_x")
	y$l_x <- 100000 # TODO - should this be cohort_size? Actually does this actually do anything in the end?
	for (i in seq_len(nrow(y))[-1L]) {
		y$l_x[i] <- y$l_x[i-1] * (1 - y$q_x[i-1])
	}
	y
}

q.male   <- l_x_est(q.male)
q.female <- l_x_est(q.female)
q.person <- merge(q.male, q.female, by="x")
colnames(q.person) <- c("x","q_male", "l_male",	"q_female", "l_female")
q.person$p.f <- with(q.person, l_female / (l_female + l_male))
q.person$l_person <- with(q.person, (p.f * l_female) + ((1 - p.f) * l_male))

q.person$bigl_x <- NA_real_
index <- seq_len(nrow(q.person) - 1)
q.person$bigl_x[index] <- with(q.person, (l_person[index] + l_person[index + 1]) / 2)
q.person$bigl_x[nrow(q.person)] <- with(q.person, l_person[nrow(q.person)]/ 2)
q.person$t_x <- rev(cumsum(rev(q.person$bigl_x)))

q.person$LE_x <- with(q.person, t_x / l_person)
dt.qol <- qol[c("low","high", "UK")]
names(dt.qol) <- c("low","high","qol_age")

setDT(q.person)
setDT(dt.qol)

qale <- q.person[dt.qol,
	on = .(x >= low, x <= high),
	nomatch = 0,
	.(x.x, l_person, bigl_x, t_x, LE_x,qol_age)]

qale[ , z_x := bigl_x * qol_age * qcm]
qale$t_adj <- rev(cumsum(rev(qale$z_x)))
qale[ , qale_x := t_adj/l_person]
qaly.calc <- qale[ , c("x.x","z_x")]

temp.q <- lapply(
	seq_len(nrow(qaly.calc)),
	function(i) qaly.calc[i:nrow(qaly.calc)]
)
temp.q <- rbindlist(temp.q, idcol = "column_label")
temp.q[ , column_label := column_label - 1]
temp.q[ , b_x := z_x / ((1+R)) ^ (x.x - column_label)] ## n.b x.x = u and column_label = x in the corresponding formulae in the CodeBook
temp.q[, index := (x.x <= time_horizon)]

total.b <- temp.q[ , .(bigb_x = sum(b_x), bigb_xfoo = sum(b_x[index])), by = column_label]

setnames(total.b, old = "column_label", new = "x.x")
qale <- merge(qale, total.b, by="x.x")
qale[ , dQALY := bigb_xfoo/l_person]

setDT(age_bands)
age_bands[ , midpoint := ceiling((low+high)/2)]
cov <- merge(qale, age_bands, by.x="x.x", by.y="midpoint", all=FALSE)
cov[,"Age Group":=paste(cov[,low],cov[,high],sep="-")]
setnames(cov, old=c("LE_x","qale_x","dQALY"),  new=c("LE","QALE","dQALY"))
agetab <- as.data.frame(cov[ , c("Age Group", "LE","QALE","dQALY")])

# -------------------------------------------------------------------------
# Save agetab -------------------------------------------------------------
# -------------------------------------------------------------------------
qs_save(agetab, outfile)
