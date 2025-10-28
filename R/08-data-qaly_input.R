# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
if (interactive()) {

	# Input
	root    <- here::here()
	infile  <- file.path(root, "data", "raw", "inputs.xlsx")
	outfile <- file.path(root, "output", "qaly_input.qs2")

	# Packages
	source(file.path(root, "R", "modify_attach.R"))
	modify_attach(qs2,    include.only = "qs_save")
	modify_attach(readxl, include.only = "read_xlsx")

} else {

	# Input
	args    <- commandArgs(trailingOnly = TRUE)
	infile  <- args[1]
	outfile <- args[2]

	# Packages
	library(qs2,    include.only = "qs_save")
	library(readxl, include.only = "read_xlsx")
}

# -------------------------------------------------------------------------
# Load data ---------------------------------------------------------------
# -------------------------------------------------------------------------
dat <- lapply(1:5, function(i) as.data.frame(read_xlsx(infile, sheet = i)))
names(dat) <- c("q.male", "q.female", "qol", "covid.age", "age_bands")

# -------------------------------------------------------------------------
# Save data ---------------------------------------------------------------
# -------------------------------------------------------------------------
qs_save(dat, outfile)
