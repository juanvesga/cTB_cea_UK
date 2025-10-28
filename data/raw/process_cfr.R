# -------------------------------------------------------------------------
# input -------------------------------------------------------------------
# -------------------------------------------------------------------------
if (interactive()) {

	# Input
	root <- here::here()
	outfile <- file.path(root,"data","raw", "parameters.qs2")

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


