#########################################################################################
# This is how i am going to combine the calcium imaging data and the transcriptome data
# together.
#########################################################################################
#load the calcium imaging software.
main_dir <- getwd()
source("Y:/Box Sync/procpharm/procPharm 170210.r")
source("Y:/Box Sync/procpharm/170103.doodles.lee.R")

# Now load the Transcriptome Information
require(xlsx)
ts_info <- read.xlsx("./FX SS library data 190402 (u0591788@utah.edu) (1).xlsx",1)
ts_info$Cell.name <- paste0("X.", ts_info$Cell.name)

# What we are interested in is Sam Espinos data
exp_rows <- grep("^Trigeminal.*", ts_info$label_experiment,T)
exp_info <- ts_info[exp_rows,]

# we also need to ensure the experiment has an rd.name
#essentially i have picked up cells for sam in the past that did 
#not get sequenced.
exp_info <- exp_info[!is.na(exp_info$rd.name),]
row.names(exp_info) <- exp_info$Cell.name
exp_info$cDNA.conc.ug.mL <- as.numeric(exp_info$cDNA.conc.ug.mL)

#Now that we have sams data i am going to load the transcriptomic data Kevin has prepared
# each collumn is a sample 
ts_data <- read.csv("./Round 5/single.cell.rft.061219.csv")
row.names(ts_data) <- make.names(ts_data$Gene.name, unique=T)

# We need to combine the transcriptome information with the transcriptome data
# To do this we creat smaller data frames in the transcriptome and transcriptome information
# data.frames
exp_ts_cols <- Reduce(c,
        lapply(as.character(exp_info$Gnomex.Label), function(x) grep(paste0(x,"$"), colnames(ts_data))) 
) 
exp_ts <- ts_data[exp_ts_cols]
exp_ts <- t(exp_ts)
exp_ts <- round(exp_ts, digits=0)

#Now we need to the load up the RD. files into our workspace to begin creating the figures
rd_to_load <- as.character( unique( exp_info$rd.name ) )
rd_paths <- lapply(rd_to_load, function(x) grep(x, list.files(all.files = T, recursive = T), value=T ) )
for( i in rd_paths){load(i)}

#lets reduce the ts_info to only what is needed
exp_info_to_add <- exp_info[c("Gnomex.Label", "Cell.name", "M.Assigned" , "rd.name", "label_cellType", "label_experiment", "cDNA.conc.ug.mL")]

#now we need to add transcriptome data to the c.dat for only the specified cells
exps <- as.character( unique(exp_info$rd.name) )
for( i in 1:length(exps) ){
    tmp_rd <- get(exps[i])

    #add the trans_info to the c.dat
    new_c.dat <- merge(tmp_rd$c.dat, exp_info_to_add[exp_info_to_add$rd.name == exps[i], ], by="row.names", all.x=T)
    row.names(new_c.dat) <- new_c.dat$Row.names
    new_c.dat <- new_c.dat[ order(as.numeric(as.character(new_c.dat$Row.names))), ]
    tmp_rd$c.dat <- new_c.dat

    #add the transcriptome to the c.dat
    exp_info <- exp_info_to_add[ exp_info_to_add$rd.name == exps[i] , ]
    exp_info_gnomex_label <- as.character(exp_info$Gnomex.Label)
    exp_ts_exp_rn <- Reduce(c, lapply( exp_info_gnomex_label, function(x) grep(x, row.names(exp_ts), value=T) ) )
    new_exp_ts <- exp_ts[exp_ts_exp_rn, ]
    row.names(new_exp_ts) <- exp_info$Cell.name
    new_c.dat <- merge(tmp_rd$c.dat, new_exp_ts, by="row.names", all.x=T)
    row.names(new_c.dat) <- new_c.dat$Row.names
    new_c.dat <- new_c.dat[ order(as.numeric(as.character(new_c.dat$Row.names))), ]
    tmp_rd$c.dat <- new_c.dat
    assign(exps[i], tmp_rd)
}

######################################################################################
######################################################################################
source("./transcriptome_pharmer.r")
source("Y:/Box Sync/procpharm/procPharm 170210.r")

graphics.off()
alarm()
exps <- as.character( unique(exp_info$rd.name) )

LinesEvery_ts(exps, exp_info)
calca
mrgprd
nefh
Cntnap2
trpa1
trpm8
trpv1
kcna1
kcna2
kcna6
cacna1h



 