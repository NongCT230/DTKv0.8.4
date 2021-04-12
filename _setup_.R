

# This is a setup environment which loads all the project libraries --------
#   and functions ---------------------------------------------------------

# shinyapps.io deployment version

# Cleanup -------------------------------------------------------------------

#clear workspace
rm(list=ls());


# Greeting ------------------------------------------------------------------

info <- "Dream.TK v0.8.4"

print(paste0("Welcome to ", info));
print(paste0("Running on ", format(Sys.Date(), "%d-%m-%Y"), " ", format(Sys.time(), "%H:%M:%S")));
print("---------------------------------------------------------------------------------------------");

# Working directory ---------------------------------------------------------

#set working dir to the one where _setup_.R resides to maintain relative directories
#trick, does not assume admin access, ONLY WORKS WHEN SOURCING
working.dir = getSrcDirectory(function(x) {x});
if(is.null(working.dir) || is.na(working.dir) || length(working.dir)==0 || working.dir == ""){
  working.dir <- ".";
}
setwd(working.dir);
print(paste0("Working directory set to ", working.dir));

# Import required libraries and functions ---------------------------------

print("Starting library imports...");
source("./_library_imports_.R");

log.level.basic <- 'DEBUG';
log.level.file <- 'DEBUG';
log.level.stdout <- 'DEBUG';
log.to.file <- TRUE;
source("./_logging_.R");
#repeat greeting in log file
loginfo(paste0("Welcome to ", info));
loginfo(paste0("Running on ", format(Sys.Date(), "%d-%m-%Y"), " ", format(Sys.time(), "%H:%M:%S")));
loginfo("-------------------------------------------------");

loginfo("Starting function imports...");
source("./_function_imports_.R");

