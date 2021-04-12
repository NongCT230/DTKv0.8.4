
# Date and time function --------------------------------------------------

#current.date.time <- function(){
#  return (cat(format(Sys.Date(), "%d-%m-%Y"), " ", format(Sys.time(), "%H:%M:%S")));
#}


# Initialize logger -------------------------------------------------------
#logger levels: basic > default configuration, file > log.log output, stdout > console output
if (!exists("log.level.basic")){  log.level.basic <- 'WARNING';  } 
if (!exists("log.level.file")){   log.level.file <- 'WARNING';  }
if(!exists("log.level.stdout")){  log.level.stdout <- 'WARNING';  }
if (!exists("log.path")) {  log.path <- paste0(working.dir, paste0("/logs/log", format(Sys.time(), "%a_%b_%d_%Y_%H_%M_%S"), ".log")); }
if (!exists("log.to.file")) { log.to.file <- FALSE }

#library("logging");
logger.initialize <- function(env=parent.frame(), log.to.file = FALSE){
  basicConfig(level=log.level.basic);
  #logReset();
  setLevel(log.level.stdout ,getHandler('basic.stdout'));
  if (log.to.file){
    close( file( log.path, open="w" ) );
    addHandler(writeToFile, file=log.path, level=log.level.file);
  }
  #addHandler(writeToConsole, level=log.level.stdout);

  loginfo("Initializing logger...");
  #logwarn("Warning test...");
  #logdebug("Debug test...");
  #logerror("Error test...");
}

logger.initialize(log.to.file = log.to.file);

#clean up
rm(logger.initialize, log.level.basic, log.level.file, log.level.stdout, log.to.file, log.path);