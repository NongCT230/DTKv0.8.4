#this script allows for running the app from R outside of RStudio: i.e. R console, RGui, etc.
#windows platform

#AUTODETECT APP DIRECTORY--------------------------------------------------------------

#set working dir to the one where R script resides to maintain relative directories
print("Attempting automatic detection of application directory...");
libraries.required <- setdiff("R.utils", rownames(installed.packages()));
if (length(libraries.required)>0){
  print("Automatic detection requires R.utils library");
  if(readline(prompt="Install from CRAN? y/n ") == 'y'){
    install.packages(libraries.required);
    library("R.utils");
  }
} else {
  #just load the library
  library("R.utils");
}

#recheck if library installed
libraries.required <- setdiff("R.utils", rownames(installed.packages()));
if (length(libraries.required)==0){
  working.dir = getSrcDirectory(function(x) {x});
  print(paste0("Detected directory: ", working.dir));
} else {
  working.dir = NULL;
}

if(is.null(working.dir) || is.na(working.dir) || length(working.dir)==0 || working.dir == ""){
  print("Automatic detection unsuccessful.");
  print("Can attempt manual detection using the file picker. (In the popup file picker choose any file in the application directory)");
  if(readline(prompt="Attempt manual detection using a file picker? y/n ") == 'y'){
    FILE <- file.choose();
    working.dir <- dirname(FILE);
  } else {
    working.dir <- ".";
  }
}
setwd(working.dir);
print(paste0("Working directory set to ", working.dir));

rm(libraries.required);

#LOAD AND LAUNCH APP-----------------------------------------------------------

source("app.R");
runApp(working.dir, launch.browser = TRUE);


