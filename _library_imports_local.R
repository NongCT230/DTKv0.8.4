
# Library imports ---------------------------------------------------------



# Obtain required library package list ------------------------------------

#check if working directory exists
if(!exists("working.dir")){
  working.dir <- "";
}

#set app.ini filepath to open, create list of import packages
source(paste0(working.dir,"/_readini_.R"));
if (!exists("ini.path")) ini.path <- paste0(working.dir,"/app.ini");
libraries.required <- readini(ini.path, "[Import library packages]");
bioc.libraries.required <- readini(ini.path, "[Import Bioconductor packages]");

# Install missing libraries -----------------------------------------------

#if not installed, ask user to install or stop program
libraries.notinstalled <- setdiff(libraries.required, rownames(installed.packages()));
if (length(libraries.notinstalled)>0){
  print(paste0("Missing libraries to run this program: ", paste0(libraries.notinstalled, collapse=", ")));
  if(readline(prompt="Install? y/n ") == 'y'){
    print("Installing libraries...");
    install.packages(libraries.notinstalled, dependencies = c("Depends", "Imports", "LinkingTo"));
  }
}
libraries.notinstalled <- setdiff(libraries.required, rownames(installed.packages()));
if (length(libraries.notinstalled)>0){
  print("Failed to install libraries. Terminating program.");
  rm(libraries.required, libraries.notinstalled);
  stop("Missing required libraries");
}

bioc.libraries.notinstalled <- setdiff(bioc.libraries.required, rownames(installed.packages()));
if (length(bioc.libraries.notinstalled)>0){
  print(paste0("Missing Bioconductor libraries to run this program: ", paste0(bioc.libraries.notinstalled, collapse=", ")));
  if(readline(prompt="Install? y/n ") == 'y'){
    print("Installing Bioconductor libraries...");
    source("https://bioconductor.org/biocLite.R");
    biocLite(bioc.libraries.notinstalled);
  }
}
bioc.libraries.notinstalled <- setdiff(bioc.libraries.required, rownames(installed.packages()));
if (length(libraries.notinstalled)>0){
  print("Failed to install Bioconductor libraries. Terminating program.");
  rm(bioc.libraries.required, bioc.libraries.notinstalled);
  stop("Missing required Bioconductor libraries");
}

# Load libraries ----------------------------------------------------------

for (loadlib in libraries.required){
  print(paste0("Importing ",loadlib,"..."));
  library(loadlib,character.only=TRUE);
}
for (loadlib in bioc.libraries.required){
  print(paste0("Importing Bioconductor ",loadlib,"..."));
  library(loadlib,character.only=TRUE);
}

#clean up
rm(libraries.required, libraries.notinstalled, bioc.libraries.required, bioc.libraries.notinstalled, loadlib, ini.path, readini);
