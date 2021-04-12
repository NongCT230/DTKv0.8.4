
# Function imports ---------------------------------------------------------

# v0.2

# imports *.R files from all of the folders that we want, as indicated in app.ini

# alternative version only imports the files indicated in "include.ini" filewithin the folders
# uncommentlines 23-27, 49; comment lines 28, 50

# List required fuctions --------------------------------------------------

#set app.ini filepath to open, create list of import packages
source(paste0(working.dir,"/_readini_.R"));
if (!exists("ini.path")) ini.path <- paste0(working.dir,"/app.ini");
function.packages.required <- readini(ini.path, "[Import script packages]");

# Load required functions into global environment -------------------------

for (package.dir in function.packages.required){
  include.dir <- str_c(working.dir, str_replace(package.dir, "./", "/"));
  # include.file <- str_c(include.dir,"/include.ini");
  # if (file.exists(include.file)){
  #   loginfo(str_c("Accessing ",include.file,"..."));
    
    # file.list <- readini(include.file,"[Script package files]");
    file.list <- list.files(include.dir,"*.R");
	
    for (file in file.list){
      file.path <- str_c(include.dir, "/", str_replace(file, "\"", ""));
      if (file.exists(file.path)){
        logdebug(str_c("Imported ",file.path));
        source(file.path);
      } else {
        logerror(str_c(file.path, " does not exist. Terminating..."));
        rm(include.files, info, include.dir, file, file.path);
        stop("Missing required function package files.");
      }
    }
    
  # } else {
  #   logerror(str_c("No include file found in ", package.dir, ". Function package not loaded. Terminating..."));
  #   rm(function.packages.required, package.dir, include.file, include.dir);
  #   stop("Missing required function package files.");
  # }
}

#cleanup memory
# rm(function.packages.required, package.dir, include.file, include.dir, ini.path, readini, file, file.list, file.path);
rm(function.packages.required, package.dir, include.dir, ini.path, readini, file, file.list, file.path);

