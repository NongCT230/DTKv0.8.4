#open file if exists, read required info
#using an inefficient method to add list items since list is small

readini = function(ini.path, subheading){
  
  output.list <- c();
  
  if (file.exists(ini.path)){
    
    ini.file <- file(ini.path, "r");
    
    while (TRUE){
      
      #keep reading lines till end of file or got all we need
      line <- readLines(ini.file, n=1);
      trimws(line, which = "both");
      
      if(length(line) > 0){

        if(line == subheading){

          while(TRUE){
            line <- readLines(ini.file, n=1);
            line <- trimws(line, which = "both");
            if(length(line) > 0 && substr(line,1,1) == "\""){
              output.list <- c(output.list, substr(line,2,nchar(line)-1));
            }else{
              break;
            }
          }
          break;
          
        }
        
      }else{ #reached EOF
        
        break;
        
      }
      
    }
    
    close(ini.file);
    
    return(output.list);

  } else {
    stop("app.ini missing in the application directoy. Unable to import required packages!");
  }
}
