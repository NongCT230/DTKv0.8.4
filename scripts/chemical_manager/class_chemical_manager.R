
# Data Class: Chemical Manager ------------------------------------------

# Class for maintaining a list of Chemical objects

# v0.7


Class.ChemicalManager <- R6Class(
  "Class.ChemicalManager",
  
  lock_class = TRUE,
  cloneable = FALSE,
  
  #self variables and functions
  private = list(
    
    numchems = NULL,
    
    #return a list of values for the varname from the chemical list
    getListOfChemicalVarsBy = function( varname ){
      clist <- lapply(self$chemicals, 
                        function(x) { return (x$chem_info[[varname]]); }
      );
      return ( clist );
    },
    
    #
    removeNullsFromList = function(list) {
      Filter(Negate(is.null), list);
    },
    
    #sort chemical list by any array of values given they are of length of self$chemlist
    sortChemicalsByValueList = function ( val_list ){
      if (!is.null(val_list)) {
        df <- tibble( names = val_list, chems = self$chemicals );
        df <- arrange(df, names );
        self$chemicals <- pull(df, chems);
      } else {
        logdebug("ChemicalManager$sortChemicalsByValueList(): Nothing to sort.");
      }
    },
    
    #
    contains = function( chemical ) {
      chemlist <- unlist ( private$getListOfChemicalVarsBy("casn") , use.names = FALSE );
      
      if ( is.null(chemlist) || !(chemical$chem_info$casn %in% chemlist) ){
        return (FALSE);
      } else {
        return (TRUE);
      }
    },
    
    #
    removeInvalidChemicalsFromList = function ( chemicals ){
      
      validchems <- unlist( 
        lapply( chemicals, function(x) self$isValidChemical(x) ),
        use.names = FALSE
        );
      
      if (length(validchems)>0) {
        chemlist <- chemicals[which(validchems)];
      } else {
        chemlist <- list();
      }
      return (chemlist);
      
    }
    
  ),
  
  #public variables and functions
  public = list(
    
    chemicals = NULL,
    
    #constructor
    initialize = function(){
      self$chemicals <- list();
      private$numchems <- 0;
    },
    
    #finalizer
    finalize = function() {},
    
    #
    isValidChemical = function( chemical ){
      if (is.null(chemical) || is.atomic( chemical ) || 
          is.null(chemical$classname) || chemical$classname != "Class.Chemical" ||
          is.null(chemical$chem_info) || is.null(chemical$chem_info$casn) || is.na(chemical$chem_info$casn)){
        return (FALSE);
      }else{
        return (TRUE);
      }
    },
    
    #
    addChemical = function ( chemical = NULL ){
     if ( self$isValidChemical( chemical ) ){

          if ( !private$contains(chemical) ){
            #self$chemicals <- append(self$chemicals, chemical); the way below way seems faster with microbenchmark()
            self$chemicals <- list(self$chemicals, chemical);
            self$chemicals <- unlist ( self$chemicals, use.names = FALSE );
            logdebug( str_c("ChemicalManager$addChemical(): Added: ", chemical$chem_info$casn) );
            private$numchems <- private$numchems + 1;
            invisible(2); #return out of bounds index
          } else {
            logwarn (str_c("ChemicalManager$addChemical(): Can't add a duplicate or NA chemical: ", chemical$chem_info$casn));
            invisible(1); #return duplicate index
          }

     } else {
       e <- "ChemicalManager$addChemical(): Can't add an object that is not a valid Chemical or missing casn!"
       logerror(e);
       stop(e);
     }
    },
    
    #a faster way to add multiple chemicals
    addListOfChemicals = function ( chemicals = NULL){
      if (!is.null(chemicals) && is.list( chemicals )){
        
        chemicals <- private$removeInvalidChemicalsFromList( chemicals );
        
        oldlength <- length(self$chemicals);
        self$chemicals <- list(self$chemicals, chemicals);
        self$chemicals <- unlist ( self$chemicals, use.names = FALSE );
        logdebug( str_c("ChemicalManager$addChemical(): Added valid Chemicals from list.") );
        private$numchems <- private$numchems + length(chemicals);
          
        duplicates <- which(duplicated( self$getListOfChemicalCasn() ));
        if (length(duplicates) > 0) {
          self$chemicals <- self$chemicals[ -duplicates ];
          logdebug( str_c("Removed ", length(duplicates), " duplicate chemical entries.") );
          private$numchems <- private$numchems - length(duplicates);
          invisible(duplicates - oldlength); #return indices of duplicate items
        } else {
          invisible(length(chemicals)+1); #return out of bounds index
        }


      } else {
        e <- "ChemicalManager$addChemical(): Can't add an object that is not a Chemical list!"
        logerror(e);
        stop(e);
      }
    },
    
    #
    removeChemicals = function ( casn = NULL ){
      #casn can be a list
      if (!is.null(casn) && is.character(casn)){
        
        if ("*" %in% casn) { #remove all
          
          logdebug(paste0( length(self$chemicals), " chemicals removed."));
          self$chemicals <- list();
          private$numchems <- 0;
          
        } else { #remove by casn list
          
          indices <- which( self$getListOfChemicalCasn() %in% casn);
          if(length(indices) > 0){
            self$chemicals <- self$chemicals[-indices];
            private$numchems <- private$numchems - length(indices);
          }
          logdebug(paste0( length(indices), " chemicals removed"));
        }
        
      } else {
        e <- "ChemicalManager$removeChemicals(): No character vector for casn provided. Nothing removed!"
        logwarn(e);
      }
    },
    
    #
    getListOfChemicals = function ( casn = "*" ){
      #casn can be a list
      if (!is.null(casn) && is.character(casn)){
        
        if ("*" %in% casn) { #return all
          return (self$chemicals);
        } else {
          clist <- lapply(self$chemicals, 
                          function(x) { 
                            if ( all(x$chem_info[["casn"]] %in% casn) ) return (x); 
                            }
          );
          return ( private$removeNullsFromList( clist ) );
        }
      } else {
        e <- "ChemicalManager$getListOfChemicals(): No character vector for casn provided. Empty list returned!"
        logwarn(e);
        return(list());
      }

    },
    
    #
    getListOfChemicalsByName = function ( names = "*" ){
      #casn can be a list
      if (!is.null(names) && is.character(names)){
        
        if ("*" %in% names) { #return all
          return (self$chemicals);
        } else {
          clist <- lapply(self$chemicals, 
                          function(x) { 
                            if ( all(x$chem_info[["name"]] %in% names) ) return (x); 
                          }
          );
          return ( private$removeNullsFromList( clist ) );
        }
      } else {
        e <- "ChemicalManager$getListOfChemicalsByName(): No character vector for names provided. Empty list returned!"
        logwarn(e);
        return(list());
      }
      
    },
    
    #
    getChemical = function (casn = NULL){
      if (length(casn) > 1) {
        logwarn("ChemicalManager$getChemical(): More than one chemical requested. Returning NULL.");
        return (NULL);
      } else {
        chemlist <- self$getListOfChemicals( casn );
        if (length(chemlist)>0){
          return ( chemlist[[1]] );
        } else {
          return (NULL);
        }
      }
    },
    
    #
    getChemicalByName = function (name = NULL){
      if (length(name) > 1) {
        logwarn("ChemicalManager$getChemicalByName(): More than one chemical requested. Returning NULL.");
        return (NULL);
      } else {
        chemlist <- self$getListOfChemicalsByName( name );
        if (length(chemlist)>0){
          return ( chemlist[[1]] );
        } else {
          return (NULL);
        }
      }
    },
    
    #
    getListOfChemicalNames = function (){
      return ( private$getListOfChemicalVarsBy("name") );
    },
    
    #
    getListOfChemicalCasn = function (){
      return ( private$getListOfChemicalVarsBy("casn") );
    },
    
    #
    sortChemicalsByName = function (){
      
      chemnames <- unlist ( self$getListOfChemicalNames(), use.names = FALSE );
      
      private$sortChemicalsByValueList ( chemnames );
      
    },
    
    #
    sortChemicalsByCasn = function(){
     
      caslist <- unlist ( self$getListOfChemicalCasn(), use.names = FALSE );
      
      private$sortChemicalsByValueList ( caslist );
       
    },
    
    #
    getNumChemicals = function(){
      return (private$numchems);
    }
    
  )
  
)
