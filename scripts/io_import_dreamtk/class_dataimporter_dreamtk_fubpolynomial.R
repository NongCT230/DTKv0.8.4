
# Data Importer class for dreamtk --------------------------------------------

# v0.7

#fub polynomial model data importer


Class.DataImporter.dreamtk.FubPolynomial <- R6Class("Class.DataImporter.dreamtk.FubPolynomial",
  
  inherit = Class.DataImporter.dreamtk,
                                   
  #self variables and functions
  private = list(
    
    model_data_rdata = NULL,
    chemicals = NULL,
    
    model_data = NULL,
    source = NULL,
    lock_data = NULL,
    
    #
    loadModelDataFromMySQL = function ( dreamtk.db = NULL, chemicals = "*"){
      
      loginfo("Importing dreamtk v1.0 FubPolynomial model data from MySQL database...");
      table <- Class.DataImporter.dreamtk.tableinfo$kowtablename;
      tablevars <- Class.DataImporter.dreamtk.tableinfo$kowtablevars;
      
      kow_data <- private$loadDataFromMySQLTable( dreamtk.db=dreamtk.db,
                                                                 filter_by = "casn",
                                                                 filter_val = chemicals,
                                                                 dbName = "dream_tk", 
                                                                 table = table,
                                                                 tablevars = tablevars);
      
      table <- Class.DataImporter.dreamtk.tableinfo$pkatablename;
      tablevars <- Class.DataImporter.dreamtk.tableinfo$pkatablevars;
      
      pka_data <- private$loadDataFromMySQLTable( dreamtk.db=dreamtk.db,
                                                  filter_by = "casn",
                                                  filter_val = chemicals,
                                                  dbName = "dream_tk", 
                                                  table = table,
                                                  tablevars = tablevars);
      
      private$model_data <- full_join(kow_data, pka_data, by = "casn");
      
      invisible(self);
    },
    
    #
    loadModelDataFromRData = function ( file.path = NULL, chemicals = "*"){
      
      if (is.null(chemicals)){
        e <- "dataimporter.dreamtk.fubpolynomial$loadModelDataFromRData(): Invalid filter list given. Use '*' to get all values.";
        logerror( e );
        stop(e);
      }
      
      private$chemicals <- chemicals;
      
      if (is.null(private$model_data_rdata)){
      
        loginfo("Importing dreamtk v1.0 FubPolynomial model data from RData file...");
        table <- Class.DataImporter.dreamtk.tableinfo$kowtablename;
        tablevars <- Class.DataImporter.dreamtk.tableinfo$kowtablevars;
        
        kow_data <- private$loadDataFromRData( file.path=file.path,
                                                              filter_by = "casn",
                                                              filter_val = "*",
                                                              table = table,
                                                              tablevars = tablevars);
        
        table <- Class.DataImporter.dreamtk.tableinfo$pkatablename;
        tablevars <- Class.DataImporter.dreamtk.tableinfo$pkatablevars;
        
        pka_data <- private$loadDataFromRData( file.path=file.path,
                                               filter_by = "casn",
                                               filter_val = "*",
                                               table = table,
                                               tablevars = tablevars);
        
        private$model_data_rdata <- full_join(kow_data, pka_data, by = "casn");
        
      }
      
      if (any("*" %in% private$chemicals)){
        private$model_data <- private$model_data_rdata;
      } else {
        private$model_data <- filter(private$model_data_rdata, casn %in% private$chemicals );
      }
      
      invisible(self);
    }

  ),
  
  #public variables and functions
  public = list(
    
    classname = NULL,
    
    #constructor
    initialize = function(){
      self$classname <- "Class.DataImporter.dreamtk.FubPolynomial";
      private$lock_data <- FALSE;
    },
    
    #finalizer
    finalize = function(){},
    
    #
    setSource = function ( source = NULL ){
      private$source = source;
    },
    
    #
    loadModelData = function ( chemicals = "*"){
      if (private$lock_data){
        w <- "dataimporter.dreamtk.fubpolynomial$loadModelData(): data not imported from source specified by $setSource() since it is locked by $setModelData(). Use $resetModelData() to unlock.";
        return ( logwarn( w ) );
      }
      
      #Removed MySQL database option #YANIC
      #if (is.character(private$source)){
        private$loadModelDataFromRData(private$source, chemicals);
      # }else{
      #   private$loadModelDataFromMySQL(private$source, chemicals);
      # }
    },
    
    #directly set model data, locks external importer functionality until unlocked by resetModelData()
    setModelData = function ( modeldata = NULL ){
      tablevars <- c(
        unlist( Class.DataImporter.dreamtk.tableinfo$kowtablevars ),
        unlist( Class.DataImporter.dreamtk.tableinfo$pkatablevars )
      ) %>% unique ();
      
      vars_exist <- tablevars %in% names(modeldata);
      if (all(vars_exist)){
        private$model_data <- modeldata;
        private$lock_data <- TRUE;
        w <- "dataimporter.dreamtk.fubpolynomial$setModelData(): data set and locked. Use $resetModelData() to unlock, to be able to $loadModelData() from source specified by $setSource().";
        logwarn( w );
      } else {
        e <- "dataimporter.dreamtk.fubpolynomial$setModelData(): Data missing required columns.";
        logerror( e );
        stop(e);
      }
    },
    
    #
    resetModelData = function (){
      private$lock_data <- FALSE;
      private$model_data <- NULL;
      w <- "dataimporter.dreamtk.fubpolynomial$resetModelData(): data reset and unlocked. $loadModelData() is now functional.";
      logwarn( w );
    },
    
    #
    getModelData = function (){
      return (private$model_data);
    }

    
  )
  
)
