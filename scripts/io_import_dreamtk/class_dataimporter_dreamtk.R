
# Data Importer class for dreamtk --------------------------------------------

# v0.7

#table variable storage container for the core importer classes
# - modify this to change variables being imported (*tablevars) from which table (*tablename)
# - modifying *filter* variable names for importer objects will require importer function modification as well
# filter vars:
# -- 'casn' and 'name' in pretty much many of the imported tables
# -- 'aeid' in the assay table
# -- 'parameter' in the physiology table
# !!!!!!!!!!!!!!!!!!!!!!!!THIS ONLY AFFECTS THE CORE IMPORTER CLASSES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Anything else downstream, such as models and GUI that uses the imported data, will need to be changed as well
Class.DataImporter.dreamtk.tableinfo <- new.env();
tableinfo <- list(
  chemtablename = "chemical_data",
  chemtablevars = list ("casn", "name", "cytotoxicity_um", "cytotoxic"),
  
  assaytablename = "assay_data",
  assaytablevars = list ("aeid", "assay_name", "assay_component_name", "assay_component_endpoint_name",
                         "organism", "tissue", "cell_short_name", "biological_process_target", 
                         "intended_target_family", "intended_target_family_sub", "gene_name", "gene_symbol"),
  
  ac50tablename = "ac50_data",
  ac50tablevars = list ("casn", "aeid", "ac50", "ac_top", "ac_cutoff", "hitc","burst_assay"),
  
  physiologytablename = "physiology_data",
  physiologytablevars = list ("parameter", "species", "value", "units", "reference"),
  
  model3csstablename = "model_3css_data",
  model3csstablevars = list ("casn", "mw", "human_funbound_plasma", "human_clint", "human_clint_unit", 
                     "human_clint_p", "human_rblood2plasma", "mouse_funbound_plasma", "rabbit_funbound_plasma",
                     "rat_funbound_plasma", "rat_clint", "rat_clint_unit", "rat_clint_p", "rat_rblood2plasma"),
  
  kowtablename = "kow_data",
  kowtablevars = list ("casn", "log_kow"),
  
  pkatablename = "pka_data",
  pkatablevars = list ("casn", "pka"),
  
  cpdattablename = "cpdat_data",
  cpdattablevars = list("casn","product_use_category","minimum_reported_weight_fraction","maximum_reported_weight_fraction","predicted_weight_fraction_mean","predicted_weight_fraction_05th_percentile", "predicted_weight_fraction_50th_percentile", "predicted_weight_fraction_95th_percentile"),
  
  shedsinfotablename = "shedsinfo_data",
  shedsinfotablevars = list("product_use_category","frequency","mass_per_use","males","duration_of_direct_use","fcont","fret","fing","direct_dermal","direct_inhalation_of_vapor","direct_inhalation_of_aerosol","direct_incidental_ingestion","indirect"),

  #YANIC
  #TODO check this out
  gentox21dataname = "gen_tox_data",
  gentox21datavars = list("casn", "compound", "assay", "endpoint", "s9_positive", "concentration", "concentration_unit", "response", "response_type"),
  
  bmddataname = "bmd_data",
  bmddatavars = list("casn", "compound", "assay", "endpoint", "s9_positive", "model_type", "aic", "log_likelihood", "bmr", "bmdl", "bmd", "bmdu")
  
)
#put the table information into an environment which can be passed around by reference
list2env(tableinfo, envir = Class.DataImporter.dreamtk.tableinfo);
rm(tableinfo);
      
#base class for data importer objects
#defines general table import functions that can be used to check and import database tables
Class.DataImporter.dreamtk <- R6Class("Class.DataImporter.dreamtk",
  
  cloneable = FALSE,
  lock_class = TRUE,
                                   
  #self variables and functions
  private = list(
    
    #
    loadDataFromMySQLTable = function (dreamtk.db = NULL, filter_by = NULL, 
                                       filter_val = NULL, dbName = NULL, 
                                       table = NULL, tablevars = NULL){
      
      if (is.null(dreamtk.db) || dreamtk.db$classname != "Database.MySQL"){
        e <- "dataimporter.dreamtk$loadDataFromMySQLTable(): Invalid database object given.";
        logerror( e );
        stop(e);
      }
      if (is.null(dbName)){
        e <- "dataimporter.dreamtk$loadDataFromMySQLTable(): Invalid database name given.";
        logerror( e );
        stop(e);
      }
      if (is.null(table)){
        e <- "dataimporter.dreamtk$loadDataFromMySQLTable(): Invalid table name given.";
        logerror( e );
        stop(e);
      }
      if (is.null(tablevars)){
        e <- "dataimporter.dreamtk$loadDataFromMySQLTable(): Invalid table field names given.";
        logerror( e );
        stop(e);
      }
      if (is.null(filter_by)){
        e <- "dataimporter.dreamtk$loadDataFromMySQLTable(): Invalid filter field name given.";
        logerror( e );
        stop(e);
      }
      if (is.null(filter_val)){
        e <- "dataimporter.dreamtk$loadDataFromMySQLTable(): Invalid filter list given. Use '*' to import all values.";
        logerror( e );
        stop(e);
      }
      
      dreamtk.db$dbConnect();
      on.exit(dreamtk.db$dbDisconnect());
      
      #ensure we are set to required database
      dreamtk.db$dbSetActiveDB(dbName);
        
      #check if database table data is valid 
      logdebug( str_c("DataImporter.dreamtk is Verifying ", dbName, " data...") );
        
      if ( !dreamtk.db$dbTablesExist( table ) ){
        e <- str_c("Missing ", table, " table for dreamtk v1.0. Check database version.");
        logerror( e );
        stop(e);
      }
          
      if ( !dreamtk.db$dbFieldsExist( table, tablevars ) ){
        e <- str_c("Missing expected table values for ", table, ". Check database version.");
        logerror( e );
        stop(e);
      }
            
      #copy the required info for given chemicals into tables if valid
      #if one of the cas numbers is *, grab everything, otherwise grab what is in filter_val
      if (all("*" %in% filter_val)){
              
        query <- dreamtk.db$dbSendQuery(
          str_c("select ", str_c(tablevars,sep="",collapse=", "), " from ", table)
        );
              
      }else{
          #Gabriel Switch case to allow diff stuff in diff tables    
        query <- dreamtk.db$dbSendQuery(
          str_c("select ", 
                str_c(tablevars,sep="",collapse=", "), 
                " from ", table,
                " where ", filter_by, 
                " in ('", str_c(filter_val,sep=" ",collapse="','"), "')")
        );
              
      }
            
      imported_data <- dreamtk.db$dbFetch(query, n = -1); #-1 is the value to retrieve all which is better than an arbitrary large value.
      dreamtk.db$dbClearResult(query);
      
      logdebug("Dreamtk chem data imported successfully.");
      return (imported_data);
            
    },
    
    #
    loadDataFromRData = function ( file.path = NULL, table = NULL, 
                                    tablevars = NULL, filter_by = NULL, 
                                    filter_val = NULL ){

      
	  
	  if (is.null(file.path) || !file.exists(file.path)){
        e <- str_c( "dataimporter.dreamtk$loadDataFromRData(): Invalid file object given.", str_c(file.path, " can not be loaded.") );
        logerror( e );
        stop(e);
      }
	  
	  #Gabriel Make this part in a separate function.
      if (is.null(table)){
        e <- "dataimporter.dreamtk$loadDataFromRData(): Invalid table name given.";
        logerror( e );
        stop(e);
      }
      if (is.null(tablevars)){
        e <- "dataimporter.dreamtk$loadDataFromRData(): Invalid table field names given.";
        logerror( e );
        stop(e);
      }
      if (is.null(filter_by)){
        e <- "dataimporter.dreamtk$loadDataFromRData(): Invalid filter field name given.";
        logerror( e );
        stop(e);
      }
      if (is.null(filter_val)){
        e <- "dataimporter.dreamtk$loadDataFromRData(): Invalid filter list given. Use '*' to import all values.";
        logerror( e );
        stop(e);
      }
      
      #load and verify tables
      for (file in file.path){
        load(file);
      }
      
      if(!exists(table)){
        e <- "dataimporter.dremtk$loadChemicalDataFromRData(): Failed to load expected RData. Check file and version.";
        logerror( e );
        stop(e);
      }
        
      #verify existence of required data columns
      vars_exist <- tablevars %in% names(eval(as.name(table)));
      if( !all(vars_exist) ) {
        e <- str_c( str_c(table, " missing expected values. Check file and version."),
                    str_c("Missing variables: ", str_c(unlist ( tablevars[which(!vars_exist)], use.names = FALSE ), collapse = " ")) 
                    );
        logerror( e );
        stop(e);
          
      }else {
        if ( all("*" %in% filter_val )) {
          #this filters data in ac50 to exclude the non hit assays.
          if(table == "ac50_data"){
            imported_data <- filter(eval(as.name(table))); 
          }else{

          imported_data <- eval(as.name(table));
          }
        } else {
			if(table == "ac50_data"){
				imported_data <- filter(eval(as.name(table))); 
			}else{
				imported_data <- filter(eval(as.name(table)), eval(as.name(filter_by)) %in% filter_val );
			}
        }
        logdebug( str_c("DataImporter.dreamtk found and loaded ", table, " data.") );
        return (imported_data);
      };

    },
    
    #importers only grab database data that is available in the databas
    #any search terms that do not have corresponding data are not imported, but may be something we want to include with NA values
    #so this function adds missing data rows, and since no data is available it is assigned NA values
    addMissing = function( data = NULL, filter_by = NULL, filter_val = NULL ){
      # we want to return lines for all values, including the ones we did not find in the import tables
      # fill the missing data with NA
      # this ensures that we can report the NA values during later stages of processing; they can always be filtered out later
      data_length <- length (  data[[filter_by]] );
      data_length_unique <- length (  unique( data[[filter_by]] ) );
      #print(data_length);
      #print(data_length_unique);
      #print(length(filter_val));
      missing <- c();
      if ( data_length == 0 ) {
        missing <- filter_val[filter_val != "*"];
      } else if ( data_length_unique < length(filter_val) ){
        present_indices <- which( (filter_val %in% unique( data[[filter_by]])) );
          if( length(present_indices) > 0 ){
            missing <- filter_val[ -present_indices ];
            missing <- missing[filter_val != "*"];
            missing <- missing[!is.na(missing)];
          }
      }
      #print(missing);
      missing_length <- length(missing);
      #print(missing_length);
      if ( missing_length == 0 ){
        return (data);
      }

      for (i in (data_length+1):(data_length+missing_length) ) {
        data <- add_row(data);
        data[[filter_by]][i] <- missing[i-data_length];
      }
      logwarn( str_c("DataImporter.dreamtk$addMissing(): Imported data missing requested values: ", paste0(missing, collapse=", "),
                     ". Missing value datarows have been added with all parameters set to NA.") );
      
      return (data);
    }
    
    
  ),
  
  #public variables and functions
  public = list(
    
    classname = NULL,
    
    #constructor
    initialize = function(){
      self$classname <- "Class.Dataimporter.dreamtk";
    },
    
    #finalizer
    finalize = function(){}
    
  )
  
)
