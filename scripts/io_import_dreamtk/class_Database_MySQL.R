

# MySQL connector class ---------------------------------------------------


# v0.9

#requires R6 class, tidyverse, logging, DBI, RMySQL packages

#creates a mysql connection object
#wrapper around RMySQL/DBI interface


Class.Database.MySQL <- R6Class("Class.Database.MySQL",
                                
  lock_class = TRUE,
  cloneable = FALSE,
  
  #private variables and functions
  private = list(
    host = NULL,
    port = NULL,
    user = NULL,
    password = NULL,
    db = NULL,
    save = NULL,
    mysql_connection = NULL, #RMySQL object
    instance_name = NULL,
    
    
    connectionObjectExists = function(){
      if (!is.null(private$mysql_connection)){
        return (TRUE);
      } else {
        return (FALSE);
      }
    },
    
    loginInfoHasBeenSaved = function(){
      if (!is.null(private$save) && private$save == "1"){
        return (TRUE);
      } else {
        return (FALSE);
      }

    },
    
    shouldSaveLoginInfo = function( login_info ){
      if(login_info$save == "1"){
        return (TRUE);
      } else {
        return (FALSE);
      }
    },
    
    saveLoginInfo = function( login_info ){
        private$host <- login_info$host;
        private$port <- login_info$port;
        private$user <- login_info$user;
        private$password <- login_info$password;
        private$db <- login_info$db;
        private$save <- login_info$save;
    },
    
    eraseLoginInfo = function(){
      private$host <- NULL;
      private$port <- NULL;
      private$user <- NULL;
      private$password <- NULL;
      private$db <- NULL;
      private$save <- NULL;
      gc(); #call garbage collector immediately
    },
    
    attemptDBConnection = function( login_info ){
      tryCatch(
        private$mysql_connection <- dbConnect(MySQL(),
                                              user = login_info$user, 
                                              password = login_info$password, 
                                              host = login_info$host, 
                                              dbname = login_info$db,
                                              port = login_info$port),
        error = function(e) { 
          logerror(paste0("MySQL Connection error in mysql$dbConnect(): ", e)); 
          private$mysql_connection <- NULL; 
          #private$eraseLoginInfo(); 
          stop("MySQL Connection error, app can not proceed.");
        } 
      );

    }
    
    
  ),
  
  #public variables and functions
  public = list(
    classname = NULL, #we want to use this to id class instead of typeof(object) which returns "environment" for all R6 classes
    
    #constructor
    initialize = function( instance_name = "MySQL" ) {
      private$instance_name <- instance_name;
      self$classname = "Database.MySQL";
    },
    
    #finalizer
    finalize = function() {
      # if(!is.null(private_mysqlconnection)){
      #   self$dbDisconnect();
      # }
    },
    
    #
    setLoginInfo = function(host = NULL, user = NULL, password = NULL, db = NULL, save = NULL, port = 3306){
        private$host <- host;
        private$port <- port;
        private$user <- user;
        private$password <- password;
        private$db <- db;
        private$save <- save;
    },
    
    #connect to DB
    dbConnect = function() {
      
      if ( private$loginInfoHasBeenSaved() ){
        
        login_info <- list(user = private$user, 
                           password = private$password, 
                           host = private$host, 
                           db = private$db,
                           port = private$port);
        private$attemptDBConnection( login_info );
        
      } else {
        
        # login_window <- Class.GUI.LoginWindow$new(str_c("MySQL Login: ", private$instance_name));
        # login_info <- login_window$getInfo(); #expected list with keys: user, password, host, db, save
        # 
        # if ( private$shouldSaveLoginInfo( login_info ) ){
        #   private$saveLoginInfo( login_info );
        # }
        # 
        # private$attemptDBConnection( login_info );
        
        #info has not been saved, we won't have it, can't connect
        logerror("MySQL Connection error, no MySQL login information specified, app can not proceed."); 
        stop("MySQL Connection error, no MySQL login information specified, app can not proceed.");
      }
      
      if( private$connectionObjectExists() ){
        logdebug("MySQL connection established.");
        return (TRUE);
      } else {
        return (FALSE);
      }

    },
    
    #disconnect from DB
    dbDisconnect = function(){
      
      if ( private$connectionObjectExists() ) {
        result <- FALSE;
        tryCatch (
          result <- dbDisconnect(private$mysql_connection),
          error = function(e) { logerror(paste0("MySQL disonnection error in mysql$dbDisconnect(): ", e)); } 
          );
        
        private$mysql_connection <- NULL;
        logdebug("MySQL connection removed.");
        return (result);
        
      }else{
        logdebug("mysql$dbDisconnect(): Attempted to remove MySQL connection which has already been removed.");
        return (FALSE);
      }
      
    },
    
    #check if connection is valid, although beware that dbIsValid from RMySQL is not super reliable...
    dbIsValid = function(){
      if ( private$connectionObjectExists() ){
        
        result <- FALSE;
        tryCatch( 
          result <- dbIsValid(private$mysql_connection), 
          error = function(e) { logerror(paste0("MySQL Connection validity error in mysql$dbIsValid(): ", e)); }
        );
        
        return(result);
        
      } else {
        logdebug("mysql$dbIsValid(): MySQL connection does not exist.");
        return (FALSE);
      }
      
    },
    
    #return connection information
    dbGetInfo = function(){
      if ( private$connectionObjectExists() ){
        tryCatch(
          dbGetInfo(private$mysql_connection),
          error = function() logerror(paste0("MySQL Connection info error in mysql$dbGetInfo(): ", e))
        );
      }else {
        logwarn("mysql$dbGetInfo(): MySQL connection does not exist. No info.");
      }
    },
    
    #return RMySQL DB connection object
    getDB = function(){
      return(private$mysql_connection);
    },
    
    #list tables
    dbListTables = function(){
      if ( private$connectionObjectExists() ){
        result <- list();
        tryCatch (
          result <- dbListTables(private$mysql_connection),
          error = function(e) { logerror(paste0("MySQL Connection table-listing error in mysql$dbListTables(): ", e, " Returning empty table list.")); } 
        );
        
        return (result);
        
      }else{
        logwarn("mysql$dbListTables(): MySQL connection does not exist. Returning empty table list.");
        return (list());
      }
    },
    
    #list fields in table belonging to the database
    dbListFields = function(table = NULL){
      if ( private$connectionObjectExists() && !is.null(table)){
        result <- list();
        tryCatch (
          result <- dbListFields(private$mysql_connection, table),
          error = function(e) { logerror(paste0("MySQL Connection field-listing error in mysql$dbListFields(): ", e, " Returning empty table field list.")); } 
        );
        
        return (result);
        
      }else{
        logwarn("mysql$dbListFields(): Either MySQL connection does not exist or no table requested. Returning empty table field list.");
        return (list());
      }
    },
    
    #read table, return dataframe
    dbReadTable = function(table = NULL){
      if ( private$connectionObjectExists() && !is.null(table)){
        result <- tibble();
        tryCatch (
          result <- as_tibble(dbReadTable(private$mysql_connection, table)),
          error = function(e) { logerror(paste0("MySQL Connection table-read error in mysql$dbReadTable(): ", e, " Returning empty tibble.")); } 
        );
        
        return (result);
        
      }else{
        logwarn("mysql$dbReadTable(): Either MySQL connection does not exist or no table requested. Returning empty tibble.");
        return (tibble());
      }
    },
    
    #send query, return query object or NULL if connection or query are invalid
    dbSendQuery = function(sql_query = ""){
      if ( private$connectionObjectExists() && !is.null(sql_query) && is.character(sql_query)){
        result <- NULL;
        tryCatch (
          result <- dbSendQuery(private$mysql_connection, sql_query),
          error = function(e) { logerror(paste0("MySQL Connection send-query error in mysql$dbSendQuery(): ", e, " Returning NULL. SQL query was: ", sql_query));
            stop("Failed to query MySQL server, which is a problem. Check connection, restart app."); 
            } 
        );
        
        return (result);
        
      }else{
        logwarn("mysql$dbSendQuery(): Either MySQL connection does not exist or sql_query is invalid. Returning NULL.");
        stop("Failed to query MySQL server, which is a problem. Check connection, restart app.");
        return (NULL);
      }
    },
    
    #fetch, return data based on query object, specify number of rows to get
    dbFetch = function(query = NULL, n = 10000){
      if ( private$connectionObjectExists() && !is.null(query) && is.numeric(n)){
        result <- tibble();
        tryCatch (
          result <- dbFetch(query, n),
          error = function(e) { logerror(paste0("MySQL Connection query-fetch error in mysql$dbFetch(): ", e, " Returning empty tibble.")); } 
        );
        
        return (result);
        
      }else{
        logwarn("mysql$dbFetch(): Either MySQL connection does not exist or query invalid. Returning empty tibble.");
        return (tibble());
      }
    },
    
    #get query, return query object or NULL if connection or query are invalid
    dbGetQuery = function(sql_query = ""){
      if ( private$connectionObjectExists() && !is.null(sql_query) && is.character(sql_query)){
        result <- tibble();
        tryCatch (
          result <- dbGetQuery(private$mysql_connection, sql_query),
          error = function(e) { logerror(paste0("MySQL Connection get-query error in mysql$dbGetQuery(): ", e, " Returning empty tibble. SQL query was: ", sql_query)); } 
        );
        
        return (result);
        
      }else{
        logwarn("mysql$dbGetQuery(): Either MySQL connection does not exist or sql_query is invalid. Returning empty tibble.");
        return (tibble());
      }
    },
    
    #clear data associated with query object
    dbClearResult = function(query = NULL){
      if ( private$connectionObjectExists() && !is.null(query)){
        result <- FALSE;
        tryCatch (
          result <- dbClearResult(query),
          error = function(e) { logerror(paste0("MySQL Connection clear-query error in mysql$dbClearResult(): ", e)); } 
        );
        
        return (result);
        
      }else{
        logwarn("mysql$dbClearResult(): Either MySQL connection does not exist or query object has not been given.");
        return (FALSE);
      }
    },
    
    #write table; give db_table = string name of table to create, and df dataframe to write
    dbWriteTable = function(db_table = NULL, df = NULL, overwrite = FALSE, append = FALSE){
      if ( private$connectionObjectExists() && 
          !is.null(db_table) && 
          is.character(db_table) && 
          !is.null(df) && 
          is.data.frame(df) &&
          !((overwrite == append) == TRUE)){
        
        result <- FALSE;
        tryCatch (
          result <- dbWriteTable(private$mysql_connection, db_table, df, overwrite = overwrite, append = append),
          error = function(e) { logerror(paste0("MySQL Connection table-write error in mysql$dbWriteTable(): ", e, " Table to be written: ", db_table)); } 
        );
        
        return (result);
        
      }else{
        logwarn("mysql$dbWriteTable(): Connection or table writing parameters invalid. Table not written.");
        return (FALSE);
      }
    },
    
    #deletes a table belonging to the database
    dbRemoveTable = function(table = NULL){
      if ( private$connectionObjectExists() && !is.null(table) && is.character(table)){
        result <- FALSE;
        tryCatch (
          result <- dbRemoveTable(private$mysql_connection, table),
          error = function(e) { logerror(paste0("MySQL Connection table-remove error in mysql$dbRemoveTable(): ", e, " Table to be removed: ", table)); } 
        );
        
        return (result);
        
      }else{
        logwarn("mysql$dbRemoveTable(): Either MySQL connection does not exist or no table requested.");
        return (FALSE);
      }
    },
    
    #check if database exists
    dbExists = function(database = NULL){
      if ( private$connectionObjectExists() && !is.null(database) && typeof(database) == "character"){
        
        #get db list
        db_list <- self$dbGetQuery("SHOW DATABASES");

        if( all(database %in% db_list$Database) ){
          return (TRUE);
        } else {
          return (FALSE);
        }
      } else {
        logwarn("mysql$dbExists(): Either MySQL connection does not exist or no database requested.");
        return (FALSE);
      }
      
    },
    
    #check if a table or list of tables exists in the currently used database
    dbTablesExist = function ( tables = NULL ){
      if ( private$connectionObjectExists() && !is.null(tables)){
        
        dbtables <- self$dbListTables();
        tables_exist <- tables %in% dbtables;
        if ( all(tables_exist) ) {
          return (TRUE);
        }else {
          logerror( str_c("Missing tables: ", str_c(unlist ( tables[which(!tables_exist)] ), collapse = " ")) );
          return (FALSE);
        };
        
      }
    },
    
    #check if table fields exist for a given table in the currently used database
    dbFieldsExist = function ( table, fields ){
      if ( private$connectionObjectExists() && !is.null(table) && !is.null(fields) ){
      
        dbtable <- self$dbListFields( table );
        fields_exist <- fields %in% dbtable;

        if( all(fields_exist) ) {
          return (TRUE);
        }else{
        
          logerror( str_c("Missing table fields: ", str_c(unlist ( fields[which(!fields_exist)] ), collapse = " ")) );
          return (FALSE);
          
        };
      
      }
    },
    
    #use a different database
    dbSetActiveDB = function(database = NULL){
      if ( private$connectionObjectExists() && !is.null(database) && typeof(database) == "character"){
        result <- TRUE;
        tryCatch(
          { q <- self$dbSendQuery( paste0("USE ", database));
            self$dbClearResult(q); },
            error = function(e) { logerror("Unable to set required database in $dbSetActiveDB(): ", e); result <- FALSE; }
          );
        return (result);
      } else {
        logwarn("mysql$dbSetActiveDB(): Either MySQL connection does not exist or no database requested.");
        return (FALSE);
      }
    },
    
    #create a new database
    dbCreate = function(database = NULL){
      if ( private$connectionObjectExists() && !is.null(database) && typeof(database) == "character"){
        
        #if DB does not exist, create, otherwise let us know
        if( !self$dbExists( database ) ) {
          
          #create DB
          #q returns nothing upon dbFetch(q)
          q <- self$dbSendQuery( paste0("CREATE DATABASE ", database));
          self$dbClearResult(q);
          
          #check successful creation of database
          if ( self$dbExists( database ) ) {
            logdebug( paste0("Databse ", database, " created.") );
            return (TRUE);
          } else {
            logerror( paste0("Failed to create database ", database, ".") );
            return (FALSE);
          }

        } else {
          
          logdebug( paste0("Databse ", database, " already exists.") );
          return (FALSE);
          
        }
      } else {
        logwarn("mysql$dbCreate(): Either MySQL connection does not exist or no database requested.");
        return (FALSE);
      }
    },
    
    #
    dbGetName = function(){
      if (is.null(private$db)) {
        return ("");
      } else {
        return (private$db);
      }
    }

  )
  
)
