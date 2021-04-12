
# Data Importer class for dreamtk --------------------------------------------

# v0.7

#model data importer for both SQL and R data.

Class.DataImporter.dreamtk.basicData <- R6Class("Class.DataImporter.dreamtk.basicData",
	inherit = Class.DataImporter.dreamtk,
	
	
	private = list(
	
		#important these are used because we load all the data at launch.
	    chemicals = NULL,
		a_chemicals = NULL,
		assays = NULL,
		params = NULL,
		
		chemical_data = NULL,
		assay_data = NULL,
		ac50_data = NULL,
		physiology_data = NULL,
		assay_info = NULL,
		
		cpdat_data = NULL,
		shedsinfo_data = NULL,
		
		#YANIC
		gentox21_data = NULL,
		bmd_data = NULL,
		
		source = NULL 
	
	),
	
	public = list(
		
		
		#constructor
		initialize = function(){
			super$initialize();
		},
		
		#finalizer
		finalize = function(){},
		
		
		#
		setSource = function ( source = NULL ){
			private$source = source;
			#Removing the MySQL option for database.. #YANIC
			#if (is.character(private$source)){
				loginfo("Importing from R data file ...");
			# }
			# else{
			# 	loginfo("Importing from mySQL database ...");
			# }
		},
		
		loadChemicalData = function ( chemicals = "*" ) {
			
			
			if (is.null(chemicals)){
				e <- "dataimporter.dreamtk.rdata$loadChemicalData(): Invalid filter list given. Use '*' to get all values.";
				logerror( e );
				stop(e);
			}

			private$chemicals <- chemicals;
			

			if ( is.null(private$chemical_data)){
				loginfo("Importing dreamtk v1.0 chemical data...");
				table <- Class.DataImporter.dreamtk.tableinfo$chemtablename;
				tablevars <- Class.DataImporter.dreamtk.tableinfo$chemtablevars;
				#if (is.character(private$source)){
					private$chemical_data <- private$loadDataFromRData ( file.path = private$source, 
																		table = table, 
																		tablevars = tablevars, 
																		filter_by = "casn", 
																		filter_val = "*" );
				# }else{
				# 	private$chemical_data <- private$loadDataFromMySQLTable( dreamtk.db=private$source,
				# 															filter_by = "casn",
				# 															filter_val = "*",
				# 															dbName = "dream_tk", 
				# 															table = table,
				# 															tablevars = tablevars);
				# }	
				private$chemical_data <- private$addMissing ( private$chemical_data, "casn", chemicals );
			}	
				invisible(self);
		},
	
	 #
		loadAssayData = function ( assays = "*" ) {
			if (is.null(assays)){
				e <- "dataimporter.dreamtk.rdata$loadAssayData(): Invalid filter list given. Use '*' to get all values.";
				logerror( e );
				stop(e);
			}
			
			private$assays <- assays;
						
			if ( is.null(private$assay_data)){
				loginfo("Importing dreamtk v1.0 assay data...");
				table <- Class.DataImporter.dreamtk.tableinfo$assaytablename;
				tablevars <- Class.DataImporter.dreamtk.tableinfo$assaytablevars;
				#if (is.character(private$source)){
					private$assay_data <- private$loadDataFromRData ( file.path = private$source, 
																	table = table, 
																	tablevars = tablevars, 
																	filter_by = "aeid", 
																	filter_val = "*" );
				# }else{
				# 	private$assay_data <- private$loadDataFromMySQLTable( dreamtk.db=private$source,
				# 													filter_by = "aeid",
				# 													filter_val = "*",
				# 													dbName = "dream_tk", 
				# 													table = table,
				# 													tablevars = tablevars);
				# }
			}
			invisible(self);
		},
		
		#
		loadAc50Data = function ( chemicals = "*" ) {
			
			
			if (is.null(chemicals)){
				e <- "dataimporter.dreamtk.rdata$loadAc50Data(): Invalid filter list given. Use '*' to get all values.";
				logerror( e );
				stop(e);
			}
			private$a_chemicals <- chemicals;

			if ( is.null(private$ac50_data)){
				loginfo("Importing dreamtk v1.0 ac50 data...");
				table <- Class.DataImporter.dreamtk.tableinfo$ac50tablename;
				tablevars <- Class.DataImporter.dreamtk.tableinfo$ac50tablevars;
				#if (is.character(private$source)){
					private$ac50_data <- private$loadDataFromRData ( file.path = private$source, 
																	table = table, 
																	tablevars = tablevars, 
																	filter_by = "casn", 
																	filter_val = "*" );
				
				# }else{
				# 	private$ac50_data <- private$loadDataFromMySQLTable( dreamtk.db=private$source,
				# 														filter_by = "casn",
				# 														filter_val = "*",
				# 														dbName = "dream_tk", 
				# 														table = table,
				# 														tablevars = tablevars);
				# }
			}
			invisible(self);
		},
		
		#
		loadPhysiologyData = function ( params ="*") {
			
			
			if (is.null(params)){
				e <- "dataimporter.dreamtk.rdata$loadPhysiologyData(): Invalid filter list given. Use '*' to get all values.";
				logerror( e );
				stop(e);
			}
			private$params <- params;
			

			if ( is.null(private$physiology_data)){
				loginfo("Importing dreamtk v1.0 physiology data...");	
				table <- Class.DataImporter.dreamtk.tableinfo$physiologytablename;
				tablevars <- Class.DataImporter.dreamtk.tableinfo$physiologytablevars;		
				#if (is.character(private$source)){
					private$physiology_data <- private$loadDataFromRData ( file.path = private$source, 
																		table = table, 
																		tablevars = tablevars, 
																		filter_by = "parameter", 
																		filter_val = "*" );
				# }else{
				# 	private$physiology_data <- private$loadDataFromMySQLTable( dreamtk.db=private$source,
				# 															filter_by = "parameter",
				# 															filter_val = "*",
				# 															dbName = "dream_tk", 
				# 															table = table,
				# 															tablevars = tablevars);
				# }
				private$physiology_data <- private$addMissing ( private$physiology_data, "parameter", params );
			}
			invisible(self);
		},
		
		#combination of ac50 data which has casn associated with aeid, and assay data which does not have casn but has aeid
		loadAssayInfo = function ( chemicals = "none" ) {

			
			if (is.null(chemicals)){
				e <- "dataimporter.dreamtk.rdata$loadAssayInfo(): Invalid filter list given. Use '*' to get all values.";
				logerror( e );
				stop(e);
			}
			

			private$a_chemicals <- chemicals;
			
			if(is.null(private$ac50_data)){
				loginfo("Importing dreamtk v1.0 assay info...");
				table <- Class.DataImporter.dreamtk.tableinfo$ac50tablename;
				tablevars <- Class.DataImporter.dreamtk.tableinfo$ac50tablevars;
				# if (is.character(private$source)){
					private$ac50_data <- private$loadDataFromRData ( file.path = private$source, 
													 table = table, 
													 tablevars = tablevars, 
													 filter_by = "casn", 
													 filter_val = "*" );
				# }else{
				# 	private$ac50_data <- private$loadDataFromMySQLTable( dreamtk.db=private$source,
				# 												filter_by = "casn",
				# 												filter_val = "*",
				# 												dbName = "dream_tk", 
				# 												table = table,
				# 												tablevars = tablevars);
				# }
			}
			
			if(is.null(private$assay_data)){
				table <- Class.DataImporter.dreamtk.tableinfo$assaytablename;
				tablevars <- Class.DataImporter.dreamtk.tableinfo$assaytablevars;
				# if (is.character(private$source)){
					private$assay_data <- private$loadDataFromRData ( file.path = private$source, 
																	table = table, 
																	tablevars = tablevars, 
																	filter_by = "aeid", 
																	filter_val = as.character(private$ac50_data$aeid) );
				# }else{
				# 	private$assay_data <- private$loadDataFromMySQLTable( dreamtk.db= private$source,
				# 											filter_by = "aeid",
				# 											filter_val = as.character(private$ac50_data$aeid),
				# 											dbName = "dream_tk", 
				# 											table = table,
				# 											tablevars = tablevars);
				# } 
				private$ac50_data$aeid <- as.numeric(private$ac50_data$aeid);
				private$assay_info <- inner_join(private$ac50_data, private$assay_data, by = "aeid");
			}
			invisible(self);
		},
		loadCPDatData = function ( chemicals = "*" ) {
			
			
			if (is.null(chemicals)){
				e <- "dataimporter.dreamtk.rdata$loadCPDatData(): Invalid filter list given. Use '*' to get all values.";
				logerror( e );
				stop(e);
			}
		

			if ( is.null(private$cpdat_data)){
				loginfo("Importing dreamtk v1.0 CPDat data...");
				table <- Class.DataImporter.dreamtk.tableinfo$cpdattablename;
				tablevars <- Class.DataImporter.dreamtk.tableinfo$cpdattablevars;
				#if (is.character(private$source)){
					private$cpdat_data <- private$loadDataFromRData ( file.path = private$source, 
																		table = table, 
																		tablevars = tablevars, 
																		filter_by = "*", 
																		filter_val = "*" );
				# }else{
				# 	private$cpdat_data <- private$loadDataFromMySQLTable( dreamtk.db=private$source,
				# 															filter_by = "*",
				# 															filter_val = "*",
				# 															dbName = "dream_tk", 
				# 															table = table,
				# 															tablevars = tablevars);
				# }	
			}	
				invisible(self);
		},
		
		loadShedData = function () {
		
			if ( is.null(private$shedsinfo_data)){
				loginfo("Importing dreamtk v1.0 Sheds data...");
				table <- Class.DataImporter.dreamtk.tableinfo$shedsinfotablename;
				tablevars <- Class.DataImporter.dreamtk.tableinfo$shedsinfotablevars;
				#if (is.character(private$source)){
					private$shedsinfo_data <- private$loadDataFromRData ( file.path = private$source, 
																		table = table, 
																		tablevars = tablevars, 
																		filter_by = "*", 
																		filter_val = "*" );
				# }else{
				# 	private$shedsinfo_data <- private$loadDataFromMySQLTable( dreamtk.db=private$source,
				# 															filter_by = "*",
				# 															filter_val = "*",
				# 															dbName = "dream_tk", 
				# 															table = table,
				# 															tablevars = tablevars);
				# }	
			}	
				invisible(self);
		},
		
		#TODO
		loadGenToxData = function() {
		  GT21Source <- "./data/GenTox21.RData";
		  if (is.null(private$gentox21_data)){
		    loginfo("Importing GenTox21 data...");
		    table <- Class.DataImporter.dreamtk.tableinfo$gentox21dataname;
		    tablevars <- Class.DataImporter.dreamtk.tableinfo$gentox21datavars;
		    
		    private$gentox21_data <- private$loadDataFromRData ( file.path = GT21Source, 
		                                                          table = table, 
		                                                          tablevars = tablevars, 
		                                                          filter_by = "*", 
		                                                          filter_val = "*" );
		  }
		  if (is.null(private$bmd_data)){
		    table <- Class.DataImporter.dreamtk.tableinfo$bmddataname;
		    tablevars <- Class.DataImporter.dreamtk.tableinfo$bmddatavars;
		    
		    private$bmd_data <- private$loadDataFromRData ( file.path = GT21Source, 
		                                                         table = table, 
		                                                         tablevars = tablevars, 
		                                                         filter_by = "*", 
		                                                         filter_val = "*" );
		  }
		  rm(GT21Source);
		},

		getChemicalData = function () {
		  if(any("*" %in% private$chemicals)){
			return( private$chemical_data );
		  }else{
			imported_data <- filter(private$chemical_data, casn %in% private$chemicals );
			imported_data <- private$addMissing ( imported_data, "casn", private$chemicals );
			return( imported_data );
		  }
		},
		
		#
		getAssayData = function () {
		  if(any("*" %in% private$assays)){
			return (private$assay_data);
		  }else {
			imported_data <- filter(private$assay_data, aeid %in% private$assays );
			return( imported_data );
		  }
		},
		
		#
		getAc50Data = function () {
		  if(any("*" %in% private$chemicals)){
			return (private$ac50_data);
		  }else{
			imported_data <- filter(private$ac50_data, casn %in% private$a_chemicals );
			return( imported_data );
		  }
		},
		
		#
		getPhysiologyData = function (){
		  if(any("*" %in% private$params)){
			return (private$physiology_data);
		  }else{
			imported_data <- filter(private$physiology_data, parameter %in% private$params );
			imported_data <- private$addMissing ( imported_data, "parameter", private$params );
			return( imported_data );
		  }
		},
		
		#
		getChemInfo = function () {
		  if(any("*" %in% private$chemicals)){
			return( private$chemical_data );
		  }else{
			imported_data <- filter(private$chemical_data, casn %in% private$chemicals );
			imported_data <- private$addMissing ( imported_data, "casn", private$chemicals );
			return( imported_data );
		  }
		},
		
		#
		getAssayInfo = function () {
		  if(any("*" %in% private$chemicals)){
			return( private$assay_info );
		  }else{
			imported_data <- filter(private$assay_info, casn %in% private$a_chemicals );
			return( imported_data );
		  }
		},
		
		getCpdatInfo = function () {
			return( private$cpdat_data );
		},
		
		getShedsInfo = function () {
			return (private$shedsinfo_data);
		},
		
		#YANIC
		getGenToxInfo = function() {
		  return (private$gentox21_data);
		},
		
		getBmdInfo = function() {
		  return (private$bmd_data);
		},
		
		#convert chem name list to casn, unknown compounds get casn = *name
		findCasnFromChemname = function( chemnamelist = "" ){
			
			chemnamelist <- tolower(chemnamelist);
			if(length(chemnamelist)==0){
				return ("");
			}
			
			# if (is.character(private$source)){
				table <- Class.DataImporter.dreamtk.tableinfo$chemtablename;
				tablevars <- Class.DataImporter.dreamtk.tableinfo$chemtablevars;
				if(is.null(private$chemical_data)){
					private$chemical_data <- private$loadDataFromRData ( file.path = private$source, 
																		table = table, 
																		tablevars = tablevars, 
																		filter_by = "name", 
																		filter_val = "*" );
				}													
				chemcastable <- select(private$chemical_data, casn, name);
				chemcastable$name <- tolower(chemcastable$name);
				chemcastable <- filter(chemcastable, name %in% chemnamelist);
			# }else{
			# 	private$database$dbConnect();
			# 	private$database$dbSetActiveDB("dream_tk");
			# 	if(is.null(private$chemical_data)){
			# 		chemcastable <- private$database$dbGetQuery( paste0("select casn, name from ", Class.DataImporter.dreamtk.tableinfo$chemtablename,
			# 						" where lower(name) in ('", str_c(chemnamelist,sep=" ",collapse="','"), "')" ) );
			# 	}else{
			# 		chemcastable <- select(private$chemical_data, casn, name);
			# 	}
			# 	chemcastable$name <- tolower(chemcastable$name);
			# 
			# 	private$database$dbDisconnect();
			# }
			chemcaslist <- chemnamelist;
			for (i in 1:length(chemnamelist)){
				tableindex <- which(chemcastable$name %in% chemnamelist[i]);
				if( length(tableindex)>0 ){
					chemcaslist[i] <- chemcastable$casn[tableindex[1]];
				}
			}
			
			return(chemcaslist);
			
		}
		

	)
)