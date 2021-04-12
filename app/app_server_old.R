# Server function

#v0.8.3

server <- function( input, output, session ){
  options(shiny.maxRequestSize = 800 * 1024 ^ 2);
  
  #run this code once user connects to app
  
  #terminate session upon disconnection
  session$onSessionEnded(stopApp);
  
  #APPLICATION TUTORIAL DATAFRAMES ---------------------------------------
  
  #help steps dataframes -> separate file due to size
  source("./app/app_tutorial.R");
  
  #APPLICATION VARIABLES AND OBJECTS ------------------------------------------
  
  #application flags
  app_flags <- reactiveValues(database_operational = NULL, enable_searching = FALSE);
  app_flags$database_operational <- FALSE;
  
  #database connection
  dtk.db <- reactiveValues(obj = NULL);
  mysql_login <- list(host = "localhost", user = "root", password = "password", db = "", port = 3306); #Gabriel ,why do we even have an option if we are gonna force it to be local host and store the data in plain text...
  
  
  #database info
  rdatabase_file <- paste0("./data/DreamTKv", app.version.major, ".", app.version.minor,".RData");
  
  
  
  #dreamtk importers
  dtk.import <- reactiveValues(obj = NULL);
  dtk.import$obj <- Class.DataImporter.dreamtk.basicData$new();
  dtk.import.3css <- Class.DataImporter.dreamtk.3compss$new();
  dtk.import.fubPoly <- Class.DataImporter.dreamtk.FubPolynomial$new();
  
  #dreamtk parser
  dtk.parser <- Class.DataParser.dreamtk$new();
  
  #dreamtk chemical manager
  chem_manager <- Class.ChemicalManager$new();
  
  #local models
  model.3css <- Class.LocalModel.Css3compss$new(); #Css from 3 compartment steady-state model
  model.fubPoly <- Class.LocalModel.FubPolynomial$new(); #Fub based on polynomial estimation
  model.oed1 <- Class.LocalModel.OED1storder$new(); #OED from 1st order model
    
  #selection and toggle variables
  chem_autocomplete <- reactiveValues(data = tribble(~casn, ~name) );
  selected_chem <- reactiveValues(casn = NULL, name = NULL, chemical = NULL);
  selected_assay <- reactiveValues(name = NULL);
  selected_assay_endpoint <- reactiveValues(name = NULL);
  options_toggle <- reactiveValues(status = "Hide", icon = "square-o");
  selected_chem_stats <- reactiveValues(selected_chems = NULL, prev_selected_chems = NULL);


  #EVENT CODE +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  #functions which get run when there is a change in inputs
  
  #TUTORIAL ELEMENTS--------------------------------------------------------
  observeEvent(c(input$link_help_overview, input$button_tour_overview), ignoreInit = TRUE,
               {
                 introjs(session,options = list(steps=steps_basics(), exitOnOverlayClick = FALSE, keyboardNavigation = TRUE),
                         events=events_basics());
               }
  );
  observeEvent(c(input$link_help_search, input$button_tour_search), ignoreInit = TRUE,
               {
                 introjs(session,options = list(steps=steps_search(), exitOnOverlayClick = FALSE, keyboardNavigation = TRUE),
                         events=events_search());
               }
  );
  observeEvent(c(input$link_help_save, input$button_tour_save), ignoreInit = TRUE,
               {
                 introjs(session,options = list(steps=steps_save(), exitOnOverlayClick = FALSE, keyboardNavigation = TRUE),
                         events=events_save());
               }
  );
  observeEvent(c(input$link_help_analysis, input$button_tour_analysis), ignoreInit = TRUE,
               {
                 introjs(session,options = list(steps=steps_analysis(), exitOnOverlayClick = FALSE, keyboardNavigation = TRUE),
                         events=events_analysis());
               }
  );
  
  #GLOBAL ELEMENTS --------------------------------------------------------
  
  #user logout panel -> for user-enabled shiny servers, does not work on free account, or offline
  output$userpanel <- renderUI({
    # session$user is non-NULL only in authenticated sessions
    if (!is.null(session$user)) {
      sidebarUserPanel(
        span("Logged in as ", session$user),
        subtitle = a(icon("sign-out"), "Logout", href="__logout__"))
    }
  })
  
  #Database: SELECT DATABASE OPTION <- set up database objects based on database selection
  observeEvent(input$radio_dbtype,

               {
                 app_flags$database_operational <- FALSE; #by default database is not operational until it has been chosen and tested
                 disableSearching( app_flags );

                 # Removed the option for the MySQL Database. #YANIC
#                  if (input$radio_dbtype == "mysql"){
#                    dtk.db$obj <- Class.Database.MySQL$new();
#                    dtk.db$obj$setLoginInfo(host = mysql_login$host, user = mysql_login$user, password = mysql_login$password, db = mysql_login$db, save = "1", port = mysql_login$port);
# 				   dtk.import$obj$setSource(dtk.db$obj)
#                    tryCatch( #attempt database connection and get common model data
#                      {
#                        setupModels(input, output, session,
#                                    dtk.db$obj, dtk.import$obj,
#                                    model.3css, dtk.import.3css,
#                                    model.fubPoly, dtk.import.fubPoly,
#                                    model.oed1);
#                        
#                        #update available search autocomplete list from database
#                        updateSearchAutocomplete(dtk.import$obj, chem_autocomplete);
#                        
#                        app_flags$database_operational <- TRUE;
#                        enableSearching( app_flags );
#                        
#                        showNotification("MySQL database connection successful.", type = "message", duration = 10);
#                        output$ui_database_status <- renderUI({
#                          h5( paste0("MySQL database connection successful."), style = "color:green; white-space:normal;" );
#                        });
#                      },
#                      warning = function(w){
#                        showNotification( "Failed DreamTK database connection test. Try a different database.", type = "error", duration = 10 );
#                        showNotification( paste0(w), type = "error", duration = 10 );
#                        output$ui_database_status <- renderUI({
#                          h5( paste0("Failed MySQL DreamTK database connection test. Try a different database."), style = "color:red; white-space:normal;" );
#                        });
#                      },
#                      error = function(e){
#                        showNotification( "Failed DreamTK database connection test. Try a different database.", type = "error", duration = 10 );
#                        showNotification( paste0(e), type = "error", duration = 10 );
#                        output$ui_database_status <- renderUI({
#                          h5( paste0("Failed MySQL DreamTK database connection test. Try a different database."), style = "color:red; white-space:normal;" );
#                        });
#                      }
#                    );
#                  } else if (input$radio_dbtype == "rdata"){
                   dtk.import$obj$setSource(rdatabase_file);
				          dtk.db$obj <- rdatabase_file;
                   tryCatch( #attempt database connection and get common model data
                     {
                       progress <- shiny::Progress$new();
                       on.exit({
                         progress$close();
                       });
                       progress$set(message = "", value = 0);
                       n <- 7;
                       #preload all rdata so that it will not need to be loaded during the first chem search
                       progress$inc(1/n, detail = "Importing chem info.");
                       dtk.import$obj$loadChemicalData();
                       progress$inc(1/n, detail = "Importing physiology info.");
                       dtk.import$obj$loadPhysiologyData();
                       progress$inc(1/n, detail = "Importing assay info.");
                       dtk.import$obj$loadAssayInfo();
					 
                       progress$inc(1/n, detail = "Setting up model imports.");
                       setupModels(input, output, session,
                                   dtk.db$obj, dtk.import$obj,
                                   model.3css, dtk.import.3css,
                                   model.fubPoly, dtk.import.fubPoly,
                                   model.oed1);
						progress$inc(1/n, detail = "Importing CPDat info.");
                       dtk.import$obj$loadCPDatData();
					   progress$inc(1/n, detail = "Importing Sheds info.");
                       dtk.import$obj$loadShedData();
                       progress$inc(1/n, detail = "Importing GenTox21 info.");
                       dtk.import$obj$loadGenToxData();
					   
                       #update available search autocomplete list from database
                       updateSearchAutocomplete(dtk.import$obj, chem_autocomplete);
                       
                       app_flags$database_operational <- TRUE;
                       enableSearching( app_flags );
                       
                       showNotification("RData database obtained.", type = "message", duration = 10);
                       output$ui_database_status <- renderUI({
                         h5( paste0("RData database obtained."), style = "color:green; white-space:normal;" );
                       });
                     },
                     warning = function(w){
                       showNotification( "Failed DreamTK file connection test. Try a different database.", type = "error", duration = 10 );
                       showNotification( paste0(w), type = "error", duration = 10 );
                       output$ui_database_status <- renderUI({
                         h5( paste0("Failed DreamTK database connection test. Try a different database."), style = "color:red; white-space:normal;" );
                       });
                     },
                     error = function(e){
                       showNotification( "Failed DreamTK file connection test. Try a different database.", type = "error", duration = 10 );
                       showNotification( paste0(e), type = "error", duration = 10 );
                       output$ui_database_status <- renderUI({
                         h5( paste0("Failed DreamTK database connection test. Try a different database."), style = "color:red; white-space:normal;" );
                       });
                     }
                   );
  #                }
  # 
              }

  );
  
  #Search: CHEM DISPLAY OPTIONS RADIO BUTTONS USED -> update chemical display based on radio button selections
  observeEvent(input$radio_listtype,
               {
                 updateChemListByRadioButtonSelection( input, output, session, selected_chem, chem_manager );
               }
  );
  
  
  #TAB: Home ---------------------------------------------------------------


  
  #TAB: Search -------------------------------------------------------------
  
  #Search: SEARCH BUTTON PRESSED -> results in chemical search, addition, and model computation cascade, update of the chemical list box
  observeEvent(input$button_search, 
               {
                 disableSearching( app_flags );
                 
                 #
                 if (!app_flags$database_operational){
                   showNotification( "Could not connect to database. Try a different database.", type = "error", duration = 10 );
                   return (NULL);
                 }

                 #progress bar
                 progress <- shiny::Progress$new();
				 progress$set(message = "", value = 0);
                 on.exit({
                   progress$close();
                 });
                 
                 #obtain and validate search field input, return CASN list
                 if (input$radio_searchtype == "casn"){
                   chems <- parseCasnSearch(input$field_search);
                 } else {
                   chems <- parseNameSearch(dtk.import$obj, input$field_search); 
                 }
                 
                 #reset search field
                 updateSelectizeInput(session, "field_search", selected = "");
                 
                 chems <- removeChemsThatExistInChemicalManager( session, input, output, chems, chem_manager );
                 
                 if (chems != ""){
                   n <- 8; #number of progress bar updates
                   progress$inc(1/n, detail = "Importing chem info.");
                   dtk.import$obj$loadChemicalData(chems);
                   progress$inc(1/n, detail = "Importing assay info.");
                   dtk.import$obj$loadAssayInfo(chems);
                   
                   #create chemical and catch a warning about chemical that may not exist
                   progress$inc(1/n, detail = "Creating chem objects.");
                   withCallingHandlers(
                     chemlist <- dtk.parser$createListOfChemicals( dtk.import$obj$getChemInfo(), dtk.import$obj$getAssayInfo()), 
                     warning=function(w) { 
                       showNotification( conditionMessage(w), type = "warning", duration = 10 );
                       invokeRestart("muffleWarning");
                     }
                   )
                   #add chemical to the chem manager
                   #if it is already in the manager, it has models computed, so do not duplicate effort
                   progress$inc(1/n, detail = "Listing chem objects.");
                   duplicates_secondpass <- chem_manager$addListOfChemicals(chemlist);
                   chemlist <- chemlist[-duplicates_secondpass];
                   
                   tryCatch({ #throw user an error in case modeling throws any
                     #Gabriel fix this once change to chemical is done
                     #compute the models

                     runModels( chemlist,
                                input, output, session,
                                model.3css, model.oed1,
                                progress, n );
                     progress$inc(1/n, detail = "Updating interface.");

                     #update chemical listbox
                     updateSelectedChemical( input, output, session, selected_chem, chem_manager);
                     updateChemListByRadioButtonSelection( input, output, session, selected_chem, chem_manager );
                     
                   }, error = function(e){
                     chem_manager$removeChemicals(chems); #error occurred processing these chemicals, remove them
                     showNotification( "error occured", type = "error", duration = 10 );
                   });
                   
                 } else {
                   showNotification( "No chemicals were entered.", type = "error", duration = 5 );
                 }
                 
                 enableSearching( app_flags );
                 
               }
  );

  #SEARCH TYPE RADIO BUTTON SELECTION <- update search autocomplete list by selected option
  observeEvent(input$radio_searchtype,
               {
                 disableSearching( app_flags );
                 
                 #update search field choices
                 if(input$radio_searchtype == "casn"){
                   choices <- select(chem_autocomplete$data, casn) %>% arrange(casn);
                 } else {
                   choices <- select(chem_autocomplete$data, name) %>% arrange(name);
                 }
                 updateSelectizeInput(session, "field_search", choices = choices );
                 
                 enableSearching( app_flags );
               }
  );
  
  #FILE REQUIREMENTS HINT LINK <- toggles table requirements text
  observeEvent(input$link_customchem_hint,
               {
                 shinyjs::toggleElement(id = "panelCustomchemHint", anim = TRUE, animType = "fade", time = 0.3);
               }
  );
  
  #DOWNLOAD CSV TEMPLATE LINK <- downloads a csv template of the custom chemical file
  output$button_customchem_template <- downloadHandler(
    filename = function() {
        paste("DreamTKv", app.version.major,"-", app.version.minor,"-customchem-template-", format(Sys.time(), "%a_%b_%d_%Y_%H_%M_%S"), ".csv",sep="" );
    },
    content = function(file) {
      showNotification("Saving custom chemical template as csv...", type = "message", duration = 5);
      logdebug("Save file dialog: Saving custom chemical template as csv...");
      
      output_table <- tibble(casn = "CAS-RN", name = "MyChemical", cytotoxicity_um = 1000, cytotoxic = "N", mw = 42.0, 
                              human_funbound_plasma = 0.42, human_clint = 42.0, human_rblood2plasma = 4.2,
                              log_kow=0.42, pka=7.4);
      
      write_excel_csv(output_table, path=file, na = "");
    }
  );
  
  #PARSE FILE BUTTON PRESS <- parses the uploaded custom chemical csv file and adds chemical objects to manager
  observeEvent(input$button_load_customchems,
               {
                 disableSearching( app_flags );
                 
                 #
                 if (!app_flags$database_operational){
                   showNotification( "Could not connect to database. Try a different database.", type = "error", duration = 10 );
                   return (NULL);
                 }
                 
                 #progress bar
                 progress <- shiny::Progress$new();
                 n <- 9; #number of progress bar updates
                 on.exit({
                   progress$close();
                 });
                 
                 progress$set(message = "", value = 0);
                 progress$inc(1/n, detail = "Verifying file.");
                 #check if file is valid
                 inFile <- input$file_customchem;
                 if(is.null(inFile)){
                   output$ui_customchems_status <- renderUI({
                     h5( paste0("No file selected."), style = "color:red" );
                   });
                   enableSearching( app_flags );
                   return (NULL);
                 }
                 
                 #obtain the file and warn if any required columns are missing
                 withCallingHandlers(
                 customchems <- read_csv(file=inFile$datapath, col_names = TRUE, 
                                          col_types= cols(casn = "c",
                                          name = "c",
                                          cytotoxicity_um = "d",
                                          cytotoxic = "c",
                                          mw = "d",
                                          human_funbound_plasma = "d",
                                          human_clint = "d",
                                          human_rblood2plasma = "d",
                                          log_kow = "d",
                                          pka = "d")),
                   warning=function(w) {
                     logwarn( paste0( w, collapse=" "));
                     showNotification( conditionMessage(w), type = "warning", duration = 10 );
                     invokeRestart("muffleWarning");
                   }
                 )
                       
                 #CHECKPOINT: stop if database not available
                 if (!app_flags$database_operational){
                   logwarn("Custom chemicals file load dialog: Selected database unavailable, unable to proceed.");
                   showNotification( "Selected database unavailable, unable to proceed", type = "error", duration = 5 );
                   output$ui_customchems_status <- renderUI({
                     h5( paste0(inFile$name, " not loaded. Selected database unavailable."), style = "color:red" );
                   });
                   return (NULL);
                 }
                       
                 #CHECKPOINT: stop if provided file does not have at least a casn column
                 if (is.null(customchems$casn) || length(customchems$casn)==0){
                   logwarn("Custom chemicals file load dialog: Selected file missing casn column or is empty, unable to proceed.");
                   showNotification( "Selected file missing casn column or is empty, unable to proceed", type = "error", duration = 5 );
                   output$ui_customchems_status <- renderUI({
                     h5( paste0(inFile$name, " not loaded. Selected file missing casn column or is empty."), style = "color:red" );
                   });
                   enableSearching( app_flags );
                   return (NULL);
                 }
                 
                 #fill in missing values from database if available  
                 progress$inc(1/n, detail = "Checking for missing values and filling them in.");
                 casn <- customchems$casn;
                 
                 #first step is to ensure all columns are present and of the correct type
                 #coalesce function used later absolutely hates non-matching types!
                 customchem_vars <- names(customchems);
                 tablevars <- c("casn", "name", "cytotoxicity_um", "cytotoxic", "mw", "human_funbound_plasma",
                                "human_clint", "human_rblood2plasma", "log_kow", "pka");
                 tablevartypes <- c("c", "c", "d", "c", "d", "d", "d", "d", "d", "d");
                 indices <- which( tablevars %in% customchem_vars );
                 missing <- tablevars[ -indices ];
                 missingtypes <- tablevartypes[ -indices ];
                 index <- 0;
                 for ( tablevar in missing ) {
                   index <- index + 1;
                   if (missingtypes[index] == "c"){
                     customchems[[tablevar]] = as.character(NA);
                   } else if (missingtypes[index] == "d"){
                     customchems[[tablevar]] = as.double(NA);
                   }
                 }

                 #second step is to import actual chemical data to grab missing values from
                 dtk.import$obj$loadChemicalData(casn);
                 dtk.import$obj$loadAssayInfo(casn);
                 dtk.import.3css$loadModelData(casn);
                 dtk.import.fubPoly$loadModelData(casn);
                         
                 cheminfo <- dtk.import$obj$getChemInfo();
                 assayinfo <- dtk.import$obj$getAssayInfo();
                 model3cssinfo <- dtk.import.3css$getModelData();
                 modelfubpolyinfo <- dtk.import.fubPoly$getModelData();
                       
                 customchems <- left_join(customchems, cheminfo, by = "casn") %>%
                   left_join(model3cssinfo, by = "casn") %>%
                   left_join(modelfubpolyinfo, by = "casn") %>%
                   mutate( name = coalesce(name.x, name.y),
                           cytotoxicity_um = coalesce(cytotoxicity_um.x, cytotoxicity_um.y),
                           cytotoxic = coalesce(cytotoxic.x, cytotoxic.y),
                           mw = coalesce(mw.x, mw.y),
                           human_funbound_plasma = coalesce(human_funbound_plasma.x, human_funbound_plasma.y),
                           human_clint = coalesce(human_clint.x, human_clint.y),
                           human_rblood2plasma = coalesce(human_rblood2plasma.x, human_rblood2plasma.y),
                           log_kow = coalesce(log_kow.x, log_kow.y),
                           pka = coalesce(pka.x, pka.y) ) %>%
                  select(casn, name, cytotoxicity_um, cytotoxic, mw, human_funbound_plasma,
                         human_clint, human_rblood2plasma, log_kow, pka);
      
                  chemtable <- tibble(casn = customchems$casn,
                                      name = customchems$name,
                                      cytotoxicity_um = customchems$cytotoxicity_um,
                                      cytotoxic = customchems$cytotoxic);

                  css3compss_data <- tibble(casn = customchems$casn,
                                            mw = customchems$mw,
                                            human_funbound_plasma = customchems$human_funbound_plasma,
                                            human_clint = customchems$human_clint,
                                            human_clint_unit = NA,
                                            human_clint_p = NA,
                                            human_rblood2plasma = customchems$human_rblood2plasma,
                                            mouse_funbound_plasma = NA,
                                            rabbit_funbound_plasma = NA,
                                            rat_funbound_plasma = NA,
                                            rat_clint = NA,
                                            rat_clint_unit = NA,
                                            rat_clint_p = NA,
                                            rat_rblood2plasma = NA);

                  fubpolynomial_data <- tibble(casn = customchems$casn,
                                               log_kow = customchems$log_kow,
                                               pka = customchems$pka);
                         
                 progress$inc(1/n, detail = "Creating and adding chem objects.");
                 #create chemical and catch a warning about chemical that may not exist
                 withCallingHandlers(
                   chemlist <- dtk.parser$createListOfChemicals( chemtable,  assayinfo), 
                   warning=function(w) { 
                     showNotification( conditionMessage(w), type = "warning", duration = 10 );
                     invokeRestart("muffleWarning");
                   }
                 )
          
                 #save casn list of chemicals we add, since if something fails we will remove them by casn
                 chemcasn <- list()
                 for (chem in chemlist){
                   chemcasn <- append(chemcasn, chem$chem_info$casn);
                 }
                 chem_manager$removeChemicals( chemcasn );
                     
                 #add chemical to the chem manager
                 chem_manager$addListOfChemicals(chemlist);

                 dtk.import.3css$setModelData( css3compss_data );
                 dtk.import.fubPoly$setModelData( fubpolynomial_data );
           
                 tryCatch({ #throw user an error in case modeling throws any
                           
                   #compute the models
                   runModels( chemlist,
                              input, output, session,
                              model.3css, model.oed1,
                              progress, n);
                   
                   dtk.import.3css$resetModelData();
                   dtk.import.fubPoly$resetModelData();
                         
                   #update chemical listbox
                   updateSelectedChemical( input, output, session, selected_chem, chem_manager);
                   updateChemListByRadioButtonSelection( input, output, session, selected_chem, chem_manager );
                         
                 }, error = function(e){
                   chem_manager$removeChemicals(chemcasn); #error occurred processing these chemicals, remove them
                   showNotification( e, type = "error", duration = 10 );
                 });

                 output$ui_customchems_status <- renderUI({
                   h5( paste0(inFile$name, " loaded."), style = "color:green");
                 });

                 updateChemListByRadioButtonSelection( input, output, session, selected_chem, chem_manager );
                 enableSearching( app_flags );
               }
  );

  observeEvent(input$select_hit,
				{
					if(!is.null(input$select_chemical)){
						CreateSelectedChemicalUI(input,output,session,selected_chem, chem_manager, selected_assay, selected_assay_endpoint, model.3css, model.oed1);
						UpdateAssayComponentList(input,output,session,selected_chem,selected_assay,selected_assay_endpoint);
					}
				}
  );

  observeEvent(input$select_background,
				{
					updatePrettySwitch(session, "analyse_background", value = input$select_background);
					if(!is.null(input$select_chemical)){
						CreateSelectedChemicalUI(input,output,session,selected_chem, chem_manager, selected_assay, selected_assay_endpoint, model.3css, model.oed1);
						UpdateAssayComponentList(input,output,session,selected_chem,selected_assay,selected_assay_endpoint);
					}

				}
  );
  
  #Search: CHEMICAL LIST SELECTION -> results in update of assay list and population of chemical info
  observeEvent(input$select_chemical,
               {
                 CreateSelectedChemicalUI(input,output,session,selected_chem, chem_manager, selected_assay, selected_assay_endpoint, model.3css, model.oed1);
				 
				 
               }
  );
  
  #Search: DELETE SELECTED BUTTON PRESS -> results in deletion of a chemical selected in the chemical list, selection of next chemical
  observeEvent( input$button_delete_chem,
                {
                  
                  chem = selected_chem$casn;
                  
                  #do nothing if nothing is selected
                  if (is.null(chem) || is.na(chem)){
                    return (NULL);
                  }
                  
                  listlength <- chem_manager$getNumChemicals();
                  
                  if (input$radio_listtype == "casn"){
                    caslist <- chem_manager$getListOfChemicalCasn();
                    selected_index <- which(caslist %in% chem);
                  }else if (input$radio_listtype == "name"){
                    chemname = selected_chem$name;
                    namelist <- chem_manager$getListOfChemicalNames();
                    selected_index <- which(namelist %in% chemname);
                  }
                  
                  if(selected_index < listlength){
                    next_index <- selected_index;
                  } else {
                    next_index <- selected_index - 1;
                  }
                  
                  chem_manager$removeChemicals( chem );

                  #if no chemicals left, reset everything
                  if (chem_manager$getNumChemicals() < 1){
                    selected_chem$casn <- NULL;
                    selected_chem$name <- NULL;
                    selected_chem$chemical <- NULL;
                    updateChemListByRadioButtonSelection( input, output, session, selected_chem, chem_manager );
                    updateSelectInput(session=session, inputId = "select_assay", selected = NULL,
                                      choices = NA
                    );
                    updateSelectInput(session=session, inputId = "select_assay_comp", selected = NULL,
                                      choices = NA
                    );
                    output$html_chemicalinfo <- renderUI({HTML("")});
                    output$html_assayinfo <- renderUI({HTML("")});
                    return (NULL);
                  }
                  
                  if (input$radio_listtype == "casn"){
                    caslist <- chem_manager$getListOfChemicalCasn();
                    selected_chem$casn <- caslist[next_index];
                    selected_chem$chemical <- chem_manager$getChemical(selected_chem$casn);
                    selected_chem$name <- selected_chem$chemical$chem_info$name;
                  }else if (input$radio_listtype == "name"){
                    namelist <- chem_manager$getListOfChemicalNames();
                    selected_chem$name <- namelist[next_index];
                    selected_chem$chemical <- chem_manager$getChemicalByName(selected_chem$name);
                    selected_chem$casn <- selected_chem$chemical$chem_info$casn;
                  }

                  updateChemListByRadioButtonSelection(input, output, session, selected_chem, chem_manager);
                  
                }
                
  );
  
  #Search: DELETE MISSING BUTTON PRESS -> results in deletion of missing chemicals in the chemical list
  observeEvent( input$button_delete_missing,
                {
                  
                  if (input$radio_listtype == "casn"){
                    chem = selected_chem$casn;
                    caslist <- unlist( chem_manager$getListOfChemicalCasn() );
                    if (length(caslist)>0){
                      missingcasn <- caslist[ startsWith(caslist, "*") ];
                    } else {
                      missingcasn <- "";
                    }
                  }else if (input$radio_listtype == "name"){
                    chem = selected_chem$name;
                    namelist <- unlist( chem_manager$getListOfChemicalNames() );
                    caslist <- unlist( chem_manager$getListOfChemicalCasn() );
                    if (length(caslist)>0){
                      missingnames <- namelist[ startsWith(caslist, "*") ];
                      missingindices <- which(namelist %in% missingnames);
                      missingcasn <- caslist[ missingindices ];
                    }else {
                      missingcasn <- "";
                    }
                  }

                  chem_manager$removeChemicals( missingcasn );

                  #if no chemicals left, reset everything
                  if (chem_manager$getNumChemicals() < 1 || (!is.null(chem) && startsWith(chem, "*")) ){
                    selected_chem$casn <- NULL;
                    selected_chem$name <- NULL;
                    selected_chem$chemical <- NULL;
                    updateChemListByRadioButtonSelection( input, output, session, selected_chem, chem_manager );
                    updateSelectInput(session=session, inputId = "select_assay", selected = NULL,
                                      choices = NA
                    );
                    updateSelectInput(session=session, inputId = "select_assay_comp", selected = NULL,
                                      choices = NA
                    );
                    output$html_chemicalinfo <- renderUI({HTML("")});
                    output$html_assayinfo <- renderUI({HTML("")});
                    return (NULL);
                  }

                  updateChemListByRadioButtonSelection(input, output, session, selected_chem, chem_manager);
                  
                }
                
  );
  
  #Search: CLEAR LIST BUTTON PRESS -> results in deletion of all listed chemicals
  observeEvent( input$button_clear_list,
                {
                  
                  chem_manager$removeChemicals( "*" );
                  
                  #reset everything
                  selected_chem$casn <- NULL;
                  selected_chem$name <- NULL;
                  selected_chem$chemical <- NULL;
                  updateChemListByRadioButtonSelection( input, output, session, selected_chem, chem_manager );
                  updateSelectInput(session=session, inputId = "select_assay", selected = NULL,
                                    choices = NA
                  );
                  updateSelectInput(session=session, inputId = "select_assay_comp", selected = NULL,
                                    choices = NA
                  );
                  output$html_chemicalinfo <- renderUI({HTML("")});
                  output$html_assayinfo <- renderUI({HTML("")});
                  return (NULL);
                  
                }
                
  );
  
  #Search: LIST MISSING BUTTON PRESS -> shows a modal dialogue with a list of missing chemicals
  observeEvent( input$button_list_missing,
                {
                  
                  if (input$radio_listtype == "casn"){
                    chem = selected_chem$casn;
                    caslist <- unlist( chem_manager$getListOfChemicalCasn() );
                    if (length(caslist)>0){
                      missingcasn <- caslist[ startsWith(caslist, "*") ];
                    } else {
                      missingcasn <- "";
                    }
                  }else if (input$radio_listtype == "name"){
                    chem = selected_chem$name;
                    namelist <- unlist( chem_manager$getListOfChemicalNames() );
                    caslist <- unlist( chem_manager$getListOfChemicalCasn() );
                    if (length(caslist)>0){
                      missingnames <- namelist[ startsWith(caslist, "*") ];
                      missingindices <- which(namelist %in% missingnames);
                      missingcasn <- caslist[ missingindices ];
                    }else {
                      missingcasn <- "";
                    }
                  }
                  
                  missingcasn <- gsub(pattern = "[*]", replacement = "", missingcasn);
                  
                  showModal(modalDialog(
                    title = "Missing chemicals",
                    paste(missingcasn, collapse=", ")
                  ))
                  
                }
                
  );
  
  #Search: CHEM LIST SAVE CSV BUTTON PRESS -> saves chemical list results to a csv
  output$button_savechems <- downloadHandler(
    filename = function() {
        paste("DreamTKv", app.version.major,"-", app.version.minor,"-chem_list-", format(Sys.time(), "%a_%b_%d_%Y_%H_%M_%S"), ".csv",sep="" );
    },
    content = function(file) {
      showNotification("Please wait. Saving chemical list as csv...", type = "message", duration = 5);
      logdebug("Save file dialog: Saving chemical list as csv...");
        
      #obtain all chemicals in the list
      chemlist <- chem_manager$getListOfChemicals("*");
      
      #obtain model names for which values will be taken
      css_model_name <- model.3css$getModelName();
      oed_model_name <- model.oed1$getModelName();
      
      output_table <- tribble(~casn, ~name, ~cytotoxic, ~cytotoxicity_um, ~css, ~css_units, ~css_model,
                              ~css_assumptions, ~min_ac50, ~min_ac50_units, ~assay_component_endpoint_name,
                              ~organism, ~tissue, ~cell_short_name, ~biological_process_target, 
                              ~intended_target_family, ~intended_target_family_sub, ~gene_name, 
                              ~min_oed, ~min_oed_units, ~oed_model);
      

	  
      for ( chem in chemlist ){
        datarow <- getBasicChemInfo( chem, css_model_name, oed_model_name, linesep="|", input$select_background, input$select_hit );
		#remove character(0) that doesn't allow for stuff to be saved.
		for(i in 1:length(datarow)){
			datarow[[i]] <- if(!(identical(datarow[[i]], character(0)))) datarow[[i]] else "";
		}
		output_table <- bind_rows(output_table, datarow);
		
      }

      write_excel_csv(output_table, path=file, na = "");
    }
  );
  
  #Search: ASSAY LIST SELECTION, CHEMICAL LIST SELECTION -> results in population of assay component list, or clearing it if  new chem selected
  observeEvent( c(input$select_assay, input$select_chemical),
                {
                 UpdateAssayComponentList(input,output,session,selected_chem,selected_assay,selected_assay_endpoint);
                  
                }
  );
  
  #Search: ASSAY COMPONENT LIST SELECTION -> results in population of assay info, reacts to changes in chemical and assay selections
  observeEvent( c(input$select_assay, input$select_chemical, input$select_assay_comp),
                {
                  selected_assay_endpoint$name <- input$select_assay_comp;
                  #get the selected chemical
                  chem <- selected_chem$chemical;
                  
                  #update assay info html text field
                  output$html_assayinfo <- renderUI(
                    {
                      
                      if( is.null(chem) ||
                          is.null(selected_assay$name) || 
                          is.null(selected_assay_endpoint$name) || 
                          selected_assay_endpoint$name == "NA" ||
                          is.null(chem$assay_info) || length(chem$assay_info)==0 ||
                          length(chem$assay_info$assay_name)==0 ){
                        
                        HTML("");
                        
                      } else {
                        #get assay info for selected endpoint
                        oed_model_name <- model.oed1$getModelName();
                        datarow <- filter(chem$assay_info, assay_component_endpoint_name == selected_assay_endpoint$name);
                        
                        #in case of duplicates, choose the first element only. this also de-lists the items
                        organism <- datarow[["organism"]][[1]];
                        tissue <- datarow[["tissue"]][[1]];
                        cell <- datarow[["cell_short_name"]][[1]];
                        bpt <- datarow[["biological_process_target"]][[1]];
                        itf <- datarow[["intended_target_family"]][[1]];
                        itfs <- datarow[["intended_target_family_sub"]][[1]];
                        gn <- datarow[["gene_name"]][[1]];
						gs <- datarow[["gene_symbol"]][[1]];
                        ac50 <- datarow[["ac50"]][[1]]; #log10 value
                        acc <- datarow[["ac_cutoff"]][[1]]; #log10 value
                        act <- log10(datarow[["ac_top"]][[1]]); #log10 value
                        aeid <- datarow[["aeid"]][[1]];
                        
                        #print(c(organism, cell, bpt, itf, ac50, acc, act, aeid));
                        
                        oed <- chem$model_results[[oed_model_name]]$oed[ chem$model_results[[oed_model_name]]$aeid == aeid ];
                        
                        if(!is.na(oed)){
                          oed_units <- chem$model_results[[oed_model_name]]$oed_units[ chem$model_results[[oed_model_name]]$aeid == aeid ];
                        } else {
                          oed_units <- "";
                        }
                        
                        HTML(
                          str_c( paste0("<strong>Assay Name: </strong>", selected_assay$name),
                                 paste0("<strong>Assay Component Endpoint Name: </strong>", selected_assay_endpoint$name),
                                 paste0("<strong>Organism: </strong>", organism),
                                 paste0("<strong>Tissue: </strong>", tissue),
                                 paste0("<strong>Cell: </strong>", cell ),
                                 paste0("<strong>Biological Process Target: </strong>", bpt),
                                 paste0("<strong>Intended Target Family: </strong>", itf),
                                 paste0("<strong>Intended Target Family Sub: </strong>", itfs),
                                 paste0("<strong>Gene: </strong>", gn),
								 paste0("<strong>Gene Symbol: </strong>",gs),
                                 paste0("<strong>Ac50 (log10): </strong>", signif(ac50, digits = 4)),
                                 paste0("<strong>Ac Cutoff (log10): </strong>", signif(acc, digits = 4)),
                                 paste0("<strong>Ac50 above cutoff: </strong>", ifelse(ac50>acc, "Yes", "No")),
                                 paste0("<strong>Ac Top (log10): </strong>", signif(act, digits = 4)),
                                 paste0("<strong>Ac Top above cutoff: </strong>", ifelse(act>acc, "Yes", "No")),
                                 paste0("<strong>OED: </strong>", signif(oed, digits = 4), " ", oed_units),
                                 paste0("<strong>OED model: </strong>", model.oed1$getModelName()),
                                 
                                 sep = '<br>')
                        );
                      }
                      
                    }
                    
                  );
                  
                }
                
  );
  
  #Search: ASSAY LIST SAVE CSV BUTTON PRESS -> saves chemical list results to a csv
  output$button_saveassays <- downloadHandler(
    filename = function() {
      if ( is.null(selected_chem$casn) ){
        paste0("No-selection.csv");
      } else {
        paste0(selected_chem$casn, ".csv");
      }
    },
    content = function(file) {
      showNotification("Please wait. Saving assay list for selected chemical as csv...", type = "message", duration = 5);
      logdebug("Save file dialog: Saving chem assay list as csv...");
      
      output_table <- tribble(~casn, ~aeid, ~ac50, ~ac_top, ~ac_cutoff, ~assay_name, ~assay_component_name,
                              ~assay_component_endpoint_name, ~organism, ~tissue, ~cell_short_name,
                              ~biological_process_target, ~intended_target_family, ~intended_target_family_sub,
                              ~gene_name, ~oed, ~oed_units, ~oed_model);
      
      if ( !is.null(selected_chem$casn) ){
        #obtain model names for which values will be taken
        css_model_name <- model.3css$getModelName();
        oed_model_name <- model.oed1$getModelName();
		tmp = selected_chem$chemical$assay_info
		if(!input$select_background){
			tmp = subset(tmp,intended_target_family != "background measurement");
		}
		if(input$select_hit){
			tmp = subset(tmp, hitc == 1);
		}
        output_table <- tmp %>% 
          left_join(selected_chem$chemical$model_results[[oed_model_name]], by = "aeid") %>%
          add_column(oed_model = oed_model_name);
        output_table$ac_top <- log10(output_table$ac_top); #convert to log10 for consistency with other values
      } else {
        output_table <- tribble(~casn, ~aeid, ~ac50, ~ac_top, ~ac_cutoff, ~assay_name, ~assay_component_name,
                                ~assay_component_endpoint_name, ~organism, ~tissue, ~cell_short_name,
                                ~biological_process_target, ~intended_target_family, ~intended_target_family_sub,
                                ~gene_name, ~oed, ~oed_units, ~oed_model);
      }
      
      write_excel_csv(output_table, path=file, na = "");
    }
  );
  
  #TAB: Analysis -------------------------------------------------------------
  
  #Stats: SELECT ALL BUTTON PRESSED <- selects all chemicals in the list
  observeEvent(input$button_stats_selectall,
               {
                 if(input$radio_listtype == "casn"){
                   selected <- chem_manager$getListOfChemicalCasn();
                 } else if (input$radio_listtype == "name"){
                   selected <- chem_manager$getListOfChemicalNames();
                 }
                 updateSelectizeInput(session, inputId = "select_chemical_stats", selected = selected);
               }
  );
  
  #Stats: DESELECT ALL BUTTON PRESSED <- deselects all chemicals in the list
  observeEvent(input$button_stats_deselectall,
               {
                 updateSelectizeInput(session, inputId = "select_chemical_stats", selected = "");
               }
  );
  #gabriel
  #Stats: SELECT OR UNSELECTS ALL THE STATS. <- selects or unselects all the statistics in the list.
  observeEvent(input$button_select_all_stats,
               {

				tempChoices =  c("Target Family Counts and ac50 values" = "tfcounts",
									"Hierarchical cluster heatmap of Target Subfamily activities" = "tfhm",
									"Hierarchical cluster heatmap of Assay Endpoint activities" = "assayhm",
									"Ac50 vs ScalarTop" = "scalartop_ac50",
									"OED vs ScalarTop" = "scalartop_oed",
									"Burst Assay vs Not Burst Assay" = "ac50_box",
									"Chemical ToxPI Plots (Individual)" = "toxpi",
									"Chemical ToxPI Plots (Cytotoxicity)" = "toxpi2",
									"Chemical ToxPI Plots (Grouped)" = "toxpigroup");
				if(is.null(input$checkbox_stats)){
			    updateCheckboxGroupInput(
					session, "checkbox_stats", 'Select statistics', 
					choices = tempChoices, selected = tempChoices)
			   }else{
					updateCheckboxGroupInput(
					session, "checkbox_stats", 'Select statistics', 
					choices = tempChoices, selected = "")
			   }
			}
                
  );
  
  #Stats: INCLUDE BACKGROUND ASSAYS BUTTON PRESSED ON ANALYSE MENU.
  
	observeEvent(input$analyse_background,
				{
					#here we just change the value of the select_background and update the chemical selected. We don't need to do anything else
					updatePrettySwitch(session, "select_background", value = input$analyse_background);
					if(!is.null(input$select_chemical)){
						CreateSelectedChemicalUI(input,output,session,selected_chem, chem_manager, selected_assay, selected_assay_endpoint, model.3css, model.oed1);
						UpdateAssayComponentList(input,output,session,selected_chem,selected_assay,selected_assay_endpoint);
					}
				}
	);
  #Stats: RUN STATS BUTTON PRESSED <- computes and displays statistics
  observeEvent(input$button_stats_run,
               {
			  
                 #clear dynamic UI stats output if exists, uses jquery selector notation with # in front of ui element id
                 uilist <- c("#ui_stats_tfcounts", "#ui_stats_tfhm", "#ui_stats_assayhm",
                             "#ui_stats_scalartop", "#ui_stats_toxpi", "#ui_stats_toxpi2", "#ui_stats_toxpi_group", "#ui_stats_Ac50box",
                             "#ui_gentox_genToxDataTable", "#ui_gentox_bmdDataTable");
							 
				clearTabUI( input, output, session, 
                                     uilist );
                 #depending on options selection the list is in casn or name format, obtain chemicals accordingly
                 if(input$radio_listtype == "casn"){
                   caslist <- input$select_chemical_stats;
                   chemicals <- chem_manager$getListOfChemicals(caslist);
                 } else if (input$radio_listtype == "name"){
                   namelist <- input$select_chemical_stats;
                   chemicals <- chem_manager$getListOfChemicalsByName(namelist);
                 }
                 numChems <- length(chemicals);
                 if( numChems>0 ){
                   
                   #progress bar
                   progress <- shiny::Progress$new();
                   on.exit({
                     progress$close();
                   });
                   
                   n <- 15; #number of progress bar updates
                   progress$set(message = "", value = 0);
                   #get analysis choices checkboxes and compute all analyses
				   
				   basicAnalysis = Class.Analysis.BasicAnalysis$new();
                   stat_choices <- input$checkbox_stats;
                   gt_choices <- c();
                   if (input$include_GenTox)
                     gt_choices <- input$checkbox_GenTox;
                   if ( length(stat_choices)>0 || length(gt_choices) > 0) {
                      computeAnalyses( input, output, session,
                                       chemicals, stat_choices, 
                                       gt_choices, basicAnalysis,
                                       model.oed1$getModelName(),
                                       progress, n );
                     
                   }
                   else {
                     showNotification("No analysis selected.", type="warning", duration=5);
                     loginfo("Stats: No analysis selected.");
                     return();
                   }

                   
                   #target family analysis chosen
                   if(any("tfcounts" %in% stat_choices)){
                     
                     progress$inc(1/n, detail = "Generating output for target family data.");
                     createTargetFamilyCountUI( input, output, session,
                                                basicAnalysis, id = "ui_stats_tfcounts" );
                   }
                   
                   
                   #target family heatmap chart
                   if(any("tfhm" %in% stat_choices)){
                     
                     progress$inc(1/n, detail = "Generating output for target family heatmap data.");
                     createTargetFamilyHeatmapUI( input, output, session,
                                                  basicAnalysis, id = "ui_stats_tfhm" );
                   }
                   
                   #assay endpoint heatmap chart
                   if(any("assayhm" %in% stat_choices)){
                     
                     progress$inc(1/n, detail = "Generating output for assay endpoint heatmap data.");
                     createAssayHeatmapUI( input, output, session,
                                           basicAnalysis, id = "ui_stats_assayhm" );
                   }


                   #ac50/oed and scalar top chart
                   if(any(c("scalartop_ac50", "scalartop_oed") %in% stat_choices)){
                     
                     progress$inc(1/n, detail = "Generating output for scalar top data.");
                     createScalarTopUI( input, output, session,
                                        basicAnalysis, stat_choices, id = "ui_stats_scalartop" );
                   }
					#ac50 box chart
                   if(any("ac50_box" %in% stat_choices)){
                     
                     progress$inc(1/n, detail = "Generating output for Burst Assay vs Not Burst Assay box chart.");
                     createAc50BoxUI( input, output, session,
                                        basicAnalysis, id = "ui_stats_Ac50box" );
                   }
                   
                   #toxpi charts
                   if(any("toxpi" %in% stat_choices)){
                     progress$inc(1/n, detail = "Generating output for ToxPI data.");
                     createToxPIUI( input, output, session,
                                    basicAnalysis, numChems, id = "ui_stats_toxpi");
                   }
                   if(any("toxpi2" %in% stat_choices)){
                     
                     progress$inc(1/n, detail = "Generating output for ToxPI data.");
                     createToxPI2UI( input, output, session,
                                    basicAnalysis, numChems, id = "ui_stats_toxpi2" );
                   }
                   if(any("toxpigroup" %in% stat_choices)){
                     
                     progress$inc(1/n, detail = "Generating output for ToxPI Grouped data.");
                     createToxPIGroupUI( input, output, session,
                                    basicAnalysis, numChems, id = "ui_stats_toxpi_group" );
                   }
                   
                   if(any("GTDataTable" %in% gt_choices)) {
                     progress$inc(1/n, detail = "Generating output for GenTox Data Table.")
                     createGenToxDataTable(input, output, session,
                                                      basicAnalysis, id = "ui_gentox_genToxDataTable");
                   }
                   
                   if(any("BMDDataTable" %in% gt_choices)) {
                     progress$inc(1/n, detail = "Generating output for BMD Data Table.")
                     createBMDDataTable(input, output, session,
                                        basicAnalysis, id = "ui_gentox_bmdDataTable");
                   }
				   
				 if(input$radio_listtype == "casn"){
                   tmp = basicAnalysis$basicData$getHitlessChemInfoByCasn();
				   if(length(tmp)>0){
						showNotification(paste0(paste(tmp, collapse=", "),  " might not appear on some plots because they lack data"), type="warning", duration=length(tmp) + 5);
				   
				   }
                 } else if (input$radio_listtype == "name"){
                   tmp = basicAnalysis$basicData$getHitlessChemInfoByName();
				   if(length(tmp)>0){
						showNotification(paste0(paste(tmp, collapse=", "), " might not appear on some plots because they lack data."), type="warning", duration= length(tmp) + 5); #displaying this for 1 extra second per chemical
				   
				   }
                 }
				
                 } else {
                   showNotification("No chemicals selected.", type="warning", duration=5);
                 }
               }
  );
  
  #TAB: MFA ---------------------------------------------------------------
  
  #MFA: SELECT ALL BUTTON PRESSED <- selects all chemicals in the list
  observeEvent(input$button_mfa_selectall,
               {
                 if(input$radio_listtype == "casn"){
                   selected <- chem_manager$getListOfChemicalCasn();
                 } else if (input$radio_listtype == "name"){
                   selected <- chem_manager$getListOfChemicalNames();
                 }
                 updateSelectizeInput(session, inputId = "select_chemical_mfa", selected = selected);
               }
  );
  
  #MFA: DESELECT ALL BUTTON PRESSED <- deselects all chemicals in the list
  observeEvent(input$button_mfa_deselectall,
               {
                 updateSelectizeInput(session, inputId = "select_chemical_mfa", selected = "");
               }
  );
  
  #MFA: RUN MFA BUTTON PRESSED <- runs mfa
  observeEvent(input$button_mfa_run,
               {
                 
                 #clear dynamic UI mfa output if exists
                 removeUI(
                   selector = "#ui_mfa_plots", immediate = TRUE
                 )
                 
                 #progress bar
                 progress <- shiny::Progress$new();
                 on.exit({
                   progress$close();
                 });
                 
                 n <- 3; #number of progress bar updates
                 progress$set(message = "", value = 0);
                 progress$inc(1/n, detail = "Building MFA table.");
                 
                 if(input$radio_listtype == "casn"){
                   caslist <- input$select_chemical_mfa;
                   chemicals <- chem_manager$getListOfChemicals(caslist);
                 } else if (input$radio_listtype == "name"){
                   namelist <- input$select_chemical_mfa;
                   chemicals <- chem_manager$getListOfChemicalsByName(namelist);
                 }
				 mfaAnalysis = Class.Analysis.MFA$new();
                 mfaAnalysis$buildMFATable( chemicals, label_by = input$radio_listtype );
				 
                 if(nrow(mfaAnalysis$getChemData())<1){
                   logwarn("Please select at least 1 valid chemical");
                   showNotification("Please select at least 1 valid chemical", type = "error", duration = 5);
                   return (NULL);
                 }
                 
                 progress$inc(1/n, detail = "Running FactoMineR MFA script.");
                 loginfo("Running FactoMineR MFA script.");

                 tryCatch({
                  mfaAnalysis$MFA(group=c(1,1,1), type=c("c","c","n"), ncp=5, name.group=c("ac50","scalar_top", "intended_target_family"));
                 },
                 error = function(e){
                   logerror("Could not run MFA. FactoMineR probably did not like something about missing data");
                   logerror(e);
                   showNotification("Could not run MFA. FactoMineR probably did not like something about missing data",
                                    type = "error", duration = 5);
                   showNotification(paste0("Error detail: ", e),
                                    type = "error", duration = 5);
                   return (NULL);
                 }
                 )
                   
                 progress$inc(1/n, detail = "Plotting MFA results.");
                 createMFAUI( input, output, session,
                              mfaAnalysis$getChemData(), mfaAnalysis$getChemDataMFA(), 
                              id = "ui_mfa_plots" );

               }
  );
  
   #TAB: BER ---------------------------------------------------------------
  
  #BER: SELECT ALL BUTTON PRESSED <- selects all chemicals in the list
  observeEvent(input$button_ber_selectall,
               {

                 if(input$radio_listtype == "casn"){
                   selected <- chem_manager$getListOfChemicalCasn();
                 } else if (input$radio_listtype == "name"){
                   selected <- chem_manager$getListOfChemicalNames();
                 }
                 updateSelectizeInput(session, inputId = "select_chemical_ber", selected = selected);
               }
  );
  
  #BER: DESELECT ALL BUTTON PRESSED <- deselects all chemicals in the list
  observeEvent(input$button_ber_deselectall,
               {
                 updateSelectizeInput(session, inputId = "select_chemical_ber", selected = "");
               }
  );
  
  #BER: RUN BER BUTTON PRESSED <- runs mfa
  observeEvent(input$button_ber_run,
               {
                 
                 #clear dynamic UI mfa output if exists
                 removeUI(
                   selector = "#ui_ber_plots", immediate = TRUE
                 )
				 
				uilist <- c("ui_stats_ber");
							 
				clearTabUI( input, output, session, 
                                     uilist );
                 
                 #progress bar
                 progress <- shiny::Progress$new();
                 on.exit({
                   progress$close();
                 });
                 
				 
				 
				 
                 n <- 2; #number of progress bar updates
                 progress$set(message = "", value = 0);
                 progress$inc(1/n, detail = "Building BER table."); 
				 
				 #depending on options selection the list is in casn or name format, obtain chemicals accordingly
                 if(input$radio_listtype == "casn"){
                   caslist <- input$select_chemical_ber;
                   chemicals <- chem_manager$getListOfChemicals(caslist);
                 } else if (input$radio_listtype == "name"){
                   namelist <- input$select_chemical_ber;
                   chemicals <- chem_manager$getListOfChemicalsByName(namelist);
                 }
				 
				 BERAnalysis = Class.Analysis.BERAnalysis$new();
				 oed_model_name <- model.oed1$getModelName();
                 BERAnalysis$BERData$buildBERStatsTable(dtk.import$obj$getCpdatInfo() ,dtk.import$obj$getShedsInfo(),input$select_chemical_ber, label_by = input$radio_listtype );
				 BERAnalysis$basicData$buildBasicStatsTable(chemicals);
				 BERAnalysis$basicData$addOEDtoBasicStatsTable (chemicals, oed_model_name );
				 
				 
				 
				 if(!BERAnalysis$BERData$BERStatsDataExists()){
				   logwarn("Please select at least 1 valid chemical");
                   showNotification("Please select at least 1 valid chemical", type = "error", duration = 5);
                   return (NULL);
				 }
				  if(input$radio_listtype == "casn"){
                   tmp = BERAnalysis$basicData$getHitlessChemInfoByCasn();
				   if(length(tmp)>0){
						showNotification(paste0(paste(tmp, collapse=", "),  " do not have any data for Ac50 values."), type="warning", duration=length(tmp) + 5);
				   
				   }
                 } else if (input$radio_listtype == "name"){
                   tmp = BERAnalysis$basicData$getHitlessChemInfoByName();
				   if(length(tmp)>0){
						showNotification(paste0(paste(tmp, collapse=", "), " do not have any data for Ac50 values."), type="warning", duration= length(tmp) + 5); #displaying this for 1 extra second per chemical
				   
				   }
                 }
				 
				 
				BERAnalysis$BERData$computeBER();
                BERAnalysis$computerMeanTableBER();  
				
                progress$inc(1/n, detail = "Plotting BER results.");
				BERAnalysis$plotBERvsAc50(label_by = input$radio_listtype);
				BERAnalysis$plotOEDvsAc50(label_by = input$radio_listtype);
				BERAnalysis$plotBER(label_by = input$radio_listtype);
				 
                 createBERUI( input, output, session, BERAnalysis,
                            id = "ui_ber_plots" );

               }
  );
  
  #TAB: GenTox Analysis -------------------------------------------------------------
  
  #GTA: SELECT ALL BUTTON PRESSED <- selects all chemicals in the list
  observeEvent(input$button_gt_selectall,
               {
                 if(input$radio_listtype == "casn"){
                   selected <- chem_manager$getListOfChemicalCasn();
                 } else if (input$radio_listtype == "name"){
                   selected <- chem_manager$getListOfChemicalNames();
                 }
                 updateSelectizeInput(session, inputId = "select_chemical_gt", selected = selected);
               }
  );
  
  #GTA: DESELECT ALL BUTTON PRESSED <- deselects all chemicals in the list
  observeEvent(input$button_gt_deselectall,
               {
                 updateSelectizeInput(session, inputId = "select_chemical_gt", selected = "");
               }
  );
  
  #GTA: Run analysis button pressed
  observeEvent(input$button_gt_createPlots,
   {
     if(input$radio_listtype == "casn"){
       caslist <- input$select_chemical_gt;
       chemicals <- chem_manager$getListOfChemicals(caslist);
     } else if (input$radio_listtype == "name"){
       namelist <- input$select_chemical_gt
       chemicals <- chem_manager$getListOfChemicalsByName(namelist);
     }
     
     if(length(chemicals) == 0)
     {
       showNotification("No chemicals selected.", type = "warning", duration = 5);
       return(NULL)
     }
     
     if (input$gt_bmr_picker == 2) #customBMR
     {
       if (input$gt_custom_choice == "")
       {
         showNotification("Please insert a custom BMR value.", type = "warning", duration = 5);
         return(NULL)
       }
       custBMR <- as.numeric(input$gt_custom_choice)
       if (is.na(custBMR))
       {
         showNotification("The custom BMR value must be a number.", type = "warning", duration = 5);
         return(NULL)
       }
       else if (custBMR <= 0 || custBMR >= 500)
       {
         showNotification("The custom BMR value must be between 0% and 500%, exclusively.", type = "warning", duration = 5);
         return(NULL)
       }
     }
     
     removeUI(
       selector = "#ui_gt_plotContainerLeft", immediate = TRUE, multiple = TRUE
     )
     removeUI(
       selector = "#ui_gt_plotContainerRight", immediate = TRUE, multiple = TRUE
     )
     removeUI(
       selector = "#ui_gt_plotContainerCenter", immediate = TRUE, multiple = TRUE
     )
     
     insertUI(
       selector = "#ui_gt",
       where = "afterBegin",
       ui = fluidRow(id = "ui_gt_plotContainerCenter")
     )
     
     insertUI(
       selector = "#ui_gt",
       where = "beforeEnd",
       ui = column(6,
                   fluidRow(id = "ui_gt_plotContainerLeft")
            )
       )
     
     insertUI(
       selector = "#ui_gt",
       where = "beforeEnd",
       ui = column(6,
                   fluidRow(id = "ui_gt_plotContainerRight")
       )
     )
     
     #progress bar
     progressB <- shiny::Progress$new();
     on.exit({
       progressB$close();
     });
     
     progressB$set(message = "", detail = "Generating GeneTox Table", value = 0);
     
     gtAnalysis = Class.Analysis.GT$new();
     gtAnalysis$buildGTTable( input$select_chemical_gt, label_by = input$radio_listtype, pValue = p, model = input$gt_modelComp, logValues = input$log_switch_gt, 
                              preCalcBMR = (input$select_graphType_gt %in% c(2,3) || input$gt_bmr_picker == 1), customBMRVal = ifelse((input$gt_modelComp %in% c(2,3) || input$gt_bmr_picker == 1), -1, as.numeric(input$gt_custom_choice)),
                              showAssayDesc = input$description_switch_gt, progress = progressB);

     if(nrow(gtAnalysis$getChemData())<1){
       logwarn("None of the selected chemicals were found with this model type");
       showNotification("None of the selected chemicals were found with this model type", type = "error", duration = 5);
       return (NULL);
     }
     
     progressB$set(1/2, detail = "Generating GeneTox Plots");
     
     if (input$select_graphType_gt == 1)
     {
      gtAnalysis$plotGeneToxSingleBMR(input, output, bmrValue = ifelse((input$gt_modelComp %in% c(2,3) || input$gt_bmr_picker == 1), input$gt_precalculated_choices, as.numeric(input$gt_custom_choice)), progress = progressB, logVal = input$log_switch_gt, showPlotData = input$plottingData_switch_gt, forceScale = input$forceScale_switch_gt, forceOrder = input$forceOrder_switch_gt);
     }
     else if (input$select_graphType_gt == 2)
     {
       gtAnalysis$plotGeneToxMultiBMR(input, output, progress = progressB, logVal = input$log_switch_gt, forceScale = input$forceScale_switch_gt);
     }
     else if (input$select_graphType_gt == 3)
     {
       gtAnalysis$plotGeneToxBMRDistribution(input, output, progress = progressB, logVal = input$log_switch_gt, showPlotData = input$plottingData_switch_gt, forceScale = input$forceScale_switch_gt);
     }
     else if (input$select_graphType_gt == 4)
     {
       gtAnalysis$plotGeneToxHeatMap(input, output, bmrValue = ifelse((input$gt_modelComp %in% c(2,3) || input$gt_bmr_picker == 1), input$gt_precalculated_choices, as.numeric(input$gt_custom_choice)), progress = progressB, logVal = input$log_switch_gt, showPlotData = input$plottingData_switch_gt)
     }
     else if (input$select_graphType_gt == 5)
     {
       gtAnalysis$plotGeneToxSingleBMRAll(input, output, bmrValue = ifelse((input$gt_modelComp %in% c(2,3) || input$gt_bmr_picker == 1), input$gt_precalculated_choices, as.numeric(input$gt_custom_choice)), progress = progressB, logVal = input$log_switch_gt, showPlotData = input$plottingData_switch_gt)
     }
     
   }
  );
  
  #TAB: Save/Load ----------------------------------------------------------
  
  #SAVE FILE
  output$button_savefile <- downloadHandler(
    filename = function() {
      if(input$radio_savetype == "cm"){
        paste("DreamTKv", app.version.major,"-", app.version.minor,"-chem_list-", format(Sys.time(), "%a_%b_%d_%Y_%H_%M_%S"), ".rds",sep="" );
      };
    },
    content = function(file) {
      if(input$radio_savetype == "cm"){
        showNotification("Saving R Chemical object list...", type = "message", duration = 5);
        logdebug("Save file dialog: Saving R Chemical object list...");
        chemlist <- chem_manager$getListOfChemicals();
        saveRDS(chemlist, file = file);
      };
    }
  );
  
  #LOAD FILE
  observeEvent(input$button_loadfile,
               {
                 inFile <- input$file_load;
                 if(is.null(inFile)){
                   showNotification("No file selected", type = "error", duration = 5);
                   output$ui_load_status <- renderUI({
                     h5( paste0("No file selected."));
                   });
                 }
                 if(input$radio_loadtype == "cm"){
                   tryCatch(
                     {
                       chemlist <- readRDS(inFile$datapath);
                       chem_manager$addListOfChemicals(chemlist);
                       logdebug("Save file dialog: R Chemical object list imported.");
                       #showNotification("R Chemical object list imported.", type = "message", duration = 5);
                       output$ui_load_status <- renderUI({
                         h5( paste0(inFile$name, " successfully loaded."), style = "color:green");
                       });
                     },
                     
                     error = function(e){
                       logwarn("Load file dialog: File not loaded. Invalid R Chemical object list.");
                       #showNotification("File not loaded. Invalid R Chemical object list.", type = "error", duration = 10);
                       output$ui_load_status <- renderUI({
                         h5( paste0(inFile$name, " not loaded. Invalid R Chemical object list."), style = "color:red" );
                       });
                     }
                   );
                   updateChemListByRadioButtonSelection( input, output, session, selected_chem, chem_manager );
                 }
               }
  );
  
  #TAB: Import Data Table ---------------------------------------------------------------
  
  #LOAD TABLE
  observeEvent(input$button_chooseFile,
               {
                 shinyjs::disable("button_chooseFile");
                 shinyjs::disable("button_chooseFolder");
                   if (!exists("tablePath")) tablePath <<- "";
                   newTable <- choose.files(tablePath, multi = FALSE);
                   if (!identical(newTable, character(0)) && !is.na(newTable) && !is.null(newTable)) {tablePath <<- newTable}
                   output$out_FilePath <- renderText({tablePath});
                   
                   if(tablePath != "")
                   {
                     tryCatch(
                       {
                         fileExtention <- sub('.*\\.', '', tablePath);
                         if (fileExtention == "txt")
                         {
                           LoadedNewTable <<- read.table(tablePath, header = TRUE, sep="\t", check.names = FALSE, stringsAsFactors = FALSE, comment.char = "", blank.lines.skip = TRUE);
                         }
                         else if (fileExtention == "csv")
                         {
                           LoadedNewTable <<- read.table(tablePath, header = TRUE, sep=",", check.names = FALSE, stringsAsFactors = FALSE, comment.char = "", blank.lines.skip = TRUE);
                         }
                         else
                         {
                           stop("Wrong file Type");
                         }
                         
                         colnames(LoadedNewTable)[is.na(colnames(LoadedNewTable)) | colnames(LoadedNewTable) == ""] <<- seq(length(colnames(LoadedNewTable)[is.na(colnames(LoadedNewTable)) | colnames(LoadedNewTable) == ""]))
                         
                         #UTF-8 Encoding can sometime leave some undesirable characters in front of the first column
                         colName <- names(LoadedNewTable)[1];
                         if (startsWith(colName, "?..")) {
                           colnames(LoadedNewTable)[1] <<- gsub("?..", "", colName, fixed = TRUE);
                         }
                         else if (startsWith(colName, "﻿")) {
                           colnames(LoadedNewTable)[1] <<- gsub("﻿", "", colName, fixed = TRUE);
                         }
                         rm(colName);
                         
                         #Remove totally empty columns and rows
                         LoadedNewTable <<- LoadedNewTable %>%
                           select_if(function(x) !all(is.na(x) | x == "")) %>%
                           filter_all(any_vars(!is.na(.) & . != ""))
                         
                         output$table_out <- renderTable({LoadedNewTable});
                         shinyjs::show("button_importTable");
                       },
                       error = function(e) {
                         output$table_out <- renderTable({matrix("The file is invalid", nrow = 1, ncol = 1, dimnames = list("Error", "Error"))});
                         cat(paste0("Loading Error: ", as.character(e)), file = stderr());
                         shinyjs::hide("button_importTable");
                       }
                     )
                   }
                   else
                   {
                     shinyjs::hide("button_importTable");
                   }
                   shinyjs::enable("button_chooseFile");
                   shinyjs::enable("button_chooseFolder");
               }
  );
  
  #Import Single Table
  observeEvent(input$button_importTable,
               {
                 load("data/GenTox21.RData");
                 
                 tableID = input$tableDropDown;
                 tryCatch(
                   {
                     if (tableID == 1)
                     {
                       newRows <- ldply(importTable(LoadedNewTable, tableID), data.frame, stringsAsFactors=FALSE);
                       gen_tox_data <- rbind(gen_tox_data, newRows);
                       trend_checker <- newRows[!duplicated(newRows[,c("compound", "assay", "endpoint", "s9_positive")]),];
                       if (nrow(trend_checker) > 0)
                       {
                         for (i in seq(nrow(trend_checker)))
                         {
                           if (nrow(trend_data) < 1)
                           {
                             trd = c();
                           }
                           else
                           {
                             trd <- which(trend_data["compound"] == trend_checker[i, "compound"] &
                                            trend_data["assay"] == trend_checker[i, "assay"] &
                                            trend_data["endpoint"] == trend_checker[i, "endpoint"] &
                                            trend_data["s9_positive"] == trend_checker[i, "s9_positive"]);
                           }
                           
                           currChemEndptData <- gen_tox_data[gen_tox_data[,"compound"] == trend_checker[i, "compound"] &
                                                        gen_tox_data["assay"] == trend_checker[i, "assay"] &
                                                        gen_tox_data["endpoint"] == trend_checker[i, "endpoint"] &
                                                        gen_tox_data["s9_positive"] == trend_checker[i, "s9_positive"]];
                           
                           logDose <- log(currChemEndptData[,"concentration"])
                           logMnResp <- log(currChemEndptData[,"response"])
                           
                           MAD <- median(abs(logMnResp-median(logMnResp)))
                           
                           paramsMn <- tcplFit(logc=logDose, resp=logMnResp, bmad=MAD, force.fit=TRUE)
                           
                           aic_hill_mn <- paramsMn$hill_aic
                           aic_gnls_mn <- paramsMn$gnls_aic
                           
                           if (length(trd) == 1)
                           {
                             trend_data[trd, "trend"] <- (!is.na(aic_hill_mn) || !is.na(aic_gnls_mn))
                           }
                           else
                           {
                             trend_data <- rbind(trend_data, data.frame(
                               casn = trend_checker[i, "casn"],
                               compound = trend_checker[i, "compound"],
                               assay = trend_checker[i, "assay"],
                               endpoint = trend_checker[i, "endpoint"],
                               s9_positive = trend_checker[i, "s9_positive"],
                               trend = (!is.na(aic_hill_mn) || !is.na(aic_gnls_mn)), stringsAsFactors = FALSE))
                           }
                         }
                       }
                     }
                     else if (tableID == 2)
                     {
                       newBMDs <- ldply(importTable(LoadedNewTable, tableID), data.frame, stringsAsFactors = FALSE);
                       if (nrow(newBMDs) == 0)
                       {
                         showNotification("The table was not imported. No valid rows were presented", type = "error", duration = 5);
                         return();
                       }
                       for (i in seq(nrow(newBMDs)))
                       {
                         if (nrow(bmd_data) < 1)
                         {
                           bmdToChange = c();
                         }
                         else
                         {
                           bmdToChange <- which(bmd_data["compound"] == newBMDs[i, "compound"] &
                                                bmd_data["assay"] == newBMDs[i, "assay"] &
                                                bmd_data["endpoint"] == newBMDs[i, "endpoint"] &
                                                bmd_data["s9_positive"] == newBMDs[i, "s9_positive"] &
                                                bmd_data["model_type"] == newBMDs[i, "model_type"] &
                                                bmd_data["bmr"] == newBMDs[i, "bmr"]);
                         }
                         
                         if (length(bmdToChange) == 1)
                         {
                           bmd_data[bmdToChange, "aic"] <- newBMDs[i, "aic"];
                           bmd_data[bmdToChange, "log_likelihood"] <- newBMDs[i, "log_likelihood"];
                           bmd_data[bmdToChange, "bmdl"] <- newBMDs[i, "bmdl"];
                           bmd_data[bmdToChange, "bmd"] <- newBMDs[i, "bmd"];
                           bmd_data[bmdToChange, "bmdu"] <- newBMDs[i, "bmdu"];
                         }
                         else
                         {
                           bmd_data <- rbind(bmd_data, newBMDs[i, ]);
                         }
                       }
                     }
                     
                     save(gen_tox_data, bmd_data, trend_data, assay_description, file="data/GenTox21.RData");
                 
                     showNotification("The table was successfully imported", type = "message", duration = 5);
                   },
                   error = function(e){
                     print(e)
                     loginfo(paste("Importer error: ", e, sep=""))
                     showNotification( "An error occured. Please check the tables and types, and try again. ", type = "error", duration = 5 );
                   }
                 )
                 shinyjs::hide("button_importTable");
                 }
             );
  
  #Load file for bulk import
  observeEvent(input$button_chooseFolder,
               {
                 shinyjs::disable("button_chooseFile");
                 shinyjs::disable("button_chooseFolder");
                   #Choose new folder to import. If a folder is not selected, it keeps the same one it had.
                   if (!exists("folderDir")) folderDir <<- "";
                   newFolder <- choose.dir("");
                   if (!is.na(newFolder) && !is.null(newFolder)) {folderDir <<- newFolder}
                   output$out_FolderPath <- renderText({folderDir});
                   
                   #list all the files to add to the project
                   filesToAdd <<- list.files(path = folderDir, all.files = TRUE, full.names = TRUE);
                   filesToAdd <<- filesToAdd[grep("\\.(txt|csv)", filesToAdd)]; #The program can only import .txt and .csv files.
                   if (is_empty(filesToAdd)) 
                   {
                     output$out_filesToLoad <- renderText({HTML("No files to import")});
                     shinyjs::hide("button_importFiles");
                   }
                   else
                   {
                     output$out_filesToLoad <- renderText({HTML(paste(filesToAdd, collapse = "<br/>"))});
                     shinyjs::show("button_importFiles");
                   }
                   shinyjs::enable("button_chooseFile");
                   shinyjs::enable("button_chooseFolder");
               }
  );
  
  #Import Multiple tables
  observeEvent(input$button_importFiles,
               {
                 tables <- list();
                 for (f in filesToAdd)
                 {
                   tab <- read.table(f, header = TRUE, sep=switch(sub('.*\\.', '', f), "txt" = "\t", "csv" = ","), stringsAsFactors = FALSE, comment.char = "", check.names = FALSE, blank.lines.skip = TRUE);
                   
                   colnames(tab)[is.na(colnames(tab)) | colnames(tab) == ""] <- seq(length(colnames(tab)[is.na(colnames(tab)) | colnames(tab) == ""]))
                   
                   #UTF-8 Encoding can sometime leave some undesirable characters in front of the first column
                   colName <- names(tab)[1];
                   if (startsWith(colName, "?..")) {
                     colnames(tab)[1] <- gsub("?..", "", colName, fixed = TRUE);
                   }
                   else if (startsWith(colName, "﻿")) {
                     colnames(tab)[1] <- gsub("﻿", "", colName, fixed = TRUE);
                   }
                   rm(colName);
                   
                   #Remove totally empty columns and rows
                   tab <- tab %>%
                     select_if(function(x) !all(is.na(x) | x == "")) %>%
                     filter_all(any_vars(!is.na(.) & . != ""))
                   
                   tables[[length(tables) + 1]] <- tab;
                 }
                 tryCatch(
                   {
                      importTables(tables, tableID = input$tableDropDown);
                   },
                   error = function(e){
                     print(e)
                     loginfo(paste("Importer error: ", e, sep=""))
                     showNotification( "An error occured. Please check the tables and types, and try again. ", type = "error", duration = 5 );
                   }
                 )
               }
           );
  
  #FILE REQUIREMENTS HINT LINK
  observeEvent(input$link_importGeneTox_hint,
               {
                 shinyjs::toggleElement(id = "panelImportGeneToxHint", anim = TRUE, animType = "fade", time = 0.3);
               }
  );
  
};
