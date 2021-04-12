# DreamTK App 

#v0.8.4
# BER updates
# GenTox updates

# Preliminary loading -----------------------------------------------------
source("./_setup_.R");
#source("./__test_environment.R");

# Main App - R Shiny ------------------------------------------------------
# Support functions -------------------------------------------------------
source("./app/app_support.R");

# Web UI ------------------------------------------------------------------
source("./app/app_ui.R");

# Server function ---------------------------------------------------------
source("./app/app_server.R");

# Start app ---------------------------------------------------------------
app <- shinyApp(ui = ui, server = server);