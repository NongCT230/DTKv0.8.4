#this script kills all database connections left open by accident!

library(RMySQL);

killMySQLDBConnections <- function () {
  
  all_cons <- dbListConnections(MySQL());
  
  #print(all_cons);
  
  for(con in all_cons) {
    dbDisconnect(con);
  };
  
  print(paste(length(all_cons), " connections killed."));
  
}
