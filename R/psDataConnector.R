#' Connect to Netezza
#'
#' Establishes a connection to the Netezza database
#'
#' This function creates a connection to the Netezza database using the supplied credentials. This connection
#' can then be used to pull data from the database into R.
#'
#' @param user.name Your user name (prid)
#' @param password Your password
#' @param driver.location Location of your netezza driver, usually "nzjdbc.jar". Can be specified as complete
#' path or the directory where the driver file is located.
#' @param db.prefix This + prid = your database
#' @param db.host URL of netezza appliance
#' @return Database connection object
#'
#' @examples
#' \dontrun{
#' ps.connect.to.netezza("user123", "password123", "/users/user123/Documents/databaseDrivers")
#' }
#' @export
ps.connect.to.netezza <- function(user.name, password, driver.location, db.prefix, db.host) {

  filesep = .Platform$file.sep # \windows /unix

  # if only folder is provided, add filename
  if (substr(driver.location, nchar(driver.location)-3, nchar(driver.location)) != ".jar"){

    if(substr(driver.location, nchar(driver.location), nchar(driver.location)) != filesep){
      driver.location <- paste(driver.location, filesep, sep="")
    }

    driver.location <- paste(driver.location, "nzjdbc.jar", sep="")
  }

  if (!file.exists(driver.location)) {
    stop(paste("Specified driver file (", driver.location, ") does not exist.", sep = ""))
  }

  # Load the driver
  drv <- RJDBC::JDBC(driverClass="org.netezza.Driver", classPath=driver.location, identifier.quote="'")

  # Configure the db connection string
  db <- paste("jdbc:netezza:/", db.host, "/", db.prefix, user.name, sep="")

  # Create the connection
  conn <- RJDBC::dbConnect(drv, db, user.name, password)
  return(conn)
}


#' Disconnect from Netezza
#'
#' Close a connection to Netezza, freeing resources
#'
#' This function closes an existing connection to the netezza appliance to free resources
#'
#' @param conn Database connection object created by ps.connect.to.netezza
#' @return NULL
#'
#' @examples
#' \dontrun{
#' ps.disconnect.from.netezza(nzconn)
#' }
#' @export
ps.disconnect.from.netezza <- function(conn) {
  RJDBC::dbDisconnect(conn)
}


#' SQL to Data.Frame
#'
#' Execute a sql query and return the results as a data.frame
#'
#' This function uses a database connection to execute a sql query string. The results returned by the query
#' are returned as an R data frame.
#'
#' @param conn Database connection object created by ps.connect.to.netezza
#' @param sql.query SQL query string, just as you might type into Aginity
#' @return Calculated accuracy of propensity scores
#'
#' @examples
#' \dontrun{
#' myData <- ps.sql2df(nzconn, "select * from database.table limit 100")
#' }
#' @export
ps.sql2df <- function(conn, sql.query) {

  res <- RJDBC::dbSendQuery(conn, sql.query)
  res.df <- RJDBC::fetch(res, n=-1)
  RJDBC::dbClearResult(res)

  return(res.df)
}
