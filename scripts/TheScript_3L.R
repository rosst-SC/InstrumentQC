# Setup in Correct Directory
setwd("/Users/r.turner/Documents/Positron_Local/InstrumentQC")
getwd()

OS <- getwd()
list.files(OS)

WorkingDirectory <- OS
setwd(WorkingDirectory)

library(stringr)
library(purrr)
library(dplyr)
library(Luciernaga)
library(lubridate)
library(git2r)

# Find out current date
Today <- Sys.Date()
Today <- as.Date(Today)

# Check for Flag Files

AnyFlags <- list.files(WorkingDirectory, pattern="Flag.csv", full.names=TRUE)

if (length(AnyFlags) == 0){

# Git Pull
RepositoryPath <- WorkingDirectory
TheRepo <- git2r::repository(RepositoryPath)
git2r::pull(TheRepo)

# Locating Archive Folder
Instrument <- "3L"
MainFolder <- file.path(WorkingDirectory, "data")
WorkingFolder <- file.path(WorkingDirectory, "data", Instrument)
StorageFolder <- file.path(WorkingFolder, "Archive")

# Gains
Gains <- list.files(StorageFolder, pattern="Archived", full.names=TRUE)
Gains <- read.csv(Gains[1], check.names = FALSE)
LastGainItem <- Gains |> dplyr::slice(1) |> dplyr::pull(DateTime)
LastGainItem <- lubridate::ymd_hms(LastGainItem)
LastGainItem <- as.Date(LastGainItem)
PotentialGainDays <- seq.Date(from = LastGainItem, to = Today, by = "day")
GainRemoveIndex <- which(PotentialGainDays == LastGainItem)
PotentialGainDays <- PotentialGainDays[-GainRemoveIndex]
  
# Usage
Apps <- list.files(StorageFolder, pattern="Application", full.names=TRUE)
Apps <- read.csv(Apps[1], check.names=FALSE)
LastAppsItem <- Apps |> dplyr::slice(1) |> dplyr::pull(DateTime)
LastAppsItem <- lubridate::ymd_hms(LastAppsItem)
LastAppsItem <- as.Date(LastAppsItem)
PotentialAppsDays <- seq.Date(from = LastAppsItem, to = Today, by = "day")
AppsRemoveIndex <- which(PotentialAppsDays == LastAppsItem)
PotentialAppsDays <- PotentialAppsDays[-AppsRemoveIndex]

if (!length(PotentialGainDays) == 0){

SetupFolder <- file.path("C:", "CytekbioExport", "Setup")
TheSetupFiles <- list.files(SetupFolder, pattern="DailyQC", full.names=TRUE)
Dates <- as.character(PotentialGainDays)
Dates <- gsub("-", "", Dates)
GainMatches <- TheSetupFiles[str_detect(TheSetupFiles, str_c(Dates, collapse = "|"))]

if (!length(GainMatches) == 0){
file.copy(GainMatches, WorkingFolder)
walk(.x=Instrument, .f=Luciernaga:::DailyQCParse, MainFolder=MainFolder)
}
} else {message("QC data has already been transferred")
  GainMatches <- NULL
  }

if (!length(PotentialAppsDays) == 0){
    SetupFolder <- file.path("C:", "CytekbioExport")
    TheSetupFiles <- list.files(SetupFolder, pattern="Application", full.names=TRUE)
    MonthStyle <- format(Today, "%Y-%m")
    MonthStyle <- sub("([0-9]{4})-([0-9]{2})", "\\2-\\1", MonthStyle)
    MonthStyle <- gsub("-", " ", MonthStyle)
    MonthStyle <- paste0(MonthStyle, ".txt")
  
    AppMatches <- TheSetupFiles[str_detect(TheSetupFiles, str_c(MonthStyle, collapse = "|"))]
    
    if (!length(AppMatches) == 0){

      if (length(GainMatches) > 0){
      file.copy(AppMatches, WorkingFolder)
      walk(.x=Instrument, .f=Luciernaga:::AppQCParse, MainFolder=MainFolder)
      }
      }
} else {message("QC data has already been transferred")
    AppMatches <- NULL
    }

if (any(length(PotentialGainDays)|length(PotentialAppsDays) > 0)){

  if (length(GainMatches) > 0){
  # Stage to Git
    git2r::add(TheRepo, "*")
  
  TheCommitMessage <- paste0("Update for ", Instrument, " on ", Today)
  git2r::commit(TheRepo, message = TheCommitMessage)
  cred <- git2r::cred_token(token = "GITHUB_PAT")
  git2r::push(TheRepo, credentials = cred)
  message("Done ", Today)
  } else {message("No files to process ", Today)}
} else {message("No files to process ", Today)}
} else {message("Automation Skipped ", Today)}