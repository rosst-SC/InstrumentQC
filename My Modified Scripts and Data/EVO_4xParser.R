# =============================================================================
# SpectroFlo 4.x DailyQC report parser for the Cytek Aurora EVO
#
# The Aurora EVO emits SpectroFlo 4.x DailyQC reports, which the Luciernaga
# parser (QC_FilePrep_DailyQC / DailyQCParse) cannot read. This script parses
# them into the SAME schema as ArchivedData5L.csv (verified column-for-column)
# and:
#   1. writes  data/EVO/Archive/ArchivedDataEVO.csv   (Gain + %rCV + flags)
#   2. appends EVO rows to data/HistoricalData.csv     (for the home summary)
#
# Re-runnable: it rebuilds the EVO archive from whatever DailyQCReport_*.CSV
# files are in data/EVO/, and replaces (not duplicates) EVO rows in
# HistoricalData.csv. A backup of HistoricalData.csv is written alongside it.
#
# NOTE: the MFI/Bead column (BeadDataEVO.csv) is produced separately by
# QCBeadParse from Before/After QC bead .fcs files and is not handled here.
#
# Usage:  Rscript "EVO_4xParser.R"   (run from this folder)
# =============================================================================

suppressMessages({
  library(dplyr); library(stringr); library(tidyr); library(lubridate); library(Luciernaga)
})

QC_FilePrep_DailyQC_4x <- function(x) {
  ReadInfo <- readLines(x, warn = FALSE)
  ReadInfo <- sub("^﻿", "", ReadInfo)            # strip BOM (line 1)

  ## DateTime from filename: DailyQCReport_<serial>_<YYYYMMDD>_<HHMMSS>
  base  <- sub("\\.CSV$", "", basename(x), ignore.case = TRUE)
  parts <- strsplit(sub("^DailyQCReport_", "", base), "_")[[1]]
  DateTime <- ymd_hms(paste(parts[length(parts) - 1], parts[length(parts)]))
  Intro <- data.frame(DateTime = DateTime)

  ## Detector table
  hdr_idx   <- grep("^Laser,Detector", ReadInfo)
  laser_idx <- grep("^Laser Settings", ReadInfo)
  det_header <- gsub(" ", "", strsplit(ReadInfo[hdr_idx], ",")[[1]])
  det_header[2] <- "Detector"                          # was "Detector(nm)"
  det_rows <- ReadInfo[(hdr_idx + 1):(laser_idx - 1)]  # no "Detector Settings" subheader in 4.x
  det_rows <- det_rows[det_rows != ""]
  TheData <- do.call(rbind, lapply(strsplit(det_rows, ","),
                                   function(r) as.data.frame(t(r), stringsAsFactors = FALSE)))
  colnames(TheData) <- det_header
  TheData$Detector <- gsub(" .*", "", TheData$Detector)   # "UV1 (373)" -> "UV1"
  TheData$Gain  <- as.integer(TheData$Gain)
  TheData$`%rCV` <- as.numeric(TheData$`%rCV`)
  TheData <- rename(TheData, DeltaGain = GainChange)       # 4.x "Gain Change"
  TheData$DeltaGain <- as.integer(TheData$DeltaGain)

  # Flags: identical logic to Luciernaga:::QC_FilePrep_DailyQC (incl. SSC precedence)
  Updated <- TheData %>%
    mutate(Baseline = Gain - DeltaGain,
           Comparison = Baseline * 2,
           GainFlag = Gain > Comparison,
           RCVFlag = `%rCV` > 6,
           RCVFlag = case_when(Detector == "SSC" | Detector == "SSC-B" & `%rCV` < 8 ~ FALSE,
                               TRUE ~ RCVFlag))

  wide <- function(df, valcol, prefix = "", suffix) {
    d <- df %>% select(Detector, all_of(valcol))
    d$Detector <- paste0(prefix, d$Detector, suffix)
    pivot_wider(d, names_from = Detector, values_from = all_of(valcol))
  }
  TheGain     <- wide(Updated, "Gain",     suffix = "-Gain")
  GainFlagged <- wide(Updated, "GainFlag", prefix = "Flag-", suffix = "-Gain")
  TheRCV      <- wide(Updated, "%rCV",     suffix = "-% rCV")
  RCVFlagged  <- wide(Updated, "RCVFlag",  prefix = "Flag-", suffix = "-% rCV")
  Assembly <- cbind(Intro, TheGain, GainFlagged, TheRCV, RCVFlagged)

  ## Laser Settings block
  lh_idx  <- laser_idx + 1
  fsc_idx <- grep("^FSC Area Scaling", ReadInfo)
  laser_rows <- ReadInfo[(lh_idx + 1):(fsc_idx - 1)]
  laser_rows <- laser_rows[laser_rows != ""]
  LaserFrame <- do.call(rbind, lapply(strsplit(laser_rows, ","),
                                      function(r) as.data.frame(t(r), stringsAsFactors = FALSE)))
  colnames(LaserFrame) <- strsplit(ReadInfo[lh_idx], ",")[[1]]

  addflags <- function(d) d %>% mutate(across(everything(), as.numeric)) %>%
                               mutate(across(everything(), ~FALSE, .names = "Flag-{.col}"))
  LaserDelay <- LaserFrame %>% select(Laser, `Laser Delay`) %>%
    pivot_wider(names_from = Laser, values_from = `Laser Delay`, names_glue = "{Laser}-Laser Delay") %>% addflags()
  LaserPower <- LaserFrame %>% select(Laser, `Laser Power`) %>%
    pivot_wider(names_from = Laser, values_from = `Laser Power`, names_glue = "{Laser}-Laser Power") %>% addflags()
  LaserArea  <- LaserFrame %>% select(Laser, `Area Scaling Factor`) %>%
    pivot_wider(names_from = Laser, values_from = `Area Scaling Factor`, names_glue = "{Laser}-Area Scaling Factor") %>% addflags()

  ## FSC area scaling, flow rate, temperature (Flow Rate + Temperature share a line in 4.x)
  fsc_line  <- strsplit(ReadInfo[fsc_idx], ",")[[1]]
  flow_line <- strsplit(ReadInfo[grep("^Flow Rate", ReadInfo)], ",")[[1]]
  FSCArea     <- data.frame(FSCAreaScalingFactor = as.numeric(fsc_line[2]),  `Flag-FSCAreaScalingFactor` = FALSE, check.names = FALSE)
  FlowRate    <- data.frame(FlowRate = as.numeric(flow_line[2]),             `Flag-FlowRate` = FALSE,            check.names = FALSE)
  Temperature <- data.frame(Temperature = as.numeric(flow_line[4]),          `Flag-Temperature` = FALSE,        check.names = FALSE)

  cbind(Assembly, LaserDelay, LaserPower, LaserArea, FSCArea, FlowRate, Temperature)
}

# ---- Driver (skipped when this file is sourced, e.g. source(..., local) -----
if (sys.nframe() == 0) {
  MainFolder <- file.path(getwd(), "data")
  reports <- list.files(file.path(MainFolder, "EVO"),
                        pattern = "^DailyQCReport_.*[.]CSV$", full.names = TRUE)
  message("EVO reports found: ", length(reports))

  Parsed <- bind_rows(lapply(reports, QC_FilePrep_DailyQC_4x))

  # 1) ArchivedDataEVO.csv (match the 5L column order if available)
  Archived <- Parsed %>% arrange(desc(DateTime))
  schemaFile <- file.path(MainFolder, "5L", "Archive", "ArchivedData5L.csv")
  if (file.exists(schemaFile)) {
    schema <- colnames(read.csv(schemaFile, check.names = FALSE, nrow = 1))
    if (setequal(colnames(Archived), schema)) Archived <- Archived[, schema]
  }
  archiveDir <- file.path(MainFolder, "EVO", "Archive")
  dir.create(archiveDir, showWarnings = FALSE, recursive = TRUE)
  write.csv(Archived, file.path(archiveDir, "ArchivedDataEVO.csv"), row.names = FALSE)
  message("Wrote ArchivedDataEVO.csv: ", nrow(Archived), " rows")

  # 2) Append EVO rows to HistoricalData.csv (idempotent; backup first)
  histPath <- file.path(MainFolder, "HistoricalData.csv")
  if (file.exists(histPath)) {
    file.copy(histPath, paste0(histPath, ".bak"), overwrite = TRUE)
    H <- read.csv(histPath, check.names = FALSE)
    H$Date <- ymd(H$Date)
    EVO <- Luciernaga:::ShinyQCSummary(x = Parsed, Instrument = "EVO")
    H <- H %>% filter(Instrument != "EVO")
    write.csv(bind_rows(EVO, H), histPath, row.names = FALSE)
    message("Updated HistoricalData.csv: +", nrow(EVO), " EVO rows")
  }
}
