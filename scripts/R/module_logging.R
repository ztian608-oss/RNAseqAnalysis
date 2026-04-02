timestamp_now <- function() {
  format(Sys.time(), "%Y-%m-%d %H:%M:%S")
}

log_info <- function(...) {
  message(sprintf("[%s] %s", timestamp_now(), paste(..., collapse = "")))
}
