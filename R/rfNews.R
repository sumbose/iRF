rfNews <- function() {
    newsfile <- file.path(system.file(package="iRF"), "NEWS")
    file.show(newsfile)
}
