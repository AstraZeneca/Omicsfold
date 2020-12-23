.package.is.installed <- function(package.name) {
  return(requireNamespace(package.name, quietly = TRUE))
}

.install.from.cran <- function(package.name, ...) {
  if (!.package.is.installed(package.name)) {
    cat(paste0("Installing package '", package.name, "' from CRAN...\n"))
    install.packages(package.name, ...)
  }
}

.install.from.bioconductor <- function(package.name, ...) {
  .install.from.cran('BiocManager', ...)

  if (!.package.is.installed(package.name)) {
    cat(paste0("Installing package '", package.name,
               "' from Bioconductor...\n"))
    BiocManager::install(package.name, ...)
  }
}

.install.specific.version <- function(package.name, version, ...) {
  if (!.package.is.installed(package.name)) {
    .install.from.cran('devtools', ...)

    cat(paste0("Installing package '", package.name,
               "' with version '", version, "'...\n"))
    devtools::install_version(package.name, version = version, ...)
  }
}

.install.from.github <- function(package.url, ref = 'master', ...) {
  .install.from.cran('devtools', ...)

  cat(paste0("Installing package with URL '", package.url,
             "' from GitHub...\n"))
  devtools::install_github(package.url, ref = ref, ...)
}

.install.from.local <- function(package.path, force = FALSE, ...) {
  .install.from.cran('devtools', ...)

  cat(paste0("Installing local package in path '", package.path, "'...\n"))
  devtools::install_local(package.path, force = force, ...)
}

install.interfold.bio = function(interfold.path = '.', force = FALSE, ...) {
  # Check the R version is high enough for InterFold.
  version = R.Version()

  major = as.numeric(version$major)
  minor = as.numeric(version$minor)

  if (major < 3 | (minor < 6.0 & major == 3)) {
    stop(paste("R versions below 3.6.0 are not supported.",
               "Please update R by visiting www.r-project.org."))
  }

  # Install the dependencies of InterFold.bio
  .install.from.github('mixOmicsTeam/mixOmics', ref = '2c22e7f', ...)
  .install.from.bioconductor('org.Mm.eg.db', ...)
  .install.from.cran('rlang', ...)
  .install.from.cran('backports', ...)

  if (major >= 4) {
    .install.from.cran('foreign', ...)
  } else {
    .install.specific.version('foreign', '0.8-76', ...)
  }

  source("http://labs.genetics.ucla.edu/horvath/htdocs/CoexpressionNetwork/GeneAnnotation/installAnRichment.R")
  installAnRichment()

  # Install InterFold itself
  .install.from.local(interfold.path, force = force, ...)
}
