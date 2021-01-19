#' Map a list of symbols to their corresponding Entrez IDs. This wraps the
#' method from the anRichmentMethods library to convert a list of symbols to
#' entrez IDs.
#'
#' @param symbols List of gene symbols.
#' @param organism The organism the symbols relate to. This can be any of the
#' supported common names, scientific names or shorthand shown when calling
#' \code{anRichmentMethods::organismLabels()}.
#'
#' @return Data frame containing pairs of symbols and Entrez IDs.
#' @export
#'
#' @examples
#' \dontrun{
#' gene.symbols.and.ids <- map.symbols.entrez.ids(gene.symbols, "mouse")
#' }
map.symbols.entrez.ids <- function(symbols, organism) {
  gene.symbols <- as.data.frame(symbols)
  colnames(gene.symbols) <- c("symbol")
  gene.symbols$entrez.id <- anRichmentMethods::convert2entrez(organism, symbols)

  return(gene.symbols)
}

#' Utility method to retrieve all gene annotated with a particular GO term, or
#' any of its children.
#'
#' @param go.id GO term.
#'
#' @return Symbol of gene annotated with this term.
#' @export
#'
getGOGenes <- function(go.id) {
  go.ids <- GOBPOFFSPRING[[go.id]]
  go.ids <- append(go.ids, go.id)
  allegs <- sapply(go.ids, function (x) try(get(x, org.Mm.egGO2ALLEGS), silent = TRUE) )
  allegs <- unlist(allegs, recursive = TRUE)
  allegs <- allegs[!str_detect(allegs, "Error", negate = FALSE)]
  genes = unlist(mget(allegs, org.Mm.egSYMBOL))
  return(genes)
}


#' Convert a list of mouse genes to their human orthologues.
#'
#' @param mouse.genes List of mouse gene symbols.
#'
#' @return Human orthologue gene symbols (not necessarily ordered).
#' @export
#'
convertMouseGeneList <- function(mouse.genes) {
  human = biomaRt::useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl", mirror = "useast")
  mouse = biomaRt::useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl", mirror = "useast")
  genesV2 = biomaRt::getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = mouse.genes , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[, 2])
  return(humanx)
}
