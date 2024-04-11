#' @param data original data, a data frame
#' @param covariates_name covariates name
#' @param text_name name of the text to be analyzed
#' @param split_by split the text by
#' @param doc.thresh lowest number of words a document should include
#' @param word.lower.thresh The lowest frequency of a word to be included
#' @param word.upper.thresh The highest frequency of a word to be included
data_preprocess <- function(data, covariates_name, text_name, time_slices, split_by, 
                            doc.thresh = 5, word.lower.thresh=40, word.upper.thresh=Inf) {
  documents <- data[,which(colnames(data)==text_name)]
  # change list to tm-type vectors
  txt <- VCorpus(VectorSource(documents))
  # incorporate covariates
  metadata <- data[, c(time_slices, covariates_name)]
  # force the time index name into "Time_slices"
  colnames(metadata)[1] <- "Time_slices"
  for (i in 1:ncol(metadata)) {
    NLP::meta(txt, colnames(metadata)[i]) <- metadata[, i]
  }
  # change it to matrix and remove documents which are shorter than 
  dtm <- DocumentTermMatrix(txt, control = list(wordLengths = c(doc.thresh, Inf),
                                                tokenize = function(x) 
                                                  gsub("^\\s+|\\s+$", "", 
                                                       unlist(strsplit(as.character(x), split_by)))))
  docindex <- unique(dtm$i)
  metadata <- meta(txt)[docindex, , drop = FALSE]
  processed <- read.slam(dtm)
  kept <- (1:length(documents) %in% unique(dtm$i))
  
  vocab <- processed$vocab
  processed <- list(
    documents=processed$documents,
    vocab=vocab, 
    meta=metadata, 
    docs.removed=which(!kept))
  
  docs.removed <- processed$docs.removed
  documents <- processed$documents
  vocab <- processed$vocab
  meta <- processed$meta
  # count the words
  wordcounts <- tabulate(dtm$j)
  # a thresh to control the frequency of words
  # to detect the words which should be removed
  toremove <- which(wordcounts <= word.lower.thresh | wordcounts >= word.upper.thresh)
  keepers <- which(wordcounts > word.lower.thresh & wordcounts < word.upper.thresh)
  droppedwords <- vocab[toremove]
  vocab <- vocab[-toremove]
  remap <- 1:length(keepers)
  
  for (i in 1:length(documents)) {
    doc <- documents[[i]]
    dockeep <- doc[1, ] %in% keepers
    doc <- doc[, dockeep, drop = FALSE]
    doc[1, ] <- remap[match(doc[1, ], keepers)]
    documents[[i]] <- doc
    if (ncol(doc) == 0) 
      docs.removed <- c(docs.removed, i)
  }
  
  documents <- lapply(documents, function(x) matrix(as.integer(x), nrow = 2))
  
  out <- list(
    documents = documents, 
    vocab = vocab, 
    meta = meta, 
    words.removed = droppedwords, 
    tokens.removed = sum(wordcounts[toremove]),
    wordcounts = wordcounts)
  
  return(out)
}


