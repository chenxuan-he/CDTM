# We cut the abstracts into phrases level, using technique from Deng et al. (2016)
## Deng, Ke, Peter K. Bol, Kate J. Li, and Jun S. Liu. 
## "On the Unsupervised Analysis of Domain-Specific Chinese Texts."
## Proceedings of the National Academy of Sciences 113, no. 22 (May 31, 2016): 6154–59.
## https://doi.org/10.1073/pnas.1516510113.
# Since the process from Deng el al. (2016) should be implemented by C++ code, 
# we store the segmented result at .txt
data <- read.csv("data_all.csv")
data$Num.author <- str_count(data$Authors, ";")+1

# Remove non-English words
documents <- str_replace_all(data$Abstract, "[^[:graph:]]", " ")
documents <- str_replace_all(data$Abstract, "[[:punct:]]", " ")
documents <- gsub("[^[:alpha:]///' ]", " ", documents)
documents <- tolower(documents)
# change plural form into single form
documents <- str_replace_all(documents, "(?<=[^aeiou])ies\\b", "y")
documents <- str_replace_all(documents, "(?<=[aeiou]y)s\\b|(?<=[sxzh])es\\b|(?<=[^sxzh])s\\b", "")
# some manual stemming

# remove stop words and custom high-frequency words
customstopwords <- readLines("HighFreqWords.txt")
documents <- lapply(documents, removeWords, unique(c(stopwords("en"), customstopwords)))

# strip white spaces
documents <- sapply(documents, stripWhitespace)
documents <- gsub("^\\s+|\\s+$", "", documents)

documents <- str_replace_all(documents, "testing", "test")
documents <- str_replace_all(documents, "hypothese", "hypothesis")
documents <- str_replace_all(documents, "hypothesi", "hypothesis")
documents <- str_replace_all(documents, "baye", "bayes")
documents <- str_replace_all(documents, "infectiou", "infectious")
documents <- str_replace_all(documents, "confidence interval ci", "ci")
documents <- str_replace_all(documents, "confidence interval", "ci")
documents <- str_replace_all(documents, "non parametric", "nonparametric")
documents <- str_replace_all(documents, "false discovery rate fdr", "fdr")
documents <- str_replace_all(documents, "false discovery rate", "fdr")
documents <- str_replace_all(documents, "gaussian process gp", "gp")
documents <- str_replace_all(documents, "gaussian process", "gp")

# stemmed abstracts are thrown into C++ code by Deng et al (2016) to generate phrase cutting

data$Abstract <- documents
writeLines(data$Abstract, "WDM package/Abstracts.txt")

#----Run WDM package .cpp file to segment the Abstract----

data$Abstract <- readLines("WDM package/WDM_SegmentedText.txt")
write.csv(data, "data_cut.csv", row.names = FALSE)

