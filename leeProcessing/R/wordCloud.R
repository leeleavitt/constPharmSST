# load("./leeProcessing/tsSuper.Rdata")

# install.packages("wordcloud")
# library(wordcloud)

# install.packages('tm')
# install.packages('RWeka')
# install.packages('quanteda')


# goTerms <- as.character(tsSuper$gene.go$GO.term.name)
# length(unique(goTerms))
# goTermSummary <- summary(as.factor(goTerms), maxsum=100000)

# ionGoTerms <- grep('voltage', names(sort(goTermSummary,TRUE)[1:300]), value=T)

# goTermSummary[ionGoTerms]

# sort(Gene.Go.Finder("ion channel"))
# sort(Gene.Go.Finder("G-protein"))


# docs <- tm::Corpus(tm::VectorSource(goTerms))

# dtm <- tm::TermDocumentMatrix(docs) 
# matrix <- as.matrix(dtm) 
# words <- sort(rowSums(matrix),decreasing=TRUE) 
# df <- data.frame(word = names(words), freq=words)

# wc <- wordcloud(
#     words = df$word, 
#     freq = df$freq, 
#     min.freq = 1,           
#     max.words=200, 
#     random.order=FALSE, 
#     rot.per=0.35,            
#     colors=brewer.pal(8, "Dark2")
# )

# sample(tsSuper$gene.go[,'GO.term.name'])[1:10]