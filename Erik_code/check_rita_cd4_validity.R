### check cd4 and rita to retain 
x <- readRDS(file = "Erik_code/tips2remove_cd4.rds")
head(t(x))
a <- readRDS(file = "Erik_code/RITASampledCloseToSequence.rds")

head(a)
sum(a$RITASampledCloseToSequence)

load("../phylo-uk/data/sub.RData")
y <- as.data.frame(cbind ( "sequenceID" = as.character(df$seqindex), "rita" = df$ritarecent ), stringsAsFactors = FALSE)
str(y)

b <- merge( a, y, 
           by = "sequenceID",
            all.x = T)
head(b)
table(b$RITASampledCloseToSequence, b$rita)
