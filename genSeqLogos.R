setwd("~/Documents/writing/MFS/MammCell paper/scripts/Ribozyme-SeqFun/")
library(ggseqlogo)


tmp=list.files(pattern = "theosubseqs_*")
myfiles=lapply(tmp,readLines)

ggseqlogo(myfiles,seq_type='RNA',ncol=1)


tmp=list.files(pattern = "./xansubseqs_*")
myfiles=lapply(tmp,readLines)
myfiles
ggseqlogo(myfiles,seq_type='RNA',ncol=1)

tmp=list.files(pattern = "FAsubseqs_*")
myfiles=lapply(tmp,readLines)
myfiles
ggseqlogo(myfiles,seq_type='RNA',ncol=1)



tmp=list.files(pattern = "cdGsubseqs_*")
myfiles=lapply(tmp,readLines)
myfiles
ggseqlogo(myfiles,seq_type='RNA',ncol=1)


tmp=list.files(pattern = "cdGIsubseqs_*")
myfiles=lapply(tmp,readLines)
myfiles
ggseqlogo(myfiles,seq_type='RNA',ncol=1)







tsnefiles=list.files(pattern = "*xantsne.txt")
myclusterfiles=lapply(tsnefiles,readLines)
ggseqlogo(myclusterfiles,seq_type='RNA',ncol=1)

ggseqlogo(myfiles[14],seq_type='RNA',ncol=1)


