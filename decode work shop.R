getwd()
setwd()
library(ggplot2)
BiocManager::install("Biostrings")
BiocManager::install("rentrez")


library("Biostrings")
library("rentrez")
organism_name = "Allobates kingsburyi[Organism]"
print(organism_name)

organism_id = entrez_search(db="nuccore",term=organism_name)
print(organism_id)

str(organism_id)
print(organism_id$ids)

sequence_data=entrez_fetch(db="nuccore",id=organism_id$ids,rettype="fasta")
print(sequence_data)


# variable
#data type: integer, float, string, character, binary
x=3
y=3.8
my_name="srisai"

print(x)
print(paste(my_name,x,y))
print(paste("Name of the candidate is:",my_name))


# data structure/container 
# vector - 1D

class_students=c("kushi","ritam","sexa","abhi","nalini","sakshi")
print(class_students)

#index - 1,2,3,4....
#fetching - []

print(class_students[1])
print(class_students[6])

print(class_students[2:5])
print(class_students[c(1,3,5)])

print(class_students[-1])
print(class_students[c(-1,-5)])

numerical_vec=c(22,44,66,88)
sum(numerical_vec)
log(numerical_vec)
mean(numerical_vec)
median(numerical_vec)
sd(numerical_vec)
min(numerical_vec)
max(numerical_vec)


# Functions - (a) Library (b) User-Defined 
# set of instructions to perform specific task
# name-of-function(arguments)
# ex: print("Avinash")
# Sometimes there is no argument 
getwd()


# Dataframe - 2D Data structure
# collection of vectors

# Lets create a Dataframe
class_students=c("kushi","ritam","sexa","abhi","nalini","sakshi")
score_students=c(65,89,76,95,93,97)
qualification_students=c("master","PhD","BSc","master","PhD","BSc")

students_data=data.frame("Name"=class_students,"Scores"=score_students,"Qualification"=qualification_students)

View(students_data)
view(students_data)
print(score_students)
print(qualification_students)

print(paste("The dimension of the student database is: ",dim(students_data)))
setwd()
dir()
getwd()

data_trail=read.csv("LR.csv")
View(data_trail)

library(dplyr)
tem_data=data_trail %>%
group_by(zipcode) %>%
summarise(Average_Price = mean(price)) %>%
arrange(desc(Average_Price))
View(tem_data)


# 1. Remove NA/Blank
dim(data_trail)[1]
num_rows=dim(data_trail)[1]
num_cols=dim(data_trail)[2]
Value=10
is.na(Value)

library(rentrez)

query = "Ameerega bilinguis[Organism] AND COI[Gene]"

query_Id=entrez_search(db="nuccore",term = query)
str(query_Id)
query_Id$ids

my_sequence=entrez_fetch(db="nuccore",id=query_Id$ids,rettype="fasta")
my_sequence

class(my_sequence)
length(my_sequence)

is.character(my_sequence)
nchar(my_sequence)

my_sequence_split = strsplit(my_sequence,split = '\n')
my_sequence_split     

class(my_sequence_split)
my_sequence_split[[1]][2:2]


my_clean_seq = paste(my_sequence_split[[1]][2:23],collapse="")

my_clean_seq
class(my_clean_seq)


big_query=c("Ameerega parvula[Organism] AND COI[Gene]", 
            "Epipedobates anthonyi[Organism] AND COI[Gene]",
            "Oophaga sylvatica[Organism] AND COI[Gene]")

big_query
my_collected_sequence=list()
for (i in 1:length(big_query)){
  query_ID=entrez_search(db="nuccore",term = big_query[i])
  my_sequence=entrez_fetch(db="nuccore",id=query_ID$ids,rettype="fasta")
  my_collected_sequence[[i]]=my_sequence
}

my_collected_sequence

length(my_collected_sequence)

my_clean_seq
pattern_looking_for="GTTTTGCATGATA"

grep(pattern=pattern_looking_for,my_clean_seq,ignore.case = TRUE)




today_seq="GTGATAATTACTCGATGATTATTTTCTACCAACCACAAAGACATCGGAACTTTATACCTAGTGTTTGGGGCATGAGCAGGCATAGTCGGCACTGCTCTTAGCCTTTTAATTCGAGCCGAATTAAGCCAGCCCGGGTCCTTACTAGGCGATGACCAGATCTACAACGTTATTGTTACCGCCCATGCTTTCGTTATAATCTTTTTTATAGTAATGCCAATTCTAATCGGTGGCTTTGGGAATTGATTAGTGCCCCTAATAATTGGAGCCCCAGACATAGCTTTTCCCCGAATAAACAATATGAGCTTTTGGCTTCTTCCCCCCTCTTTCCTACTACTCCTAGCATCCGCAGGCGTTGAAGCAGGCGCCGGTACTGGCTGAACTGTGTACCCTCCCCTTGCAGGCAACCTAGCTCATGCTGGCCCATCAGTTGATTTAACTATTTTTTCACTTCATCTCGCCGGTGTTTCTTCTATTCTAGGGGCAATTAACTTTATTACAACAACCTTAAACATAAAACCCCCTTCATTAACACAATATCAAACCCCATTATTTGTCTGATCTGTATTAATTACTGCAGTCCTTCTTCTTCTCTCCCTCCCAGTTCTGGCTGCCGGAATCACTATACTCTTGACTGACCGAAACCTAAACACCACCTTCTTTGACCCAGCAGGTGGAGGCGACCCTGTCCTGTACCAACACCTGTTCTGATTCTTTGGTCACCCCGAAGTCTACATCCTTATCCTGCCTGGATTTGGTATCATCTCCCATGTTGTCACATTCTACTCTAGCAAAAAAGAACCCTTCGGCTATATAGGAATAGTCTGAGCTATAATATCGATTGGTCTCCTAGGTTTCATTGTTTGAGCTCACCACATATTCACAACAGACCTTAATGTAGACACTCGAGCCTACTTTACCTCAGCTACTATAATCATCGCTATCCCAACAGGTGTCAAAGTCTTTAGCTGACTTGCCACCATGCACGGAGGAATTATTAAATGAGACGCCGCCATATTATGGGCTCTCGGATTCATCTTTTTATTTACAGTTGGAGGACTAACTGGAATCGTTTTAGCCAACTCCTCTTTAGACATTGTTTTGCATGATACATATTATGTAGTAGCCCACTTTCACTACGTTCTTTCTATGGGGGCAGTATTTGCCATTATAGCCGGCTTCGTACACTGATTTCCTCTCTTTTCCGGATTTACCCTTCATGAAGCCTGAACAAAAATTCAATTTGGCGTCATATTTACCGGCGTAAATTTAACATTCTTCCCCCAGCATTTCTTAGGTCTCGCAGGCATGCCTCGACGTTATTCAGACTACCCTGACGCCTACACATTATGAAACACCGTTTCATCAATCGGCTCTTTAATCTCTCTAGTTGCAGTAATCATTATGATGTTTATCATTTGAGAAGCTTTCTCTTCCAAACGCCTACCTCTACCTGCAGAAATAACCCCAACTAATGTAGAATGATTATACGGATCCCCCCCACCTTACCACACTTTTGAGGAAGCCGTTTACTCCAAAATT"

today_seq
pattern_looking_for="GTTTTGCATGATA"

grep(pattern=pattern_looking_for,my_clean_seq,ignore.case = TRUE)

sequence_vec=c()
sequence_vec[1]=today_seq

today_seq_2="AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
sequence_vec[2]=today_seq_2

sequence_vec[3]="GTTTTGCATGATAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAGTTTTGCATGAT"

grep(pattern=pattern_looking_for,sequence_vec,ignore.case = TRUE)

new_today_seq = sub(pattern_looking_for,"QQQQQ",today_seq)
new_today_seq
pattern_looking_for="QQQQQ"
grep(pattern=pattern_looking_for,new_today_seq,ignore.case = TRUE)

P1="QQQQQ"
combined_pattern=paste0(pattern_looking_for,"|",P1)
combined_pattern
new_today_seq
today_seq
grep(pattern=pattern_looking_for,new_today_seq,ignore.case = TRUE)

#comparing, ordering, and simplifying vectors

my_vector_1 = c("Allobates", "Epipedobates", "Dendrobates", "Hyloxalus", "Mannophryne", "Hyloxalus", "Silverstoneia")
my_vector_2 = c("Epipedobates", "Epipedobates","Dendrobates", "Hyloxalus", "Mannophryne", "Hyloxalus", "Rheobates")

unique(my_vector_1)
unique(my_vector_2)

union(my_vector_1,my_vector_2)
intersect(my_vector_1,my_vector_2)

setdiff(my_vector_1,my_vector_2)
setdiff(my_vector_2,my_vector_1)

setequal(my_vector_1,my_vector_2)
is.element("Allobates",my_vector_1)
is.element("Allobates",my_vector_2)

library(stringr)
today_seq

adenine=str_count(today_seq,pattern="A")
cytosine=str_count(today_seq,pattern="C")
guanine=str_count(today_seq,pattern="G")
thymidine=str_count(today_seq,pattern="T")

adenine
cytosine
guanine
thymidine

count_nucleotide_DB=data.frame("Adenine"=adenine,
                               "Cytosine"=cytosine,
                               "Guanine"=guanine,
                               "Thymidine"=thymidine)
View(count_nucleotide_DB)

adenine=str_count(sequence_vec,pattern="A")
adenine

gc=str_count(sequence_vec,pattern="GC")

count_nucleotide_DB=data.frame("Adenine"=adenine,
                               "Cytosine"=cytosine,
                               "Guanine"=guanine,
                               "Thymidine"=thymidine,
                               "GC Content"=gc)
View(count_nucleotide_DB)

restiction_site=str_locate_all(today_seq,pattern="CTAG")

restiction_site


#uniprot(https://www.uniprot.org/)
library(UniprotR)

uniprot_id= c("P17787", "Q9ERK7", "P09484", "A0A3Q1MJN8", "E7F4S7", "H3B5V5", "H9GMJ0", "A4IIS6")

prot_Seq=UniprotR::GetSequences(uniprot_id)

class(prot_Seq)
View(prot_Seq)

prot_Seq$Sequence
prot_Seq[1]



View(iris)
unique(iris$Species)

obs=nrow(iris)
train=sample(obs,0.5*obs)

library(dplyr)
test=obs%>%
  seq_len()%>%
  setdiff(train)

length(train)
length(test)

traindata=iris[train,]
testdata=iris[test,]

table(traindata$Species)
table(testdata$Species)

dim(traindata)
dim(testdata)

model_tree=rpart(Species~.,data=traindata,method= 'class')
model_tree=rpart(Species~.,data=traindata,method = 'class')
prediction_tree=predict(model_tree,testdata)
prediction_tree=predict(model_tree,testdata,type='class')
prediction_tree
observed_species=testdata$Species

table(prediction_tree,observed_species)










