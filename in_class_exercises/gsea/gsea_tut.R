#install required R and bioconductor packages
tryCatch(expr = { library("RCurl")}, 
         error = function(e) {  
           install.packages("RCurl")}, 
         finally = library("RCurl"))

# define file paths

gsea_jar <- "/home/rstudio/GSEA_4.3.3/gsea-cli.sh"
rnk_file <- "MesenvsImmuno_RNASeq_ranks.rnk"
working_dir <- "./projects/gsea"
output_dir <- "./projects/gsea"
analysis_name <- "mesen-vs-immuno"
gmt_url <- "https://download.baderlab.org/EM_Genesets/September_01_2025/Human/symbol/"

# Get file

#list all the files on the server
filenames = getURL(gmt_url)
tc = textConnection(filenames)
contents = readLines(tc)
close(tc)

#get the gmt that has all the pathways and does not include terms 
# inferred from electronic annotations(IEA)
#start with gmt file that has pathways only and GO Biological Process only.
rx = gregexpr(
  "(?<=<a href=\")(.*.GOBP_AllPathways_noPFOCR_no_GO_iea.*.)(.gmt)(?=\">)",
              contents, perl = TRUE)
gmt_file = unlist(regmatches(contents, rx))

dest_gmt_file <- file.path(output_dir,gmt_file)

#check if this gmt file already exists
if(!file.exists(dest_gmt_file)){
  download.file(
    paste(gmt_url,gmt_file,sep=""),
    destfile=dest_gmt_file
  )
}

# run GSEA program through command line.
tic()
command <- paste(gsea_jar,  
                 "GSEAPreRanked -gmx", dest_gmt_file, 
                 "-rnk" ,file.path(working_dir,rnk_file), 
                 "-collapse false -nperm 1000 -scoring_scheme weighted", 
                 "-rpt_label",analysis_name,
                 "-plot_top_x 20 -rnd_seed 12345  -set_max 200",  
                 "-set_min 15 -zip_report false",
                 "-out" ,output_dir, 
                 "> ./projects/gsea/gsea_output.txt",sep=" ")

system(command)
toc() #print to see just how long this takes ...
