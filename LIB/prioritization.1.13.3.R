# Generalized file for variant filtering and prioritization

try(sink())
try(sink())
try(sink())

# Section 0 - Libraries ---------------
library(biomaRt)
library(GenomicRanges)
library(rtracklayer)
library(jsonlite)
library(knitr)
library(dplyr)

# questa funzione è per annotare le varianti con le informazioni sul fenotipo,
#  secondo ensembl
annotate.vars<-function(var_ids) {
  library(biomaRt)
  variation = useEnsembl(biomart="snp", dataset="hsapiens_snp")
  res<-getBM(attributes=c('clinical_significance','phenotype_description','variation_source','p_value','pmid'), 
             filters = 'snp_filter', values =var_ids, mart = variation)
  return(res)
}

# questa funzione serve per imputare il genotipo
#  NA si può considerare HOMREF
# Accetta come argomento la string GT
do.zyg.lax<-function(GT) {
  res<-unlist(strsplit(x=GT,split="[|/]"))
  if (length(res)!=2) { return("HOMREF") } 
  if (res[1]==res[2]) {
    if(res[1] %in% c("0",".")) { return("HOMREF") } 
    else return("HOM") }
  else {return ("HET")} 
}

do.zyg<-function(GT) {
  res<-unlist(strsplit(x=GT,split="[|/]"))
  if (length(res)!=2) { return(NA) } 
  if (res[1]==res[2]) {
    if(res[1]=="0") { return("HOMREF") } 
    else if (res[1]==".") { return(NA) } 
    else return("HOM") }
  else {return ("HET")} 
}


# legge l'elenco dei geni candidati e lo trasforma in una lista di vettori di stringhe
#  l'input è una cartella con i file dei pannelli. 
#  Questi files sono composti da un elenco di gene symbol
# l'output è una lista di vettori, ogni vettore è un pannello, in nome del vettore è il nome del file
read.panel<-function(candidate_gene_file) {
  panels<-list()
  for (panel in candidate_gene_file) {
    cat(" . reading "); cat(panel); cat("\n")
    panels[[panel]]<-readLines(panel)
  }
  return(panels)
}

read.variants<-function(candidate_variant_file) {
  variants<-list()
  for (vars in candidate_variant_file) {
    cat(" . reading "); cat(vars); cat("\n")
    variants[[vars]]<-import.bed(con=vars)
  }
  return(variants)
}

get.file.name<-function(x) {s=unlist(strsplit(x,"/")); return(s[length(s)])}




# Questa funzione classifica i casi in base 
#  ai modelli di trasmissione e alla segregazione nel pedigree
get.cases_by_segregation<-function(ped) {
  ped[,6]<-as.numeric(ped[,6])
  if (min(ped[,6],na.rm=T)==0 & max(ped[,6],na.rm=T)==1) {
    cat("\nped range of status column between 0 and 1, increased to +1/n")
    ped[,6]<-ped[,6]+1
  }
  dominant.cases<-ped[ped[,6]==2,2]
  dominant.controls<-ped[ped[,6]==1,2]
  
  matcases<-ped[ped[,2] %in% ped[ped[,6]==2,4],2]
  patcases<-ped[ped[,2] %in% ped[ped[,6]==2,3],2]
  
  recessive.cases<-ped[ped[,6]==2,2]
  recessive.controls<-ped[ped[,6]==1,2]
  recessive.carriers<-ped[ped[,6]==1 & ped[,2] %in% c(matcases,patcases),2]
  recessive.carriers.paternal<-ped[ped[,6]==1 & ped[,2] %in% c(patcases),2]
  recessive.carriers.maternal<-ped[ped[,6]==1 & ped[,2] %in% c(matcases),2]
  
  xlinked.cases<-ped[ped[,6]==2,2]
  # i carrier sono femmine genitori di un caso
  xlinked.carriers<-ped[ped[,6]==1 & ped[,2] %in% matcases & ped[,5]==2,2]
  # i controlli sono maschi non affetti
  xlinked.controls<-ped[ped[,6]==1 & ped[,5]==1 ,2]
  
  return(list(samples=ped[,2],
              dominant.cases=dominant.cases,
              dominant.controls=dominant.controls,
              recessive.cases=recessive.cases,
              recessive.controls=recessive.controls,
              recessive.carriers=recessive.carriers,
              recessive.carriers.paternal=recessive.carriers.paternal,
              recessive.carriers.maternal=recessive.carriers.maternal,
              xlinked.cases=xlinked.cases,
              xlinked.controls=xlinked.controls,
              xlinked.carriers=xlinked.carriers,
              mothers.of.affected=matcases,
              fathers.of.affected=patcases,
              parents.of.affected=c(matcases,patcases)))
}


# questa funzione classifica le varianti in base alla presenza in determinati casi
#  equivale ad un AND
do.segregation<-function(samples,genotypes) {
  if (is.null(dim(samples))) {
    if (!is.na(length(samples))) {
      return (samples %in% genotypes)
    }
      return (NA)
  } 
  if (ncol(samples)>0) {
    return (apply(samples,1,function(x) sum(x %in% genotypes)==length(x))) 
    } 
  else {return(rep(FALSE,nrow(samples)))}
}

# questa funzione classifica le varianti in base alla presenza in determinati casi
#  equivale ad un AND
# questa versione 2 accetta penetranza incompleta, 
#  cioè un numero di controlli pari a PenetranceTollerance con la variante 
#  PenetranceTollerance è solo per i controlli!
do.segregation2<-function(samplesGenotype,queryGenotypes,PenetranceTollerance=0) {
  # se è un solo individuo
  if (is.vector(samplesGenotype)) {
    # se è ammessa penetranza incompleta mette tutti FALSE (è un falso controllo)
    if (PenetranceTollerance>0) {
      return(rep(FALSE,length(samplesGenotype))) 
    }
    # controlla se i samplesGenotype siano i queryGenotypes
    return (samplesGenotype %in% queryGenotypes)
  } 
  # se è nessun o più individui
  else if (is.data.frame(samplesGenotype)) {               
    # se è più individui  
    if (ncol(samplesGenotype)>0) {
      # conta il nuimero di samplesGenotype che sono queryGenotypes
      #  e devono essere tutti TRUE meno il numero di PenetranceTollerance
      return (apply(samplesGenotype,1,function(x) sum(x %in% queryGenotypes)>= (length(x)-PenetranceTollerance)  )) 
    } 
    # se è nessun individuo
    else {
      return(rep(FALSE,nrow(samplesGenotype)))
    }
  }
}



# function per trovare i geni comp_het
# gtab è la tabella di UN SOLO GENE, con le varianti che sono state filtrate
#  questa tabella è filtrata precedentemente con il vettore AR (vedi dopo)
# gmat è il numero della colonna con il genotipo materno
# gpat è il numero della colonna con il genotipo paterno
do_comp_het_index<-function(gtab,gmat,cpat) {
  m=sum(gtab[cmat] == "HET")>0
  p=sum(gtab[cpat] == "HET")>0
  return (m & p)
}


# Questa funzione restituisce i geni che hanno varianti comp-het
#   accetta in input la lista delle varianti che già 
#   sono state filtrate per essere eterorizoti nel probabo e in un solo genitore
#   le varianti possono anche essere filtrate per altri fattori
# EXAMPLE: tab<-outxt[ outxt$AR & outxt$rare & outxt$CompTrio & outxt[,cfilt]=="PASS",]	
find_comp_het_genes<-function(tab) {
  comphet_gene<-c()
  for (g in unique(tab[,cgene])) {
    gtab<-tab[ tab[,cgene]==g, ]
    if (sum(gtab$HTm)>1 & sum(gtab$HTp)>1) comphet_gene<-c(comphet_gene,g)
  }
  results<-list()
  results[["variants"]]<-tab[tab[,cgene] %in% comphet_gene,]
  results[["genes"]]<-comphet_gene
  return (results)
}


# Report
# EXAMPLE: tab<-outxt[outxt$known & outxt$XL & outxt[,cfreq]<0.5,]
do_report<-function(tab){
  wtab<-tab[order(tab[,1],tab[,2]),]
  
  #for (c in 1:nrow(tab)) cat(do_row(tab,c))
  genes<-unique(tab[,cgene])
  
  results<-list()
  results[["genes"]]<-genes
  results[["tab"]]<-tab
  
  return(results)
}

# conta il numero di campioni nella tabella "variants"
#  che hanno il genotipo "genotype"
#  per i soggetti dei membri di "segregation" del tipo "fam_members"
call.cases<-function(variants,segregt,fam_memb,genotype) {
  c=which(names(segregt)==fam_memb)
  tab<-data.frame(variants[,segregt[[c]]])
  res=as.vector(apply(X = tab,1,function(x) sum(x==genotype)))
  return(res)
}

# crea la variante genomica in versione HGVS
make.HGVS<-function(tab,cchrom,cpos,cref,calt) {
  return(paste0("g.",tab[cchrom],":",as.numeric(tab[cpos]),tab[cref],">",tab[calt]))
}


# per print.geneinfo
do_row<-function(tab,c,segregation) {
  g_prob<-paste("Aff=",paste(tab[c,segregation$dominant.cases],collapse=","),sep="")
  g_F<-paste("F=",paste(tab[c,segregation$fathers.of.affected],collapse=","),sep="")
  g_M<-paste("M=",paste(tab[c,segregation$mothers.of.affected],collapse=","),sep="")
  #  g_con<-paste("Unaff=",paste(tab[c,segregation$recessive.controls],collapse=","),sep="")
  g_con=paste0("Nctrls=",as.character(sum(tab[,segregation$dominant.controls] %in% c("HET","HOM"))))
  g_var<-paste0("g.",tab[c,cchrom],":",tab[c,cpos],tab[c,cref],">",tab[c,calt],"\t",tab[c,caa],"\t",tab[c,ccds])
  g_var2<-paste0(tab[c,cdbSNP],"\t",tab[c,cfreq],"\t",tab[c,varianteffect_cname])
  g_var3<-paste(tab[c,cdiseasevar],collpase="\t")
  g_var4<-paste(tab[c,cphenotypes],collapse = "\t")
  g_var5<-paste0(tab$candidate[c])
  
  this_row<-paste(g_prob,g_F,g_M,g_con,g_var,g_var2,g_var3,g_var4,"\n",sep="\t")
  return(this_row)
}



do.row2<-function(tab,c,segregation) {
  print (c)
  g_prob<-paste(tab[c,segregation$dominant.cases],collapse=",")
  g_F   <-paste(tab[c,segregation$fathers.of.affected],collapse=",")
  g_M   <-paste(tab[c,segregation$mothers.of.affected],collapse=",")
  g_con <-sum(tab[c,segregation$dominant.controls] %in% c("HET","HOM"))
  g_var <-paste0("g.",tab[c,cchrom],":",tab[c,cpos],tab[c,cref],">",tab[c,calt])
  g_var1<-tab[c,caa]
  g_var2<-tab[c,ccds]
  g_var3<-tab[c,cdbSNP]
  g_var4<-tab[c,cfreq]
  g_var5<-tab[c,varianteffect_cname]
  g_var6<-paste(tab[c,cdiseasevar],collapse=";")
  g_var7<-paste(tab[c,cphenotypes],collapse = ";")
  g_var8<-paste(tab$candidate[c],collapse = ";")
  
  this_row<-c(g_prob,g_F,g_M,g_con,g_var,g_var1,g_var2,g_var3,g_var4,g_var5,g_var6,g_var7,g_var8)
  #names(this_row)<-c("Affected","Fathers","Mothers","Num_Controls_var","variant","amino_acid","cds","RS","Freq","Effect","Clinical_relevance","Phenotype","Panels")
  return(this_row)
}

do.row3<-function(tab,segregation) {
 g<-data.frame(Affected=character(),Fathers=character(),Mothers=character(),
                  Unaffected=numeric(),variant=character(),amino_acid=character(),
                  cds=character(),RS=character(),Freq=numeric(),
                  Effect=character(),Clinical_relevance=character(),Phenotype=character(),
                  Panels=character())
  
  g[1,1]<-paste0("HOM=",sum(tab[segregation$dominant.cases] %in% c("HOM"))," HET=",sum(tab[segregation$dominant.cases] %in% c("HET")))
  g[1,2]<-paste0("HOM=",sum(tab[segregation$fathers.of.affected] %in% c("HOM"))," HET=",sum(tab[segregation$fathers.of.affected] %in% c("HET")))
  g[1,3]<-paste0("HOM=",sum(tab[segregation$mothers.of.affected] %in% c("HOM"))," HET=",sum(tab[segregation$mothers.of.affected] %in% c("HET")))
  g[1,4]<-paste0("HOM=",sum(tab[segregation$dominant.controls] %in% c("HOM"))," HET=",sum(tab[segregation$dominant.controls] %in% c("HET")))
  g[1,5]<-paste0("g.",tab[cchrom],":",tab[cpos],tab[cref],">",tab[calt])
  g[1,6]<-tab[caa] 
  g[1,7]<-tab[ccds]
  g[1,8]<-tab[cdbSNP]
  g[1,9]<-tab[cfreq]
  g[1,10]<-tab[varianteffect_cname]
  g[1,11]<-paste(tab[cdiseasevar],collapse=";")
  g[1,12]<-paste(tab[cphenotypes],collapse = ";")
  g[1,13]<-paste(tab[names(tab)=="candidate"],collapse = ";")
  
  #this_row<-c(g_prob,g_F,g_M,g_con,g_var,g_var1,g_var2,g_var3,g_var4,g_var5,g_var6,g_var7,g_var8)
  #names(this_row)<-c("Affected","Fathers","Mothers","Num_Controls_var","variant","amino_acid","cds","RS","Freq","Effect","Clinical_relevance","Phenotype","Panels")
  return(g)
}


do.row4<-function(tab,segregation,printcolumn) {
 g<-data.frame(Affected=character(),Fathers=character(),Mothers=character(),
                  Unaffected=numeric(),variant=character(),amino_acid=character(),
                  cds=character(),RS=character(),Freq=numeric(),
                  Effect=character(),Clinical_relevance=character(),Phenotype=character(),
                  Panels=character())
  
  
  g[1,1]<-paste0("HOM=",sum(tab[segregation$dominant.cases] %in% c("HOM"))," HET=",sum(tab[segregation$dominant.cases] %in% c("HET")))
  g[1,2]<-paste0("HOM=",sum(tab[segregation$fathers.of.affected] %in% c("HOM"))," HET=",sum(tab[segregation$fathers.of.affected] %in% c("HET")))
  g[1,3]<-paste0("HOM=",sum(tab[segregation$mothers.of.affected] %in% c("HOM"))," HET=",sum(tab[segregation$mothers.of.affected] %in% c("HET")))
  g[1,4]<-paste0("HOM=",sum(tab[segregation$dominant.controls] %in% c("HOM"))," HET=",sum(tab[segregation$dominant.controls] %in% c("HET")))
  g[1,5]<-paste0("g.",tab[cchrom],":",tab[cpos],tab[cref],">",tab[calt])
  g[1,6]<-tab[caa] 
  g[1,7]<-tab[ccds]
  g[1,8]<-tab[cdbSNP]
  g[1,9]<-tab[cfreq]
  g[1,10]<-tab[varianteffect_cname]
  g[1,11]<-paste(tab[cdiseasevar],collapse=";")
  g[1,12]<-paste(tab[cphenotypes],collapse = ";")
  g[1,13]<-paste(tab[names(tab)=="candidate"],collapse = ";")
  
  #this_row<-c(g_prob,g_F,g_M,g_con,g_var,g_var1,g_var2,g_var3,g_var4,g_var5,g_var6,g_var7,g_var8)
  #names(this_row)<-c("Affected","Fathers","Mothers","Num_Controls_var","variant","amino_acid","cds","RS","Freq","Effect","Clinical_relevance","Phenotype","Panels")
  return(g)
}


do.rows<-function(tab,segregation) {
  res<-data.frame(Affected=character(),Fathers=character(),Mothers=character(),
                  Unaffected=numeric(),variant=character(),amino_acid=character(),
                  cds=character(),RS=character(),Freq=numeric(),
                  Effect=character(),Clinical_relevance=character(),Phenotype=character(),
                  Panels=character())
  for (r in 1:nrow(tab)) {
#    if (length(res)<13) print(res)
    try(res[r,]<-do.row2(tab,r,segregation))
  }
  return(res)
}

do.rows2<-function(tab,segregation) {
  res<-apply(tab,1,function(x) do.row3(x,segregation)) %>% bind_rows()
  report<-c(Affected=TRUE,
                  Fathers=TRUE,
                  Mothers=TRUE,
                  Unaffected=TRUE,
                  variant=TRUE,
                  amino_acid=TRUE,
                  cds=FALSE,
                  RS=TRUE,
                  Freq=TRUE,
                  Effect=TRUE,
                  Clinical_relevance=TRUE,
                  Phenotype=TRUE,
                  Panels=TRUE)
  return(res[report])
}

do.rows3<-function(tab,summarycolumns) {
	res<-tab[,summatycolumns]

}
	

# stampa un report dei geni con le informzioni del gene e le varianti selezionate per quel gene
print.geneinfo<-function(gene_list, gene_full_table, tab, 
                         cgene, atab, print_other=FALSE, cfreq, panels, 
                         check=FALSE, panelApp=FALSE,segregation,pubmed=FALSE, genetab=FALSE) {
  for (g in unique(gene_list)) {
    gtab<-tab[tab[,cgene]==g,]
    if (check) {
     if (sum(gtab$AD)<1 & sum(gtab$DN)<1 & sum(gtab$HM)<1 & sum(gtab$XL)<1 & sum(gtab$HTm | gtab$HTp)<2 ) next
    } else {
     if (sum(gtab$AD)<1 & sum(gtab$DN)<1 & sum(gtab$HM)<1 & sum(gtab$XL)<1 & sum(gtab$HTm | gtab$HTp)<2 ) gtab<-gtab[0,]
    }
    cat("\n------------------------------------\n")
    cat("-**********************************-\n")
    cat("------------------------------------\n")
    for (k in c("Symbol","description","type_of_gene","OMIM_Title","ClinVar_PhenotypeList","HGMD_disease","HPO_term-name","pLI","pRec","pNull"))
    {
      cat("- ")
      cat(k)
      cat(":\t\t")
      cat(unique(gene_full_table[gene_full_table$Symbol==g,k]))
      cat("\n")
    }
    cat("Panels:\n")
    for (p in 1:length(panels)) if (g %in% panels[[p]]) print(get.file.name(names(panels)[p]))
    cat("\n")
    if (genetab) try(print(kable(gene_full_table[ gene_full_table$Symbol==g,])))
    
    cat("\n")
    if (panelApp) query.panelApp(g)
    if (panelApp) try(print(kable(query.gene.omim(g))))
    cat("\n")
    cat("\n")
#    try(if (sum(atab[,cgene]==g)>20) aatab<-atab[atab[,cgene]==g & atab[,cfreq]<0.001 ,])
    aatab<-atab[atab[,cgene]==g,]
    cat ("Numero di varianti selezionate ="); cat(nrow(gtab)); cat("\n")
    cat ("Numero di tutte le varianti ="); cat(nrow(atab[atab[,cgene]==g,])); cat("\n")

    print(kable(gtab[,summarycolumns]))
    if (print_other) {
      cat(">  more variant in the gene\n")
      print(kable(aatab[,summarycolumns]))
    }
    cat("\n")
    if (pubmed) {
      cat(g,"\n")
      pubres<-do.pubmed.summary(g)
      cat(pubres[["outrow"]])
      for (pubmedterms in pubmed_terms) {
        cat("Terms: ",paste(pubmedterms, collapse=" AND "),"\n")
        pubres<-do.pubmed.summary(c(g,pubmedterms))
        cat(pubres[["outrow"]])
      }
    }
    
  }
}


# esegue la query in panelApp e formatta il risultato
query.panelApp<-function(genename,if.cat=TRUE) {
  if (nchar(genename)<1) return()
  querytext=paste0("https://bioinfo.extge.co.uk/WebServices/search_genes/",genename,"/")
  res <- fromJSON(querytext)
  if (res$meta$numOfResults<1) return()
  
  # questi sono i campi dei risultati di panelapp che riporta
  fields=c("SpecificDiseaseName","LevelOfConfidence",
           "ModeOfInheritance","Penetrance","Publications","Evidences")
  
  # qui scrive i risultati per esportarli dalla funzione
  qres=list()
  
  for (l in 1:nrow(res$results)) {

    
    # aggiunge un numero di elementi alla lista quanti sono i campi 
    qres[[l]]=rep(NA,length(fields))
    names(qres[[l]])<-fields
    
    for (ff in 1:length(fields)) {
      # questo per far in modo che abbia un indice numerico da 1 a N di campi
      #  e in f mette il valore di fields
      f=fields[ff]
      # se if.cat=FALSE non restituisce output a schermo ma solo come return
      if (if.cat) {
        cat(f)
        cat(":  ")
        cat(paste(unlist(res$results[l,f]),collapse="; "))
        cat("\n")
      }
      names(qres[[l]][ff])=f
      qres[[l]][f]=paste(unlist(res$results[l,f]),collapse="; ")
      
    }
    cat("\n")
  }
  return(qres)
}

do.pubmed.summary<-function(pubmedterms) {
  library(RISmed)
  
  querytext<-paste(pubmedterms,collapse = " ")

  res <- EUtilsSummary(querytext, type="esearch", db="pubmed", datetype='pdat', retmax=50)
  
  m<-QueryCount(res)
  if (m>20) {n=20} else {n=m}
  
  ti<-ArticleTitle(EUtilsGet(res))
  ab<-AbstractText(EUtilsGet(res))
  ye<-YearPubmed(EUtilsGet(res))
  pmid<-PMID(EUtilsGet(res))
  
  outrow<-""
  
  for (i in 1:n) outrwo<-paste0(outrow,paste0(ye,"\t",ti,"\t","https://www.ncbi.nlm.nih.gov/pubmed/",pmid,"\n"))
  
  return(list(outrow=outrwo,titl=ti,abst=ab,year=ye))

}

query.gene.omim<-function(genes) {
  library(biomaRt)
  # ora usa GRCh38
  
  mart = useMart(biomart="ENSEMBL_MART_ENSEMBL", host="www.ensembl.org", path="/biomart/martservice", dataset="hsapiens_gene_ensembl")
  attr<-c('hgnc_symbol','description','band','mim_morbid_description','mim_gene_accession')
  annotations<-getBM(attributes = attr, uniqueRows = T, filters = 'hgnc_symbol', values = genes, mart = mart, bmHeader = T)
  for (a in 1:length(names(annotations))) {while(grepl(" ",names(annotations)[a])) names(annotations)[a] <- sub(" ", "_", names(annotations)[a])}
  return(annotations)
}


# section 1 - Loading ---------- 

# default variables
if (!exists("ZYGsuff")) {cat("Default ZYG suffix: _ZYG");ZYGsuff='_ZYG'}
if (!exists("if_bestHit_gene_report")) {cat("No Best hits gene report");if_bestHit_gene_report=FALSE}
if (!exists("if_gene_report")) {cat("No Gene reports");if_gene_report=FALSE}
if (!exists("if_search_panel")) {cat("No searching for variant in gene panels");if_search_panel=FALSE}
if (!exists("if_write_00")) {cat("No searching for variant in gene panels");if_write_00=FALSE}
if (!exists("if_write_R")) {cat("No searching for variant in gene panels");if_write_R=FALSE}
if (!exists("penetranceTollerance_dominantControls")) {cat("No variant penetranceTollerance_dominantControls");penetranceTollerance_dominantControls=0}
if (!exists("penetranceTollerance_xlinkedControls")) {cat("No variant penetranceTollerance_xlinkedControls");penetranceTollerance_xlinkedControls=0}
if (!exists("penetranceTollerance_recessiveControls")) {cat("No variant penetranceTollerance_recessiveControls");penetranceTollerance_recessiveControls=0}

#FRANK: I need to send only environment variables. This will pick my previous arguments from the PERL script
#args = commandArgs(trailingOnly=TRUE)
#if (length(args)>0) {
#	cat("reading parameter file...\n")
#	paramfile=args[1]
#	cat("Reading parameter file:",paramfile,"\n")
#	source(paramfile)
	
#	cat("reading column file...\n")
#	colfiles=args[2]
#	cat("Reading column file:",colfiles,"\n")
#	source(colfiles)	
#	} else cat ("Using environment variables...\n")

cat("Loading data...\n")

# data files
setwd(wdir)
options(stringsAsFactors = FALSE)

#load(datafile)
# guessing filetype (tsv VS Rdata)
cat("Loading exome data file:\n")
cat("- ",datafile," in ",wdir,"\n")
if (grepl(".tsv$",datafile)) {filetype="table"} else if (grepl(".[Rr][Dd]at[a]?$",datafile)) {filetype="rdata"}
if (filetype=="table") {
  cat("reading tsv file: ",datafile,"\n")
  outxt<-read.table(datafile,sep="\t",heade=T,quote = "")
  } else if (filetype=="rdata") {
  cat("reading Rdata file: ",datafile,"\n")
  load(datafile)
 } else {
 cat("Unable to read exome data\n")
 break
}

# output file name
filename=paste0(wdir,"/",analysis_name,"_OUT")

cat("Loading panel data files...\n")
# panels

cat("candidate_gene_dir\n")
candidate_gene_file=dir(candidate_gene_dir,full.names=T)
panels<-read.panel(candidate_gene_file)

cat("annotation_gene_dir\n")
annotation_gene_file=dir(annotation_gene_dir,full.names=T)
gene_annotations<-read.panel(annotation_gene_file)

cat("candidate_variants_dir\n")
candidate_variant_file=dir(candidate_variants_dir,full.names=T)
variants<-read.variants(candidate_variant_file)

cat("annotation_variants_dir\n")
annotation_variant_file=dir(annotation_variants_dir,full.names=T)
variant_annotation<-read.variants(annotation_variant_file)

panel_annotation<-c(panels,gene_annotations)

# annotation file
# questo dataset si chiama "gene_full_table"
load(gene_table_file)

# Section 1.1 Loading PED data -----

cat("Deriving segregation...\n")
cat("Loading ped data file...\n")
# read pedigree data
ped<-read.table(pedfile,sep="\t",header=F,comment.char = "#")
# remove samples with status=NA
ped<-ped[!is.na(ped[,6]),]
# create dummy samples
ped<-rbind(ped,c("dummy","dummy1","","","1","0"))
ped<-rbind(ped,c("dummy","dummy2","","","2","0"))
# derive segregation
segregation<-get.cases_by_segregation(ped)
# print ped info
print(segregation)
cat("Number of affected males   = ",sum(ped$V5==1 &  ped$V6==2),"\n")
cat("Number of affected females = ",sum(ped$V5==2 &  ped$V6==2),"\n")
# adding dummy-control columns
cat("Adding dummy controls\n")
outxt[,paste0("dummy1",GTsuff)]<-rep("0/0",nrow(outxt))
outxt[,paste0("dummy2",GTsuff)]<-rep("0/0",nrow(outxt))
cat("Refining genotype call (lax/strict)...\n")
if (if_lax_segregation) {cat("Lax segregation applied\n")}
for (s in segregation[["samples"]]) {
  incolname=paste0(s,GTsuff)
  print(incolname)
  outcolname=paste0(s,ZYGsuff)
  print(outcolname)
  if (if_lax_segregation) {
    outxt[,outcolname]<-unlist(lapply(outxt[,incolname],do.zyg.lax)) 
  } else {
    outxt[,outcolname]<-unlist(lapply(outxt[,incolname],do.zyg))  
  }
}

# reassign the names to samples
cat("Reassign sample sames with ZYG suffix")
ped$V2<-paste0(ped$V2,ZYGsuff)
ped$V3<-paste0(ped$V3,ZYGsuff)
ped$V4<-paste0(ped$V4,ZYGsuff)
segregation<-get.cases_by_segregation(ped)
print(segregation)

for (cas in segregation$samples) {
  cat ("Sample ", cas,"has ", 
       sum(do.segregation2(outxt[,cas],c("HET","HOM"))), " variants, i.e. ",
       sum(do.segregation2(outxt[,cas],c("HET","HOM")))/nrow(outxt), "\n")
}

# Numero di varianti nel 
#  membri familiari delle diverse classi
for (n in names(segregation)) {
  for (g in c("HET","HOM")) {
    colName=paste0(n,"_",g)
    print(colName)
    outxt[,colName]<-call.cases(outxt,segregation,n,g)
  }
}



# Combina le colonne in base al genotipo dei tipi di familiari
printcolumn<-list()
printcolumn[["affected"]]<-c(HOM="dominant.cases_HOM",HET="dominant.cases_HET")
printcolumn[["fathers"]]<-c(HOM="fathers.of.affected_HOM",HET="fathers.of.affected_HET")
printcolumn[["mothers"]]<-c(HOM="mothers.of.affected_HOM",HET="mothers.of.affected_HET")
printcolumn[["controls"]]<-c(HOM="recessive.controls_HOM",HET="recessive.controls_HET")

for (k in names(printcolumn)) {
  if (length(printcolumn[[k]])==2) {
    outxt[,k]=apply(data.frame(outxt[,printcolumn[[k]]]),MARGIN = 1, 
          function(x) paste0(names(printcolumn[[k]][1]),"=",x[1],"; ",
                             names(printcolumn[[k]][2]),"=",x[2]))
  }
}

# Section 1.2 Search for colnames ----

cat("Assigning column names...\n")
# columns
cid=which(names(outxt)==id_cname)
cfreq=which(names(outxt)==maxfrequency_cname)
#cfunc=which(names(outxt)==varianteffect_cname)
cref=which(names(outxt)==ref_cname)
calt=which(names(outxt)==alt_cname)
cchrom=which(names(outxt)==chromosome_cname)
cpos=which(names(outxt)==position_cname)
cgene=which(names(outxt)==gene_cname)
cdbSNP<-which(names(outxt)==rs_cname)
ccds<-which(names(outxt)==cds)
caa<-which(names(outxt)==aa)
cphenotypes=which(names(outxt)==phenotype)

# Nomenclatura HGVS
outxt$HGVS<-apply(outxt,1,function(x) make.HGVS(x,cchrom,cpos,cref,calt))


# Section 2 - Segregation -----
#
# o---------------------------o
# | 2) SEGREGAZIONE           |
# o---------------------------o
# 



cat("Models of variant segregation..\n")

# Varianti con possibile effetto AD, per trasmissione
cat("Searching for AD variants..\n")
outxt$AD=rep(FALSE,nrow(outxt))
outxt$AD<-do.segregation2(outxt[,segregation$dominant.cases],c("HOM","HET")) &
  do.segregation2(outxt[,segregation$dominant.controls],c("HOMREF"),penetranceTollerance_dominantControls) &
  !do.segregation2(outxt[,segregation$parents.of.affected],c("HOMREF")) 
if (sum(outxt$AD)==0) cat("* Something wrong? No variant in the case(s)\n")

# Varianti con possibile effetto XL, per trasmissione
cat("Searching for XL variants..\n")
outxt$XL=rep(FALSE,nrow(outxt))
outxt$XL<-do.segregation2(outxt[,segregation$xlinked.cases],c("HOM","HET")) &
  do.segregation2(outxt[,segregation$xlinked.carriers],c("HET")) &
  do.segregation2(outxt[,segregation$xlinked.controls],c("HOMREF"),penetranceTollerance_xlinkedControls) &
  outxt[,cchrom]=="chrX"
if (sum(outxt$XL)==0) cat("* Something wrong? No variant in the case(s)\n")

# Varianti Omozigoti
cat("Searching for HM variants..\n")
outxt$HM=rep(FALSE,nrow(outxt))
outxt$HM<-do.segregation2(outxt[,segregation$recessive.cases],c("HOM")) &
  do.segregation2(outxt[,segregation$recessive.carriers],c("HET")) &
  do.segregation2(outxt[,segregation$recessive.controls],c("HET","HOMREF"),penetranceTollerance_recessiveControls) 
if (sum(outxt$HM)==0) cat("* Something wrong? No variant in the case(s)\n")

# Varianti Eterozigoti trasmesse da un solo genitore
outxt$HTm=rep(FALSE,nrow(outxt))
outxt$HTp=rep(FALSE,nrow(outxt))
cat("Searching for Maternal variants..\n")
outxt$HTm<-
   do.segregation2(outxt[, segregation$recessive.cases            ],c("HET")) &
  !do.segregation2(outxt[, segregation$recessive.carriers.paternal],c("HET")) &
   do.segregation2(outxt[, segregation$recessive.carriers.maternal],c("HET")) & 
  !do.segregation2(outxt[, segregation$recessive.controls         ],c("HOM"),penetranceTollerance_recessiveControls) 
if (sum(outxt$HTm)==0) cat("* Something wrong? No variant in the case(s)\n")
cat("Searching for Paternal variants..\n")
outxt$HTp<-
  do.segregation2(outxt[, segregation$recessive.cases            ],c("HET")) &
  do.segregation2(outxt[, segregation$recessive.carriers.paternal],c("HET")) &
  !do.segregation2(outxt[, segregation$recessive.carriers.maternal],c("HET")) & 
  !do.segregation2(outxt[, segregation$recessive.controls         ],c("HOM"),penetranceTollerance_recessiveControls) 
if (sum(outxt$HTp)==0) cat("* Something wrong? No variant in the case(s)\n")

# Varianti de novo
cat("Searching for DN variants..\n")
outxt$DN=rep(FALSE,nrow(outxt))
outxt$DN<-do.segregation2(outxt[,segregation$dominant.cases],c("HET","HOM")) &
  do.segregation2(outxt[,segregation$dominant.controls],c("HOMREF"),penetranceTollerance_dominantControls) & 
  do.segregation2(outxt[,segregation$parents.of.affected],c("HOMREF"))
if (sum(outxt$DN)==0) cat("* Something wrong? No variant in the case(s)\n")

# Variants in all cases
outxt$no_all_cases_HOMREF<-!do.segregation2(outxt[,segregation$dominant.cases],c("HOMREF"))

# Varianti note per malattia. -known-
#  Mette true a tutte le varianti delle colonne in cdiseasevar che hanno grep 
#  per le stringhe in diseasevar_greps
cat("Looking for known disease-causing variants...\n")
outxt$known=rep(FALSE,nrow(outxt))
for (c in cdiseasevar) {
	for (g in diseasevar_grep[[c]]) {
		cat(" searching for ",g," in ",c,"\n")
		outxt$known[grepl(g, outxt[,c])]<-TRUE 
	}
}
if (sum(outxt$known)==0) cat("* Something wrong? No variant in the case(s)\n")

# Varianti rare
outxt$rare<-outxt[,cfreq]<Fmax

# Variants annotation
outxt$candidate=rep(NA,nrow(outxt))
outxt.gr<-GRanges(seqnames=outxt[,cchrom],ranges=IRanges(outxt[,cpos],end=outxt[,cpos]+nchar(outxt[,cref])-1))
for (v in 1:length(variant_annotation)) {
  ov <- findOverlaps(outxt.gr,variant_annotation[[v]])@from
  outxt$candidate[ov]<-paste0(outxt$candidate[ov]," | ",get.file.name(names(variant_annotation)[v]))
}
outxt$candidate_var<-!is.na(outxt$candidate)

# Gene annotation to-do
# panel_annotation
if (length(panel_annotation)>0) {
  for (p in 1:length(panel_annotation)) {
    
    panel=panel_annotation[[p]]
    panel_name=names(panel_annotation)[p]
    selected_var<-outxt[,cgene] %in% panel & outxt[,cfreq]<Fmax & (outxt$HTm | outxt$HTp | outxt$DN | outxt$HM | outxt$XL | outxt$AD)
    # make a column for each gene-panel and report if the variant fall in the gene
    outxt[,get.file.name(panel_name)]<-selected_var
    
  }
}

# Functionally relevant
outxt$functionally_relevant<-rep(FALSE,nrow(outxt))
outxt$functionally_relevant<-outxt[,varianteffect_cname] %in% effects
outxt$functionally_relevant[!is.na(outxt$candidate)]<-TRUE

outxt[is.na(outxt[,cgene]),cgene]<-""

# Models
# gmodels=c("known_xl","known_ad","known_ar","new_ar","new_dn","new_xl","new_ad","all_ad","all_ar","all_xl","all_dn")
# outxt$gmodel<-apply(outxt[,which(names(outxt) %in% gmodels)], 1, function(x) paste0(names(which(x)),collapse=","))
# 
cmodels=c("AD","XL","HM","HTm","HTp","DN","known","rare","functionally_relevant","candidate_var")
outxt$vmodel<-apply(outxt[,which(names(outxt) %in% cmodels)], 1, function(x) paste0(names(which(x)),collapse=","))
outxt$vmodel[outxt$vmodel==""]<-"."

# Section 3 - Filtering ----

cat("Variant filering ...\n")

tab_list=list()
gene_list=""

# ** known XL 
cat(" . filering known XL...\n")
tab<-outxt[outxt$known & outxt$XL & outxt[,cfreq]<0.05,]
known_xl<-do_report(tab)
tab_list[["known_xl"]]<-known_xl[["tab"]]
# aggiunge una colonna nel outxt con le varianti selezionate tenendo conto della segregazione nella famiglia
outxt$known_xl<- outxt[,cid] %in% known_xl[["tab"]][,cid]

# ** known AD/DN
cat(" . filering known AD/DN...\n")
tab<-outxt[outxt$known & outxt$AD & outxt[,cfreq]<0.05,]
known_ad<-do_report(tab)
tab_list[["known_ad"]]<-known_ad[["tab"]]
# aggiunge una colonna nel outxt con le varianti selezionate tenendo conto della segregazione nella famiglia
outxt$known_ad<- outxt[,cid] %in% known_ad[["tab"]][,cid]

# ** known AR
cat(" . filering known AR...\n")
comphet_gene<-find_comp_het_genes(outxt[ (outxt$HTm | outxt$HTp) & outxt$known & outxt[,cfreq]<0.05,])[["genes"]]
homozygous_gene<-outxt[ outxt$HM & outxt$known & outxt[,cfreq]<0.05,cgene]
tab<-rbind(outxt[outxt[,cgene] %in% comphet_gene & outxt$known & outxt[,cfreq]<0.05 & (outxt$HTm | outxt$HTp),],
           outxt[outxt[,cgene] %in% homozygous_gene & outxt$HM & outxt$known & outxt[,cfreq]<0.05,])
known_ar<-do_report(tab)
tab_list[["known_ar"]]<-known_ar[["tab"]]
# aggiunge una colonna nel outxt con le varianti selezionate tenendo conto della segregazione nella famiglia
outxt$known_ar<- outxt[,cid] %in% known_ar[["tab"]][,cid]

# ** Novel AR
cat(" . filering novel AR...\n") 
comphet_gene<-find_comp_het_genes(outxt[ (outxt$HTm | outxt$HTp)  & outxt$functionally_relevant & outxt[,cfreq]<Fmax,])[["genes"]]
homozygous_gene<-outxt[outxt$functionally_relevant & outxt$HM &  outxt[,cfreq]<Fmax,cgene]
tab<-rbind(outxt[outxt[,cgene] %in% comphet_gene & (outxt$HTm | outxt$HTp) & outxt$functionally_relevant & outxt[,cfreq]<Fmax,],
           outxt[outxt[,cgene] %in% homozygous_gene & outxt$functionally_relevant & outxt$HM & outxt[,cfreq]<Fmax,])
new_ar<-do_report(tab)
tab_list[["new_ar"]]<-new_ar[["tab"]]
# aggiunge una colonna nel outxt con le varianti selezionate tenendo conto della segregazione nella famiglia
outxt$new_ar<- outxt[,cid] %in% new_ar[["tab"]][,cid]

# ** Novel De Novo
cat(" . filering novel DN...\n") 
tab<-outxt[outxt$DN & outxt$rare & outxt$functionally_relevant,]
new_dn<-do_report(tab)
tab_list[["new_dn"]]<-new_dn[["tab"]]
# aggiunge una colonna nel outxt con le varianti selezionate tenendo conto della segregazione nella famiglia
outxt$new_dn<- outxt[,cid] %in% new_dn[["tab"]][,cid]

# ** Novel XL
cat(" . filering novel XL...\n") 
tab<-outxt[outxt$XL & outxt$rare & outxt$functionally_relevant,]
new_xl<-do_report(tab)
tab_list[["new_xl"]]<-new_xl[["tab"]]
# aggiunge una colonna nel outxt con le varianti selezionate tenendo conto della segregazione nella famiglia
outxt$new_xl<- outxt[,cid] %in% new_xl[["tab"]][,cid]


# ** Novel AD
cat(" . filering novel AD...\n") 
tab<-outxt[outxt$AD & outxt$rare & outxt$functionally_relevant,]
new_ad<-do_report(tab)
tab_list[["new_ad"]]<-new_ad[["tab"]]
# aggiunge una colonna nel outxt con le varianti selezionate tenendo conto della segregazione nella famiglia
outxt$new_ad<- outxt[,cid] %in% new_ad[["tab"]][,cid]

# ** All AD
cat(" . filering all AD...\n") 
tab<-outxt[outxt$AD & outxt$rare,]
all_ad<-do_report(tab)
tab_list[["all_ad"]]<-all_ad[["tab"]]
# aggiunge una colonna nel outxt con le varianti selezionate tenendo conto della segregazione nella famiglia
outxt$all_ad<- outxt[,cid] %in% all_ad[["tab"]][,cid]

# ** All AR
cat(" . filering all AR...\n") 
comphet_gene<-find_comp_het_genes(outxt[ (outxt$HTm | outxt$HTp) & outxt$rare,])[["genes"]]
homozygous_gene<-outxt[outxt$functionally_relevant & outxt$HM & outxt$rare,cgene]
tab<-rbind(outxt[outxt[,cgene] %in% comphet_gene & (outxt$HTm | outxt$HTp) & outxt$rare,],
           outxt[outxt[,cgene] %in% homozygous_gene & outxt$HM & outxt$rare,])
all_ar<-do_report(tab)
tab_list[["all_ar"]]<-all_ar[["tab"]]
# aggiunge una colonna nel outxt con le varianti selezionate tenendo conto della segregazione nella famiglia
outxt$all_ar<- outxt[,cid] %in% all_ar[["tab"]][,cid]

# ** All XL
cat(" . filering all XL...\n") 
tab<-outxt[outxt$XL & outxt$rare,]
all_xl<-do_report(tab)
tab_list[["all_xl"]]<-all_xl[["tab"]]
# aggiunge una colonna nel outxt con le varianti selezionate tenendo conto della segregazione nella famiglia
outxt$all_xl<- outxt[,cid] %in% all_xl[["tab"]][,cid]

# ** All De Novo
cat(" . filering all DN...\n") 
tab<-outxt[outxt$DN & outxt$rare,]
all_dn<-do_report(tab)
tab_list[["all_dn"]]<-all_dn[["tab"]]
# aggiunge una colonna nel outxt con le varianti selezionate tenendo conto della segregazione nella famiglia
outxt$all_dn<- outxt[,cid] %in% all_dn[["tab"]][,cid]

# ** Gene Panel
cat(" . filering panels...\n") 
if (length(panels)>0 & if_search_panel) {
  for (p in 1:length(panels)) {
    
    panel=panels[[p]]
    panel_name=names(panels)[p]
    cat(" .. filtering panel ",panel_name,"\n")
    selected_var<-outxt[,cgene] %in% panel & outxt[,cfreq]<Fmax & (outxt$HTm | outxt$HTp | outxt$DN | outxt$HM | outxt$XL | outxt$AD)
    # make a column for each gene-panel and report if the variant fall in the gene
    outxt[,get.file.name(panel_name)]<-selected_var
    
    tab<-outxt[selected_var,]
    all_panel<-do_report(tab)
    genes<-all_panel[["tab"]][,cgene]
    tab_list[[panel_name]]<-all_panel[["tab"]]
  }
}

# Section 4 - write files -----

generaloutfilename<-paste0(filename,"_GeneralOut.txt")
cat("#Segregation File\n",file=generaloutfilename)
for (l in names(segregation)) cat(l,paste(segregation[[l]],collapse=", "),"\n",file = generaloutfilename,append = T)

# scrive lista di geni
gene_list<-c()
for (t in 1:length(tab_list)) {
  gene_list<-c(gene_list,unlist(tab_list[[t]][cgene]))
}
gene_list<-unique(gene_list)

# scrive report PER GENE
#   (aggiugendo i modelli di trasmissione come pannelli)

# elenco di modelli e pannelli per l'annotazione dei geni
models<-list()
for (t in 1:7) {
  models[[names(tab_list)[t]]]<-unlist(tab_list[[t]][cgene])
}
panel_annotation<-c(models,panel_annotation)

## 4.1 Best hits, gene report -----

cat("Writing BestHits GeneReport\n")

best_models<-list()
for (t in 1:7) {
  best_models[[names(tab_list)[t]]]<-unlist(tab_list[[t]][cgene])
}
glist<-unique(unlist(best_models))
cat("Number of Best Hit genes; ")
cat(length(glist))
cat("\n")
tab<-outxt[ outxt$rare & (outxt$HM | outxt$HTm | outxt$HTp | outxt$DN | outxt$AD) & outxt$functionally_relevant,]

print(summarycolumns)
print(summarycolumns %in% names(outxt))
sinkfile=paste0(filename,"_BestHits_GeneReport.txt")
sink(sinkfile,split=F)
if (if_bestHit_gene_report) print.geneinfo(glist, gene_full_table, tab, cgene, outxt, TRUE, cfreq, panel_annotation, TRUE, if_panel_app,segregation,if_pubmed_besthits,TRUE)
sink()

# selezione le variabili di best-hits e le salve in un bed file
sinkfile=paste0(filename,"_BestHits.bed")
bestVar<-data.frame(chrom=character(),pos=numeric(),vmodel=character())
for (t in 1:7) {
  bestVar<-rbind(bestVar,cbind(tab_list[[t]][,c(cchrom,cpos)],tab_list[[t]]$vmodel))
}
write.table(data.frame(bestVar[,1],bestVar[,2]-1,bestVar[,2],bestVar[,3]),sinkfile,sep="\t",col.names = F,row.names = F,quote = FALSE)

sinkfile=paste0(filename,"_BestHits_genes.txt")
writeLines(glist,sinkfile)
sinkfile=paste0(filename,"_BestHits.tsv")
write.table(data.frame(tab),sinkfile,sep="\t",col.names = T,row.names = F,quote=FALSE)

## 4.2 Other genes, gene report -----

cat("Writing GeneReport\n")
#sinkfile=paste0(filename,"_GeneReport.txt")
#sink(sinkfile,split=F)
# 
# models<-list()
# for (t in 1:7) {
#   models[[names(tab_list)[t]]]<-unlist(tab_list[[t]][cgene])
# }
# glist<-unique(unlist(panels))
# glist<-glist[!glist %in% gene_list]
# 
# panel_annotation<-c(models,panel_annotation)
# 
# tab<-outxt[ outxt[,cfreq]<Fmax  & (outxt$AD | outxt$HTm | outxt$HTp | outxt$DN | outxt$HM) & outxt$functionally_relevant,]
# 
# if (if_gene_report) print.geneinfo(gene_list, gene_full_table, tab, cgene, outxt, TRUE, cfreq, panel_annotation, TRUE,FALSE,segregation,if_pubmed_others,FALSE)
# cat("\n\n-----------------------------------------------------------------------------------------\n\n")
# if (if_gene_report) print.geneinfo(glist, gene_full_table, tab, cgene, outxt, TRUE, cfreq, panel_annotation, TRUE, FALSE,segregation,if_pubmed_others,FALSE)
# sink()

cat("Writing gene list\n")
 writeLines(unique(gene_list),paste0(filename,"_gene_list.txt"))
 
if (if_write_vartable) {
  cat("Writing Variant tables\n")
  best_hits_variants=""
  for (t in 1:length(tab_list)) {
    this_filename=paste0(filename,"_",formatC(t, width=2, flag="0"),"_",get.file.name(names(tab_list)[t]),".tsv")
    write.table(file=paste0(this_filename),x=tab_list[[t]],sep="\t",col.names=T,row.names=F,quote=F)
    cat("Writing gene report for ")
    cat(get.file.name(names(tab_list)[t]))
    cat("\n")
    if (if_gene_report) {
      sinkfile=paste0(filename,"_",get.file.name(names(tab_list)[t]),"_GeneReport.txt")
      sink(sinkfile,split=F)
      glist<-unique(tab_list[[t]][,cgene])
      tab<-tab_list[[t]]
      print.geneinfo(gene_list = glist, 
                     gene_full_table = gene_full_table, 
                     tab = tab, 
                     cgene = cgene, 
                     atab = outxt, 
                     print_other = TRUE, 
                     cfreq = cfreq, 
                     panels = panel_annotation, 
                     check = TRUE, 
                     panelApp = if_panel_app, 
                     segregation = segregation,
                     pubmed = if_pubmed_others,
                     genetab = FALSE)
      cat("\n\n-----------------------------------------------------------------------------------------\n\n")
        print.geneinfo(gene_list = glist,
                       gene_full_table =  gene_full_table,
                       tab =  tab,
                       cgene =  cgene,
                       atab =  outxt,
                       print_other =  TRUE,
                       cfreq = cfreq,
                       panels =  panel_annotation,
                       check = TRUE,
                       panelApp =  if_panel_app,
                       segregation = segregation,
                       pubmed = if_pubmed_others,
                       genetab = FALSE)
      sink()
    }
  }
}
 
cat("Writing variant_lists and gene_lists\n")
for (t in 1:length(tab_list)) {
  print (t)
  tname=get.file.name(names(tab_list)[t])
  cat(". writing ",tname,"\n")
  tab=tab_list[[t]]
  sinkfile=paste0(filename,"_",tname,".bed")
  write.table(data.frame(tab[,cchrom],tab[,cpos]-1,tab[,cpos],tab$vmodel),sinkfile,sep="\t",col.names = F,row.names = F,quote = FALSE)
  sinkfile=paste0(filename,"_",tname,"_genes.txt")
  writeLines(unique(tab[,cgene]),sinkfile)
  tname=get.file.name(names(tab_list)[t])
  if (if_gene_report) {
    
  }
}


# Variants in all cases
cat ("Writing variants in all cases .bed\n")
this_filename=paste0(filename,"_variant_list.bed")
tab<-outxt[outxt$no_all_cases_HOMREF,]
tab2<-data.frame(tab[,cchrom],tab[,cpos]-1,tab[,cpos],tab$vmodel)
write.table(file=this_filename,tab2,sep="\t",col.names=F,row.names=F,quote=F)

cat("Writing 00_preprocessed\n")
this_filename=paste0(filename,"_00_preprocessed.tsv")
if (if_write_00) write.table(file=paste0(this_filename),x=outxt,sep="\t",col.names=T,row.names=F,quote=F)

cat("Writing R data\n")
if (if_write_R) save.image(file=paste0(filename,".Rdata"))

alarm()
