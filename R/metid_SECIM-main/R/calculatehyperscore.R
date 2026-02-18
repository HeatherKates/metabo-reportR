CalculateHyperScore <- function(mSetObj=NA){
  mSetObj <- mSet
  
  # make a clean dataSet$cmpd data based on name mapping
  # only valid hmdb name will be used
  
  # SSP input for enrich is different
  if(mSetObj$analSet$type=="msetssp"){
    ora.vec <- mSetObj$dataSet$cmpd;
  }else{
    nm.map <- GetFinalNameMap(mSetObj);
    valid.inx <- !(is.na(nm.map$hmdb)| duplicated(nm.map$hmdb));
    ora.vec <- nm.map$hmdb[valid.inx];
  }
  
  q.size<-length(ora.vec);
  
  if(all(is.na(ora.vec)) || q.size==0) {
    AddErrMsg("No valid HMDB compound names found!");
    return(0);
  }
  
  # move to api only if R package + KEGG msets
  if(!.on.public.web & grepl("kegg", mSetObj$analSet$msetlibname)){
    
    # make this lazy load
    #if(!exists("my.hyperscore.kegg")){ # public web on same user dir
      source("/blue/timgarrett/hkates/SECIM_Reporting/R/util_api.R")
    #}
    
    mSetObj$api$oraVec <- ora.vec; 
    
    if(mSetObj$api$filter){
      mSetObj$api$filterData <- mSetObj$dataSet$metabo.filter.kegg
      toSend <- list(mSet = mSetObj, libNm = mSetObj$api$libname, filter = mSetObj$api$filter,
                     oraVec = mSetObj$api$oraVec, filterData = mSetObj$api$filterData,
                     excludeNum = mSetObj$api$excludeNum)
    }else{
      toSend <- list(mSet = mSetObj,libNm = mSetObj$api$libname, 
                     filter = mSetObj$api$filter, oraVec = mSetObj$api$oraVec, excludeNum = mSetObj$api$excludeNum)
    }
    saveRDS(toSend, "tosend.rds")
    return(my.hyperscore.kegg());
  }
  
  current.mset <- current.msetlib$member;
  
  # make a clean metabilite set based on reference metabolome filtering
  # also need to update ora.vec to the updated mset
  if(mSetObj$dataSet$use.metabo.filter && !is.null(mSetObj$dataSet$metabo.filter.hmdb)){
    current.mset <- lapply(current.mset, function(x){x[x %in% mSetObj$dataSet$metabo.filter.hmdb]})
    mSetObj$dataSet$filtered.mset <- current.mset;
  }
  
  set.size<-length(current.mset);
  if(set.size ==1){
    AddErrMsg("Cannot perform enrichment analysis on a single metabolite set!");
    return(0);
  }
  
  # now perform enrichment analysis
  
  # update data & parameters for ORA stats, based on suggestion
  # https://github.com/xia-lab/MetaboAnalystR/issues/168
  # https://github.com/xia-lab/MetaboAnalystR/issues/96
  # https://github.com/xia-lab/MetaboAnalystR/issues/34
  
  # the universe based on reference metabolome (should not matter, because it is already filtered based on the reference previously)
  
  my.univ <- unique(unlist(current.mset, use.names=FALSE));
  if(!is.null(mSetObj$dataSet$metabo.filter.hmdb)){
    my.univ <- unique(mSetObj$dataSet$metabo.filter.hmdb);
  }
  uniq.count <- length(my.univ);
  ora.vec <- ora.vec[ora.vec %in% my.univ];
  q.size <- length(ora.vec); 
  
  hits<-lapply(current.mset, function(x){x[x %in% ora.vec]});  
  hit.num<-unlist(lapply(hits, function(x) length(x)), use.names = FALSE);
  
  if(sum(hit.num>0)==0){
    AddErrMsg("No match was found to the selected metabolite set library!");
    return(0);
  }
  
  set.num<-unlist(lapply(current.mset, length), use.names = FALSE);
  
  # prepare for the result table
  res.mat<-matrix(NA, nrow=set.size, ncol=6);        
  rownames(res.mat)<-names(current.mset);
  colnames(res.mat)<-c("total", "expected", "hits", "Raw p", "Holm p", "FDR");
  for(i in 1:set.size){
    res.mat[i,1]<-set.num[i];
    res.mat[i,2]<-q.size*(set.num[i]/uniq.count);
    res.mat[i,3]<-hit.num[i];
    
    # use lower.tail = F for P(X>x)
    # phyper("# of white balls drawn", "# of white balls in the urn", "# of black balls in the urn", "# of balls drawn")
    res.mat[i,4]<-phyper(hit.num[i]-1, set.num[i], uniq.count-set.num[i], q.size, lower.tail=F);
  }
  
  # adjust for multiple testing problems
  res.mat[,5] <- p.adjust(res.mat[,4], "holm");
  res.mat[,6] <- p.adjust(res.mat[,4], "fdr");
  
  res.mat <- res.mat[hit.num>0,];
  
  ord.inx<-order(res.mat[,4]);
  mSetObj$analSet$ora.mat <- signif(res.mat[ord.inx,],3);
  mSetObj$analSet$ora.hits <- hits;
  
  fast.write.csv(mSetObj$analSet$ora.mat, file="msea_ora_result.csv");
  return(.set.mSet(mSetObj));
}

