cai_computation<-function(transcr,maf_file,pos_tr,chrm_maf,pos_chrm,pos_mut_maf,pos_type_mut,seq_al,cod_tab,limit_min=NA,limit_max=NA)
{
  library(stringr)
  library(GenomicRanges)
  library(phastCons100way.UCSC.hg19)
  gsco <- phastCons100way.UCSC.hg19
  tr_maf<-maf_file[which(maf_file[,pos_tr]==transcr),]
  if(!is.na(limit_min))
    tr_maf<-tr_maf[which(as.numeric(tr_maf[,pos_chrm])>=as.numeric(limit_min)),]
  if(!is.na(limit_max))
    tr_maf<-tr_maf[which(as.numeric(tr_maf[,pos_chrm])<=as.numeric(limit_max)),]
  #if(!is.na(type_mut))
  #tr_maf<-tr_maf[which(tr_maf[,pos_type_mut]==type_mut),]
  #Now I have the maf file
  seq_tr<-seqinr::getSequence(seq_al, as.string = FALSE)
  length(seq_tr)
  orig_cai<-as.numeric(cai(seq_tr,w=as.numeric(cod_tab)))
  wt_rscu<-uco(seq_tr,index="rscu")
  #Now I have to introduce the mutations in the sequence.
  #I will add the most frequent mutation
  #I create the a vector with the position of each codon.
  pos_codon<-seq(from=1,to=length(seq_tr),by=3)
  freq_mut<-table(tr_maf[,pos_mut_maf])
  pos_rm_mut<-grep("\\NA|\\?",names(freq_mut))
  pos_rm_mut<-c(pos_rm_mut,grep("-",names(freq_mut)))
  pos_null<-which(names(freq_mut)=="")
  pos_rm_mut<-c(pos_rm_mut,pos_null)
  if(length(pos_rm_mut)>0)
    freq_mut<-freq_mut[-pos_rm_mut]
  most_freq<-freq_mut[which(freq_mut==max(freq_mut))[1]]
  chrm_pos<-tr_maf[which(tr_maf[,pos_mut_maf]==names(most_freq))[1],pos_chrm]
  type_mut<-tr_maf[which(tr_maf[,pos_mut_maf]==names(most_freq))[1],pos_type_mut]
  chrm<-tr_maf[which(tr_maf[,pos_mut_maf]==names(most_freq))[1],chrm_maf]
  chrm<-paste("chr",chrm,sep="")
  #res<-score(gsco, GRanges(seqnames=chrm, IRanges(start=chrm_pos, width=1)))
  score_mut<-score(gsco, GRanges(seqnames=chrm, IRanges(start=as.numeric(chrm_pos), width=1)))
  mut_changed<-names(most_freq)
  pos_mut<-gsub("[A-Z]","",mut_changed)
  pos_mut<-as.numeric(gsub(">","",pos_mut))
  if(pos_mut>length(seq_tr))
    return(res<-NA) else{
      mut_from<-strsplit(gsub("[0-9]","",mut_changed),">")[[1]][1]
      mut_to<-strsplit(gsub("[0-9]","",mut_changed),">")[[1]][2]
      new_seq<-seq_tr
      new_seq[pos_mut]<-mut_to
      mut_rscu<-uco(new_seq,index="rscu")
      #Now I want also to know the frequency of the new codons compared to the old one
      #Identification of the codon.
      diff_pos<-abs(pos_codon-pos_mut)
      dist_min3<-which((diff_pos<3)==T)
      dist_cod<-which((pos_codon<=pos_mut)==T)
      min_dist<-intersect(dist_min3,dist_cod)
      orig_cod<-paste(seq_tr[pos_codon[min_dist]:(pos_codon[min_dist]+2)],collapse="")
      freq_orig_cod<-cod_tab[which(names(cod_tab)==orig_cod)]
      mut_cod<-paste(new_seq[pos_codon[min_dist]:(pos_codon[min_dist]+2)],collapse="")
      freq_mut_cod<-cod_tab[which(names(cod_tab)==mut_cod)]
      cai_mut<-as.numeric(cai(new_seq,w=cod_tab))
      ratio_mut_norm<-cai_mut/orig_cai
      rscu_orig<-wt_rscu[which(names(wt_rscu)==str_to_lower(orig_cod))]
      rscu_mut<-mut_rscu[which(names(mut_rscu)==str_to_lower(mut_cod))]
      rscu_ratio<-rscu_mut/rscu_orig
      if(rscu_ratio<1)
        rscu_ratio=1/rscu_ratio
      #gene_name<-ens2gene(transcr,ens_type="ensembl_transcript_id")
      res<-cbind(transcr,orig_cai,cai_mut,ratio_mut_norm,mut_changed,most_freq,limit_min,limit_max,type_mut,orig_cod,mut_cod,
                 freq_orig_cod,freq_mut_cod,rscu_orig,rscu_mut,rscu_ratio,score_mut,chrm_pos)
      #nona_pos<-which(!is.na(res[1,]))
      #res<-res[,nona_pos]
      #name_to_assign<-c("seq","original_cai","mutated_cai","ratio_Mut_N","mutations","freq_mut","pos_min","pos_max","mutation type","wt codon",
      #             "mutated codon","freq wt codon","freq mut codon","rscu_orig","rscu_mut","conserved_score","chrm_pos")
      #names(res)<-name_to_assign[nona_pos]
      #names(res)<-name_to_assign
      return(res)}
  #Now I introduce the mutation in the sequence
}