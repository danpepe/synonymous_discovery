window_mut2_cai<-function(maf_info,pos_chrm,pos_chr,pos_gene,wind_size,pos_barcode,pos_variant,name_silent,col_tr_name,thr_samples,tot_samples,thr_mut,tr_leng,nperm,cancer_gene,ens_type,cod_tb,aa_pos)
{
  #I will run before an analysis for the removal of transcripts 
  #represented by less than thr_samples.
  #Genes in the maf file 
  genes_name<-unique(maf_info[,pos_gene])
  #which(genes_name=="ENST00000407796")
  gene_list_mut<-list()
  n_tr<-NULL
  s<-1
  for(i in 1:length(genes_name))
  {
    #print(i)
    #I extract the information about that gene
    gene_info<-maf_info[which(maf_info[,pos_gene]==genes_name[i]),]
    wind_to_save<-list()
    n_mut<-nrow(gene_info)
    c<-1
    if(thr_samples<=length(unique(gene_info[,pos_barcode])))
    {
      n_tr<-c(n_tr,i)
      #Unique position in the chromosome
      chr_pos<-unique(gene_info[,pos_chrm])
      j<-1
      for(j in 1:length(chr_pos))
      {
        #print(j)
        #I look for the rows where the first chromosome position appears. 
        pos_to_cons<-which(gene_info[,pos_chrm]==chr_pos[j])[1]
        wind<-subset(gene_info, gene_info[,pos_chrm]>=gene_info[pos_to_cons,pos_chrm] & gene_info[,pos_chrm]<=(as.numeric(gene_info[pos_to_cons,pos_chrm])+wind_size))
        #Computation of the window
        wind_lim<-max(as.numeric(wind[,pos_chrm]))-min(as.numeric(wind[,pos_chrm]))
        #I want to save the frequency of the sample and the mutations
        freq_m<-nrow(wind)/n_mut*100
        freq_m_abs<-nrow(wind)
        freq_s_tr_abs<-length(unique(gene_info[,pos_barcode]))
        freq_s_tr_perc<-length(unique(gene_info[,pos_barcode]))/tot_samples*100
        freq_s_abs<-length(unique(wind[,pos_barcode]))
        freq_s<-freq_s_abs/tot_samples*100
        perc_syn<-table(wind[,pos_variant])
        pos_sil<-which(names(perc_syn)==name_silent)
        if(length(pos_sil)==0)
          sil_freq<-0
        else
          sil_freq<-perc_syn[pos_sil]/sum(perc_syn)
        if(freq_m>thr_mut)
        {
          #I don't have to check if it is the first window to save. 
          if(c==1)
          {
            wind<-wind[,c(pos_chr,pos_chrm,pos_gene,pos_variant)]
            wind<-wind[!duplicated(wind),]
            #wind_hot<-rep(paste("Hotspot",c,sep=""),nrow(wind))
            wind<-cbind(wind,freq_s,freq_m,freq_m_abs,freq_s_abs,freq_s_tr_abs,freq_s_tr_perc,n_mut,wind_lim,sil_freq)
            wind_to_save[[c]]<-wind
            c<-c+1
          }
          else
          {
            trial<-merge(wind_to_save[[c-1]],wind)
            if(nrow(trial)!=nrow(wind))
            {
              wind<-wind[,c(pos_chr,pos_chrm,pos_gene,pos_variant)]
              wind<-wind[!duplicated(wind),]
              #wind_hot<-rep(paste("Hotspot",c,sep=""),nrow(wind))
              wind<-cbind(wind,freq_s,freq_m,freq_m_abs,freq_s_abs,freq_s_tr_abs,freq_s_tr_perc,n_mut,wind_lim,sil_freq)
              wind_to_save[[c]]<-wind
              c<-c+1
            }
          }
          #I should see if the one I am considering is a complete subset of the previous one.
        }
      }
      if(length(wind_to_save)!=0)
      {
        gene_list_mut[[s]]<-wind_to_save
        #names(gene_list_mut)[s]<-genes_name[i]
        s<-s+1
      }
    }
  }
  #Ok. Now I can work on this dataset.
  data_fr_res<-data.frame()
  if(length(gene_list_mut)==0)
    return(data_fr_res)
  else
  {
    list1 <- unlist(gene_list_mut, recursive = FALSE)
    data_fr_res<-do.call(rbind,list1)
    #Here I have the dataframe on which I can apply the permutation test
    #to see if the hotspots are random or not.
    #First step, detection of the column with the freq_mut
    freq_m_col<-which(colnames(data_fr_res)=="freq_m")
    tr_col<-which(colnames(data_fr_res)==col_tr_name)
    n_mut_col<-which(colnames(data_fr_res)=="n_mut")
    wind_lim<-which(colnames(data_fr_res)=="wind_lim")
    sil_freq_n<-which(colnames(data_fr_res)=="sil_freq")
    pos_sm_tr_abs<-which(colnames(data_fr_res)=="freq_s_tr_abs")
    pos_sm_tr_perc<-which(colnames(data_fr_res)=="freq_s_tr_perc")
    datafr_1<- data_fr_res[,c(freq_m_col,tr_col,n_mut_col,wind_lim,sil_freq_n,pos_sm_tr_abs,pos_sm_tr_perc)]
    data_fr2<-unique(datafr_1)
    #unique(data_fr2[,2])=="ENST00000349310"
    #Now. I want to get a unique value for each transcript.
    #I can perform the permutation test for each hotspot.
    #i<-11
    p_value_hot<-vector()
    #data_fr2[1:10,]
    for(i in 1:nrow(data_fr2))
    {
      print(paste("Permutation hotspot",i," out of",nrow(data_fr2)))
      tr_l<-tr_leng[which(tr_leng[,1]==data_fr2[i,2]),3]
      p_value_hot[i]<-window_randistr(data_fr2[i,3],tr_l,data_fr2[i,1],data_fr2[i,4],nperm)
    }
    adj_pvalue<-p.adjust(p_value_hot, method = "BH", n = length(p_value_hot))
    datafr_2<-cbind(data_fr2,adj_pvalue)
    hot_procdata<-hotspot_processing(list(data_fr_res,datafr_2),cancer_gene,3,6,ens_type)
    hot_agg<-hotspot_unit_results(hot_procdata)
    colnames(hot_agg)
    hot_agg[1:10,]
    colnames(hot_agg)<- c("Ensembl_tr","chrm","min_chrm","max_chrm","%samples_hotspot","%Mutations_hotspot", "#mutations_hotspot",
                          "#samples_hotspot", "#samples_transcript","%samples_transcript","#Mutations_transcript","size_window","%Reference mutation","Cancer_gene","Symbol","pvalue")
    # hot_agg<-mc3SMs_hot[[3]]
    cai_df<-data.frame()
    i<-33
    for(i in 1:nrow(hot_agg))
    {
      print(i)
      seq_x_cai<-tr_leng[which(tr_leng[,1]==hot_agg[i,1]),2]
      res_cai<-cai_computation(hot_agg[i,1],maf_info,pos_gene,pos_chr,pos_chrm,aa_pos,pos_variant,seq_x_cai,cod_tb,
                               hot_agg[i,3],hot_agg[i,4])
      cai_df<-rbind(cai_df,res_cai,stringsAsFactors=F)
      #colnames(cai_conc_nSMs)<-names(res_cai)
    } 
    hot_agg<-cbind(hot_agg,cai_df[,c(5,6,9:18)])
    return(list(data_fr_res,datafr_2,hot_agg))
  }
}

window_randistr<-function(n_mut,gene_length,thrs_hot,wind_size,nperm)
{
  c<-1
  hot_pos<-vector()
  hot_perm<-vector()
  i<-1
  for(i in 1:nperm)
  {
    pos_mut<-sample(1:gene_length,n_mut,replace=T)
    pos_mut<-pos_mut[order(pos_mut,decreasing = F)]
    tab_mut<-table(pos_mut)
    tab_mut<-cbind(tab_mut,as.numeric(names(tab_mut)))
    #I create a table where the first column is the frequency and the last
    #represents the positions. 
    j<-1
    for(j in 1:nrow(tab_mut))
    {
      #I look for the rows where the first chromosome position appears. 
      wind<-subset(tab_mut[,2], tab_mut[,2]>=tab_mut[j,2] & tab_mut[,2]<=(tab_mut[j,2]+wind_size))
      #Computation of the window
      #I want to save the frequency of the sample and the mutations
      freq_m<-sum(tab_mut[tab_mut[,2]%in%wind,1])/n_mut*100
      if(freq_m>=thrs_hot)
      {
        #I don't have to check if it is the first window to save. 
        hot_pos[j]<-1
      }
      else
      {
        hot_pos[j]<-0
      }
    }
    if(sum(hot_pos)>0)
      hot_perm[i]=1
    else
      hot_perm[i]=0
  }
  num_extr1<-(length(which(hot_perm==1))+1)/(nperm+1)
  return(num_extr1)
}