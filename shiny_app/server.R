#source("https://bioconductor.org/biocLite.R")
#biocLite("IMPCdata")

library(IMPCdata)
library(stringr)
library(data.table)
#library(MASS)
library(sfsmisc)
library(rlm)
library(magrittr)

library(shiny)
library(DT)
library(plotly)
library(crosstalk)

# TO DO:
# maybe add O2 smooth but very problematic
# CHANGE shortCenter and METABOLISM REPORT PHENOCENTER NAME !!!

phenoCenters = c('bcm', 'hmgu','ics','kmpc','marc','mrc','rbrc','tcp','ucd')

impc_folder = '/IMPC/' ### TODO CHANGE FOLDER

rdsPath = paste0( impc_folder, 'calorimetryData/iterationKO/shiny_table/')
shortCenter = 'hmgu'
mgiList <- readRDS(paste0(rdsPath,shortCenter,'_mgiList.rds'))

#getIMPC arg: phenCenterName, pipelineId, procedureId, parameterId, alleleId, strainId

####################################
# download necessary IMPCtable
####################################

### get HMGU mgi list
# printPhenCenters()
# pPhen = 'HMGU'
# printPipelines(pPhen)
# pPip = 'HMGU_001'
# printProcedures(pPhen,pPip)
# pProc = 'IMPC_CAL_003'
# printParameters(pPhen,pPip,pProc)
# pPar = 'IMPC_CAL_001_001'
# printStrains(pPhen,pPip,pProc,pPar,n=5)
# 

# printPhenCenters()
# pPhen = 'KMPC'
# printPipelines(pPhen)
# pPip = 'IMPC_001'
# printProcedures(pPhen,pPip)
# pProc = 'IMPC_CAL_003'
# printParameters(pPhen,pPip,pProc)
# pPar = 'IMPC_CAL_001_001'
# printStrains(pPhen,pPip,pProc,pPar,n=5)

### get list of all genotypes, last column: column to download it, only need list of strains
outputTableIC = paste0(impc_folder, 'calorimetryData/iterationKO/',shortCenter,'_impcTable_IC')
# getIMPCTable(outputTableIC,pPhen,pPip,pProc,pPar)
##getIMPCTable(outputTableIC,'HMGU','HMGU_001','IMPC_CAL_003','IMPC_CAL_001_001')  # calorimetry strains


####################################
# dl mgi data from impc and make mgi list
####################################

### handle data and extract only necessary things
# bwSfilterData <- function(dataF){
#   df = dataF
#   df = df[c('external_sample_id','pipeline_name','Strain','Center','date_of_birth','age_in_weeks','Sex','biological_sample_group','MetadataGroup','Weight')]
#   df$Weight = as.numeric(as.character(df$Weight))
#   df$external_sample_id = as.character(df$external_sample_id)
#   df = unique(df)
#   return(df)   }
# 
# bwEfilterData <- function(dataF){
#   df = dataF
#   #df = df[c('external_sample_id','Weight','Value','discrete_point')]
#   df = df[c('external_sample_id','age_in_weeks','Weight')]
#   df$Weight = as.numeric(as.character(df$Weight))
#   #df$Value = as.numeric(as.character(df$Value))
#   #df$discrete_point = as.numeric(as.character(df$discrete_point))
#   df$external_sample_id = as.character(df$external_sample_id)
#   df = unique(df)
#   return(df)   }
# 
# ### merge both dataframes
# bwSEmerge <- function(dfS,dfE) {
#   names(dfS)[names(dfS) == 'Weight'] <- 'bw_start'
#   names(dfE)[names(dfE) == 'Weight'] <- 'bw_end'
#   names(dfE)[names(dfE) == 'age_in_weeks'] <- 'age_in_weeks_end'
#   bwDf = merge(dfS,dfE)
#   bwDf$weightGain = bwDf$bw_end - bwDf$bw_start
#   bwDf$weightGainProc = (bwDf$bw_end-bwDf$bw_start)/bwDf$bw_start
#   return(bwDf)
#   }
# 
# 
# 
# ### extract weight gain table for external_sample_id, returns mgiList with all mgiNO: list(icData,heartData,merged)
# mgiTable = fread(paste0(outputTableIC,"_1.csv"))
# regexMgi = str_match(mgiTable$Allele, "MGI:\\d+")[,1]
# mgiStrains = na.omit(regexMgi)
# ### get MGIs from report
# mgiTable = fread(paste0( impc_folder, 'release-8.0/reports/metabolismCalorimetryReport.csv') )
# # mgiTable = mgiTable[mgiTable$`Phenotyping Center`=='MRC Harwell',]   # <- CHANGE HERE: MARC, HMGU, UC Davis, BCM, WTSI, ICS, TCP, MRC Harwell !!!!!!!!!!!!!!!!!!!!!!!
# mgiTable = mgiTable[mgiTable$`Phenotyping Center`==toupper(shortCenter),]
# mgiStrains2 = unique(mgiTable$`MGI Gene Id`)
# mgiStrains = c(mgiStrains,mgiStrains2)
# 
# #mgiStrains = c('MGI:5501061', 'MGI:5501065','MGI:5501071')
# mgiList = list()
# for(mgi in mgiStrains) {
#   ### specify both experiments pipeline to weight gain, bwS: IC data, bwE: bw at heart weight
#   # bwS = getIMPCDataset('HMGU','HMGU_001','IMPC_CAL_003','IMPC_CAL_001_001',mgi)
#   # bwE = getIMPCDataset('HMGU','HMGU_001','IMPC_HWT_001','IMPC_HWT_007_001',mgi)
#   # bwS = getIMPCDataset('UC Davis','UCD_001','IMPC_CAL_003','IMPC_CAL_001_001',mgi)
#   # bwE = getIMPCDataset('UC Davis','UCD_001','IMPC_HWT_001','IMPC_HWT_007_001',mgi)
#   # bwS = getIMPCDataset('BCM','IMPC_001','IMPC_CAL_003','IMPC_CAL_001_001',mgi)
#   # bwE = getIMPCDataset('BCM','IMPC_001','IMPC_HWT_001','IMPC_HWT_007_001',mgi)
#   # bwS = getIMPCDataset('ICS','ICS_001','IMPC_CAL_003','IMPC_CAL_001_001',mgi)
#   # bwE = getIMPCDataset('ICS','ICS_001','IMPC_HWT_001','IMPC_HWT_007_001',mgi)
#   # bwS = getIMPCDataset('MRC Harwell','HRWL_001','IMPC_CAL_002','IMPC_CAL_001_001',mgi)
#   # bwE = getIMPCDataset('MRC Harwell','HRWL_001','HRWL_OWT_001','IMPC_HWT_007_001',mgi)
#   # bwS = getIMPCDataset('RBRC','IMPC_001','IMPC_CAL_002','IMPC_CAL_001_001',mgi)
#   # bwE = getIMPCDataset('RBRC','IMPC_001','IMPC_HWT_001','IMPC_HWT_007_001',mgi)
#   # bwS = getIMPCDataset('TCP','TCP_001','IMPC_CAL_002','IMPC_CAL_001_001',mgi)
#   # bwE = getIMPCDataset('TCP','TCP_001','IMPC_HWT_001','IMPC_HWT_007_001',mgi)
#   # bwS = getIMPCDataset('MARC','IMPC_001','IMPC_CAL_002','IMPC_CAL_001_001',mgi)
#   # bwE = getIMPCDataset('MARC','IMPC_001','IMPC_HWT_001','IMPC_HWT_007_001',mgi)
#   bwS = getIMPCDataset('KMPC','IMPC_001','IMPC_CAL_002','IMPC_CAL_001_001',mgi)
#   bwE = getIMPCDataset('KMPC','IMPC_001','IMPC_HWT_001','IMPC_HWT_007_001',mgi)
#   
#   if(is.data.frame(bwS) & is.data.frame(bwE)) {   # error at dl data
#     bwS = bwSfilterData(bwS)
#     bwE = bwEfilterData(bwE)
#     bwSE = bwSEmerge(bwS,bwE)
#     mgiList[[gsub(':','',mgi)]] = list('bwStart'=bwS,'bwEnd'=bwE, 'bwMerged' = bwSE)
#   }
#   print(paste0('### ',match(mgi,mgiStrains) ) )
#   closeAllConnections()
# }
# 
# ### make mgi lists for other institutes
# saveRDS(mgiList, paste0(rdsPath,shortCenter,'_mgiList.rds'))


### RUN TILL HERE !!!

####################################
# import pvalue data for mgis 
####################################

getAsterix <- function(pvalue){
  symnum(pvalue, corr = FALSE, na = FALSE, legend = FALSE,
         cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
         symbols = c("***", "**", "*", ".", ""))
}

### import significance p-value report
pvalDf = fread(paste0( impc_folder, 'release-8.0/reports/impcPValuesReport.csv') )
pvalDf = pvalDf[,c('Genotype','Gene Symbol','MGI Gene Id','Center','Respiratory Exchange Ratio(IMPC_CAL_017_001)','Activity (body position)(IMPC_CSD_029_001)','Area under glucose response curve(IMPC_IPG_012_001)','Lean/Body weight(IMPC_DXA_008_001)','Fat/Body weight(IMPC_DXA_009_001)','BMC/Body weight(IMPC_DXA_007_001)')]
colnames(pvalDf) = c('genotype','gene_symbol','MGI','center','RER_diff','activity','AUC_GTT','lean_bw','fat_bw','BMC_bw')
pvalDf$MGI = gsub(':','',pvalDf$MGI)

### get sig value output
checkPval <- function(mgiNo) {
  r = pvalDf[pvalDf$MGI == mgiNo,]
  sigText = '-'
  if(nrow(r)!=0) {
    sigText=''
    a = getAsterix(min(r$RER_diff, na.rm = T))
    if(a!='') {
      sigText = paste0(sigText,'RER_diff',a,' ')
    }
    a = getAsterix(min(r$activity, na.rm = T))
    if(a!='') {
      sigText = paste0(sigText,'act',a,' ')
    }
    a = getAsterix(min(r$AUC_GTT, na.rm = T))
    if(a!='') {
      sigText = paste0(sigText,'AUC_GTT',a,' ')
    }
    a = getAsterix(min(r$lean_bw, na.rm = T))
    if(a!='') {
      sigText = paste0(sigText,'lean/bw',a,' ')
    }
    a = getAsterix(min(r$fat_bw, na.rm = T))
    if(a!='') {
      sigText = paste0(sigText,'fat/bw',a,' ')
    }
    a = getAsterix(min(r$BMC_bw, na.rm = T))
    if(a!='') {
      sigText = paste0(sigText,'BMC/bw',a,' ')
    }   }
  return(trimws(sigText))
}


####################################
# import calorimetry data + smoothed ones
####################################

### import calorimetry data of each single mouse
calDf = fread(paste0( impc_folder, 'release-8.0/reports/metabolismCalorimetryReport.csv') )
colnames(calDf) <- c("Mouse_Id","Sample_Type","Gene_Symbol","MGI_Gene_Id","Allele","Zygosity","Sex","Colony_Id","Phenotyping_Center","Metadata_Group","bw_before","bw_after","O2_min","O2_max","O2_mean","O2_count","CO2_min","CO2_max","CO2_mean","CO2_count","heat_min","heat_max","heat_mean","heat_count","activity_min","activity_max","activity_mean","activity_count","total_activity_min","total_activity_max","total_activity_mean","total_activity_count","total_food_intake","cumulative_food_intake_min","cumulative_food_intake_max","cumulative_food_intake_mean","cumulative_food_intake_count","RER","total_water","cumulative_water_intake_min","cumulative_water_intake_max","cumulative_water_intake_mean","cumulative_water_intake_count","RER_min","RER_max","RER_mean","RER_count","RER_derived_min","RER_derived_max","RER_derived_mean","RER_derived_count","Metadata")
calDf[calDf=='No information available'] <-NA
calDf = cbind(calDf[,c(1:10,52)], as.data.frame(lapply(calDf[,11:51], as.numeric)))
calDf$RER_derived_amp = calDf$RER_derived_max - calDf$RER_derived_min
calDf$O2_amp = calDf$O2_max - calDf$O2_min
calDf$CO2_amp = calDf$CO2_max - calDf$CO2_min
calDf$heat_amp = calDf$heat_max - calDf$heat_min
calDf$total_activity_amp = calDf$total_activity_max - calDf$total_activity_min
calDf$cumulative_food_intake_amp = calDf$cumulative_foot_intake_max - calDf$cumulative_food_intake_min

### add smoothed data to IC report
smoothDf = data.frame()
for(ce in phenoCenters) {
  smoothDf1 = fread(paste0(impc_folder, 'calorimetryData/getRawData/',ce,'_RER_featureTable.csv'))
  smoothDf = rbind(smoothDf,smoothDf1, fill=T)
}
names(smoothDf) = c('Mouse_Id','RER_smooth_min','RER_smooth_max','RER_smooth_mean','RER_smooth_amp','maxPos', 'minPos', 'maxMinDis', 'gradientMinMax', 'AUCall')
calDf = merge(calDf, smoothDf[,1:5], by = 'Mouse_Id', all.x = T)


### get difference between finished and smoothed data
# diffDf = merge(c, calDf, by="Mouse_Id")
# diffRERmin = round(diffDf$RER_derived_min.x - diffDf$RER_derived_min.y, 3)
# diffRERmax = round(diffDf$RER_derived_max.x - diffDf$RER_derived_max.y,3 )
# diffRERmean = round(diffDf$RER_derived_mean.x - diffDf$RER_derived_mean.y,3)
# diffRERamp = round(diffDf$RER_derived_amp.x - diffDf$RER_derived_amp.y,3)
# differenceDf = data.frame('Mouse_Id' = diffDf$Mouse_Id, 'rer_min_diff'=diffRERmin, 'rer_max_diff' = diffRERmax, 'rer_mean_diff' = diffRERmean, 'rer_amp_diff' = diffRERamp)
# summary(differenceDf)
# write.csv(differenceDf, paste0( impc_folder, 'calorimetryData/getRawData/rerdiff.csv', row.names = F, quote = F) )


####################################
# get df of sig associations
####################################

### return summary dataframe of all possible values
getSummaryDf = function(correlationAttr) {
    #correlationAttr = 'RER_derived_max'
    summaryDf = data.frame('mgi'=character(),'sex'=character(),'numMice'=integer() ,'KO_RER_pval'=double(),'WT_RER_pval'=double(), 'KOWT_RER_pval'=double(),'sigWeightDiff'=character(),'supportOurResults'=character(), 'ancova_pval'=double())
    #summaryDf = data.frame('mgi'=character(),'sex'=character(),'numMice'=integer() , 'KO_RER_slope'=double(),'KO_RER_pval'=double(),'WT_RER_slope'=double(),'WT_RER_pval'=double(),'KOWT_RER_slope'=double(),'KOWT_RER_pval'=double(),'sigWeightDiff'=character(),'supportOurResults'=character(), 'ancova_pval'=double())
    
    ### iterate over all mgi and check
    for(mgiNo in names(mgiList)) {
      for(sVec in c('male','female')){
        bwDf = subset(mgiList[[mgiNo]]$bwMerged, Sex==sVec)
        matchDf = merge(calDf, bwDf, by.x='Mouse_Id', by.y = 'external_sample_id')
        
        wtRLM = NULL
        wtPval = NULL
        koRLM = NULL
        koPval = NULL
        kowtRLM = NULL
        koPval = NULL
    
        ### get linear regression models and check sig towards our results
        supRes=""
    
        tryCatch({
          #wtRLM = rlm(matchDf[matchDf$Sample_Type=='control'][[correlationAttr]] ~ matchDf[matchDf$Sample_Type=='control'][['weightGainProc']])
          #wtPval = f.robftest(wtRLM)$p.value
          wtRLM = lm(matchDf[matchDf$Sample_Type=='control'][[correlationAttr]] ~ matchDf[matchDf$Sample_Type=='control'][['weightGainProc']],na.action=na.omit)
          wtPval = cor.test(matchDf[matchDf$Sample_Type=='control'][[correlationAttr]],matchDf[matchDf$Sample_Type=='control'][['weightGainProc']],method="pearson")$p.value
          if(wtRLM$coefficients[[2]]<(-0.1) & wtPval < 0.1) { supRes = paste0('WT') }
        },  error = function(e) {
          wtRLM = NULL
          wtPval = NULL
        })
        tryCatch({
          #koRLM = rlm(matchDf[matchDf$Sample_Type!='control'][[correlationAttr]] ~ matchDf[matchDf$Sample_Type!='control'][['weightGainProc']])
          #koPval = f.robftest(koRLM)$p.value
          koRLM = lm(matchDf[matchDf$Sample_Type!='control'][[correlationAttr]] ~ matchDf[matchDf$Sample_Type!='control'][['weightGainProc']], na.action=na.omit)
          koPval = cor.test(matchDf[matchDf$Sample_Type!='control'][[correlationAttr]],matchDf[matchDf$Sample_Type!='control'][['weightGainProc']],method="pearson")$p.value
        if(koRLM$coefficients[[2]]<(-0.1) & koPval < 0.1) {  supRes = paste0(supRes,' KO ') }
        },  error = function(e) {
          koRLM = NULL
          koPval = NULL
        })
        tryCatch({
          # kowtRLM = rlm(matchDf[[correlationAttr]] ~ matchDf[['weightGainProc']])
          # kowtPval =  f.robftest(kowtRLM)$p.value
          kowtRLM = lm(matchDf[[correlationAttr]] ~ matchDf[['weightGainProc']], na.action=na.omit)
          kowtPval =  cor.test(matchDf[[correlationAttr]],matchDf[['weightGainProc']], method='pearson')$p.value
          if(kowtRLM$coefficients[[2]]<(-0.1) & kowtPval < 0.1) { supRes = paste0(supRes, ' KO+WT')  }
        },  error = function(e) {
          kowtRLM = NULL
          kowtPval = NULL
        })
    
        supRes = trimws(supRes)
        tryCatch( {
          ko_slope = ifelse(is.null(koRLM), 0, koRLM$coefficients[[2]])
          wt_slope = ifelse(is.null(wtRLM), 0, wtRLM$coefficients[[2]])
          kowt_slope = ifelse(is.null(kowtRLM), 0, kowtRLM$coefficients[[2]])
          ko_pval = ifelse(is.null(koPval), NA, koPval)
          wt_pval = ifelse(is.null(wtPval), NA, wtPval)
          kowt_pval = ifelse(is.null(kowtPval), NA, kowtPval)
    
          matchDf$Sample_Type = as.factor(matchDf$Sample_Type)
          formu = as.formula(paste("weightGainProc~",correlationAttr,'*Sample_Type' ))
          ancova = aov(formu,data = matchDf, na.action=na.omit)
          ancova_pval = summary(ancova)[[1]]$'Pr(>F)'[3]
    
        addDf = data.frame('mgi'=mgiNo,'sex'=sVec,'numMice'=dim(matchDf)[1],'KO_RER_pval'=ko_pval,'WT_RER_pval'=wt_pval,'KOWT_RER_pval'=kowt_pval,'sigWeightDiff'=checkPval(mgiNo),'supportOurResults'=supRes, 'ancova_pval'=ancova_pval)
        #addDf = data.frame('mgi'=mgiNo,'sex'=sVec,'numMice'=dim(matchDf)[1],'KO_RER_slope'=ko_slope,'KO_RER_pval'=ko_pval,'WT_RER_slope'=wt_slope,'WT_RER_pval'=wt_pval,'KOWT_RER_slope'=kowt_slope,'KOWT_RER_pval'=kowt_pval,'sigWeightDiff'=checkPval(mgiNo),'supportOurResults'=supRes, 'ancova_pval'=ancova_pval)
        
        summaryDf = rbind(summaryDf, addDf)
        },  error = function(e) {
          print(paste0(mgiNo,' - ', sVec))
        })

        } # end sex
    } # end mgiNo
    summaryDf$ancova_pval_adj = p.adjust(summaryDf$ancova_pval, method = 'fdr')
    summaryDfRounded = data.frame(lapply(summaryDf, function(y) if(is.numeric(y)) round(y, 4) else y)) 
    return(summaryDfRounded)
}


### iterate over all attributes to generate summary tables
# for(ce in phenoCenters) {
#   mgiList <- readRDS(paste0(rdsPath,ce,'_mgiList.rds'))
#   summaryMeas = c('RER_derived_', 'O2_', 'CO2_', 'heat_','RER_smooth_')
#   summaryFeat = c('max','min','mean','amp')
#   for(m in summaryMeas) {
#     for(f in summaryFeat){
#       ta = getSummaryDf(paste0(m,f))
#       write.csv(ta, paste0(rdsPath,'summary_table/',ce,'_',m,f,'.csv') )
#       print(paste0("#### saved ",m,f,' ####'))
#     }
#   }
#   print(paste0('#### ', ce,' ####') )
# }
# system("say the script is finished")

# write.csv(summaryDf, paste0( impc_folder, 'calorimetryData/iterationKO/hmgu_summary_table.csv') )



### make mgi lists for other institutes
#saveRDS(mgiList, paste0(rdsPath,'hmgu_mgiList.rds'))
# mgiList <- readRDS(paste0(rdsPath,'hmgu_mgiList.rds'))

# redefine mgilist 


####################################
# raw data outlier removal - removed
####################################
### read in IC data and look for big variance, call them and later remove
# e = read.csv(paste0( impc_folder, 'calorimetryData/getRawData/hmgu_RER.csv'), header = T)
# e$time=NULL
# v = apply(e,2,var, na.rm=T)
# # hist(v, breaks=20)   # plot histogram
# # lines(density(v))
# # abline(v=0.015, col='red')
# # var(e$X30274104, na.rm=T)
# varLimit = 0.015 
# outlierRawData = gsub('X', '', names(v[v>varLimit]) )


####################################
# shiny server
####################################
server = function(input,output, session) {
  
  #sharedDf <<- SharedData$new(read.csv(paste0(rdsPath,'summary_table/hmgu_RER_derived_max.csv')))
  sharedDf <<- NULL
  mgiList <<- NULL
  allDf <<- fread(paste0(impc_folder, 'calorimetryData/compareRawData/allCentersDf.csv'))
  
  observe({
    mgiList <<- readRDS(paste0(rdsPath,input$selectInst,'mgiList.rds'))
    ta = read.csv(paste0(rdsPath,'summary_table/',input$selectInst,input$selectMeas,input$selectFeature,'.csv'))
    #sharedDf <<- SharedData$new(ta)
    ta[,1] = NULL  # remove X
    sharedDf <<- ta
    output$tableKOmice = DT::renderDataTable(DT::datatable(sharedDf,selection = 'single', rownames=F, extensions='Buttons', filter='bottom', options = list(searching=T, pageLength=10, dom='Blfrtip',lengthMenu=list( c(10,20,-1), c('10','20','all'))))  ,server=FALSE)
    updateSelectInput(session, "selectIC", selected = input$selectMeas  )
    })
  
  
  ### function to plot attributes for each mgi
  plotAttr = function(yAttr, xAttr) {
    selectedRow = input$tableKOmice_rows_selected
    if(length(selectedRow)) {
      mgiRow = sharedDf[selectedRow[length(selectedRow)],]
      #print(str(mgiRow))
      bwDf = subset(mgiList[[as.character(mgiRow$mgi)]]$bwMerged, Sex==as.character(mgiRow$sex))
      matchDf = merge(calDf, bwDf, by.x='Mouse_Id', by.y = 'external_sample_id')
      matchDf$biological_sample_group= factor(matchDf$biological_sample_group, levels=c('control','experimental'))
      #print(matchDf)
      matchDf = matchDf[!is.na(matchDf[[xAttr]]) & !is.na(matchDf[[yAttr]]),]
      #matchDf = matchDf[!matchDf$Mouse_Id %in% outlierRawData,] # REMOVE OUTLIERS !!!
      yLimit = c(min(matchDf[[yAttr]],na.rm = T)-0.1, max(matchDf[[yAttr]], na.rm = T)+0.1)
      xLimit = c(min(matchDf[[xAttr]],na.rm = T)-0.1, max(matchDf[[xAttr]],na.rm = T)+0.1)
      if(input$selectIC =='RER_derived_') { xLimit[2] = min(xLimit[2], 1.4) }  # only show values < 1.4, do not get removed for statistics !!
      
      plot(matchDf[matchDf$biological_sample_group=='control'][[xAttr]],matchDf[matchDf$biological_sample_group=='control'][[yAttr]],pch=20, cex=0.5, col='black', xlab='',ylab='', xlim=xLimit, ylim=yLimit)
      wtReg = lm(matchDf[matchDf$biological_sample_group=='control'][[yAttr]] ~ matchDf[matchDf$biological_sample_group=='control'][[xAttr]])
      abline(wtReg, col='black')
      
      try({   # no experimental values available
        points(matchDf[matchDf$biological_sample_group!='control'][[xAttr]],matchDf[matchDf$biological_sample_group!='control'][[yAttr]],pch=20, cex=0.8, col='red')
        koReg = lm(matchDf[matchDf$biological_sample_group!='control'][[yAttr]] ~ matchDf[matchDf$biological_sample_group!='control'][[xAttr]])
        if(nrow(matchDf[matchDf$biological_sample_group!='control'])>1){ abline(koReg, col='red') }
        frm<-paste0(yAttr,'~',xAttr,'*Sample_Type')
        ancova = aov(formula(frm),data = matchDf, na.action=na.omit)
        ancova_pval = summary(ancova)[[1]]$'Pr(>F)'[3]
      })
      
      title(xAttr, line = 1.3)
      try({  if(ancova_pval<0.05) { mtext( paste0("[sig ancova ", round(ancova_pval,4),' ]'), line=2.3, cex=0.8, col='red') } else {  mtext( paste0("[ancova ", round(ancova_pval,4),' ]'), line=2.3, cex=0.8) }    })
      title(ylab = yAttr, xlab= xAttr, line=2)
      try({
      if( (anova(wtReg)[,5][1] >0.05) || (anova(koReg)[,5][1] >0.05 ) ) {
          mtext(paste0('p-val: WT:',round(anova(wtReg)[,5][1],4),' KO: ',round(anova(koReg)[,5][1],4)), side=3, line=0.2, cex = 0.8)  }
      else {  mtext(paste0('p-val: WT:',round(anova(wtReg)[,5][1],4),' KO: ',round(anova(koReg)[,5][1],4)), side=3, line=0.2, cex = 0.8, col='red') }
      })
      
      #pp = matchDf %>% plot_ly() %>% add_trace(x=~weightGainProc, y=~RER_derived_max, color=~factor(biological_sample_group), colors=c('black','red'), text=matchDf$Mouse_Id, hoverinfo='text') %>% layout(title=paste0(mgiRow$mgi,': RER - weight gain[%]'))
    }
  }


output$plotKOmice_min = renderPlot({ plotAttr('weightGainProc', paste0(input$selectIC,'min') )  })
output$plotKOmice_max = renderPlot({ plotAttr('weightGainProc',paste0(input$selectIC,'max') )  })
output$plotKOmice_mean = renderPlot({ plotAttr('weightGainProc',paste0(input$selectIC,'mean') ) })
output$plotKOmice_amp = renderPlot({ plotAttr('weightGainProc',paste0(input$selectIC,'amp') ) })


  ### interactive plotly for experimental outlier info, for smoothed night data only
  plotlyScatterAttr = function(sexM,attr, onlyControl) {
    plotDf=allDf
    if(sexM != 'both'){  plotDf = plotDf[plotDf$sex==sexM,] }
    if(onlyControl=='ko'){ plotDf = plotDf[plotDf$group!='control',] }
    else{   if(onlyControl=='wt') {  plotDf = plotDf[plotDf$group=='control',] }   }
    
    p = plotDf %>% plot_ly(mode='markers', type='scattergl', showlegend=T, source='plotAll')  %>%
      add_trace(y=~weightGainProc, x=plotDf[[attr]], color=~factor(center), text=paste0(plotDf$sample,'<br>',plotDf$mgi),key=plotDf$sample, hoverinfo='text') %>%
      layout(title=paste0('RER ',sexM,' ',attr, ' against weightGain'), xaxis = list(title=paste0('RER ',attr), range=c(min(plotDf[[attr]])-0.1,max(plotDf[[attr]])+0.1) ), yaxis=list(title="weightGain [*100=%]", range=c(min(plotDf[['weightGainProc']])-0.02,max(plotDf[['weightGainProc']])+0.02)) )

    # mark in table clicked mgi     
    selectedRow = input$tableKOmice_rows_selected
    if(length(selectedRow)) {
      mgiNo = sharedDf[selectedRow[length(selectedRow)],]$mgi
      markedDf = plotDf[grep(mgiNo,plotDf$mgi),]
      p = p %>% add_trace(data = markedDf, y=~weightGainProc, x=markedDf[[attr]],color=~factor(group), name='selected mgi', text=paste0(markedDf$sample,'<br>',markedDf$mgi,'<br>',markedDf$center),key=markedDf$sample, hoverinfo='text')              
      # p = p %>% add_trace(data = markedDf, y=~weightGainProc, x=markedDf[[attr]],color=~factor(group), marker=list(size=12, color='rgba(255,255,0,0.5)', line=list(color='rgba(255,0,0,0.8', width=2)), colors=I('blue'), name='selected mgi', text=paste0(markedDf$sample,'<br>',markedDf$mgi,'<br>',markedDf$center),key=markedDf$sample, hoverinfo='text')              
    }
    return(p)
  }
  
  output$plotlyAll = renderPlotly( plotlyScatterAttr(input$selectAllSex, input$selectAllFeature, input$selectAllControl) )

  
  ### plot single sampleId raw data curve with smoothing
  plotSingleDatacurve = function(sampleId) {
    getSample = allDf[allDf$sample==sampleId,]
    getSampleCe = getSample$center
    curveDf = fread(paste0(impc_folder,'calorimetryData/getRawData/',getSampleCe,'_RER_night.csv'))
    curveDf = curveDf[[sampleId]]
    curveTime = 1:length(curveDf)
    plot(curveTime,curveDf,pch=20, col='grey', ylab='RER', xlab='time steps', type='o', ylim = c((min(curveDf,na.rm = T)-0.1),(max(curveDf,na.rm = T)+0.1)) )
    smooth = loess(curveDf ~ curveTime )
    lines( 1:length(predict(smooth)), predict(smooth), col = 'red', lwd=2)
    title(paste0('ID ',sampleId,' from ',getSampleCe ), line = 1)
  }
  
  output$plot_singleMouse = renderPlot( {
    sa = event_data("plotly_click", source='plotAll')$key
    # print(sa)
    if(!is.null(sa)){  plotSingleDatacurve(sa)  }
  } )
  
  
  ### gene info difference
  getGenePvalues = function(geneId){
    dfPvalues = fread(paste0( impc_folder, 'release-8.0/reports/impcPValuesReport.csv') )
    # df = dfPvalues[mapply(grepl, toupper(geneId), toupper(dfPvalues$`Gene Symbol`)),]
    df = dfPvalues[toupper(dfPvalues$`Gene Symbol`)==toupper(geneId),]
    outStr = ''
    for(ro in nrow(df)){
      outStr = paste0(outStr,'Gene:',df[ro,3],'  ',df[ro,4],'  Center:',df[ro,5], '\n' )
      for(sigCol in 6:dim(df)[2]){  # iterate all columns to find significant ones
        pval = df[ro,][[sigCol]]
        if(!is.na(pval) && pval<0.05) {
          sigName = strsplit(names(df)[sigCol],"\\(")[[1]][1]
          outStr = paste0(outStr,sigName,':  ',round(pval,5), '\n')
        }
      }
    } # end for ro
    return(outStr)
  }
  
  ### gene phenotype info
  getGenePhenotype = function(geneId) {
    dfPheno = fread(paste0( impc_folder, 'release-8.0/reports/phenotypeOverviewPerGeneReport.csv') )
    df = dfPheno[mapply(grepl, toupper(geneId), toupper(dfPheno$`Gene Symbol`)),]
    outStr = 'Phenotypes: '
    if(nrow(df)!=0){
      hits = gsub('::','\n',df$`Phenotype Hits`)
      outStr = paste0(outStr,df$`MGI Gene Id`,'  Gene:',df$`Gene Symbol`,'\n',hits,'\n\n')
    }
    return(outStr)
  }
  
  observeEvent( input$buttonShowGene, {
    selectedRow = input$tableKOmice_rows_selected
    if(length(selectedRow)) {
      mgiNo = as.character(sharedDf[selectedRow[length(selectedRow)],]$mgi)  # get gene symbol of mgi
      mgiNo = paste0(substr(mgiNo,1,3),':',substr(mgiNo,4,nchar(mgiNo))) # insert :
      mgi_gene_list = fread(paste0( impc_folder, 'calorimetryData/iterationKO/shiny_table/mgi_gene_list.csv') )
      # print(mgiNo)
      # mgiNo = 'MGI:5665212'
      # mgiNo = 'MGI:5548328'
      geneId = mgi_gene_list[mgi_gene_list$mgi==mgiNo,]$gene
      withProgress(message = 'searching gene phenotypes', value = 0, {
        if(length(geneId)!=0) {
            renderStr = paste0(getGenePvalues(geneId),'\n',getGenePhenotype(geneId) )
            incProgress(80)
            output$textShowGene = renderText({ renderStr })
          } else {
          ### get gene symbol for unknown mgi from JAX
          mgiCmd = paste0('perl /Users/stefanloipfinger/Documents/IMPC/calorimetryData/iterationKO/shiny_table/mgiGetGene.pl ',mgiNo)
          geneId = system(mgiCmd, intern = T)
          incProgress(30)
          if(geneId!='no_match'){
            renderStr = paste0(getGenePvalues(geneId),'\n',getGenePhenotype(geneId) )
            incProgress(80)
            output$textShowGene = renderText({ renderStr })
          }
          else{ output$textShowGene = renderText({ 'not found' }) }
        }
      }) # progress end
      
      
    }
  })
  
    
  
  
} # end server
  










