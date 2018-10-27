# /lustre/software/target/R-3.3.2/bin/R
#
# Date：       2018-10-18
# Function：   对<指定基因>在<指定类型数据库>中<指定癌症类型>的数据绘制生存图
# Description: 1. 对 mRNA、methylation、rnaseq 三类数据库，用表达量的高低画生存图
#              2. matations 数据库根据有无突变绘制生存图
#              3. 相应数据表（csv格式）、相应的生存图（pdf格式）输出
#
#------------------------------------------------------------------------------------
# Problems：
#          1. only cancer samples ？？？
#          2. matations 生存图正确？
#          3. 生存时间有负数，什么含义？
#          4. pdf 空白页？（更改为 png）
#          5. 如何在网上查找对应数据库的数据？
#
#------------------------------------------------------------------------------------

#----------------- Time
cat(rep('-',times=70), '\nStart:', sep='')
print(Sys.time())


#----------------- parameters
OutPath     <- '.'
Gene        <- 'ELMO2'
Cancer_type <- 'OV'
DB_type     <- 'mRNA'

Gene_OK <- sub('\\|\\d+', '', Gene)
SurvPNG <- paste(OutPath, paste(Cancer_type,DB_type,Gene_OK,'Survival_Plot.png', sep='_'), sep='/')
OutFile <- paste(OutPath, paste(Cancer_type,DB_type,Gene_OK,'Survival_Data.csv', sep='_'), sep='/')


#----------------- 确认癌症数据版本和癌症表达数据
if(class(try(library(package=paste('RTCGA', DB_type, '20160128', sep='.'),character.only=T),silent=T)) == 'try-error'){
    if(class(try(library(package=paste('RTCGA', DB_type, sep='.'),character.only=T),silent=T)) == 'try-error'){
        stop('!!!!!!!!!! Error database input:', DB_type)
    }else{
        Exp_DB <- try(eval(parse(text=paste(Cancer_type, DB_type, sep='.'))), silent=T)
    }
}else{
    Exp_DB <- try(eval(parse(text=paste(Cancer_type, DB_type, '20160128', sep='.'))), silent=T)
}
library(survminer)
if( DB_type == 'mutations' ){
    Expressions <- mutationsTCGA(Exp_DB) %>%
        dplyr::filter(Hugo_Symbol == Gene_OK) %>%
        dplyr::filter(substr(bcr_patient_barcode, 14, 15) == '01') %>%
        dplyr::mutate(bcr_patient_barcode = substr(bcr_patient_barcode, 1, 12))
}else{
    if(length(Exp_DB) != 1 && Gene %in% colnames(Exp_DB)){
        Expressions <- expressionsTCGA(Exp_DB, extract.cols=Gene) %>%
                        dplyr::filter(substr(bcr_patient_barcode,14,15)=='01') %>% #-------- only cancer samples ??????
                        dplyr::mutate(bcr_patient_barcode=substr(bcr_patient_barcode,1,12))
    }else{
        stop("!!!!!!!!!! Error gene name: ", Gene)
    }
}


#----------------- 确认生存数据版本
if(class(try(library(RTCGA.clinical.20160128),silent=T)) == 'try-error'){
    if(class(try(library(RTCGA.clinical),silent=T)) == 'try-error'){
        stop('!!!!!!!!!! No clinical database !!!!!!')
    }else{
        SurvDB <- try(eval(parse(text=paste(Cancer_type, 'clinical', sep='.'))), silent=T)
    }
}else{
    SurvDB <- try(eval(parse(text=paste(Cancer_type, 'clinical.20160128', sep='.'))), silent=T)
}
if( class(SurvDB) == 'try-error' ) stop('!!!!!!!! Error cancer type in clinical database.')


#---------------- 合并信息、分类、绘图、保存数据文件
library(survival)
if( DB_type == 'mutations' ){
    Surv_Mut <- survivalTCGA(SurvDB) %>% dplyr::left_join(Expressions, by='bcr_patient_barcode') %>% 
                    dplyr::mutate(Classification = ifelse(!is.na(Variant_Classification), 'Mutated', 'Unmutated')) %>%
                    dplyr::select(times, bcr_patient_barcode, patient.vital_status, Classification)

    fit <- survfit(Surv(times, patient.vital_status)~Classification, data=Surv_Mut)

    #- survival plot
    png(SurvPNG, width = 7 ,height = 7,units = "in",res = 300)
    print(ggsurvplot(fit,
                surv.median.line = "hv",
                pval=TRUE,
                conf.int=TRUE,
                conf.int.style = 'ribbon',
                legend=c(0.85,0.9),
                legend.title=Gene_OK,
                legend.labs=c('Mutated', 'Unmutated'),
                risk.table='abs_pct',
                risk.table.fontsize=3.2,
                tables.height=0.2,
                tables.theme=theme_RTCGA(),
                ggtheme=theme_bw(),
    ))
    dev.off()

    #- data
    write.csv(Surv_Mut, file=OutFile, row.names=F)

}else{
    CutFile <- paste(OutPath, paste(Cancer_type,DB_type,Gene_OK,'Cut_value.txt', sep='_'), sep='/')
    CutPNG <- paste(OutPath, paste(Cancer_type,DB_type,Gene_OK,'Expression_Cut_Plot.PNG', sep='_'), sep='/')

    Surv_Exp <- survivalTCGA(SurvDB) %>% dplyr::left_join(Expressions, by='bcr_patient_barcode') %>% dplyr::filter(!is.na(dataset))
    colnames(Surv_Exp)[4:5] <- c('Regulation', Gene_OK)
    Surv_Exp_cut <- surv_cutpoint(Surv_Exp, time='times', event='patient.vital_status', variables=Gene_OK)
    Surv_Exp_cat <- surv_categorize(Surv_Exp_cut)
    colnames(Surv_Exp_cat)[3] <- 'Gene'
    Surv_Exp[,4] <- Surv_Exp_cat[,'Gene']
    fit <- survfit(Surv(times, patient.vital_status)~Gene, data=Surv_Exp_cat)

    #- cut plot
    png(CutPNG, width = 7 ,height = 7,units = "in",res = 300)
    # bins :Number of bins for histogram. Defaults to 30.
    print(plot(Surv_Exp_cut, Gene_OK, palette='npg'), bins=60)
    dev.off()

    #- survival plot
    png(SurvPNG, width = 7 ,height = 7,units = "in",res = 300)
    print(ggsurvplot(fit,
                surv.median.line = "hv",
                pval=TRUE,
                conf.int=TRUE,
                conf.int.style = 'ribbon',
                legend=c(0.85,0.9),
                legend.title=Gene_OK,
                legend.labs=c('High', 'Low'),
                risk.table='abs_pct',
                risk.table.fontsize=3.2,
                tables.height=0.2,
                tables.theme=theme_RTCGA(),
                ggtheme=theme_bw(),
    ))
    dev.off()

    #- data
    write.csv(Surv_Exp, file=OutFile, row.names=F)
    write.table(summary(Surv_Exp_cut), file=CutFile, quote=F, sep='\t', row.names=F)
}


#----------------- Time
cat('End:', sep='')
print(Sys.time())
cat(rep('-',times=70), '\n\n', sep='')


#------------------------------------------- END -----------------------------------------------
