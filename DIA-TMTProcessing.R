##Read and modify DDA files from Hyperplex containing Reporter Ion intensity Values
quansdda<- function (x) {
  data<-read.table(x, sep='\t', header=T, stringsAsFactors = F)
  data<-data[!duplicated(data$First.Scan),]
  data<-data[data$RT.in.min < 65,] #Remove scans from wash peak
  data<-data.frame(data[,grepl("Reporter.Ion.Intensity.",names(data))],
                    rtmin=data$RT.in.min,
                    prec=data$Precursor.mz.in.Da,
                    scan=data$First.Scan,
                    scan.id=paste(data$First.Scan, file, sep="_"),
                    na=apply(is.na(data[,grepl("Reporter.Ion.Intensity.",names(data))]), 1, sum))
  rtindex<-seq(from = 0, to = 80, by=0.2); mzindex<-rep(375:1200)
  data$mz<-match.closest(data$prec, mzindex, tolerance = 1, nomatch = NA_integer_) #Match for isolation window
  data$rt<-match.closest(data$rtmin, rtindex, tolerance = 0.2, nomatch = NA_integer_) #Match for RT
  }
##Read and modify DIA data from Hyperplex containing Reporter Ion intensity Values
quansdia<- function (x) {data<-read.table(x, sep='\t', header=T, stringsAsFactors = F)
data<-data[!duplicated(data$First.Scan),]
data <- data[order(data$First.Scan), ]
data<-data[data$RT.in.min < 65,] #Remove scans from wash peak
data<-data.frame(data[,grepl("Reporter.Ion.Intensity.",names(data))],
                 mz=rep(1:80, each=1, length.out=nrow(data)),
                 rt=rep(1:620, each=80, length.out=nrow(data)),
                 First.Scan=data$First.Scan, RT.in.min=data$RT.in.min, prec=data$Precursor.mz.in.Da, 
                 scan=paste(data$First.Scan, file, sep="_"))
                        }
##Normalizing Reporter Ion Intensity according to their loadings from DDA files
ndda1<-function(df1){
  target <- mean(c(colSums(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], na.rm = T)))
  norm_facs <- target / colSums(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], na.rm = T)
  mix1_sl <- data.frame(sweep(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], 2, norm_facs, FUN = "*"),
                        mz=df1$mz, rt=df1$rt, na=df1$na, scan=df1$scan)
  }
ndda3<-function(df1, df2, df3){
  target <- mean(c(colSums(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], na.rm = T),
                   colSums(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], na.rm = T),
                   colSums(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], na.rm = T)))
  norm_facs <- target / colSums(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], na.rm = T)
  mix1_sl <- data.frame(sweep(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], 2, norm_facs, FUN = "*"),
                        mz=df1$mz, rt=df1$rt, na=df1$na, scan=df1$scan)
  norm_facs <- target / colSums(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], na.rm = T)
  mix2_sl <- data.frame(sweep(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], 2, norm_facs, FUN = "*"),
                        mz=df2$mz, rt=df2$rt, na=df2$na, scan=df2$scan)
  norm_facs <- target / colSums(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], na.rm = T)
  mix3_sl <- data.frame(sweep(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], 2, norm_facs, FUN = "*"),
                        mz=df3$mz, rt=df3$rt, na=df3$na, scan=df3$scan)
  join<-inner_join(mix1_sl, mix2_sl, by=c('mz', 'rt'))
  join<-join[order(join$na.x),]; join<-join[!duplicated(join$scan.x),]
  join<-join[order(join$na.y),]; join<-join[!duplicated(join$scan.y),]
  join<-inner_join(join,mix3_sl, by=c('mz', 'rt'))
  join<-join[order(join$na.x),];join<-join[!duplicated(join$scan.x),]
  }
ndda4<-function(df1, df2, df3, df4){
  target <- mean(c(colSums(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], na.rm = T),
                   colSums(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], na.rm = T),
                   colSums(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], na.rm = T),
                   colSums(df4[,grepl("Reporter.Ion.Intensity.",names(df4))], na.rm = T)))
  norm_facs <- target / colSums(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], na.rm = T)
  mix1_sl <- data.frame(sweep(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], 2, norm_facs, FUN = "*"),
                        mz=df1$mz, rt=df1$rt, na=df1$na, scan=df1$scan)
  norm_facs <- target / colSums(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], na.rm = T)
  mix2_sl <- data.frame(sweep(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], 2, norm_facs, FUN = "*"),
                        mz=df2$mz, rt=df2$rt, na=df2$na, scan=df2$scan)
  norm_facs <- target / colSums(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], na.rm = T)
  mix3_sl <- data.frame(sweep(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], 2, norm_facs, FUN = "*"),
                        mz=df3$mz, rt=df3$rt, na=df3$na, scan=df3$scan)
  norm_facs <- target / colSums(df4[,grepl("Reporter.Ion.Intensity.",names(df4))], na.rm = T)
  mix4_sl <- data.frame(sweep(df4[,grepl("Reporter.Ion.Intensity.",names(df4))], 2, norm_facs, FUN = "*"),
                        mz=df4$mz, rt=df4$rt, na=df4$na, scan=df4$scan)
  join<-inner_join(mix1_sl, mix2_sl, by=c('mz', 'rt'))
  join<-join[order(join$na.x),]; join<-join[!duplicated(join$scan.x),]
  join<-join[order(join$na.y),]; join<-join[!duplicated(join$scan.y),]
  join<-inner_join(join,mix3_sl, by=c('mz', 'rt'))
  join<-join[order(join$na.x),];join<-join[!duplicated(join$scan.x),]
  join<-inner_join(join, mix4_sl, by=c('mz', 'rt'))
  join<-join[order(join$na.x),];join<-join[!duplicated(join$scan.x),]
  }
ndda12<-function(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12){
  target <- mean(c(colSums(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], na.rm = T),
                   colSums(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], na.rm = T),
                   colSums(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], na.rm = T),
                   colSums(df4[,grepl("Reporter.Ion.Intensity.",names(df4))], na.rm = T),
                   colSums(df5[,grepl("Reporter.Ion.Intensity.",names(df5))], na.rm = T),
                   colSums(df6[,grepl("Reporter.Ion.Intensity.",names(df6))], na.rm = T),
                   colSums(df7[,grepl("Reporter.Ion.Intensity.",names(df7))], na.rm = T),
                   colSums(df8[,grepl("Reporter.Ion.Intensity.",names(df8))], na.rm = T),
                   colSums(df9[,grepl("Reporter.Ion.Intensity.",names(df9))], na.rm = T),
                   colSums(df10[,grepl("Reporter.Ion.Intensity.",names(df10))], na.rm = T),
                   colSums(df11[,grepl("Reporter.Ion.Intensity.",names(df11))], na.rm = T),
                   colSums(df12[,grepl("Reporter.Ion.Intensity.",names(df12))], na.rm = T)))
  norm_facs <- target / colSums(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], na.rm = T)
  mix1_sl <- data.frame(sweep(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], 2, norm_facs, FUN = "*"),
                        mz=df1$mz, rt=df1$rt, na=df1$na, scan=df1$scan)
  norm_facs <- target / colSums(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], na.rm = T)
  mix2_sl <- data.frame(sweep(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], 2, norm_facs, FUN = "*"),
                        mz=df2$mz, rt=df2$rt, na=df2$na, scan=df2$scan)
  norm_facs <- target / colSums(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], na.rm = T)
  mix3_sl <- data.frame(sweep(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], 2, norm_facs, FUN = "*"),
                        mz=df3$mz, rt=df3$rt, na=df3$na, scan=df3$scan)
  norm_facs <- target / colSums(df4[,grepl("Reporter.Ion.Intensity.",names(df4))], na.rm = T)
  mix4_sl <- data.frame(sweep(df4[,grepl("Reporter.Ion.Intensity.",names(df4))], 2, norm_facs, FUN = "*"),
                        mz=df4$mz, rt=df4$rt, na=df4$na, scan=df4$scan)
  norm_facs <- target / colSums(df5[,grepl("Reporter.Ion.Intensity.",names(df5))], na.rm = T)
  mix5_sl <- data.frame(sweep(df5[,grepl("Reporter.Ion.Intensity.",names(df5))], 2, norm_facs, FUN = "*"),
                        mz=df5$mz, rt=df5$rt, na=df5$na, scan=df5$scan)
  norm_facs <- target / colSums(df6[,grepl("Reporter.Ion.Intensity.",names(df6))], na.rm = T)
  mix6_sl <- data.frame(sweep(df6[,grepl("Reporter.Ion.Intensity.",names(df6))], 2, norm_facs, FUN = "*"),
                        mz=df6$mz, rt=df6$rt, na=df6$na, scan=df6$scan)
  norm_facs <- target / colSums(df7[,grepl("Reporter.Ion.Intensity.",names(df7))], na.rm = T)
  mix7_sl <- data.frame(sweep(df7[,grepl("Reporter.Ion.Intensity.",names(df7))], 2, norm_facs, FUN = "*"),
                        mz=df7$mz, rt=df7$rt, na=df7$na, scan=df7$scan)
  norm_facs <- target / colSums(df8[,grepl("Reporter.Ion.Intensity.",names(df8))], na.rm = T)
  mix8_sl <- data.frame(sweep(df8[,grepl("Reporter.Ion.Intensity.",names(df8))], 2, norm_facs, FUN = "*"),
                        mz=df8$mz, rt=df8$rt, na=df8$na, scan=df8$scan)
  norm_facs <- target / colSums(df9[,grepl("Reporter.Ion.Intensity.",names(df9))], na.rm = T)
  mix9_sl <- data.frame(sweep(df9[,grepl("Reporter.Ion.Intensity.",names(df9))], 2, norm_facs, FUN = "*"),
                        mz=df9$mz, rt=df9$rt, na=df9$na, scan=df9$scan)
  norm_facs <- target / colSums(df10[,grepl("Reporter.Ion.Intensity.",names(df10))], na.rm = T)
  mix10_sl <- data.frame(sweep(df10[,grepl("Reporter.Ion.Intensity.",names(df10))], 2, norm_facs, FUN = "*"),
                         mz=df10$mz, rt=df10$rt, na=df10$na, scan=df10$scan)
  norm_facs <- target / colSums(df11[,grepl("Reporter.Ion.Intensity.",names(df11))], na.rm = T)
  mix11_sl <- data.frame(sweep(df11[,grepl("Reporter.Ion.Intensity.",names(df11))], 2, norm_facs, FUN = "*"),
                         mz=df11$mz, rt=df11$rt, na=df11$na, scan=df11$scan)
  norm_facs <- target / colSums(df12[,grepl("Reporter.Ion.Intensity.",names(df12))], na.rm = T)
  mix12_sl <- data.frame(sweep(df12[,grepl("Reporter.Ion.Intensity.",names(df12))], 2, norm_facs, FUN = "*"),
                         mz=df12$mz, rt=df12$rt, na=df12$na, scan=df12$scan)
  join<-inner_join(mix1_sl, mix2_sl, by=c('mz', 'rt'))
  join<-join[order(join$na.x),]; join<-join[!duplicated(join$scan.x),]
  join<-join[order(join$na.y),]; join<-join[!duplicated(join$scan.y),]
  join1<-inner_join(join,mix3_sl, by=c('mz', 'rt'))
  join1<-join1[order(join1$na.x),];join1<-join1[!duplicated(join1$scan.x),]
  
  join<-inner_join(mix4_sl, mix5_sl, by=c('mz', 'rt'))
  join<-joina[order(join$na.x),]; join<-join[!duplicated(join$scan.x),]
  join<-joina[order(join$na.y),]; join<-join[!duplicated(join$scan.y),]
  join2<-inner_join(join,mix6_sl, by=c('mz', 'rt'))
  join2<-join2[order(join2$na.x),];join2<-join2[!duplicated(join2$scan.x),]
  
  join<-inner_join(mix7_sl, mix8_sl, by=c('mz', 'rt'))
  join<-join[order(join$na.x),]; join<-join[!duplicated(join$scan.x),]
  join<-join[order(join$na.y),]; join<-join[!duplicated(join$scan.y),]
  join3<-inner_join(join,mix9_sl, by=c('mz', 'rt'))
  join3<-join3[order(join3$na.x),];join3<-join3[!duplicated(join3$scan.x),]
  
  join<-inner_join(mix10_sl, mix11_sl, by=c('mz', 'rt'))
  join<-join[order(join$na.x),]; join<-join[!duplicated(join$scan.x),]
  join<-join[order(join$na.y),]; join<-join[!duplicated(join$scan.y),]
  join4<-inner_join(join,mix12_sl, by=c('mz', 'rt'))
  join4<-join4[order(join4$na.x),];join4<-join4[!duplicated(join4$scan.x),]
  
  join<-inner_join(join1, join2, by=c('mz', 'rt'))
  join<-inner_join(join, join3, by=c('mz', 'rt'))
  join<-inner_join(join, join4, by=c('mz', 'rt'))
  }
ndda16<-function(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15, df16){
  target <- mean(c(colSums(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], na.rm = T),
                   colSums(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], na.rm = T),
                   colSums(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], na.rm = T),
                   colSums(df4[,grepl("Reporter.Ion.Intensity.",names(df4))], na.rm = T),
                   colSums(df5[,grepl("Reporter.Ion.Intensity.",names(df5))], na.rm = T),
                   colSums(df6[,grepl("Reporter.Ion.Intensity.",names(df6))], na.rm = T),
                   colSums(df7[,grepl("Reporter.Ion.Intensity.",names(df7))], na.rm = T),
                   colSums(df8[,grepl("Reporter.Ion.Intensity.",names(df8))], na.rm = T),
                   colSums(df9[,grepl("Reporter.Ion.Intensity.",names(df9))], na.rm = T),
                   colSums(df10[,grepl("Reporter.Ion.Intensity.",names(df10))], na.rm = T),
                   colSums(df11[,grepl("Reporter.Ion.Intensity.",names(df11))], na.rm = T),
                   colSums(df12[,grepl("Reporter.Ion.Intensity.",names(df12))], na.rm = T),
                   colSums(df13[,grepl("Reporter.Ion.Intensity.",names(df13))], na.rm = T),
                   colSums(df14[,grepl("Reporter.Ion.Intensity.",names(df14))], na.rm = T),
                   colSums(df15[,grepl("Reporter.Ion.Intensity.",names(df15))], na.rm = T),
                   colSums(df16[,grepl("Reporter.Ion.Intensity.",names(df16))], na.rm = T)))
  norm_facs <- target / colSums(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], na.rm = T)
  mix1_sl <- data.frame(sweep(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], 2, norm_facs, FUN = "*"),
                        mz=df1$mz, rt=df1$rt, na=df1$na, scan=df1$scan)
  norm_facs <- target / colSums(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], na.rm = T)
  mix2_sl <- data.frame(sweep(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], 2, norm_facs, FUN = "*"),
                        mz=df2$mz, rt=df2$rt, na=df2$na, scan=df2$scan)
  norm_facs <- target / colSums(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], na.rm = T)
  mix3_sl <- data.frame(sweep(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], 2, norm_facs, FUN = "*"),
                        mz=df3$mz, rt=df3$rt, na=df3$na, scan=df3$scan)
  norm_facs <- target / colSums(df4[,grepl("Reporter.Ion.Intensity.",names(df4))], na.rm = T)
  mix4_sl <- data.frame(sweep(df4[,grepl("Reporter.Ion.Intensity.",names(df4))], 2, norm_facs, FUN = "*"),
                        mz=df4$mz, rt=df4$rt, na=df4$na, scan=df4$scan)
  norm_facs <- target / colSums(df5[,grepl("Reporter.Ion.Intensity.",names(df5))], na.rm = T)
  mix5_sl <- data.frame(sweep(df5[,grepl("Reporter.Ion.Intensity.",names(df5))], 2, norm_facs, FUN = "*"),
                        mz=df5$mz, rt=df5$rt, na=df5$na, scan=df5$scan)
  norm_facs <- target / colSums(df6[,grepl("Reporter.Ion.Intensity.",names(df6))], na.rm = T)
  mix6_sl <- data.frame(sweep(df6[,grepl("Reporter.Ion.Intensity.",names(df6))], 2, norm_facs, FUN = "*"),
                        mz=df6$mz, rt=df6$rt, na=df6$na, scan=df6$scan)
  norm_facs <- target / colSums(df7[,grepl("Reporter.Ion.Intensity.",names(df7))], na.rm = T)
  mix7_sl <- data.frame(sweep(df7[,grepl("Reporter.Ion.Intensity.",names(df7))], 2, norm_facs, FUN = "*"),
                        mz=df7$mz, rt=df7$rt, na=df7$na, scan=df7$scan)
  norm_facs <- target / colSums(df8[,grepl("Reporter.Ion.Intensity.",names(df8))], na.rm = T)
  mix8_sl <- data.frame(sweep(df8[,grepl("Reporter.Ion.Intensity.",names(df8))], 2, norm_facs, FUN = "*"),
                        mz=df8$mz, rt=df8$rt, na=df8$na, scan=df8$scan)
  norm_facs <- target / colSums(df9[,grepl("Reporter.Ion.Intensity.",names(df9))], na.rm = T)
  mix9_sl <- data.frame(sweep(df9[,grepl("Reporter.Ion.Intensity.",names(df9))], 2, norm_facs, FUN = "*"),
                        mz=df9$mz, rt=df9$rt, na=df9$na, scan=df9$scan)
  norm_facs <- target / colSums(df10[,grepl("Reporter.Ion.Intensity.",names(df10))], na.rm = T)
  mix10_sl <- data.frame(sweep(df10[,grepl("Reporter.Ion.Intensity.",names(df10))], 2, norm_facs, FUN = "*"),
                         mz=df10$mz, rt=df10$rt, na=df10$na, scan=df10$scan)
  norm_facs <- target / colSums(df11[,grepl("Reporter.Ion.Intensity.",names(df11))], na.rm = T)
  mix11_sl <- data.frame(sweep(df11[,grepl("Reporter.Ion.Intensity.",names(df11))], 2, norm_facs, FUN = "*"),
                         mz=df11$mz, rt=df11$rt, na=df11$na, scan=df11$scan)
  norm_facs <- target / colSums(df12[,grepl("Reporter.Ion.Intensity.",names(df12))], na.rm = T)
  mix12_sl <- data.frame(sweep(df12[,grepl("Reporter.Ion.Intensity.",names(df12))], 2, norm_facs, FUN = "*"),
                         mz=df12$mz, rt=df12$rt, na=df12$na, scan=df12$scan)
  norm_facs <- target / colSums(df13[,grepl("Reporter.Ion.Intensity.",names(df13))], na.rm = T)
  mix13_sl <- data.frame(sweep(df13[,grepl("Reporter.Ion.Intensity.",names(df13))], 2, norm_facs, FUN = "*"),
                         mz=df13$mz, rt=df13$rt, na=df13$na, scan=df13$scan)
  norm_facs <- target / colSums(df14[,grepl("Reporter.Ion.Intensity.",names(df14))], na.rm = T)
  mix14_sl <- data.frame(sweep(df14[,grepl("Reporter.Ion.Intensity.",names(df14))], 2, norm_facs, FUN = "*"),
                         mz=df14$mz, rt=df14$rt, na=df14$na, scan=df14$scan)
  norm_facs <- target / colSums(df15[,grepl("Reporter.Ion.Intensity.",names(df15))], na.rm = T)
  mix15_sl <- data.frame(sweep(df15[,grepl("Reporter.Ion.Intensity.",names(df15))], 2, norm_facs, FUN = "*"),
                         mz=df15$mz, rt=df15$rt, na=df15$na, scan=df15$scan)
  norm_facs <- target / colSums(df16[,grepl("Reporter.Ion.Intensity.",names(df16))], na.rm = T)
  mix16_sl <- data.frame(sweep(df16[,grepl("Reporter.Ion.Intensity.",names(df16))], 2, norm_facs, FUN = "*"),
                         mz=df16$mz, rt=df16$rt, na=df16$na, scan=df16$scan)
  
  join<-inner_join(mix1_sl, mix2_sl, by=c('mz', 'rt'))
  join<-join[order(joina$na.x),]; join<-join[!duplicated(join$scan.x),]
  join<-join[order(join$na.y),]; join<-join[!duplicated(join$scan.y),]
  join<-inner_join(join,mix3_sl, by=c('mz', 'rt'))
  join<-join[order(join$na.x),];join<-join[!duplicated(join$scan.x),]
  join1t<-inner_join(join, mix4_sl, by=c('mz', 'rt'))
  join1t<-join1t[order(join1t$na.x),];join1t<-join1t[!duplicated(join1t$scan.x),]
  
  join<-inner_join(mix5_sl, mix6_sl, by=c('mz', 'rt'))
  join<-join[order(join$na.x),]; join<-join[!duplicated(join$scan.x),]
  join<-join[order(join$na.y),]; join<-join[!duplicated(join$scan.y),]
  join<-inner_join(join,mix7_sl, by=c('mz', 'rt'))
  join<-join[order(join$na.x),];join<-join[!duplicated(join$scan.x),]
  join2t<-inner_join(join, mix8_sl, by=c('mz', 'rt'))
  join2t<-join2t[order(join2t$na.x),];join2t<-join2t[!duplicated(join2t$scan.x),]
  
  join<-inner_join(mix9_sl, mix10_sl, by=c('mz', 'rt'))
  join<-join[order(join$na.x),]; join<-join[!duplicated(join$scan.x),]
  join<-join[order(join$na.y),]; join<-join[!duplicated(join$scan.y),]
  join<-inner_join(join,mix11_sl, by=c('mz', 'rt'))
  join<-join[order(join$na.x),];join<-join[!duplicated(join$scan.x),]
  join3t<-inner_join(join, mix12_sl, by=c('mz', 'rt'))
  join3t<-join3t[order(join3t$na.x),];join3t<-join3t[!duplicated(join3t$scan.x),]
  
  join<-inner_join(mix13_sl, mix14_sl, by=c('mz', 'rt'))
  join<-join[order(join$na.x),]; join<-join[!duplicated(join$scan.x),]
  join<-join[order(join$na.y),]; join<-join[!duplicated(join$scan.y),]
  join<-inner_join(join,mix15_sl, by=c('mz', 'rt'))
  join<-join[order(join$na.x),];join<-join[!duplicated(join$scan.x),]
  join4t<-inner_join(join, mix16_sl, by=c('mz', 'rt'))
  join4t<-join4t[order(join4t$na.x),];join4t<-join4t[!duplicated(join4t$scan.x),]
  
  join<-inner_join(join1t, join2t, by=c('mz', 'rt'))
  join<-inner_join(join, join3t, by=c('mz', 'rt'))
  join<-inner_join(join, join4t, by=c('mz', 'rt'))
  }
##Normalizing Reporter Ion Intensity according to their loadings from DIA files
ndia1<-function(df1){
  target <- mean(c(colSums(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], na.rm = T)))
  norm_facs <- target / colSums(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], na.rm = T)
  mix1_sl <- data.frame(sweep(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], 2, norm_facs, FUN = "*"),
                        mz=df1$mz, rt=df1$rt)
}
ndia3<-function(df1, df2, df3){
  target <- mean(c(colSums(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], na.rm = T),
                   colSums(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], na.rm = T),
                   colSums(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], na.rm = T)))
  norm_facs <- target / colSums(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], na.rm = T)
  mix1_sl <- data.frame(sweep(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], 2, norm_facs, FUN = "*"),
                        mz=df1$mz, rt=df1$rt)
 norm_facs <- target / colSums(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], na.rm = T)
  mix2_sl <- data.frame(sweep(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], 2, norm_facs, FUN = "*"),
                        mz=df2$mz, rt=df2$rt)
 norm_facs <- target / colSums(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], na.rm = T)
  mix3_sl <- data.frame(sweep(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], 2, norm_facs, FUN = "*"),
                        mz=df3$mz, rt=df3$rt)
  join<-inner_join(mix1_sl, mix2_sl, by=c('mz', 'rt'))
  join<-inner_join(join, mix3_sl, by=c('mz', 'rt'))
  }
ndia4<-function(df1, df2, df3, df4){
  target <- mean(c(colSums(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], na.rm = T),
                   colSums(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], na.rm = T),
                   colSums(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], na.rm = T),
                   colSums(df4[,grepl("Reporter.Ion.Intensity.",names(df4))], na.rm = T)))
  norm_facs <- target / colSums(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], na.rm = T)
  mix1_sl <- data.frame(sweep(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], 2, norm_facs, FUN = "*"),
                        mz=df1$mz, rt=df1$rt)
  norm_facs <- target / colSums(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], na.rm = T)
  mix2_sl <- data.frame(sweep(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], 2, norm_facs, FUN = "*"),
                        mz=df2$mz, rt=df2$rt)
  norm_facs <- target / colSums(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], na.rm = T)
  mix3_sl <- data.frame(sweep(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], 2, norm_facs, FUN = "*"),
                        mz=df3$mz, rt=df3$rt)
  norm_facs <- target / colSums(df4[,grepl("Reporter.Ion.Intensity.",names(df4))], na.rm = T)
  mix4_sl <- data.frame(sweep(df4[,grepl("Reporter.Ion.Intensity.",names(df4))], 2, norm_facs, FUN = "*"),
                        mz=df4$mz, rt=df4$rt)
  
  join<-inner_join(mix1_sl, mix2_sl, by=c('mz', 'rt'))
  join<-inner_join(join, mix3_sl, by=c('mz', 'rt'))
  join<-inner_join(join, mix4_sl, by=c('mz', 'rt'))
  }
ndia12<-function(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12){
  target <- mean(c(colSums(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], na.rm = T),
                   colSums(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], na.rm = T),
                   colSums(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], na.rm = T),
                   colSums(df4[,grepl("Reporter.Ion.Intensity.",names(df4))], na.rm = T),
                   colSums(df5[,grepl("Reporter.Ion.Intensity.",names(df5))], na.rm = T),
                   colSums(df6[,grepl("Reporter.Ion.Intensity.",names(df6))], na.rm = T),
                   colSums(df7[,grepl("Reporter.Ion.Intensity.",names(df7))], na.rm = T),
                   colSums(df8[,grepl("Reporter.Ion.Intensity.",names(df8))], na.rm = T),
                   colSums(df9[,grepl("Reporter.Ion.Intensity.",names(df9))], na.rm = T),
                   colSums(df10[,grepl("Reporter.Ion.Intensity.",names(df10))], na.rm = T),
                   colSums(df11[,grepl("Reporter.Ion.Intensity.",names(df11))], na.rm = T),
                   colSums(df12[,grepl("Reporter.Ion.Intensity.",names(df12))], na.rm = T)))
  norm_facs <- target / colSums(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], na.rm = T)
  mix1_sl <- data.frame(sweep(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], 2, norm_facs, FUN = "*"),
                        mz=df1$mz, rt=df1$rt)
  norm_facs <- target / colSums(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], na.rm = T)
  mix2_sl <- data.frame(sweep(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], 2, norm_facs, FUN = "*"),
                        mz=df2$mz, rt=df2$rt)
  norm_facs <- target / colSums(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], na.rm = T)
  mix3_sl <- data.frame(sweep(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], 2, norm_facs, FUN = "*"),
                        mz=df3$mz, rt=df3$rt)
  norm_facs <- target / colSums(df4[,grepl("Reporter.Ion.Intensity.",names(df4))], na.rm = T)
  mix4_sl <- data.frame(sweep(df4[,grepl("Reporter.Ion.Intensity.",names(df4))], 2, norm_facs, FUN = "*"),
                        mz=df4$mz, rt=df4$rt)
  norm_facs <- target / colSums(df5[,grepl("Reporter.Ion.Intensity.",names(df5))], na.rm = T)
  mix5_sl <- data.frame(sweep(df5[,grepl("Reporter.Ion.Intensity.",names(df5))], 2, norm_facs, FUN = "*"),
                        mz=df5$mz, rt=df5$rt)
  norm_facs <- target / colSums(df6[,grepl("Reporter.Ion.Intensity.",names(df6))], na.rm = T)
  mix6_sl <- data.frame(sweep(df6[,grepl("Reporter.Ion.Intensity.",names(df6))], 2, norm_facs, FUN = "*"),
                        mz=df6$mz, rt=df6$rt)
  norm_facs <- target / colSums(df7[,grepl("Reporter.Ion.Intensity.",names(df7))], na.rm = T)
  mix7_sl <- data.frame(sweep(df7[,grepl("Reporter.Ion.Intensity.",names(df7))], 2, norm_facs, FUN = "*"),
                        mz=df7$mz, rt=df7$rt)
  norm_facs <- target / colSums(df8[,grepl("Reporter.Ion.Intensity.",names(df8))], na.rm = T)
  mix8_sl <- data.frame(sweep(df8[,grepl("Reporter.Ion.Intensity.",names(df8))], 2, norm_facs, FUN = "*"),
                        mz=df8$mz, rt=df8$rt)
  norm_facs <- target / colSums(df9[,grepl("Reporter.Ion.Intensity.",names(df9))], na.rm = T)
  mix9_sl <- data.frame(sweep(df9[,grepl("Reporter.Ion.Intensity.",names(df9))], 2, norm_facs, FUN = "*"),
                        mz=df9$mz, rt=df9$rt)
  norm_facs <- target / colSums(df10[,grepl("Reporter.Ion.Intensity.",names(df10))], na.rm = T)
  mix10_sl <- data.frame(sweep(df10[,grepl("Reporter.Ion.Intensity.",names(df10))], 2, norm_facs, FUN = "*"),
                         mz=df10$mz, rt=df10$rt)
  norm_facs <- target / colSums(df11[,grepl("Reporter.Ion.Intensity.",names(df11))], na.rm = T)
  mix11_sl <- data.frame(sweep(df11[,grepl("Reporter.Ion.Intensity.",names(df11))], 2, norm_facs, FUN = "*"),
                         mz=df11$mz, rt=df11$rt)
  norm_facs <- target / colSums(df12[,grepl("Reporter.Ion.Intensity.",names(df12))], na.rm = T)
  mix12_sl <- data.frame(sweep(df12[,grepl("Reporter.Ion.Intensity.",names(df12))], 2, norm_facs, FUN = "*"),
                         mz=df12$mz, rt=df12$rt)
  
  join<-inner_join(mix1_sl, mix2_sl, by=c('mz', 'rt'))
  join1<-inner_join(join, mix3_sl, by=c('mz', 'rt'))
  join<-inner_join(mix4_sl, mix5_sl, by=c('mz', 'rt'))
  join2<-inner_join(join, mix6_sl, by=c('mz', 'rt'))
  join<-inner_join(mix7_sl, mix8_sl, by=c('mz', 'rt'))
  join3<-inner_join(join, mix9_sl, by=c('mz', 'rt'))
  join<-inner_join(mix10_sl, mix11_sl, by=c('mz', 'rt'))
  join4<-inner_join(join, mix12_sl, by=c('mz', 'rt'))
  
  join<-inner_join(join1, join2, by=c('mz', 'rt'))
  join<-inner_join(join, join3, by=c('mz', 'rt'))
  join<-inner_join(join, join4, by=c('mz', 'rt'))
  }
ndia16<-function(df1, df2, df3, df4, df5, df6, df7, df8, df9, df10, df11, df12, df13, df14, df15, df16){
  target <- mean(c(colSums(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], na.rm = T),
                   colSums(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], na.rm = T),
                   colSums(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], na.rm = T),
                   colSums(df4[,grepl("Reporter.Ion.Intensity.",names(df4))], na.rm = T),
                   colSums(df5[,grepl("Reporter.Ion.Intensity.",names(df5))], na.rm = T),
                   colSums(df6[,grepl("Reporter.Ion.Intensity.",names(df6))], na.rm = T),
                   colSums(df7[,grepl("Reporter.Ion.Intensity.",names(df7))], na.rm = T),
                   colSums(df8[,grepl("Reporter.Ion.Intensity.",names(df8))], na.rm = T),
                   colSums(df9[,grepl("Reporter.Ion.Intensity.",names(df9))], na.rm = T),
                   colSums(df10[,grepl("Reporter.Ion.Intensity.",names(df10))], na.rm = T),
                   colSums(df11[,grepl("Reporter.Ion.Intensity.",names(df11))], na.rm = T),
                   colSums(df12[,grepl("Reporter.Ion.Intensity.",names(df12))], na.rm = T),
                   colSums(df13[,grepl("Reporter.Ion.Intensity.",names(df13))], na.rm = T),
                   colSums(df14[,grepl("Reporter.Ion.Intensity.",names(df14))], na.rm = T),
                   colSums(df15[,grepl("Reporter.Ion.Intensity.",names(df15))], na.rm = T),
                   colSums(df16[,grepl("Reporter.Ion.Intensity.",names(df16))], na.rm = T)))
  norm_facs <- target / colSums(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], na.rm = T)
  mix1_sl <- data.frame(sweep(df1[,grepl("Reporter.Ion.Intensity.",names(df1))], 2, norm_facs, FUN = "*"),
                        mz=df1$mz, rt=df1$rt)
  norm_facs <- target / colSums(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], na.rm = T)
  mix2_sl <- data.frame(sweep(df2[,grepl("Reporter.Ion.Intensity.",names(df2))], 2, norm_facs, FUN = "*"),
                        mz=df2$mz, rt=df2$rt)
  norm_facs <- target / colSums(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], na.rm = T)
  mix3_sl <- data.frame(sweep(df3[,grepl("Reporter.Ion.Intensity.",names(df3))], 2, norm_facs, FUN = "*"),
                        mz=df3$mz, rt=df3$rt)
  norm_facs <- target / colSums(df4[,grepl("Reporter.Ion.Intensity.",names(df4))], na.rm = T)
  mix4_sl <- data.frame(sweep(df4[,grepl("Reporter.Ion.Intensity.",names(df4))], 2, norm_facs, FUN = "*"),
                        mz=df4$mz, rt=df4$rt)
  norm_facs <- target / colSums(df5[,grepl("Reporter.Ion.Intensity.",names(df5))], na.rm = T)
  mix5_sl <- data.frame(sweep(df5[,grepl("Reporter.Ion.Intensity.",names(df5))], 2, norm_facs, FUN = "*"),
                        mz=df5$mz, rt=df5$rt)
  norm_facs <- target / colSums(df6[,grepl("Reporter.Ion.Intensity.",names(df6))], na.rm = T)
  mix6_sl <- data.frame(sweep(df6[,grepl("Reporter.Ion.Intensity.",names(df6))], 2, norm_facs, FUN = "*"),
                        mz=df6$mz, rt=df6$rt)
  norm_facs <- target / colSums(df7[,grepl("Reporter.Ion.Intensity.",names(df7))], na.rm = T)
  mix7_sl <- data.frame(sweep(df7[,grepl("Reporter.Ion.Intensity.",names(df7))], 2, norm_facs, FUN = "*"),
                        mz=df7$mz, rt=df7$rt)
  norm_facs <- target / colSums(df8[,grepl("Reporter.Ion.Intensity.",names(df8))], na.rm = T)
  mix8_sl <- data.frame(sweep(df8[,grepl("Reporter.Ion.Intensity.",names(df8))], 2, norm_facs, FUN = "*"),
                        mz=df8$mz, rt=df8$rt)
  norm_facs <- target / colSums(df9[,grepl("Reporter.Ion.Intensity.",names(df9))], na.rm = T)
  mix9_sl <- data.frame(sweep(df9[,grepl("Reporter.Ion.Intensity.",names(df9))], 2, norm_facs, FUN = "*"),
                        mz=df9$mz, rt=df9$rt)
  norm_facs <- target / colSums(df10[,grepl("Reporter.Ion.Intensity.",names(df10))], na.rm = T)
  mix10_sl <- data.frame(sweep(df10[,grepl("Reporter.Ion.Intensity.",names(df10))], 2, norm_facs, FUN = "*"),
                         mz=df10$mz, rt=df10$rt)
  norm_facs <- target / colSums(df11[,grepl("Reporter.Ion.Intensity.",names(df11))], na.rm = T)
  mix11_sl <- data.frame(sweep(df11[,grepl("Reporter.Ion.Intensity.",names(df11))], 2, norm_facs, FUN = "*"),
                         mz=df11$mz, rt=df11$rt)
  norm_facs <- target / colSums(df12[,grepl("Reporter.Ion.Intensity.",names(df12))], na.rm = T)
  mix12_sl <- data.frame(sweep(df12[,grepl("Reporter.Ion.Intensity.",names(df12))], 2, norm_facs, FUN = "*"),
                         mz=df12$mz, rt=df12$rt)
  norm_facs <- target / colSums(df13[,grepl("Reporter.Ion.Intensity.",names(df13))], na.rm = T)
  mix13_sl <- data.frame(sweep(df13[,grepl("Reporter.Ion.Intensity.",names(df13))], 2, norm_facs, FUN = "*"),
                         mz=df13$mz, rt=df13$rt)
  norm_facs <- target / colSums(df14[,grepl("Reporter.Ion.Intensity.",names(df14))], na.rm = T)
  mix14_sl <- data.frame(sweep(df14[,grepl("Reporter.Ion.Intensity.",names(df14))], 2, norm_facs, FUN = "*"),
                         mz=df14$mz, rt=df14$rt)
  norm_facs <- target / colSums(df15[,grepl("Reporter.Ion.Intensity.",names(df15))], na.rm = T)
  mix15_sl <- data.frame(sweep(df15[,grepl("Reporter.Ion.Intensity.",names(df15))], 2, norm_facs, FUN = "*"),
                         mz=df15$mz, rt=df15$rt)
  norm_facs <- target / colSums(df16[,grepl("Reporter.Ion.Intensity.",names(df16))], na.rm = T)
  mix16_sl <- data.frame(sweep(df16[,grepl("Reporter.Ion.Intensity.",names(df16))], 2, norm_facs, FUN = "*"),
                         mz=df16$mz, rt=df16$rt)
  join<-inner_join(mix1_sl, mix2_sl, by=c('mz', 'rt'))
  join<-inner_join(join, mix3_sl, by=c('mz', 'rt'))
  join1t<-inner_join(join, mix4_sl, by=c('mz', 'rt'))
  
  join<-inner_join(mix5_sl, mix6_sl, by=c('mz', 'rt'))
  join<-inner_join(join, mix7_sl, by=c('mz', 'rt'))
  join2t<-inner_join(join, mix8_sl, by=c('mz', 'rt'))
  
  join<-inner_join(mix9_sl, mix10_sl, by=c('mz', 'rt'))
  join<-inner_join(join, mix11_sl, by=c('mz', 'rt'))
  join3t<-inner_join(join, mix12_sl, by=c('mz', 'rt'))
  
  join<-inner_join(mix13_sl, mix14_sl, by=c('mz', 'rt'))
  join<-inner_join(join, mix15_sl, by=c('mz', 'rt'))
  join4t<-inner_join(join, mix16_sl, by=c('mz', 'rt'))
  
  join<-inner_join(join1t, join2t, by=c('mz', 'rt'))
  join<-inner_join(join, join3t, by=c('mz', 'rt'))
  join<-inner_join(join, join4t, by=c('mz', 'rt'))
  }
impute<-function(dfc) {
  dfc[is.na(dfc)] <- 10^rnorm(sum(is.na(dfc)), mean(log10(dfc)-1.8, na.rm=TRUE), 0.3)
  dfc
}

x<-"20200312_Exploris_RSLC9_Waters_CC_DDA_M1_0c5ng_QuanSpectra.txt";file<-"M1_0c5"; dda1p1<-quansdda(x)
x<-"20200312_Exploris_RSLC9_Waters_CC_DDA_M1_1ng_QuanSpectra.txt";file<-"M1_1"; dda1p2<-quansdda(x)
x<-"20200312_Exploris_RSLC9_Waters_CC_DDA_M1_5ng_QuanSpectra.txt";file<-"M1_5"; dda1p3<-quansdda(x)
x<-"20200312_Exploris_RSLC9_Waters_CC_DDA_M1_10ng_QuanSpectra.txt";file<-"M1_10"; dda1p4<-quansdda(x)

x<-"20200312_Exploris_RSLC9_Waters_CC_DDA_M2_0c5ng_QuanSpectra.txt";file<-"M2_0c5"; dda2p1<-quansdda(x)
x<-"20200311_Exploris_RSLC9_Waters_CC_DDA_M2_1ng_QuanSpectra.txt";file<-"M2_1"; dda2p2<-quansdda(x)
x<-"20200311_Exploris_RSLC9_Waters_CC_DDA_M2_5ng_QuanSpectra.txt";file<-"M2_5"; dda2p3<-quansdda(x)
x<-"20200312_Exploris_RSLC9_Waters_CC_DDA_M2_10ng_QuanSpectra.txt";file<-"M2_10"; dda2p4<-quansdda(x)

x<-"20200314_Exploris_RSLC9_Waters_CC_DDA_M3_0c5ng_QuanSpectra.txt";file<-"M3_0c5"; dda3p1<-quansdda(x)
x<-"20200316_Exploris_RSLC9_Waters_CC_DDA_M3_1ng_QuanSpectra.txt";file<-"M3_1"; dda3p2<-quansdda(x)
x<-"20200314_Exploris_RSLC9_Waters_CC_DDA_M3_5ng_QuanSpectra.txt";file<-"M3_5"; dda3p3<-quansdda(x)
x<-"20200314_Exploris_RSLC9_Waters_CC_DDA_M3_10ng_QuanSpectra.txt";file<-"M3_10"; dda3p4<-quansdda(x)

x<-"20200314_Exploris_RSLC9_Waters_CC_DDA_TKO_0c5ng_QuanSpectra.txt";file<-"TKO_0c5"; dda4p1<-quansdda(x)
x<-"20200314_Exploris_RSLC9_Waters_CC_DDA_TKO_1ng_QuanSpectra.txt";file<-"TKO_1"; dda4p2<-quansdda(x)
x<-"20200314_Exploris_RSLC9_Waters_CC_DDA_TKO_5ng-(1)_QuanSpectra.txt";file<-"TKO_5"; dda4p3<-quansdda(x)
x<-"20200314_Exploris_RSLC9_Waters_CC_DDA_TKO_10ng-(1)_QuanSpectra.txt";file<-"TKO_10"; dda4p4<-quansdda(x)

##Combine DDA data by their input
jdda1<-ndda3(dda1p1, dda2p1, dda3p1); jdda2<-ndda3(dda1p2, dda2p2, dda3p2); jdda3<-ndda3(dda1p3, dda2p3, dda3p3); jdda4<-ndda3(dda1p4, dda2p4, dda3p4)
##Combine DDA HeLa HEK data
jdda<-ndda12(dda1p1, dda2p1, dda3p1, dda1p2, dda2p2, dda3p2,dda1p3, dda2p3, dda3p3,dda1p4, dda2p4, dda3p4)
##Combine all DDA data
jddat<-ndda16(dda1p1, dda2p1, dda3p1, dda4p1, dda1p2, dda2p2, dda3p2, dda4p2, dda1p3, dda2p3, dda3p3, dda4p3,dda1p4, dda2p4, dda3p4, dda4p4)

x<-"20200312_Exploris_RSLC9_Waters_CC_DIA_M1_0c5ng_QuanSpectra.txt";file<-"M1_0c5"; dia1p1<-quansdia(x)
x<-"20200312_Exploris_RSLC9_Waters_CC_DIA_M1_1ng_QuanSpectra.txt";file<-"M1_1"; dia1p2<-quansdia(x)
x<-"20200312_Exploris_RSLC9_Waters_CC_DIA_M1_5ng_QuanSpectra.txt";file<-"M1_5"; dia1p3<-quansdia(x)
x<-"20200312_Exploris_RSLC9_Waters_CC_DIA_M1_10ng_QuanSpectra.txt";file<-"M1_10"; dia1p4<-quansdia(x)

x<-"20200312_Exploris_RSLC9_Waters_CC_DIA_M2_0c5ng_QuanSpectra.txt";file<-"M2_0c5"; dia2p1<-quansdia(x)
x<-"20200311_Exploris_RSLC9_Waters_CC_DIA_M2_1ng_QuanSpectra.txt";file<-"M2_1"; dia2p2<-quansdia_m2(x)
x<-"20200311_Exploris_RSLC9_Waters_CC_DIA_M2_5ng_QuanSpectra.txt";file<-"M2_5"; dia2p3<-quansdia_m2(x)
x<-"20200312_Exploris_RSLC9_Waters_CC_DIA_M2_10ng_QuanSpectra.txt";file<-"M2_10"; dia2p4<-quansdia(x)

x<-"20200314_Exploris_RSLC9_Waters_CC_DIA_M3_0c5ng_QuanSpectra.txt";file<-"M3_0c5"; dia3p1<-quansdia(x)
x<-"20200316_Exploris_RSLC9_Waters_CC_DIA_M3_1ng_QuanSpectra.txt";file<-"M3_1"; dia3p2<-quansdia(x)
x<-"20200314_Exploris_RSLC9_Waters_CC_DIA_M3_5ng_QuanSpectra.txt";file<-"M3_5"; dia3p3<-quansdia(x)
x<-"20200314_Exploris_RSLC9_Waters_CC_DIA_M3_10ng_QuanSpectra.txt";file<-"M3_10"; dia3p4<-quansdia(x)

x<-"20200314_Exploris_RSLC9_Waters_CC_DIA_TKO_0c5ng-(1)_QuanSpectra.txt";file<-"TKO_0c5"; dia4p1<-quansdia(x)
x<-"20200314_Exploris_RSLC9_Waters_CC_DIA_TKO_1ng-(1)_QuanSpectra.txt";file<-"TKO_1"; dia4p2<-quansdia(x)
x<-"20200314_Exploris_RSLC9_Waters_CC_DIA_TKO_5ng-(1)_QuanSpectra.txt";file<-"TKO_5"; dia4p3<-quansdia(x)
x<-"20200314_Exploris_RSLC9_Waters_CC_DIA_TKO_10ng-(1)_QuanSpectra.txt";file<-"TKO_10"; dia4p4<-quansdia(x)

##Combine DIA data by their input
jdia1<-ndia3(dia1p1, dia2p1, dia3p1); jdia2<-ndia3(dia1p2, dia2p2, dia3p2); jdia3<-ndia3(dia1p3, dia2p3, dia3p3); jdia4<-ndia3(dia1p4, dia2p4, dia3p4)
##Combine DIA HeLa HEK data
jdia<-ndia12(dia1p1, dia2p1, dia3p1,dia1p2, dia2p2, dia3p2,dia1p3, dia2p3, dia3p3,dia1p4, dia2p4, dia3p4)
##Combine all DIA data
jdiat<-ndia16(dia1p1, dia2p1, dia3p1, dia4p1,dia1p2, dia2p2, dia3p2, dia4p2,dia1p3, dia2p3, dia3p3, dia4p3,dia1p4, dia2p4, dia3p4, dia4p4)

###Spectromine Search
data<-read.table("20200403_223208_DDA_Search_Peptide precursor Report_20200403_231810.tsv", sep="\t", header=T, stringsAsFactors = F); data<-data.frame(data, ID=gsub(":.*", "\\1", data$R.FileName), stringsAsFactors = F); data[10:19]<-sapply(data[10:19], as.numeric)
sdda1p1<-data[grep("20200312_Exploris_RSLC9_Waters_CC_DDA_M1_0c5ng.raw", data$ID),]; sdda1p1 <- sdda1p1[ !duplicated(sdda1p1$PEP.StrippedSequence), ]; sdda1p1$RT<-paste('M1_0c5', sdda1p1$PP.EmpiricalRT, sep='_')
sdda1p2<-data[grep("20200312_Exploris_RSLC9_Waters_CC_DDA_M1_1ng.raw", data$ID),];sdda1p2 <- sdda1p2[ !duplicated(sdda1p2$PEP.StrippedSequence), ]; sdda1p2$RT<-paste('M1_1', sdda1p2$PP.EmpiricalRT, sep='_')
sdda1p3<-data[grep("20200312_Exploris_RSLC9_Waters_CC_DDA_M1_5ng.raw", data$ID),];sdda1p3 <- sdda1p3[ !duplicated(sdda1p3$PEP.StrippedSequence), ]; sdda1p3$RT<-paste('M1_5', sdda1p3$PP.EmpiricalRT, sep='_')
sdda1p4<-data[grep("20200312_Exploris_RSLC9_Waters_CC_DDA_M1_10ng.raw", data$ID),];sdda1p4 <- sdda1p4[ !duplicated(sdda1p4$PEP.StrippedSequence), ]; sdda1p4$RT<-paste('M1_10', sdda1p4$PP.EmpiricalRT, sep='_')

sdda2p1<-data[grep("20200312_Exploris_RSLC9_Waters_CC_DDA_M2_0c5ng.raw", data$ID),];sdda2p1 <- sdda2p1[ !duplicated(sdda2p1$PEP.StrippedSequence), ]; sdda2p1$RT<-paste('M2_0c5', sdda2p1$PP.EmpiricalRT, sep='_')
sdda2p2<-data[grep("20200311_Exploris_RSLC9_Waters_CC_DDA_M2_1ng.raw", data$ID),];sdda2p2 <- sdda2p2[ !duplicated(sdda2p2$PEP.StrippedSequence), ]; sdda2p2$RT<-paste('M2_1', sdda2p2$PP.EmpiricalRT, sep='_')
sdda2p3<-data[grep("20200311_Exploris_RSLC9_Waters_CC_DDA_M2_5ng.raw", data$ID),];sdda2p3 <- sdda2p3[ !duplicated(sdda2p3$PEP.StrippedSequence), ]; sdda2p3$RT<-paste('M2_5', sdda2p3$PP.EmpiricalRT, sep='_')
sdda2p4<-data[grep("20200312_Exploris_RSLC9_Waters_CC_DDA_M2_10ng.raw", data$ID),];sdda2p4 <- sdda2p4[ !duplicated(sdda2p4$PEP.StrippedSequence), ]; sdda2p4$RT<-paste('M2_10', sdda2p4$PP.EmpiricalRT, sep='_')

sdda3p1<-data[grep("20200314_Exploris_RSLC9_Waters_CC_DDA_M3_0c5ng.raw", data$ID),];sdda3p1 <- sdda3p1[ !duplicated(sdda3p1$PEP.StrippedSequence), ]; sdda3p1$RT<-paste('M3_0c5', sdda3p1$PP.EmpiricalRT, sep='_')
sdda3p2<-data[grep("20200316_Exploris_RSLC9_Waters_CC_DDA_M3_1ng.raw", data$ID),];sdda3p2 <- sdda3p2[ !duplicated(sdda3p2$PEP.StrippedSequence), ]; sdda3p2$RT<-paste('M3_1', sdda3p2$PP.EmpiricalRT, sep='_')
sdda3p3<-data[grep("20200314_Exploris_RSLC9_Waters_CC_DDA_M3_5ng.raw", data$ID),];sdda3p3 <- sdda3p3[ !duplicated(sdda3p3$PEP.StrippedSequence), ]; sdda3p3$RT<-paste('M3_5', sdda3p3$PP.EmpiricalRT, sep='_')
sdda3p4<-data[grep("20200314_Exploris_RSLC9_Waters_CC_DDA_M3_10ng.raw", data$ID),];sdda3p4 <- sdda3p4[ !duplicated(sdda3p4$PEP.StrippedSequence), ]; sdda3p4$RT<-paste('M3_10', sdda3p4$PP.EmpiricalRT, sep='_')

data<-read.table("D:/Claudia/DIA-test/Analysis - Spectronaut/20200416_181548_CC_DDA_TKO_Peptide precursor Report_20200416_182524.tsv", sep="\t", header=T, stringsAsFactors = F);data<-data.frame(data, ID=gsub(":.*", "\\1", data$R.FileName), stringsAsFactors = F);data[10:19]<-sapply(data[10:19], as.numeric)
sdda4p1<-data[grep("20200314_Exploris_RSLC9_Waters_CC_DDA_TKO_0c5ng", data$ID),]; sdda4p1 <- sdda4p1[ !duplicated(sdda4p1$PEP.StrippedSequence), ]; sdda4p1$RT<-paste('TKO_0c5', sdda4p1$PP.EmpiricalRT, sep='_')
sdda4p2<-data[grep("20200314_Exploris_RSLC9_Waters_CC_DDA_TKO_1ng", data$ID),]; sdda4p2 <- sdda4p2[ !duplicated(sdda4p2$PEP.StrippedSequence), ]; sdda4p2$RT<-paste('TKO_1', sdda4p2$PP.EmpiricalRT, sep='_')
sdda4p3<-data[grep("20200314_Exploris_RSLC9_Waters_CC_DDA_TKO_5ng", data$ID),]; sdda4p3 <- sdda4p3[ !duplicated(sdda4p3$PEP.StrippedSequence), ]; sdda4p3$RT<-paste('TKO_5', sdda4p3$PP.EmpiricalRT, sep='_')
sdda4p4<-data[grep("20200314_Exploris_RSLC9_Waters_CC_DDA_TKO_10ng", data$ID),]; sdda4p4 <- sdda4p4[ !duplicated(sdda4p4$PEP.StrippedSequence), ]; sdda4p4$RT<-paste('TKO_10', sdda4p4$PP.EmpiricalRT, sep='_')

#Impute and normalize ID-based data for DDA files with imputation
join<-full_join(sdda1p1, sdda2p1, by=c("PEP.StrippedSequence")); isdda1<-full_join(join, sdda3p1, by=c("PEP.StrippedSequence"))
isdda1<-data.frame(apply(isdda1[,grepl("PG.TMT",names(isdda1))], 2, impute), isdda1[,!grepl("PG.TMT",names(isdda1))])
isdda1t<-full_join(isdda1, sdda4p1, by=c("PEP.StrippedSequence")); isdda1t<-data.frame(apply(isdda1t[,grepl("PG.TMT",names(isdda1t))], 2, impute), isdda1t[,!grepl("PG.TMT",names(isdda1t))])
join<-full_join(sdda1p2, sdda2p2, by=c("PEP.StrippedSequence")); isdda2<-full_join(join, sdda3p2, by=c("PEP.StrippedSequence"))
isdda2<-data.frame(apply(isdda2[,grepl("PG.TMT",names(isdda2))], 2, impute), isdda2[,!grepl("PG.TMT",names(isdda2))])
isdda2t<-full_join(isdda2, sdda4p2, by=c("PEP.StrippedSequence")); isdda2t<-data.frame(apply(isdda2t[,grepl("PG.TMT",names(isdda2t))], 2, impute), isdda2t[,!grepl("PG.TMT",names(isdda2t))])
join<-full_join(sdda1p3, sdda2p3, by=c("PEP.StrippedSequence")); isdda3<-full_join(join, sdda3p3, by=c("PEP.StrippedSequence"))
isdda3<-data.frame(apply(isdda3[,grepl("PG.TMT",names(isdda3))], 2, impute), isdda3[,!grepl("PG.TMT",names(isdda3))])
isdda3t<-full_join(isdda3, sdda4p3, by=c("PEP.StrippedSequence")); isdda3t<-data.frame(apply(isdda3t[,grepl("PG.TMT",names(isdda3t))], 2, impute), isdda3t[,!grepl("PG.TMT",names(isdda3t))])
join<-full_join(sdda1p4, sdda2p4, by=c("PEP.StrippedSequence")); isdda4<-full_join(join, sdda3p4, by=c("PEP.StrippedSequence"))
isdda4<-data.frame(apply(isdda4[,grepl("PG.TMT",names(isdda4))], 2, impute), isdda4[,!grepl("PG.TMT",names(isdda4))])
isdda4t<-full_join(isdda4, sdda4p4, by=c("PEP.StrippedSequence")); isdda4t<-data.frame(apply(isdda4t[,grepl("PG.TMT",names(isdda4t))], 2, impute), isdda4t[,!grepl("PG.TMT",names(isdda4t))])
join<-full_join(isdda1, isdda2, by=c("PEP.StrippedSequence")); join<-full_join(join, isdda3, by=c("PEP.StrippedSequence")); isdda<-full_join(join, isdda4, by=c("PEP.StrippedSequence"))
isdd<-data.frame(apply(isdda[,grepl("PG.TMT",names(isdda))], 2, impute), isdda[,!grepl("PG.TMT",names(isdda))])
join<-full_join(isdda1t, isdda2t, by=c("PEP.StrippedSequence")); join<-full_join(join, isdda3t, by=c("PEP.StrippedSequence")); isddat<-full_join(join, isdda4t, by=c("PEP.StrippedSequence"))
isddat<-data.frame(apply(isddat[,grepl("PG.TMT",names(isddat))], 2, impute), isddat[,!grepl("PG.TMT",names(isddat))])

##Impute and normalize ID-based data for DDA files without imputation
join<-inner_join(sdda1p1, sdda2p1, by=c("PEP.StrippedSequence")); jsdda1<-inner_join(join, sdda3p1, by=c("PEP.StrippedSequence"))
jsdda1t<-inner_join(jsdda1, sdda4p1, by=c("PEP.StrippedSequence"))
join<-inner_join(sdda1p2, sdda2p2, by=c("PEP.StrippedSequence")); jsdda2<-inner_join(join, sdda3p2, by=c("PEP.StrippedSequence"))
jsdda2t<-inner_join(jsdda2, sdda4p2, by=c("PEP.StrippedSequence"))
join<-inner_join(sdda1p3, sdda2p3, by=c("PEP.StrippedSequence")); jsdda3<-inner_join(join, sdda3p3, by=c("PEP.StrippedSequence"))
jsdda3t<-inner_join(jsdda3, sdda4p3, by=c("PEP.StrippedSequence"))
join<-inner_join(sdda1p4, sdda2p4, by=c("PEP.StrippedSequence")); jsdda4<-inner_join(join, sdda3p4, by=c("PEP.StrippedSequence"))
jsdda4t<-inner_join(jsdda4, sdda4p4, by=c("PEP.StrippedSequence"))
join<-inner_join(jsdda1, jsdda2, by=c("PEP.StrippedSequence")); join6<-inner_join(join, jsdda3, by=c("PEP.StrippedSequence")); jsdda<-inner_join(join, jsdda3, by=c("PEP.StrippedSequence"))
join<-inner_join(jsdda1t, jsdda2t, by=c("PEP.StrippedSequence")); join6t<-inner_join(join, jsdda3t, by=c("PEP.StrippedSequence")); jsddat<-inner_join(join, jsdda4t, by=c("PEP.StrippedSequence"))

##Spectronaut search
data<-read.table("DIA_PiPro_0c5to10ng_M1to3_10ngLib.tsv", sep="\t", header=T); data<-data.frame(data, ID=gsub(":.*", "\\1", data$R.FileName), stringsAsFactors = F)
sdia1p1<-data[grep("20200312_Exploris_RSLC9_Waters_CC_DIA_M1_0c5ng", data$ID),]; sdia1p1 <- sdia1p1[ !duplicated(sdia1p1$PEP.StrippedSequence), ]; sdia1p1$RT<-paste('M1_0c5', sdia1p1$PP.EmpiricalRT, sep='_')
sdia1p2<-data[grep("20200312_Exploris_RSLC9_Waters_CC_DIA_M1_1ng", data$ID),]; sdia1p2 <- sdia1p2[ !duplicated(sdia1p2$PEP.StrippedSequence), ]; sdia1p2$RT<-paste('M1_1', sdia1p2$PP.EmpiricalRT, sep='_')
sdia1p3<-data[grep("20200312_Exploris_RSLC9_Waters_CC_DIA_M1_5ng", data$ID),]; sdia1p3 <- sdia1p3[ !duplicated(sdia1p3$PEP.StrippedSequence), ]; sdia1p3$RT<-paste('M1_5', sdia1p3$PP.EmpiricalRT, sep='_')
sdia1p4<-data[grep("20200312_Exploris_RSLC9_Waters_CC_DIA_M1_10ng", data$ID),]; sdia1p4 <- sdia1p4[ !duplicated(sdia1p4$PEP.StrippedSequence), ]; sdia1p4$RT<-paste('M1_10', sdia1p4$PP.EmpiricalRT, sep='_')

sdia2p1<-data[grep("20200312_Exploris_RSLC9_Waters_CC_DIA_M2_0c5ng", data$ID),]; sdia2p1 <- sdia2p1[ !duplicated(sdia2p1$PEP.StrippedSequence), ]; sdia2p1$RT<-paste('M2_0c5', sdia2p1$PP.EmpiricalRT, sep='_')
sdia2p2<-data[grep("20200311_Exploris_RSLC9_Waters_CC_DIA_M2_1ng", data$ID),]; sdia2p2 <- sdia2p2[ !duplicated(sdia2p2$PEP.StrippedSequence), ]; sdia2p2$RT<-paste('M2_1', sdia2p2$PP.EmpiricalRT, sep='_')
sdia2p3<-data[grep("20200311_Exploris_RSLC9_Waters_CC_DIA_M2_5ng", data$ID),]; sdia2p3 <- sdia2p3[ !duplicated(sdia2p3$PEP.StrippedSequence), ]; sdia2p3$RT<-paste('M2_5', sdia2p3$PP.EmpiricalRT, sep='_')
sdia2p4<-data[grep("20200312_Exploris_RSLC9_Waters_CC_DIA_M2_10ng", data$ID),]; sdia2p4 <- sdia2p4[ !duplicated(sdia2p4$PEP.StrippedSequence), ]; sdia2p4$RT<-paste('M2_10', sdia2p4$PP.EmpiricalRT, sep='_')

sdia3p1<-data[grep("20200314_Exploris_RSLC9_Waters_CC_DIA_M3_0c5ng", data$ID),]; sdia3p1 <- sdia3p1[ !duplicated(sdia3p1$PEP.StrippedSequence), ]; sdia3p1$RT<-paste('M3_0c5', sdia3p1$PP.EmpiricalRT, sep='_')
sdia3p2<-data[grep("20200316_Exploris_RSLC9_Waters_CC_DIA_M3_1ng", data$ID),]; sdia3p2 <- sdia3p2[ !duplicated(sdia3p2$PEP.StrippedSequence), ]; sdia3p1$RT<-paste('M3_0c5', sdia3p1$PP.EmpiricalRT, sep='_')
sdia3p3<-data[grep("20200314_Exploris_RSLC9_Waters_CC_DIA_M3_5ng", data$ID),]; sdia3p3 <- sdia3p3[ !duplicated(sdia3p3$PEP.StrippedSequence), ]; sdia3p1$RT<-paste('M3_0c5', sdia3p1$PP.EmpiricalRT, sep='_')
sdia3p4<-data[grep("20200314_Exploris_RSLC9_Waters_CC_DIA_M3_10ng", data$ID),]; sdia3p4 <- sdia3p4[ !duplicated(sdia3p4$PEP.StrippedSequence), ]; sdia3p1$RT<-paste('M3_0c5', sdia3p1$PP.EmpiricalRT, sep='_')

data<-read.table("20200416_161110_20200314_CC_TKO_10ngLib_Report.tsv", sep="\t", header=T, stringsAsFactors = F);data<-data.frame(data, ID=gsub(":.*", "\\1", data$R.FileName), stringsAsFactors = F)
sdia4p1<-data[grep("20200314_Exploris_RSLC9_Waters_CC_DIA_TKO_0c5ng", data$ID),]; sdia4p1 <- sdia4p1[ !duplicated(sdia4p1$PEP.StrippedSequence), ]; sdia4p1$RT<-paste('TKO_0c5', sdia4p1$PP.EmpiricalRT, sep='_')
sdia4p2<-data[grep("20200314_Exploris_RSLC9_Waters_CC_DIA_TKO_1ng", data$ID),]; sdia4p2 <- sdia4p2[ !duplicated(sdia4p2$PEP.StrippedSequence), ]; sdia4p2$RT<-paste('TKO_1', sdia4p2$PP.EmpiricalRT, sep='_')
sdia4p3<-data[grep("20200314_Exploris_RSLC9_Waters_CC_DIA_TKO_5ng", data$ID),]; sdia4p3 <- sdia4p3[ !duplicated(sdia4p3$PEP.StrippedSequence), ]; sdia4p3$RT<-paste('TKO_5', sdia4p3$PP.EmpiricalRT, sep='_')
sdia4p4<-data[grep("20200314_Exploris_RSLC9_Waters_CC_DIA_TKO_10ng", data$ID),]; sdia4p4 <- sdia4p4[ !duplicated(sdia4p4$PEP.StrippedSequence), ]; sdia4p4$RT<-paste('TKO_10', sdia4p4$PP.EmpiricalRT, sep='_')
#Impute and normalize ID-based data for DIA files with imputation
join<-full_join(sdia1p1, sdia2p1, by=c("PEP.StrippedSequence")); isdia1<-full_join(join, sdia3p1, by=c("PEP.StrippedSequence"))
isdia1<-data.frame(apply(isdia1[,grepl("RI",names(isdia1))], 2, impute), isdia1[,!grepl("RI",names(isdia1))])
isdia1t<-full_join(isdia1, sdia4p1, by=c("PEP.StrippedSequence"))
isdia1t<-data.frame(apply(isdia1t[,grepl("RI",names(isdia1t))], 2, impute), isdia1t[,!grepl("RI",names(isdia1t))])
join<-full_join(sdia1p2, sdia2p2, by=c("PEP.StrippedSequence")); isdia2<-full_join(join, sdia3p2, by=c("PEP.StrippedSequence"))
isdia2<-data.frame(apply(isdia2[,grepl("RI",names(isdia2))], 2, impute), isdia2[,!grepl("RI",names(isdia2))])
isdia2t<-full_join(isdia2, sdia4p2, by=c("PEP.StrippedSequence"))
isdia2t<-data.frame(apply(isdia2t[,grepl("RI",names(isdia2t))], 2, impute), isdia2t[,!grepl("RI",names(isdia2t))])
join<-full_join(sdia1p3, sdia2p3, by=c("PEP.StrippedSequence")); isdia3<-full_join(join, sdia3p3, by=c("PEP.StrippedSequence"))
isdia3<-data.frame(apply(isdia3[,grepl("RI",names(isdia3))], 2, impute), isdia3[,!grepl("RI",names(isdia3))])
isdia3t<-full_join(isdia3, sdia4p3, by=c("PEP.StrippedSequence"))
isdia3t<-data.frame(apply(isdia3t[,grepl("RI",names(isdia3t))], 2, impute), isdia3t[,!grepl("RI",names(isdia3t))])
join<-full_join(sdia1p4, sdia2p4, by=c("PEP.StrippedSequence")); isdia4<-full_join(join, sdia3p4, by=c("PEP.StrippedSequence"))
isdia4<-data.frame(apply(isdia4[,grepl("RI",names(isdia4))], 2, impute), isdia4[,!grepl("RI",names(isdia4))])
isdia4t<-full_join(isdia4, sdia4p4, by=c("PEP.StrippedSequence"))
isdia4t<-data.frame(apply(isdia4t[,grepl("RI",names(isdia4t))], 2, impute), isdia4t[,!grepl("RI",names(isdia4t))])
join<-full_join(isdia1, isdia2, by=c("PEP.StrippedSequence")); join<-full_join(join, isdia3, by=c("PEP.StrippedSequence")); isdia<-full_join(join, isdia4, by=c("PEP.StrippedSequence"))
isdia<-data.frame(apply(isdia[,grepl("RI",names(isdia))], 2, impute), isdia[,!grepl("RI",names(isdia))])
join<-full_join(isdia1t, isdia2t, by=c("PEP.StrippedSequence")); join<-full_join(join, isdia3t, by=c("PEP.StrippedSequence")); isdiat<-full_join(join, isdia4t, by=c("PEP.StrippedSequence"))
isdiat<-data.frame(apply(isdiat[,grepl("RI",names(isdiat))], 2, impute), isdiat[,!grepl("RI",names(isdiat))])
#Impute and normalize ID-based data for DIA files without imputation
join<-inner_join(sdia1p1, sdia2p1, by=c("PEP.StrippedSequence"));jsdia1<-inner_join(join, sdia3p1, by=c("PEP.StrippedSequence"))
jsdia1t<-inner_join(jsdia1, sdia4p1, by=c("PEP.StrippedSequence"))
join<-inner_join(sdia1p2, sdia2p2, by=c("PEP.StrippedSequence"));jsdia2<-inner_join(join, sdia3p2, by=c("PEP.StrippedSequence"))
jsdia2t<-inner_join(jsdia2, sdia4p2, by=c("PEP.StrippedSequence"))
join<-inner_join(sdia1p3, sdia2p3, by=c("PEP.StrippedSequence"));jsdia3<-inner_join(join, sdia3p3, by=c("PEP.StrippedSequence"))
jsdia3t<-inner_join(jsdia3, sdia4p3, by=c("PEP.StrippedSequence"))
join<-inner_join(sdia1p4, sdia2p4, by=c("PEP.StrippedSequence"));jsdia4<-inner_join(join, sdia3p4, by=c("PEP.StrippedSequence"))
jsdia4t<-inner_join(jsdia4, sdia4p4, by=c("PEP.StrippedSequence"))
join<-inner_join(jsdia1, jsdia2, by=c("PEP.StrippedSequence"));join<-inner_join(join, jsdia3, by=c("PEP.StrippedSequence"));jsdia<-inner_join(join, jsdia4, by=c("PEP.StrippedSequence"))
join<-inner_join(jsdia1t, jsdia2t, by=c("PEP.StrippedSequence"));join6t<-inner_join(join, jsdia3t, by=c("PEP.StrippedSequence")); jsdiat<-inner_join(join, jsdia4t, by=c("PEP.StrippedSequence"))
