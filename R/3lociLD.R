

sa= 0.25
sb=0.7

pa=0.5
pb=0.5

Dab=0.05

ld3loci<-function(sa=0.2,sb=0.2,sc=0,pa=0.5,pb=0.5,pc=0,Dab=0.05,Dac=0,Dbc=0,
                whatreturn='Rchange' ){

genos<-expand.grid(c('A','a'),c('B','b'),c('C','c'))
genos$g<-paste(genos$Var1,genos$Var2,genos$Var3,sep='')
genos<-data.frame(genos)

genos$x<-1
genos$w<-1

genos[which(genos$g =='ABC'),'w']<- (1+sa)*(1+sb)*(1+sc)
genos[which(genos$g =='aBC'),'w']<- (1+sb)*(1+sc)
genos[which(genos$g =='AbC'),'w']<- (1+sa)*(1+sc)
genos[which(genos$g =='abC'),'w']<- (1+sc)
genos[which(genos$g =='ABc'),'w']<- (1+sa)*(1+sb)
genos[which(genos$g =='aBc'),'w']<- (1+sb)
genos[which(genos$g =='Abc'),'w']<- (1+sa)
genos[which(genos$g =='abc'),'w']<- 1


Dab_ori=Dab
Dac_ori=Dac
Dbc_ori=Dbc

Dab=Dab/2
Dac=Dac/2
Dbc=Dbc/2

genos[which(genos$g =='ABC'),'x']<- pa*pb*pc *(1+ Dab+Dac+Dbc)
genos[which(genos$g =='aBC'),'x']<- (1-pa)*pb*pc * (1-Dab-Dac+Dbc)
genos[which(genos$g =='AbC'),'x']<- pa*(1-pb)*pc  * (1- Dab+Dac-Dbc)
genos[which(genos$g =='abC'),'x']<- (1-pa)*(1-pb)*pc * (1+ Dab-Dac-Dbc)
genos[which(genos$g =='ABc'),'x']<- pa*pb*(1-pc)  * (1+ Dab-Dac-Dbc)
genos[which(genos$g =='aBc'),'x']<- (1-pa)*pb*(1-pc) * (1-Dab+Dac-Dbc)
genos[which(genos$g =='Abc'),'x']<- pa*(1-pb)*(1-pc) * (1-Dab-Dac+Dbc)
genos[which(genos$g =='abc'),'x']<- (1-pa)*(1-pb)*(1-pc) * (1+Dab+Dac+Dbc)
genos

if(!all(genos$x<0)){
print(genos)
genos$x[genos$x<0]<-0
genos$x<-genos$x/sum(genos$x) # two tricks
}

genos$xprime<-(genos$w *genos$x)
genos$xprime<-genos$xprime/ sum(genos$xprime)

if(!sum(genos$xprime)==1){
# stopifnot(sum(genos$xprime)==1)
print(genos)
genos$xprime<-genos$xprime/ sum(genos$xprime)
}

pa1= sum(dplyr::filter(genos,Var1 =='A')$xprime)
pb1 =sum(dplyr::filter(genos,Var2 =='B')$xprime)
pc1 =sum(dplyr::filter(genos,Var3 =='C')$xprime)


Dac_infer<-(
  sum(dplyr::filter(genos,Var1 =='A' & Var2 =='C')$x) * sum(dplyr::filter(genos,Var1 =='a' & Var2 =='c')$x)
)-(
  sum(dplyr::filter(genos,Var1 =='A' & Var2 =='c')$x) * sum(dplyr::filter(genos,Var1 =='a' & Var2 =='C')$x)
)
Dbc_infer<-(
  sum(dplyr::filter(genos,Var1 =='C' & Var2 =='B')$x) * sum(dplyr::filter(genos,Var1 =='c' & Var2 =='b')$x)
)-(
  sum(dplyr::filter(genos,Var1 =='C' & Var2 =='b')$x) * sum(dplyr::filter(genos,Var1 =='c' & Var2 =='B')$x)
)
Dab_infer<-(
  sum(dplyr::filter(genos,Var1 =='A' & Var2 =='B')$x) * sum(dplyr::filter(genos,Var1 =='a' & Var2 =='b')$x)
)-(
  sum(dplyr::filter(genos,Var1 =='A' & Var2 =='b')$x) * sum(dplyr::filter(genos,Var1 =='a' & Var2 =='B')$x)
)

if(Dac_infer != Dac_ori | Dab_infer != Dab_ori | Dbc_infer != Dbc_ori){
  message('The input and D values do not coincide!')
  message(Dab_ori, " =? ", Dab_infer)
  message(Dac_ori, " =? ", Dac_infer)
  message(Dbc_ori, " =? ", Dbc_infer)
  message('Results might not hold')
}

# (
#   sum(dplyr::filter(genos,Var1 =='A' & Var2 =='B')$xprime) * sum(dplyr::filter(genos,Var1 =='a' & Var2 =='b')$xprime)
# )-(
#   sum(dplyr::filter(genos,Var1 =='A' & Var2 =='b')$xprime) * sum(dplyr::filter(genos,Var1 =='a' & Var2 =='B')$xprime)
# )
# (
#   sum(dplyr::filter(genos,Var1 =='A' & Var2 =='B')$xprime)
# )-(
#   pa1 * pb1
# )


if(whatreturn=='Rchange')
(
  sum(dplyr::filter(genos,Var1 =='A' & Var2 =='B')$x) * sum(dplyr::filter(genos,Var1 =='a' & Var2 =='b')$x)
)/(
  sum(dplyr::filter(genos,Var1 =='A' & Var2 =='b')$x) * sum(dplyr::filter(genos,Var1 =='a' & Var2 =='B')$x)
) %>% message('g = 0 -> R = ',.)
(
  sum(dplyr::filter(genos,Var1 =='A' & Var2 =='B')$xprime) * sum(dplyr::filter(genos,Var1 =='a' & Var2 =='b')$xprime)
)/(
  sum(dplyr::filter(genos,Var1 =='A' & Var2 =='b')$xprime) * sum(dplyr::filter(genos,Var1 =='a' & Var2 =='B')$xprime)
)%>% message('g = 1 -> R = ',.)

}



# # r2 change
# ((
#   (( sum(dplyr::filter(genos,Var1 =='A' & Var2 =='B')$xprime) * sum(dplyr::filter(genos,Var1 =='a' & Var2 =='b')$xprime))-( sum(dplyr::filter(genos,Var1 =='A' & Var2 =='b')$xprime) * sum(dplyr::filter(genos,Var1 =='a' & Var2 =='B')$xprime))
#   )/(sqrt(pa*pb*(1-pa)*(1-pb)))
# ))-((
#   (( sum(dplyr::filter(genos,Var1 =='A' & Var2 =='B')$x) * sum(dplyr::filter(genos,Var1 =='a' & Var2 =='b')$x))-( sum(dplyr::filter(genos,Var1 =='A' & Var2 =='b')$x) * sum(dplyr::filter(genos,Var1 =='a' & Var2 =='B')$x)
#   ))/(sqrt(pa1*pb1*(1-pa1)*(1-pb1)))
# ))
