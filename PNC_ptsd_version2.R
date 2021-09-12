#dim: 973 by 933
data <- read.csv("phs_brain_img.csv")
#compute height based on existing variables
height<-data$MED001*12+data$MED002
#compute age based on existing variables
age<-2009-data$Med_birth_year
#compute bmi based on existing variables
bmi<-data$MED003 * 703 / (height*height)
data<-cbind(data, age, bmi)
dim(data)
#[1] 973 935
#select participants with ptsd 
ind1<-which(data$smry_ptd==4)
data_ptsd<-data[ind1, ]
dim(data_ptsd)
#[1] 104 935
#exclude participants with medical ratings higher than 1
ind2<-which(data_ptsd$Med_Rating<=1)
data_ptsd2<-data_ptsd[ind2, ]
dim(data_ptsd2)
#[1]  79 935
#Note that the youngest one is 8 years old, the highest bmi is 40, for many participants the bmi value is NA
# should we exclude obese participants? should we exclude young kids?
cbind(data_ptsd2$Sex, data_ptsd2$age, data_ptsd2$bmi, data_ptsd2$Med_Rating)
#all summary variables on psycopathology disorders
sum<-rowSums(data[, 4:30])
length(sum)
#only 20 subjects have all zeros
length(which(sum==0))
#take all these psychiatric summary variables
subdata<-data[, 4:30]
indicator<-rep(0, nrow(subdata) )
indicator<-matrix(indicator, nrow(subdata), 1)
for (i in 1: 27){
	#see whether the score is >=2
	ind<-which(subdata[, i]>1)
	indicator[ind]<-indicator[ind]+1
}
# we get 104 subjects with all phychiatric summary variables <=1
length(which(indicator==0))
health_ind<-which(indicator==0)
health_data<- data[health_ind, ]
dim(health_data)
#among these 104 subjects, 68 of them have Med_rating <=1
length(which(health_data$Med_Rating<=1))
#[1]  68 935
health_group<-health_data[which(health_data$Med_Rating<=1),   ] 

#select participants with ptsd summary variable==1 and other major disorders <=2 including Anxiety, Mood, #Behavior, Psychosis, Eating 
ind1<-which(data$smry_ptd==1 & data$smry_anx_cat<=2 & data$smry_mood_cat<=2 & data$smry_eat_cat<=2 & data$smry_psy_cat<=2 & data$smry_beh_cat<=2 )
data_ptsd1<-data[ind1, ]
dim(data_ptsd1)
#[1] 83 935
#exclude participants with medical ratings higher than 1
ind2<-which(data_ptsd1$Med_Rating<=1)
data_ptsd11<-data_ptsd1[ind2, ]
dim(data_ptsd11)
#59 935

#combine the ptsd group and the healthy group, data_ptsd2, health_group

group_ind<-c(  rep(0, nrow(health_group)),rep(1, nrow(data_ptsd11) ),rep(2, nrow(data_ptsd2) ) )
group_ind<-as.factor(group_ind)
# combine data from the two groups
data_com<- rbind(health_group,data_ptsd11, data_ptsd2)
dim(data_com)
#[1] 206 935
data_com<- cbind(data_com, group_ind)
dim(data_com)
#[1] 206 935
table(data_com$group_ind)
# 0   1     2 
#68  59   79

fit<-aov(data_com$bmi ~ data_com$group_ind)
summary(fit) # p value 0.158
fit<-aov(data_com$age ~ data_com$group_ind)
summary(fit) # p value 0.164
library(MASS)
tbl<- table(data_com$Sex, data_com$group_ind)
tbl
#      0   1    2
#  F 33 15  56
# M 35  44  23
chisq.test(tbl) #X-squared = 28.081, df = 2, p-value = 7.985e-07

write.table(data_com, file="PTSD_TraumaWithoutSympton_Control.csv", row.names=FALSE, sep=",")
#instead of matching, we could just use Sex as a covariate in later analyses
#install.packages("MatchIt")
library(MatchIt)
data_com_out = matchit(group_ind ~ Sex, 
                   data = data_com, method = "nearest",
                   ratio = 1)

summary(data_com_out)
data_com_matched <- match.data(data_com_out)










