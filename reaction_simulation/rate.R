require(gdata)
# intrinsic rate
ksurf=c("k_1kD", "k_10kD", "k_100kD", "k_inf")
# your output directory
root="OUTPUT_DIRECTORY"
OUTPUT_DIRECTORY="OUTPUT_DIRECTORY"
# the adsorption interactions between species and the nanoparticle
ep = c("0.1", "1.0", "2.0", "3.0", "4.0", "5.0", "6.0", "7.0", "8.0", "9.0", "10.0")
# initial number of A or B
N = c(500, 1000, 2000)


kd = 34.23855
rate_temp = data.frame(entry = seq(1,10,1), rate = 0)
ad_temp = data.frame(entry = seq(1,10,1), A_mean = 0, B_mean = 0)
rate = data.frame(N = N, rate = 0, se = 0, NANB = 0, se_NANB = 0, norm_rate = 0, se_norm_rate = 0)

e_count = 1

for(n in N){
directory = paste(root, "/N_", n, sep = "")

file = paste(directory, "/particle_num_k_inf_10.0.txt", sep = "")
rxn_event = read.table(file, header = FALSE)
start = length(rxn_event$V3) / 10 * 2
rt = data.frame(t = seq(1,length(rxn_event$V3),1), num = rxn_event$V3)
l = (length(rxn_event$V3) - start) %/% 10
for(count in seq(1,10,1)){
head = start + 1 + (count - 1) * l
tail = start + count * l
lf = summary(lm(rt$num[head:tail] ~ rt$t[head:tail]))
rate_temp$rate[count] = lf$coefficients[2] * 10^4
}
rate$rate[e_count] = mean(rate_temp$rate) / (kd * n)
rate$se[e_count] = sd(rate_temp$rate) / sqrt(10) / (kd * n)
rate$NANB[e_count] = mean(rxn_event$V1 * rxn_event$V2)
for(i in c(1:10)){
head = (i - 1) * (length(rxn_event$V1) / 10) + 1
tail = i * length(rxn_event$V1) / 10
ad_temp$A_mean[i] = mean(rxn_event$V1[head:tail])
ad_temp$B_mean[i] = mean(rxn_event$V2[head:tail])
}
rate$se_NANB[e_count] = mean(ad_temp$A_mean) * mean(ad_temp$B_mean) * sqrt(sd(ad_temp$A_mean)^2 / mean(ad_temp$A_mean)^2 + sd(ad_temp$B_mean)^2 / mean(ad_temp$B_mean)^2) / sqrt(10)

e_count = e_count + 1
}
write.fwf(rate, paste("OUTPUT_DIRECTORY", sep = ""), sep = " ")

