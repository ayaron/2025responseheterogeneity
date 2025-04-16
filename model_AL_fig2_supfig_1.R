library(tidyverse)
library(ggridges)
library(pracma)


myligands = colorRampPalette(c("midnightblue","skyblue3","green", "green4")) # ligand colors
# model AL
model_AL = function(A0 = 1, kmin = 1, kmax = 1, kres = 1, cmin = -3, cmax = 3, e = 1, cres = 101){ #for a concentration or affinity of one, input zero to the min and max 
  C0 = logspace(cmin,cmax,cres) #the ligand concentration in log space
  K = logspace(kmin,kmax,kres) #the ligand affinity in log space
  F = matrix(nrow = kres, ncol = cres) # the amount of the full complex, F
  S = matrix(nrow = kres, ncol = cres) # the sensitivity to teh receptor, S
  for (c in seq_along(C0)) {
    for (k in seq_along(K)) {
      F[k, c] =  K[k]*C0[c]*A0/(1+K[k]*C0[c])
      S[k, c] = 1 # instead of writing the whole calculation, see supplementary data
    } 
   
  }

  dynamics = tibble(expand_grid(C0 = C0, K = K), F = c(F), S = c(S), E = c(e*F))
  dynamics
}

#a data frame holding the data for figures 2A,B and sup fig 1A
AL = tibble(F = model_AL(kmin = -3, kmax = 3, kres = 5, cmin = -5, cmax = 5, cres = 101 )$F,
            S = model_AL(kmin = -3, kmax = 3, kres = 5, cmin = -5, cmax = 5, cres = 101 )$S,
            C0 = model_AL(kmin = -3, kmax = 3, kres = 5, cmin = -5, cmax = 5, cres = 101 )$C0,
            K = model_AL(kmin = -3, kmax = 3, kres = 5, cmin = -5, cmax = 5, cres = 101 )$K)

#figure 2B left
AL %>% 
  mutate(labs = paste0("10^", log10(K))) %>% 
  ggplot(aes(x = K, y = as.factor(K), color = as.factor(K))) +
  geom_point(pch = 4, size = 2)+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), limits = c(10^-3,10^3),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  theme_bw(base_size = 20,base_rect_size = 2,
           base_line_size = 0)+
  annotation_logticks(sides = "b",outside = T)+
  coord_cartesian(clip  = "off") +
  scale_color_manual(name = "K",values = myligands(length(unique(AL$K))))+
  labs(title = "Ligands", x = "Affinity (KL)", y = "") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.ticks = element_line(size = 1),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank(),
        legend.position = "none") 

#supplementary figure 1A
AL %>% 
  ggplot(aes(C0)) +
  geom_point(aes(y = F), color = "black", size = 2) +
  geom_point(aes(y = F, color = as.factor(K)), size = 1.5) +
  theme_bw(base_size = 30,base_rect_size = 2,
           base_line_size = 0)+
  annotation_logticks(size = 1, sides = "b")+
  scale_color_manual(values = myligands(5))+
  scale_x_log10(name = "C0L",breaks = scales::trans_breaks("log10", function(x) 10^x), limits = c(10^-5,10^5),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_continuous("FL") +
  coord_fixed(5)+
  theme(legend.position = "none")

#figure 2E
AL %>% 
  ggplot(aes(C0)) +
  geom_point(aes(y = S), color = "black", size = 2) +
  geom_point(aes(y = S, color = as.factor(K)), size = 1.5) +
  theme_bw(base_size = 30,base_rect_size = 2,
           base_line_size = 0)+
  annotation_logticks(size = 1, sides = "b")+
  scale_color_manual(values = myligands(5))+
  scale_x_log10(name = "C0L",breaks = scales::trans_breaks("log10", function(x) 10^x), limits = c(10^-5,10^5),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  scale_y_continuous("S", limits = c(0,2)) +
  coord_fixed(5)+
  theme(legend.position = "none")

set.seed(100000)  
Fdist =  tibble(a0s = rgamma(n = 100000, shape = 1/(0.5^2), scale = 0.5^2), 
                f1 = NA, f2 = NA, f3 = NA, f4 = NA, f5 = NA)

#figure 2B right
Fdist %>% 
  ggplot(aes(a0s)) +
  geom_histogram(bins = 100,color = "#FFF4C7", size = 2.5) +
  geom_histogram(bins = 100,fill = "#FFD400") +
  theme_bw(base_size = 30,base_rect_size = 2,
           base_line_size = 0)+labs(x = "log10(A0) distribution") +
  scale_x_log10(name = "receptors (A0)",
                breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) +
  annotation_logticks(side = "b", size = 1)


for (a0 in seq_along(Fdist$a0s)) {
  Fdist$f1[a0] = model_AL(A0 = Fdist$a0s[[a0]], cres = 1, kres = 5, kmin = -3, kmax = 3, cmin = 0, cmax = 0)$F[1]
  Fdist$f2[a0] = model_AL(A0 = Fdist$a0s[[a0]], cres = 1, kres = 5, kmin = -3, kmax = 3, cmin = 0, cmax = 0)$F[2]
  Fdist$f3[a0] = model_AL(A0 = Fdist$a0s[[a0]], cres = 1, kres = 5, kmin = -3, kmax = 3, cmin = 0, cmax = 0)$F[3]
  Fdist$f4[a0] = model_AL(A0 = Fdist$a0s[[a0]], cres = 1, kres = 5, kmin = -3, kmax = 3, cmin = 0, cmax = 0)$F[4]
  Fdist$f5[a0] = model_AL(A0 = Fdist$a0s[[a0]], cres = 1, kres = 5, kmin = -3, kmax = 3, cmin = 0, cmax = 0)$F[5]
}

colnames(Fdist) = c("a0s", (round(logspace(-3,3,5),3)))

Fdist = Fdist %>% pivot_longer(-c(a0s), names_to = "K", values_to = "F") %>% mutate(K = factor(K, levels = unique(K)))

#figure 2C
Fdist %>% 
  ggplot(aes(F,K, fill = K))+
  geom_density_ridges(alpha = 0.7)+labs(y = "K")+
  theme_bw(base_size = 30,base_rect_size = 2,
           base_line_size = 0)+coord_fixed(1)+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), limits = c(10^-4,10^0.75),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(size = 1, sides = "b")+
  scale_fill_manual(name = "K",values = myligands(length(unique(Fdist$K))))+
  theme(legend.position = "none")

#figure 2D
Fdist %>% 
  group_by(K) %>% 
  summarise(cv = sd(F)/mean(F)) %>%
  ggplot(aes(cv,K, fill = K))+
  geom_col(position = "identity",alpha = 0.7)+labs(x = "CV", y = "K")+
  theme_bw(base_size = 30,base_rect_size = 2,
           base_line_size = 0)+
  scale_x_continuous(breaks = c(0,0.5,1))+
  coord_fixed(ratio = 0.5,xlim = c(0,1))+
  scale_fill_manual(name = "K",values = myligands(length(unique(Fdist$K))))+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.ticks = element_line(size = 1),
        legend.position = "none")



