library(tidyverse)
library(ggridges)
library(pracma)

myligands = colorRampPalette(c("midnightblue","skyblue3","green", "green4")) # ligand colors
# model ALA
model_ALA = function(kfmin = -3,kfmax = 3,cmin = 0, cmax = 0, cres = 1, e = 1, kpmax = 3, kpmin = -3, A0 = 1, kfres = 101, kpres = 101){
  #preparing the different matrices for the different parameters
  C0 = logspace(cmin, cmax, cres)
  Kf = logspace(kfmin, kfmax, kfres)
  Kp = logspace(kpmin, kpmax, kpres)
  F = array(dim = c(kpres, kfres, cres)) #full complex
  P =  array(dim = c(kpres, kfres, cres)) #partial complex A bound to ligand
  A =  array(dim = c(kpres, kfres, cres)) #receptor A
  S = array(dim = c(kpres, kfres, cres)) #sensitivity
  
  #filling the different variables based on the different parameters
  for (c in seq_along(C0)) {
    for (kp in seq_along(Kp)) {
      for (kf in seq_along(Kf)) {
        x = 8*A0*Kf[kf]*Kp[kp]*C0[c] + (1 + Kp[kp]*C0[c])^2
        F[kf, kp, c] = (4*A0*Kf[kf]*Kp[kp]*C0[c] + (1 + Kp[kp]*C0[c])^2 - (1 + Kp[kp]*C0[c])*sqrt(x))/
          (8*Kf[kf]*Kp[kp]*C0[c])
        P[kf, kp, c] = Kp[kp]*C0[c]*(A0 - 2*F[kf, kp, c])/(1 + Kp[kp]*C0[c])
        A[kf, kp, c] = A0 - P[kf, kp, c] - 2*F[kf, kp, c]
        S[kf, kp, c] = (4*A0*Kf[kf]*Kp[kp]*C0[c])/
          (4*A0*Kf[kf]*Kp[kp]*C0[c] + (1 + Kp[kp]*C0[c])^2 - (1 + Kp[kp]*C0[c])*sqrt(x))*
          (1 - (1 + Kp[kp]*C0[c])/(sqrt(x)))
      }
    }
  }
  
  dynamics = tibble(expand.grid(C0 = C0, Kf = Kf, Kp = Kp), F = c(F), S = c(S), E = c(e*F),
                    P = c(P), A = c(A))
  dynamics
}  

#a data frame holding the data for figures 3B, C and sup fig 2A
ALA = tibble(F = model_ALA()$F,
            S = model_ALA()$S,
            C0 = model_ALA()$C0,
            Kf = model_ALA()$Kf,
            Kp = model_ALA()$Kp,
            A = model_ALA()$A)

# figure 3B
ALA %>% 
  ggplot(aes(Kp, Kf, fill = S)) + 
  geom_tile() + 
  scale_fill_viridis_c(limits = c(min(ALA$S), max(ALA$S)),
                       breaks = c(min(ALA$S), max(ALA$S)),
                       guide =
                         guide_colorbar(label = TRUE,frame.linewidth = 2,
                                        ticks = F,frame.colour = "black")
  )+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), limits = c(10^-3,10^3),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                expand = c(0,0))+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), limits = c(10^-3,10^3),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                expand = c(0,0))+
  annotation_logticks(outside = T)+
  coord_cartesian(clip = "off") +
  labs(title = "Sensitivity") + 
  theme_bw(base_line_size = 0, base_size = 20,base_rect_size = 2)+
  theme(text=element_text(size=20))

#figure 3C
ALA %>% 
  ggplot(aes(F/max(F), S)) + # the x axis gives the percerntage of bound receptor
  geom_point(size = 0.5) + 
  theme_bw()+
  labs(x = "Bound receptor (%)")

#supplementary figure 2A
ALA %>% 
  ggplot(aes(Kp, Kf, fill = F)) + 
  geom_tile() + 
  scale_fill_viridis_c(limits = c(min(ALA$F), max(ALA$F)),
                       breaks = c(min(ALA$F), max(ALA$F)),
                       guide =
                         guide_colorbar(label = TRUE,frame.linewidth = 2,
                                        ticks = F,frame.colour = "black")
  )+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), limits = c(10^-3,10^3),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                expand = c(0,0))+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), limits = c(10^-3,10^3),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                expand = c(0,0))+
  annotation_logticks(outside = T)+
  coord_cartesian(clip = "off") +
  labs(title = "Full Complexes") + 
  theme_bw(base_line_size = 0, base_size = 20,base_rect_size = 2)+
  theme(text=element_text(size=20))

#preparing figures 2D-F

set.seed(10^6)  
Edist =  tibble(a0s = rgamma(n = 100000, shape = 1/(0.5^2), scale = 0.5^2), 
                f1 = NA, f2 = NA, f3 = NA, f4 = NA, f5 = NA)
Ks = tibble(Kf = c(-2,-1,0,1,2.7), Kp = c(-2.7,1,0,-1,2), e = NA)

for (a0 in seq_along(Edist$a0s)) {
  for (k in seq_along(Ks$Kf)) {
    e = 1/model_ALA(A0 = 1, kfres = 1, kpres = 1,
                  kfmin = Ks$Kf[k], kfmax = Ks$Kf[k], 
                  kpmax = Ks$Kp[k], kpmin = Ks$Kp[k])$F #to get all the vectors to average at one
    Edist[a0,k+1] = model_ALA(A0 = Edist$a0s[[a0]], kfres = 1, kpres = 1,
                             kfmin = Ks$Kf[k], kfmax = Ks$Kf[k], 
                             kpmax = Ks$Kp[k], kpmin = Ks$Kp[k], e = e)$E
    if(is.na(Ks$e[k])){
      Ks$e[k] = e
    }
  }
}

#figure 3D left
Ks %>% 
  ggplot(aes(10^Kp, 10^Kf)) +
  geom_point(aes(color = as.factor(Kf))) +
  scale_color_manual(values = myligands(length(Ks$Kf)))+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), limits = c(10^-3,10^3),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                expand = c(0,0))+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), limits = c(10^-3,10^3),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                expand = c(0,0))+
  labs(x = "KP", y = "KF")+
  annotation_logticks(outside = T)+
  coord_cartesian(clip = "off") +
  theme_bw(base_line_size = 0, base_size = 20,base_rect_size = 2)+
  theme(text=element_text(size=20))
  
#figure 3D right
Ks %>% 
  ggplot(aes(0, e)) +
  geom_point(aes(color = as.factor(Kf)), size = 2.5)+
  scale_color_manual(values = myligands(length(Ks$Kf)))+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                limits = c(10^0, 10^5),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                expand = c(0,0))+
  labs(x = "eL", y = "")+
  annotation_logticks(sides = "l", size = 0.5)+
  coord_fixed() +
  theme_bw(base_line_size = 0, base_size = 20,base_rect_size = 0)+
  theme(text=element_text(size=20), axis.text.x = element_blank())
  
  
colnames(Edist) = c("a0s", paste(round(10^Ks$Kf,3), round(10^Ks$Kp,3), round(Ks$e), sep = ";"))


# figure 2E
Edist %>% 
  pivot_longer(cols = c(2:6), names_to = "Ligand (KF;KP)") %>% 
  ggplot(aes(value, `Ligand (KF;KP)`, fill = `Ligand (KF;KP)`))+
  geom_density_ridges(alpha = 0.7)+
  theme_bw(base_size = 30,base_rect_size = 2, base_line_size = 0)+
  coord_fixed()+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(size = 1, sides = "b")+
  scale_fill_manual(values = myligands(length(Ks$Kf)))+
  theme(legend.position = "none")


Ks$Ecv = unlist(lapply(Edist[2:6], std))

#figure 3F
Ks %>% 
  ggplot(aes(Ecv, as.factor(Kf))) +
  geom_col(aes(fill = as.factor(Kf)), alpha = 0.7) +
  theme_bw(base_size = 30,base_rect_size = 2,
           base_line_size = 0)+
  scale_x_continuous(limits = c(0, max(Ks$Ecv)), expand = c(0,0)) +
  scale_fill_manual(name = "K111",values = myligands(length(unique(Ks$Kf))))+
  labs(y = "Ligand", x = expression(paste(sigma, "(E)"))) +
  theme(axis.ticks = element_line(size = 1),
    legend.position = "none")  

# preparing supplementary figure 2B
Ks$S = NA 
for (s in seq_along(Ks$S)) {
  Ks$S[s] = model_ALA(A0 = 1, kfres = 1, kpres = 1,
                      kfmin = Ks$Kf[s], kfmax = Ks$Kf[s], 
                      kpmax = Ks$Kp[s], kpmin = Ks$Kp[s])$S
}

#supplementary figure 2B
Ks %>% 
  ggplot(aes(S, Ecv)) + 
  geom_smooth(color = "black",se = F,method = "lm", formula = y~x)+
  geom_point(aes(color = as.factor(round(Kf, 3))))+
  scale_color_manual(values = myligands(5))+
  labs(x = "Sensitivity (S)", y = expression(paste("Response Variability (", sigma, "(E))"))) +
  coord_fixed() + 
  theme_bw(base_size = 30,base_rect_size = 2,
           base_line_size = 0)+
  theme(legend.position = "none")
