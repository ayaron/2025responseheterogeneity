library(pracma)
library(tidyverse)
library(ggridges)
library(plot3D)

myligands = colorRampPalette(c("midnightblue","skyblue3","green", "green4")) # ligand colors
#model ALB with ordered ligand binding sequence
model_ALB_ordered = function(A0 = 0.5, B0 = 0.5, kpmin = -3, kpmax = 3,
                     kfmin = -3, kfmax = 3, cmin = 0,cmax = 0,
                     kpres = 101, kfres = 101, cres = 1, e = 1){
  #preparing the different matrices for the different parameters
  C0 = logspace(cmin, cmax, cres)
  Kp = logspace(kpmin, kpmax, kpres)
  Kf = logspace(kfmin, kfmax, kfres)
  F = array(dim = c(kpres, kfres, cres)) #full complex
  P =  array(dim = c(kpres, kfres, cres)) #partial complex A bound to ligand
  A =  array(dim = c(kpres, kfres, cres)) #receptor A
  B =  array(dim = c(kpres, kfres, cres)) #receptor B
  Sa = array(dim = c(kpres, kfres, cres)) #sensitvity to receptor A
  Sb =  array(dim = c(kpres, kfres, cres)) #sensitvity to receptor A
  
  #filling the different variables based on the different parameters
  for (c in seq_along(C0)) {
    for (kp in seq_along(Kp)) {
      for (kf in seq_along(Kf)) {
        y = ((Kf[kf]*Kp[kp]*C0[c])^2)*((A0-B0)^2) +
          2*Kf[kf]*Kp[kp]*C0[c]*(A0 + B0)*(1 + Kp[kp]*C0[c]) + (1 + Kp[kp]*C0[c])^2
        F[kf, kp, c] = (Kf[kf]*Kp[kp]*C0[c]*(A0+B0) + 1 + Kp[kp]*C0[c] - sqrt(y))/(2*Kf[kf]*Kp[kp]*C0[c])
        P[kf, kp, c] = Kp[kp]*C0[c]*(A0 - F[kf,kp, c])/(1 + Kp[kp]*C0[c]) 
        A[kf, kp, c] = A0 - P[kf, kp, c] - F[kf, kp, c]
        B[kf, kp, c] = B0 - F[kf, kp, c]
        Sa[kf, kp, c] = (A0*Kf[kf]*Kp[kp]*C0[c])/
          (Kf[kf]*Kp[kp]*C0[c]*(A0 + B0) + 1 + Kp[kp]*C0[c] - sqrt(y))*
          (1 - (Kf[kf]*Kp[kp]*C0[c]*(A0 - B0) + 1 + Kp[kp]*C0[c])/(sqrt(y)))
        Sb[kf, kp, c] = (B0*Kf[kf]*Kp[kp]*C0[c])/
          (Kf[kf]*Kp[kp]*C0[c]*(A0 + B0) + 1 + Kp[kp]*C0[c] - sqrt(y))*
          (1 - (Kf[kf]*Kp[kp]*C0[c]*(B0 - A0) + 1 + Kp[kp]*C0[c])/(sqrt(y)))
      }
    }
  }

  dynamics = tibble(expand.grid(C0 = C0, Kf = Kf, Kp = Kp), F = c(F), Sa = c(Sa),
                    Sb = c(Sb), E = c(e*F), P = c(P), A = c(A), B = c(B))
  dynamics
  
}
#model ALB with unordered ligand binding sequence
model_ALB_unordered = function(A0 = 0.5, B0 = 0.5, kpamin = -2, kpamax = 2,
                               kpbmin = -2, kpbmax = 2, kfmin = -2.8, kfmax = 2.8, cmin = 0,cmax = 0,
                               kpares = 101, kpbres = 101, kfres = 101, cres = 1, e = 1){
  #preparing the different matrices for the different parameters
  C0 = logspace(cmin, cmax, cres)
  Kpa = logspace(kpamin, kpamax, kpares)
  Kpb = logspace(kpbmin, kpbmax, kpbres)
  Kf = logspace(kfmin, kfmax, kfres)
  F = array(dim = c(kpares, kpbres, kfres, cres)) #full complex
  Pa =  array(dim = c(kpares, kpbres, kfres, cres)) #partial complex A bound to ligand
  Pb =  array(dim = c(kpares, kpbres, kfres, cres)) #partial complex B bound to ligand
  A =  array(dim = c(kpares, kpbres, kfres, cres)) #receptor A
  B =  array(dim = c(kpares, kpbres, kfres, cres)) #receptor B
  Sa = array(dim = c(kpares, kpbres, kfres, cres)) #sensitvity to receptor A
  Sb = array(dim = c(kpares, kpbres, kfres, cres)) #sensitvity to receptor A
  
  #filling the different variables based on the different parameters
  for (c in seq_along(C0)) {
    for (kpb in seq_along(Kpb)) {
      for (kpa in seq_along(Kpa)) {
        for (kf in seq_along(Kf)) {
          z = ((Kf[kf]*Kpb[kpb]*C0[c])^2)*(A0-B0)^2 +
            2*Kf[kf]*Kpb[kpb]*C0[c]*(A0 + B0)*(1 + Kpa[kpa]*C0[c])*(1 + Kpb[kpb]*C0[c]) + ((1 + Kpa[kpa]*C0[c])*(1 + Kpb[kpb]*C0[c]))^2
          F[kf, kpa, kpb, c] = (Kf[kf]*Kpb[kpb]*C0[c]*(A0+B0) + (1 + Kpa[kpa]*C0[c])*(1 + Kpb[kpb]*C0[c]) - sqrt(z))/
            (2*Kf[kf]*Kpb[kpb]*C0[c])
          Pa[kf, kpa, kpb, c] = Kpa[kpa]*C0[c]*(A0 - F[kf, kpa, kpb, c])/(1 + Kpa[kpa]*C0[c]) 
          Pb[kf, kpa, kpb, c] = Kpb[kpb]*C0[c]*(B0 - F[kf, kpa, kpb, c])/(1 + Kpb[kpb]*C0[c]) 
          A[kf, kpa, kpb, c] = A0 - Pa[kf, kpa, kpb, c] - F[kf, kpa, kpb, c]
          B[kf, kpa, kpb, c] = B0 - Pb[kf, kpa, kpb, c] - F[kf, kpa, kpb, c]
          Sa[kf, kpa, kpb, c] = (A0*Kf[kf]*Kpb[kpb]*C0[c])/
            (Kf[kf]*Kpb[kpb]*C0[c]*(A0 + B0) + (1 + Kpa[kpa]*C0[c])*(1 + Kpb[kpb]*C0[c]) - sqrt(z))*
            (1 - (Kf[kf]*Kpb[kpb]*C0[c]*(A0 - B0) + (1 + Kpa[kpa]*C0[c])*(1 + Kpb[kpb]*C0[c]))/(sqrt(z)))
          Sb[kf, kpa, kpb, c] = (B0*Kf[kf]*Kpb[kpb]*C0[c])/
            (Kf[kf]*Kpb[kpb]*C0[c]*(A0 + B0) + (1 + Kpa[kpa]*C0[c])*(1 + Kpb[kpb]*C0[c]) - sqrt(z))*
            (1 - (Kf[kf]*Kpb[kpb]*C0[c]*(B0 - A0) + (1 + Kpa[kpa]*C0[c])*(1 + Kpb[kpb]*C0[c]))/(sqrt(z)))
        }
      }
    }
  }
  
  dynamics = tibble(expand.grid(C0 = C0, Kf = Kf, Kpa = Kpa, Kpb = Kpb), F = c(F), Sa = c(Sa),
                    Sb = c(Sb), E = c(e*F), Pa = c(Pa), Pb = c(Pb), A = c(A), B = c(B))
  dynamics
  
}

#preperation for supplementary figure 3A-D, G, figure 4B ii-iii
A0 = c(0.001, 0.25, 0.5, 0.75, 0.999) #different amounts of A0
ALB_o = tibble(.rows = 1)
for (a0 in seq_along(A0)) {
  temp = tibble(F = model_ALB_ordered(A0 = A0[a0], B0 = 1 - A0[a0])$F,
    Sa = model_ALB_ordered(A0 = A0[a0], B0 = 1 - A0[a0])$Sa,
    Sb = model_ALB_ordered(A0 = A0[a0], B0 = 1 - A0[a0])$Sb,
    C0 = model_ALB_ordered(A0 = A0[a0], B0 = 1 - A0[a0])$C0,
    Kf = model_ALB_ordered(A0 = A0[a0], B0 = 1 - A0[a0])$Kf,
    Kp = model_ALB_ordered(A0 = A0[a0], B0 = 1 - A0[a0])$Kp,
    A = model_ALB_ordered(A0 = A0[a0], B0 = 1 - A0[a0])$A,
    B = model_ALB_ordered(A0 = A0[a0], B0 = 1 - A0[a0])$B,
    A0 = A0[a0],
    B0 = 1-A0[a0]
  )
  ALB_o = bind_rows(ALB_o, temp)
  rm(temp)
}
ALB_o = ALB_o%>%slice(-1)

#figure 4B i
ALB_o %>% 
  dplyr::filter(A0 == 0.001) %>%
  pivot_longer(cols = c(Sa, Sb), names_to = "S") %>% 
  ggplot(aes(Kf, Kp, fill = value)) + 
  geom_tile()+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), limits = c(10^-3,10^3),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                expand = c(0,0))+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), limits = c(10^-3,10^3),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                expand = c(0,0))+
  annotation_logticks(outside = T)+
  theme_bw(base_line_size = NA, base_size = 20,base_rect_size = 1)+
  scale_fill_viridis_c(name = "SL",
                       guide =
                         guide_colorbar(label = TRUE,frame.linewidth = 2,
                                        ticks = F,frame.colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 18),
  )+
  coord_cartesian(clip = "off") + 
  facet_wrap(vars(S), ncol = 1)

#figure 4B iii
ALB_o %>% 
  dplyr::filter(A0 == 0.999) %>%
  pivot_longer(cols = c(Sa, Sb), names_to = "S") %>% 
  ggplot(aes(Kf, Kp, fill = value)) + 
  geom_tile()+
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), limits = c(10^-3,10^3),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                expand = c(0,0))+
  scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x), limits = c(10^-3,10^3),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                expand = c(0,0))+
  annotation_logticks(outside = T)+
  theme_bw(base_line_size = NA, base_size = 20,base_rect_size = 1)+
  scale_fill_viridis_c(name = "SL",
                       guide =
                         guide_colorbar(label = TRUE,frame.linewidth = 2,
                                        ticks = F,frame.colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 18),
  )+
  coord_cartesian(clip = "off") + 
  facet_wrap(vars(S), ncol = 1)


#supplementary figure 3A-D
scatter3D(x = log10(ALB_o$Kf), y = log10(ALB_o$Kp), z = ALB_o$A0, pch = 20, bty = 'b2',
          xlab = "KF(A0+B0)", ylab = "KPC0", zlab = "A0 relative abbundance", main = "Full Complex", theta = 130,
          colvar = ALB_o$F, col = viridisLite::cividis(256),  alpha = 0.35, ticktype = "detailed") #A
scatter3D(x = log10(ALB_o$Kf), y = log10(ALB_o$Kp), z = ALB_o$A0, pch = 20, bty = 'b2',
          xlab = "KF(A0+B0)", ylab = "KPC0", zlab = "A0 relative abbundance", main = "Sensitivities towards A0", theta = 130,
          colvar = ALB_o$Sa, col = viridisLite::cividis(256),  alpha = 0.35, ticktype = "detailed") #B
scatter3D(x = log10(ALB_o$Kf), y = log10(ALB_o$Kp), z = ALB_o$A0, pch = 20, bty = 'b2',
          xlab = "KF(A0+B0)", ylab = "KPC0", zlab = "A0 relative abbundance", main = "Sensitivities towards B0",theta = 130,
          colvar = ALB_o$Sb, col = viridisLite::cividis(256),  alpha = 0.35, ticktype = "detailed")# C
scatter3D(x = log10(ALB_o$Kf), y = log10(ALB_o$Kp), z = ALB_o$A0, pch = 20, bty = 'b2',
          xlab = "KF(A0+B0)", ylab = "KPC0", zlab = "A0 relative abbundance", main = "SB + SA", theta = 130,
          colvar = ALB_o$Sa + ALB_o$Sb, clim =  c(0, max(ALB_o$Sa + ALB_o$Sb)), col = viridisLite::cividis(256),  alpha = 0.35, ticktype = "detailed") #D

#supplementary figure 3G
ALB_o %>% 
  pivot_longer(cols = c(Sa, Sb), names_to = "Stype", values_to = "Sval") %>% 
  dplyr::filter(A0 != .001 & A0 !=.999) %>%
  dplyr::filter(Kp %in% logspace(-3,3,21), Kf %in% logspace(-3,3,21))%>%
  group_by(A0) %>%
  mutate(nF = F/max(F), 
         title = case_when(A0 > B0 ~"A0>B0",
                           A0 == B0 ~"A0=B0",
                           A0 < B0 ~"A0<B0"))%>% 
  ggplot(aes(nF))+
  geom_point(aes(y = Sval, color = Stype))+
  facet_wrap(~title, nrow = 1, scales = "free_x")+
  theme_bw() +
  labs(y = "SL", x = "Bound resceptors (%)") +
  scale_color_manual(name = "", values = c("#FFD400","#92278B"))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1), text = element_text(size = 18),
        strip.background = element_blank())

#preperation for figure 4B ii and supplementary figure 3E
delta_percent = 0.01 #how many points along the receptor range will be sampled
receptor_a_percent = seq(from = 0.0, to = 1, by = delta_percent) #A0 percentage of the total receptor amount, from 0 to 100%
receptor_sum = logseq(0.5, 2, 3) #the total amount of receptors sampled
Sa_range = zeros(2,length(receptor_a_percent)) #vector for Sa range along the changes in percentage of A0
Sb_range = zeros(2,length(receptor_a_percent))#vector for Sb range along the changes in percentage of B0
S_range = tibble(.rows = 1)
Sa_range_u = zeros(2,length(receptor_a_percent))
Sb_range_u = zeros(2,length(receptor_a_percent))
S_range_u = tibble(.rows = 1)
for (rs in seq_along(receptor_sum)) {
  for (rap in seq_along(receptor_a_percent)) {
    a0 = receptor_sum[rs]*receptor_a_percent[rap]
    b0 = receptor_sum[rs]-a0
    Sa_range[1,rap] = range(model_ALB_ordered(A0 = a0, B0 = b0, kpres = 3, kfres = 3)$Sa)[2]
    Sa_range[2,rap] = range(model_ALB_ordered(A0 = a0, B0 = b0, kpres = 3, kfres = 3)$Sa)[1]
    Sb_range[1,rap] = range(model_ALB_ordered(A0 = a0, B0 = b0, kpres = 3, kfres = 3)$Sb)[2]
    Sb_range[2,rap] = range(model_ALB_ordered(A0 = a0, B0 = b0, kpres = 3, kfres = 3)$Sb)[1]
    Sa_range_u[1,rap] = range(model_ALB_unordered(A0 = a0, B0 = b0, kpares = 3, kpbres = 3, kfres = 3)$Sa)[2]
    Sa_range_u[2,rap] = range(model_ALB_unordered(A0 = a0, B0 = b0, kpares = 3, kpbres = 3, kfres = 3)$Sa)[1]
    Sb_range_u[1,rap] = range(model_ALB_unordered(A0 = a0, B0 = b0, kpares = 3, kpbres = 3, kfres = 3)$Sb)[2]
    Sb_range_u[2,rap] = range(model_ALB_unordered(A0 = a0, B0 = b0, kpares = 3, kpbres = 3, kfres = 3)$Sb)[1]

  }
  temp = tibble(receptor_a_percent = receptor_a_percent, Sa_min = Sa_range[2,], Sa_max = Sa_range[1,],
               Sb_min = Sb_range[2,], Sb_max = Sb_range[1,], receptor_sum = receptor_sum[rs])
  temp_u = tibble(receptor_a_percent = receptor_a_percent, Sa_min = Sa_range_u[2,], Sa_max = Sa_range_u[1,],
                  Sb_min = Sb_range_u[2,], Sb_max = Sb_range_u[1,], receptor_sum = receptor_sum[rs])
  S_range = bind_rows(S_range,temp)
  S_range_u = bind_rows(S_range_u, temp_u)
  rm(temp, temp_u)
}
S_range = S_range %>% slice(-1)
S_range_u = S_range_u %>% slice(-1)

#supplementary figure 3E
S_range%>% 
  pivot_longer(cols = c(Sa_min, Sa_max, Sb_min, Sb_max), names_to = "S_range", values_to = "S_val") %>% 
  separate(col = S_range, into = c("S", "range"), sep = "_") %>%
  pivot_wider(names_from = range, values_from = S_val) %>%
  ggplot(aes(x = receptor_a_percent, alpha = as.factor(receptor_sum)))+
  geom_ribbon(aes(ymin = min, ymax = max, fill = as.character(S))) +
  scale_y_continuous(name = "SL", breaks = seq(0,1,0.2)) +
  scale_fill_manual(name = "", values = c("orange", "purple1")) +
  theme_bw(base_line_size = 0, base_rect_size = 1, base_size = 20) +
  labs(x = "A0/A0+B0", title = "SL range")+
  scale_alpha_discrete(range = c(0.55,0.25), name = "A0+B0") +
  coord_fixed(0.5)

#figure 4B ii
S_range%>% 
  dplyr::filter(receptor_sum == 1) %>% 
  pivot_longer(cols = c(Sa_min, Sa_max, Sb_min, Sb_max), names_to = "S_range", values_to = "S_val") %>% 
  separate(col = S_range, into = c("S", "range"), sep = "_") %>%
  pivot_wider(names_from = range, values_from = S_val) %>%
  ggplot(aes(x = receptor_a_percent))+
  geom_ribbon(aes(ymin = min, ymax = max, fill = as.character(S)), alpha = 0.7) +
  geom_rect(aes(xmin = 0.98, xmax = 1, ymin = -Inf, ymax = Inf), color = "orange3", size = 1, fill = NA)+
  geom_rect(aes(xmin = 0.0, xmax = 0.02, ymin = -Inf, ymax = Inf), color = "purple4", size = 1, fill = NA)+
  scale_y_continuous(name = "SL", breaks = seq(0,1,0.2)) +
  scale_fill_manual(name = "", values = c("orange", "purple1")) +
  theme_bw(base_line_size = 0, base_rect_size = 1, base_size = 20) +
  labs(x = "A0/A0+B0", title = "SL range")+
  coord_fixed(0.5)

#supplementary figure 3F
S_range_u %>% 
  dplyr::filter(receptor_sum == 1) %>% 
  pivot_longer(cols = c(Sa_min, Sa_max, Sb_min, Sb_max), names_to = "S_range", values_to = "S_val") %>% 
  separate(col = S_range, into = c("S", "range"), sep = "_") %>%
  pivot_wider(names_from = range, values_from = S_val)%>%
  ggplot(aes(x = receptor_a_percent))+
  geom_ribbon(aes(ymin = min, ymax = max, fill = as.character(S)), alpha = 0.7) +
  scale_y_continuous(name = "SL", breaks = seq(0,1,0.2)) +
  scale_fill_manual(name = "", values = c("orange", "purple1")) +
  theme_bw(base_line_size = 0, base_rect_size = 1, base_size = 20) +
  labs(x = "A0/A0+B0", title = "SL range - nonsequential receptor binding")+
  coord_fixed(0.5)

#preparation for figure 4E

Cells = tibble(a0mean = c(0.25, 0.5, 0.75, 0.9), b0mean = c(0.75, 0.5, 0.25, 0.1),
            a0V = c(0.001, 0.25, 0.1, 0.5), b0V = c(0.75, 0.25, 0.5, 0.01), 
            a0_q25 = NA, b0_q25 = NA, a0_q75 = NA, b0_q75 = NA)
Ks = tibble(Kf = c(-3,0,3), Kp = c(3,3,0), e = NA)
Edist = tibble(.rows = 1)

set.seed(1000000)
for (cell in seq_along(Cells$a0mean)) {
  a = rgamma(n = 100000, shape = Cells$a0mean[[cell]]/(Cells$a0V[[cell]]^2),
             scale = Cells$a0V[[cell]]^2)
  Cells$a0_q25[cell] = quantile(a, prob = c(.25,.75))[1]
  Cells$a0_q75[cell] = quantile(a, prob = c(.25,.75))[2]
  b = rgamma(n = 100000, shape = Cells$b0mean[[cell]]/(Cells$b0V[[cell]]^2),
             scale = Cells$b0V[[cell]]^2)
  Cells$b0_q25[cell] = quantile(b, prob = c(.25,.75))[1]
  Cells$b0_q75[cell] = quantile(b, prob = c(.25,.75))[2]
  fdist = tibble(A0 = a, B0 = b, k1 = NA, k2 = NA, k3 = NA,
                 a0_mean_sd = paste(Cells$a0mean[cell], Cells$a0V[cell], sep = "_"),
                 b0_mean_sd = paste(Cells$b0mean[cell], Cells$b0V[cell], sep = "_"))
  for (k in seq_along(Ks$Kf)) {
    for (ab in seq_along(fdist$A0)) {
      fdist[ab, k + 2] = model_ALB_ordered(A0 = fdist$A0[[ab]], B0 = fdist$B0[[ab]], kfres = 1, kpres = 1,
                                 kfmin = Ks$Kf[k], kfmax = Ks$Kf[k], 
                                 kpmax = Ks$Kp[k], kpmin = Ks$Kp[k])$F
      
    }
    if (is.na(Ks$e[k])) {
      e1 =  unlist(fdist[, k+2])/mean(unlist(fdist[, k+2]), na.rm = T)
      Ks$e[k] = mean(e1/unlist(fdist[,k+2]))
      rm(e1)
    }
  }
  Edist = bind_rows(Edist, fdist)
  rm(fdist)
}
Edist = Edist %>% slice(-1)

#figure 4E middle
Cells%>%
  ggplot(aes(b0mean, a0mean, color = as.factor(a0mean)))+
  geom_point(size = 5)+
  geom_errorbarh(aes(xmax = b0mean+b0_q75, xmin = b0mean-b0_q25),size = 1, height = 0.05)+
  geom_errorbar(aes(ymax = a0mean+a0_q75, ymin = a0mean-a0_q25), size = 1, width = 0.05)+
  coord_fixed()+
  xlim(0,2.25)+
  ylim(0,2.25)+
  theme_bw(base_size = 20, base_rect_size = 2)+
  scale_color_manual(values = colorRampPalette(c("purple", "pink2","orange","red"))(4))+
  labs(x = "B0 mean", y = "A0 mean", title = "Receptor Space") +
  theme(legend.position = "none")
  
#figure 4E density plots
Edist[,3:5] = t(t(matrix(unlist(Edist[,3:5]),ncol = 3))*unlist(Ks$e))
colnames(Edist)[3:5] = c(paste(round(10^Ks$Kf,3), round(10^Ks$Kp,3), round(Ks$e), sep = ";"))


Edist %>% 
  pivot_longer(cols = 3:5, names_to = "ligand", values_to = "E") %>% 
  unite("population", a0_mean_sd:b0_mean_sd, sep = ";") %>% 
  ggplot(aes(E, ligand, fill = ligand)) +
  geom_density_ridges(aes(height= ..ndensity..),scale = 0.9,alpha = 0.5, size = 1) +
  scale_x_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x))) + 
  annotation_logticks(size = 1, sides = "b")+
  theme_ridges()+
  scale_fill_manual(values = myligands(3))+
  facet_wrap(~population, nrow = 1)+
  theme_bw(base_size = 30,base_rect_size = 2,
           base_line_size = 1)+coord_fixed(0.5)+
  theme(legend.position = "none")

#figure 4E ligand legend
Ks %>% 
  mutate(Ligand = c("ligand1", "ligand2", "ligand3"), e= round(log10(e), 2)) %>% 
  pivot_longer(c(Kf, Kp, e), names_to = "par") %>% 
  mutate(labs = paste0("10^", value)) %>% 
  ggplot(aes(par, Ligand)) +
  geom_tile(aes( fill = Ligand), alpha = 0.7)+
  geom_text(aes(label = labs), parse = T)+
  scale_fill_manual(values = myligands(3))+
  theme_bw(base_line_size = NA, base_size = 20) +
  labs(x = "", y = "") + 
  scale_x_discrete(expand = c(0,0), position = "top")+
  scale_y_discrete(expand = c(0,0)) +
  theme(legend.position = "none")





#preparation for supplementary figure 4A
Edist_x = tibble(expand.grid(
  A0mean = c(0.001, 0.25, 0.5, 0.75, 0.999), 
  A0_Var = c(0.001,0.5,1), B0_Var = c(0.001,0.5,1),
  Kf = c(-3,0,3), Kp = c(-3,0,3)),B0mean = 1-A0mean, E_Var = NA, e = NA)


set.seed(10^6)
for (E in seq_along(Edist_x$A0mean)) {
  Fs = rep(NA,10^5)
  for (f in seq_along(Fs)) {
    a = rgamma(n = 1, shape = Edist_x$A0mean[[E]]/(Edist_x$A0_Var[[E]]^2),
               scale = Edist_x$A0_Var[[E]]^2)
    b = rgamma(n = 1, shape = Edist_x$B0mean[[E]]/(Edist_x$B0_Var[[E]]^2),
               scale = Edist_x$B0_Var[[E]])
    Fs[f] = model_ALB_ordered(A0 = a, B0 = b, kpmin = Edist_x$Kp[E],
                              kpmax = Edist_x$Kp[E], kfmin = Edist_x$Kf[E],
                              kfmax = Edist_x$Kf[E], kpres = 1, kfres = 1)$F
  }
  Es = Fs/mean(Fs, na.rm = T) #making the mean of all response distributions 1
  Edist_x$E_Var[[E]] = std(Es)
  Edist_x$e[[E]] = mean(Es/Fs)
}

#supplemetary figure 4A, main. For the zoom on the middle add filter(A0mean>0.01, B0mean>0.01) 
Edist_x %>% 
  group_by(A0mean, B0_Var, A0_Var) %>%
  mutate(nE_var = E_Var/min(E_Var, na.rm = T))%>%
  ggplot(aes(interaction(Kf,A0_Var), interaction(Kp,B0_Var))) + 
  geom_tile(color = "black", size = 1.1)+
  geom_tile(aes(fill = nE_var))+
  geom_hline(yintercept = 3.5, color = "black")+
  geom_hline(yintercept = 6.5, color = "black")+
  geom_vline(xintercept = 3.5, color = "black")+
  geom_vline(xintercept = 6.5, color = "black")+
  labs(x = "log(KF);A0SD", y = "log(KP);B0SD")+
  coord_fixed()+
  facet_wrap(.~A0mean, nrow = 1, labeller = "label_both")+
  scale_fill_viridis_c(na.value = "white",
                       name = "Normalized\nE CV",
                       guide =
                         guide_colorbar(label = TRUE,frame.linewidth = 2,
                                        ticks = F,frame.colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust =1), text = element_text(size = 18),,
        strip.background = element_blank(),)


#preperation for supplementary figure 4B
vars = unique(Edist_x %>% select(c(A0_Var, B0_Var)))
S_comb = tibble(.rows = 1)
A0_1 = A0 #as there is another A0 in the ALB_o function which seems to interfere with the function here
for (a0 in seq_along(A0_1)) {
  temp1 = ALB_o %>% dplyr::filter(Kp %in% logspace(3,-3,3)& Kf %in% logspace(-3,3,3)& A0 == A0_1[a0])
  temp2 = tibble(Kp = 1:9)
  for (s in seq_along(vars$A0_Var)) {
    temp2$Kp = temp1$Kp
    temp2$Kf = temp1$Kf
    temp2$A0_Var = vars$A0_Var[[s]]
    temp2$B0_Var = vars$B0_Var[[s]]
    temp2$A0mean = temp1$A0[[a0]]
    temp2$B0mean = temp1$B0[[a0]]
    temp2$S_tot = temp1$Sa*vars$A0_Var[[s]] + temp1$Sb*vars$B0_Var[[s]]
    S_comb = bind_rows(S_comb, temp2)
  }
  
}
S_comb = S_comb %>% slice(-1)

#supplemetary figure 4B, main. For the zoom on the middle add filter(A0mean>0.01, B0mean>0.01) 
S_comb%>%
  group_by(A0mean, A0_Var, B0_Var) %>%
  mutate(nS_tot = S_tot/min(S_tot, na.rm = T),
         Kf = paste0("10^", log10(Kf)))%>%
  ggplot(aes(interaction(Kf,A0_Var), interaction(Kp,B0_Var))) + 
  geom_tile(color = "black", size = 1.1)+
  geom_tile(aes(fill = nS_tot))+
  geom_hline(yintercept = 3.5, color = "black")+
  geom_hline(yintercept = 6.5, color = "black")+
  geom_vline(xintercept = 3.5, color = "black")+
  geom_vline(xintercept = 6.5, color = "black")+
  labs(x = "log(KF);A0SD", y = "log(KP);B0SD")+
  coord_fixed()+
  facet_wrap(.~A0mean, nrow = 1, labeller = "label_both")+
  scale_fill_viridis_c(na.value = "white",
                       name = "Normalized\ntotal S",
                       guide =
                         guide_colorbar(label = TRUE,frame.linewidth = 2,
                                        ticks = F,frame.colour = "black"))+
  theme(axis.text.x = element_text(angle = 45, hjust =1), text = element_text(size = 18),,
        strip.background = element_blank(),)




#suppementary figure 4C
left_join(Edist_x%>% mutate(Kf = 10^(Kf), Kp = 10^(Kp)),S_comb)%>% 
  group_by(A0mean, A0_Var, B0_Var) %>%
  mutate(nS_tot = S_tot/min(S_tot, na.rm = T),
         nE_var = E_Var/min(E_Var, na.rm = T))%>%
  ggplot(aes(nE_var, sqrt(nS_tot)))+
  geom_point(aes(color = interaction(Kp,Kf), size = interaction(A0_Var,B0_Var)), alpha = 0.5)+
  scale_size_discrete(name = "A0SD;B0SD") + 
  scale_color_manual(name = "Kp;Kf",values = myligands(9))+
  labs(x = "Normalized CV(E)", y = expression(sqrt("S*SD"))) +
  theme_bw(base_size = 30,base_rect_size = 2,
           base_line_size = 0)


