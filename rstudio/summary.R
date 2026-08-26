library(tidyverse)
library(magrittr)
library(ggthemes)
library(colorspace)
library(latex2exp)
library(scales)
library(directlabels)
library(geomtextpath)
library(plotly)
library(showtext)
library(stringr)
library(sysfonts)
library(ggnewscale)
library(patchwork)
font_add_google("Inter")

source('theme-mcm.R')

# DeCelles 2015
arizaro_tectonic_subs_t_Ma = c(20.54475457170356, 18.806544754571703, 17.836381135707413, 15.626564003849856, 14.171318575553418, 14.063522617901828, 13.753609239653514, 13.659287776708373, 8.45813282001925, 8.350336862367662)
arizaro_tectonic_subs_depth_m = c(0, 84.96732026143809, 166.66666666666652, 1477.124183006536, 1568.6274509803923, 1568.6274509803923, 1627.450980392157, 1630.7189542483661, 1633.986928104575, 1630.7189542483661)

# Garzione 2008?
altiplano_elev_t_Ma = rev(c(4.977426568691192, 6.580135541647188, 8.65688502209836, 10.090293260238296, 11.602708817957016, 12.957110976509998, 14.130925684592516, 14.706546012981415))
altiplano_elev_km = rev(c(3.986783739607607, 4, 1.845814893223917, 1.845814893223917, 1.872246688039992, 1.6607930554801018, 1.4889867521599691, 1.3964758332880616))

# He & Kapp
bidahochi_total_subs_t_Ma = c(20.047036055670603, 18.269049276987495, 15.964249993217514, 17.017875003391236, 15.470365982691101, 13.824082330279706, 13.725304081245078, 8.325493990721558, 6.942614782191913, 5.921918266578645, 5.197552880745889, 4.176856365132615)
bidahochi_total_subs_depth_m = c(-9.80389939322069, -11.764690042776678, 3.9215812991120202, -19.60785264100066, 25.49022458966867, 70.5882479657786, 82.35293800855527, 94.11768190589127, 101.96079064955599, 131.3725965383366, 76.47061991444666, 7.843162598223984)

yr = 3.154e7 # s per year
pi = 3.14159265359
e = 2.71828
FIG_WIDTH = 6.5
FIG_HEIGHT = 6.5

log10_breaks <- c(0.001,0.003,0.01,0.03,0.1,0.3,1,3,10,30,100,300,1000)

parse_array <- function(s) {
  list(array(as.numeric(unlist(strsplit(str_sub(s,2,-2), split = ",")))))
}

extract_value <- function(version, key) {
  rx = sprintf("(?<=(_|^)%s)(.+?)(?=(_|$))", key)
  val = str_extract(version, rx) %>% gsub("p", ".", .) %>% as.numeric
  return(val)
}

drip_type_colors <- c(
  "st" = "#F8766D", 
  "md" = "#7CAE00",
  "hy" = "#00BFC4",
  "wk" = "#C77CFF", 
  "?" = "gray25"
)

dat <- read_csv('../plots/pp/summary_2026.csv',show_col_types = FALSE) %>%
  rowwise() %>%
  mutate(
    depths=parse_array(depths),
    times=parse_array(times),
    topos=parse_array(topos),
    t_diff2 = times[4]-times[2],
    t0=times[1],
  ) %>%
  mutate(
    label=ifelse(version=="A1_G50_ch2_cs0p9_dmg1_ecx0p5_fr15_fs0p9_k2p5_lcx1_lin0_mlx0p2_mtz0_rec250_ruc2750_u0p8_ucx0p8", "A1",
      ifelse(version=="A1_G50_ch2_cs0p9_dmg1_ecx0p5_fr15_fs0p9_k2p5_lcx0_lin0_mlx0p5_mtz0_rec70_ruc2750_u0p8_ucx0p2", "A2",
      ifelse(version=="A1_G50_ch2_cs0p9_dmg1_ecx0p5_fr15_fs0p9_k2p5_lcx0p5_lin0_mlx0p5_mtz0_rec250_ruc2750_u0p8_ucx0p2", "A3", 
      ifelse(version=="A1_G50_ch2_cs0p9_dmg1_ecx0p5_fr15_fs0p9_k2p5_lcx1_lin0_mlx0p2_mtz0_rec70_ruc2750_u0p8_ucx0","A4","")))),
    # setup
    well_behaved_time=t_detach<=80,
    well_behaved_size=drip_width >= 60 & drip_width <= 150,
    well_behaved=well_behaved_size & well_behaved_time,
    m_extending = ifelse(m_extending==0, 15, m_extending),
    ver = str_split_1(version,"(?=_)") %>% .[str_detect(.,"_ml") | str_detect(.,"_lc") | str_detect(.,"_uc") | str_detect(.,"_rec") | str_detect(.,"_u") | str_detect(.,"_ec")] %>% paste(collapse="") %>% sub('.','',.),
    drip_type = ifelse(well_behaved, drip_type, "?"),
    drip_type = factor(drip_type, ordered = FALSE, levels=c("st","wk","hy","md","?")),
    suite = ifelse(str_detect(version,"G1e11"),"E",ifelse(str_detect(version,"u0p4"),"B",ifelse(str_detect(version,"ecx1"),"C",ifelse(str_detect(version,"ecx0p2"),"D","A")))),
    # Strength ratios
    str_ratio_LSE = (str_uc_km_MPa + str_lc_km_MPa)/(str_ml_km_MPa + str_uc_km_MPa + str_lc_km_MPa),
    str_ratio_yse = (str_uc_yse + str_lc_yse)/(str_ml_yse + str_uc_yse + str_lc_yse),
    # Density contrasts
    delta_rho_uc_ml = 3300-2750,
    delta_rho_lc_ml = 3300-2900,
    delta_rho_uc_ec = delta_rho_uc_ml+delta_rho,
    delta_rho_lc_ec = delta_rho_lc_ml+delta_rho,
    # Timescales
    # RT (Molnar & Houseman)
    wavelength=drip_width*3*1e3,
    h_L=90e3, # thickness of negatively buoyant part of lithosphere (LC + ML)
    eta_L = 10^(eta_Ec_avg/4 + eta_LC_avg/4 + eta_ML_avg/2), # avg visc. of lower part of lithosphere
    kh = 2*pi/wavelength*h_L,
    qp = (sinh(kh)*cosh(kh)-kh)/(kh*cosh(kh)*cosh(kh)+kh^3), # q prime
    T_RT = (2*eta_L)/((delta_rho+10)*9.81*h_L), # 10 added to account for thermal buoyancy
    tau_RT = T_RT/qp,
    tau_RT_Myr = tau_RT/yr/1e6,
    tau_RT_measured = -t_diff2/log(1/4),
    # McKenzie
    tau_crust = 4*pi*(10^(eta_UC_avg)/2 + 10^(eta_LC_avg)/2 + 10^(eta_ML_avg))/((delta_rho_lc_ec/2 + delta_rho_uc_ec/2)*9.81*wavelength)/yr/1e6, # McKenzie
  )

print_info <- function(d) {
  print("Number of each type of drip:")
  print(c("Strong",sum(d$drip_type=="st")))
  print(c("Weak",sum(d$drip_type=="wk")))
  print(c("Hybrid",sum(d$drip_type=="hy")))
  print(c("Intermediate",sum(d$drip_type=="md")))
  print(c("None",sum(d$drip_type=="?")))
  print(c("Total",sum(d$drip_type!="")))
  print(c("A",sum(d$suite=="A")))
  print(c("B",sum(d$suite=="B")))
  print(c("C",sum(d$suite=="C")))
  print(c("D",sum(d$suite=="D")))
  print(c("E",sum(d$suite=="E")))
}
dat %>% print_info


RT_plot <- function(d) {
  p <- ggplot(d) + 
    scale_x_log10(breaks=log10_breaks) + 
    scale_y_log10(breaks=log10_breaks) +
    scale_color_manual(values=drip_type_colors) +
    scale_fill_manual(values=drip_type_colors) +
    scale_shape_manual(values=c(19,3)) +
    coord_fixed(ratio=1) +
    scale_alpha_identity() +
    guides(fill="none", colour="none", size="none", shape="none") +
    geom_abline(slope=1)+
    geom_point(aes(x=tau_RT_measured, y=tau_RT_Myr, colour=drip_type, fill=drip_type, group=version, shape=!well_behaved, size=2, alpha=0.75)) +
    geom_point(aes(x=tau_RT_measured, y=tau_RT_Myr, alpha=ifelse(label=="",0,1), fill=drip_type), colour="black", size=5, stroke=1, shape=21) 
  return(p)
}
dat %>% filter(suite!="E") %>% RT_plot %T>% ggsave(file="tb_vs_tdiff.png", width=FIG_WIDTH, height=FIG_HEIGHT)


# CLASSIFICATION PLOTS

# All data
timescale_plot <- function(d, time_var_x, time_var_y, colour, factor=0.25, intercept=0,xlims=c(3,40), ylims=c(0.2,40)) {
  p <- ggplot(d) +
    scale_shape_manual(values=c(19,3)) +
    scale_x_log10(breaks=log10_breaks ) +
    scale_y_log10(breaks=log10_breaks) +
    coord_cartesian(xlim=xlims, ylim=ylims) +
    geom_function(fun=function(x) factor*x+intercept, color="black", size=0.75,xlim=c(0,100)) +
    scale_alpha_identity() + 
    scale_colour_manual(values = drip_type_colors) + 
    scale_fill_manual(values=alpha(drip_type_colors,0.75)) + 
    guides(fill="none", colour="none", size="none", shape="none") +
    geom_point(aes(x={{time_var_x}}, y={{time_var_y}}, colour={{colour}}, fill={{colour}}, group=ver, shape=!well_behaved, size=2, alpha=0.75)) +
    geom_point(aes(x={{time_var_x}}, y={{time_var_y}}, alpha=ifelse(label=="",0,1), fill={{colour}}), colour="black", size=5, stroke=1, shape=21) 
  return(p)
}
dat %>% filter(suite!="E") %>% timescale_plot(tau_crust, tau_RT_Myr, drip_type) %T>% ggsave(file="tcrust_vs_tb.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(suite!="E") %>% timescale_plot(tau_crust, tau_RT_measured, drip_type, ylims=c(0.05,20)) %T>% ggsave(file="tcrust_vs_tdiff.png", width=FIG_WIDTH, height=FIG_HEIGHT)

dat %>% filter(suite!="E") %>% timescale_plot(tau_P, tau_RT_Myr, drip_type, xlims=c(1,40)) %T>% ggsave(file="tP_vs_tb.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(suite!="E") %>% timescale_plot(tau_P, tau_RT_measured, drip_type,ylims=c(.05,20),xlims=c(1,40))  %T>% ggsave(file="tP_vs_tdiff.png", width=FIG_WIDTH, height=FIG_HEIGHT)

dat %>% filter(suite=="A") %>% timescale_plot(tau_crust, tau_RT_Myr, drip_type) %T>% ggsave(file="tcrust_vs_tb_A.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(suite=="A") %>% timescale_plot(tau_crust, tau_RT_measured, drip_type,ylims=c(.05,20))  %T>% ggsave(file="tcrust_vs_tdiff_A.png", width=FIG_WIDTH, height=FIG_HEIGHT)

dat %>% filter(suite=="B") %>% timescale_plot(tau_crust, tau_RT_Myr, drip_type) %T>% ggsave(file="tcrust_vs_tb_B.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(suite=="B") %>% timescale_plot(tau_crust, tau_RT_measured, drip_type,ylims=c(.05,20)) %T>% ggsave(file="tcrust_vs_tdiff_B.png", width=FIG_WIDTH, height=FIG_HEIGHT)

dat %>% filter(suite=="C") %>% timescale_plot(tau_crust, tau_RT_Myr, drip_type) %T>% ggsave(file="tcrust_vs_tb_C.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(suite=="C") %>% timescale_plot(tau_crust, tau_RT_measured, drip_type,ylims=c(.05,20)) %T>% ggsave(file="tcrust_vs_tdiff_C.png", width=FIG_WIDTH, height=FIG_HEIGHT)

dat %>% filter(suite=="D") %>% timescale_plot(tau_crust, tau_RT_Myr, drip_type) %T>% ggsave(file="tcrust_vs_tb_D.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(suite=="D") %>% timescale_plot(tau_crust, tau_RT_measured, drip_type,ylims=c(.05,20)) %T>% ggsave(file="tcrust_vs_tdiff_D.png", width=FIG_WIDTH, height=FIG_HEIGHT)

dat %>% filter(suite=="E") %>% timescale_plot(tau_crust, tau_RT_Myr, drip_type, 0.08, xlims=c(7,90)) %T>% ggsave(file="tcrust_vs_tb_E.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(suite=="E") %>% timescale_plot(tau_crust, tau_RT_measured, drip_type, 0.08, ylims=c(.05,10),xlims=c(7,90)) %T>% ggsave(file="tcrust_vs_tdiff_E.png", width=FIG_WIDTH, height=FIG_HEIGHT)


dat %>% filter(suite!="E"&delta_rho==70) %>% timescale_plot(tau_crust, tau_RT_Myr, drip_type) %T>% ggsave(file="tcrust_vs_tb_rec70.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(suite!="E"&delta_rho==70) %>% timescale_plot(tau_crust, tau_RT_measured, drip_type, ylims=c(0.05,20))  %T>% ggsave(file="tcrust_vs_tdiff_rec70.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(suite!="E"&delta_rho==160) %>% timescale_plot(tau_crust, tau_RT_Myr, drip_type) %T>% ggsave(file="tcrust_vs_tb_rec160.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(suite!="E"&delta_rho==160) %>% timescale_plot(tau_crust, tau_RT_measured, drip_type, ylims=c(0.05,20))  %T>% ggsave(file="tcrust_vs_tdiff_rec160.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(suite!="E"&delta_rho==250) %>% timescale_plot(tau_crust, tau_RT_Myr, drip_type) %T>% ggsave(file="tcrust_vs_tb_rec250.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(suite!="E"&delta_rho==250) %>% timescale_plot(tau_crust, tau_RT_measured, drip_type, ylims=c(0.05,20))  %T>% ggsave(file="tcrust_vs_tdiff_rec250.png", width=FIG_WIDTH, height=FIG_HEIGHT)


classify_plot <- function(d, time_var, str_var, colour=drip_type,xlim=c(0.07,30),ylim=c(0.2,1),const_x=Inf, const_y=Inf) {
  p <- ggplot(d) +
    scale_x_log10(breaks=log10_breaks) +
    scale_y_continuous(breaks=seq(0,1,0.2)) +
    coord_cartesian(xlim=xlim, ylim=ylim) +
    scale_alpha_identity() + 
    scale_colour_manual(values=drip_type_colors) + 
    scale_fill_manual(values=alpha(drip_type_colors, 0.75)) +
    scale_shape_manual(values=c(19,3)) +
    geom_abline(slope=0, intercept=const_y, color="black", size=0.75,xlim=c(0,100),linetype="dashed") +
    geom_vline(xintercept=const_x, color="black", size=0.75,xlim=c(0,100), linetype="dashed") +
    guides(fill="none", colour="none", size="none", shape="none") +
    geom_point(aes(x={{time_var}}, y={{str_var}}, colour={{colour}}, fill={{colour}}, group=ver, shape=!well_behaved, size=2, alpha=0.75)) +
    geom_point(aes(x={{time_var}}, y={{str_var}}, alpha=ifelse(label=="",0,1),fill={{colour}}), colour="black", size=5, stroke=1, shape=21) 
  return(p)
}

dat %>% filter(suite!="E") %>% classify_plot(tau_RT_Myr, str_ratio_LSE) %T>% ggsave(file="tb_vs_str_LSE.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(suite!="E") %>% classify_plot(tau_RT_measured, str_ratio_LSE) %T>% ggsave(file="tdiff_vs_str_LSE.png", width=FIG_WIDTH,height=FIG_HEIGHT)
dat %>% filter(suite!="E") %>% classify_plot(tau_RT_Myr, str_ratio_yse, ylim=c(0.35,1),xlim=c(0.2,33),const_x=1.2, const_y=0.85) %T>% ggsave(file="tb_vs_str_yse.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(suite!="E") %>% classify_plot(tau_RT_measured, str_ratio_yse,ylim=c(0.35,1),xlim=c(0.05,30),const_x=1.2, const_y=0.85) %T>% ggsave(file="tdiff_vs_str_yse.png", width=FIG_WIDTH, height=FIG_HEIGHT)

dat %>% filter(suite=="A") %>% classify_plot(tau_RT_Myr, str_ratio_yse) %T>% ggsave(file="tb_vs_str_yse_A.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(suite=="A") %>% classify_plot(tau_RT_measured, str_ratio_yse)  %T>% ggsave(file="tdiff_vs_str_yse_A.png", width=FIG_WIDTH, height=FIG_HEIGHT)

dat %>% filter(suite=="B") %>% classify_plot(tau_RT_Myr,str_ratio_yse) %T>% ggsave(file="tb_vs_str_yse_B.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(suite=="B") %>% classify_plot(tau_RT_measured,str_ratio_yse) %T>% ggsave(file="tdiff_vs_str_yse_B.png", width=FIG_WIDTH, height=FIG_HEIGHT)

dat %>% filter(suite=="C") %>% classify_plot(tau_RT_Myr,str_ratio_yse) %T>% ggsave(file="tb_vs_str_yse_C.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(suite=="C") %>% classify_plot(tau_RT_measured,str_ratio_yse) %T>% ggsave(file="tdiff_vs_str_yse_C.png", width=FIG_WIDTH, height=FIG_HEIGHT)

dat %>% filter(suite=="D") %>% classify_plot(tau_RT_Myr,str_ratio_yse) %T>% ggsave(file="tb_vs_str_yse_D.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(suite=="D") %>% classify_plot(tau_RT_measured,str_ratio_yse) %T>% ggsave(file="tdiff_vs_str_yse_D.png", width=FIG_WIDTH, height=FIG_HEIGHT)

dat %>% filter(suite=="E") %>% classify_plot(tau_RT_Myr,str_ratio_yse) %T>% ggsave(file="tb_vs_str_yse_E.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(suite=="E") %>% classify_plot(tau_RT_measured,str_ratio_yse) %T>% ggsave(file="tdiff_vs_str_yse_E.png", width=FIG_WIDTH, height=FIG_HEIGHT)

# VIOLIN/HISTOGRAM PLOTS
violin_uplift <- function(d) {
  p <- ggplot(d, aes(x=drip_type, y=net_uplift, fill=drip_type, colour=drip_type, alpha=0.75)) +
    scale_fill_manual(values=alpha(drip_type_colors, 0.75)) +
    scale_colour_manual(values=drip_type_colors) + 
    scale_alpha_identity() +
    guides(fill="none", colour="none", size="none", shape="none") +
    geom_violin()
  return(p)
}
dat %>% filter(well_behaved && suite!="E") %>% violin_uplift %T>% ggsave(file="violin_uplift.png", width=FIG_WIDTH, height=FIG_HEIGHT)

violin_size <- function(d) {
  p <- ggplot(d, aes(x=drip_type, y=drip_width, fill=drip_type, colour=drip_type)) +
    scale_fill_manual(values=alpha(drip_type_colors, 0.75)) +
    scale_colour_manual(values=drip_type_colors) + 
    scale_alpha_identity() +
    guides(fill="none", colour="none", size="none", shape="none") +
    geom_violin()
  return(p)
}
dat %>% filter(well_behaved && suite!="E") %>% violin_size %T>% ggsave(file="violin_size.png", width=FIG_WIDTH, height=FIG_HEIGHT)

histogram_size <- function(d) {
  p <- ggplot(d, aes(x=drip_width, fill=drip_type, colour=drip_type)) +
    scale_fill_manual(values=alpha(drip_type_colors, 0.75)) +
    scale_colour_manual(values=drip_type_colors) + 
    scale_alpha_identity() +
    guides(fill="none", colour="none", size="none", shape="none") +
    geom_histogram(binwidth=15)
  return(p)
}
dat %>% filter(well_behaved && suite!="E") %>% histogram_size %T>% ggsave(file="histogram_size.png", width=FIG_WIDTH, height=FIG_HEIGHT)

histogram_eta_ML <- function(d) {
  p <- ggplot(d, aes(x=eta_ML_top, fill=drip_type, colour=drip_type)) +
    scale_fill_manual(values=alpha(drip_type_colors, 0.85)) +
    scale_colour_manual(values=drip_type_colors) + 
    scale_alpha_identity() +
    guides(fill="none", colour="none", size="none", shape="none") +
    coord_cartesian(ylim=c(0,60),expand=c(0,0)) +
    geom_histogram()
  return(p)
}
dat %>% filter(well_behaved && suite!="E") %>% histogram_eta_ML %T>% ggsave(file="histogram_eta_ML.png", width=FIG_WIDTH, height=FIG_HEIGHT)

histogram_eta_ratio <- function(d) {
  p <- ggplot(d, aes(x=10^(eta_UC_avg/2 + eta_LC_avg/2)/(10^eta_ML_avg), fill=drip_type, colour=drip_type)) +
    scale_fill_manual(values=alpha(drip_type_colors, 0.75)) +
    scale_colour_manual(values=drip_type_colors) + 
    scale_alpha_identity() +
    guides(fill="none", colour="none", size="none", shape="none") +
    geom_histogram(binwidth=10)
  return(p)
}
dat %>% filter(well_behaved && suite!="E") %>% histogram_eta_ratio %T>% ggsave(file="histogram_eta_ratio.png", width=FIG_WIDTH, height=FIG_HEIGHT)

histogram_delta_rho <- function(d) {
  p <- ggplot(d, aes(x=delta_rho, fill=drip_type, colour=drip_type)) +
    scale_fill_manual(values=alpha(drip_type_colors, 0.75)) +
    scale_colour_manual(values=drip_type_colors) + 
    scale_alpha_identity() +
    guides(fill="none", colour="none", size="none", shape="none") +
    geom_histogram(binwidth=50) 
  return(p)
}
dat %>% filter(well_behaved && suite!="E") %>% histogram_delta_rho %T>% ggsave(file="histogram_drho.png", width=FIG_WIDTH, height=FIG_HEIGHT)

# TOPGRAPHIC CHANGE
topo_plot <- function(d) {
  p <- ggplot() + scale_x_log10(breaks=c(1.3, 1.5, 2.1, 3.0, 4.5, 6.5, 9.0, 12.0))
  for(i in 1:nrow(d)) {
    d_row <- d[i,]
    times <- d_row$times[[1]]
    depths <- d_row$depths[[1]]
    depths[depths==13] <- NaN
    topos <- d_row$t[[1]]
    groups <- d_row$drip_type
    # use aes_ instead of aes
    p <- p + geom_line(aes_(x=depths, y=topos, color=groups, alpha=0.65)) +
      scale_color_manual(values=drip_type_colors) +
      guides(fill="none", colour="none", size="none", shape="none") +
      coord_cartesian(expand=c(0,0)) +
      scale_alpha_identity()
  }
  return (p + xlab(latex2exp::TeX(r"($z'$)")) + ylab(latex2exp::TeX(r"($\Delta z$)")))
}

dat %>% filter(well_behaved & suite!="E") %>% topo_plot %T>% ggsave(file="topo.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(well_behaved & suite!="E" & (drip_type=="wk" | drip_type=="st")) %>% topo_plot %T>% ggsave(file="topo_wk_st.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(well_behaved & suite=="A") %>% topo_plot %T>% ggsave(file="topo_A.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(well_behaved & suite=="B") %>% topo_plot %T>% ggsave(file="topo_B.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(well_behaved & suite=="C") %>% topo_plot %T>% ggsave(file="topo_C.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(well_behaved & suite=="D") %>% topo_plot %T>% ggsave(file="topo_D.png", width=FIG_WIDTH, height=FIG_HEIGHT)
dat %>% filter(well_behaved & suite=="E") %>% topo_plot %T>% ggsave(file="topo_E.png", width=FIG_WIDTH, height=FIG_HEIGHT)

# version, time shift, absolute elevation
comparisons <- list(
  # Arizaro comp
  list("A1_G50_ch2_cs0p9_dmg1_ecx0p5_fr15_fs0p9_k2p5_lcx1_lin0_mlx0p2_mtz0_rec250_ruc2750_u0p8_ucx0p2", -37, FALSE),
  # Bidahochi comp
  list("A1_G50_ch2_cs0p9_dmg1_ecx0p5_fr15_fs0p9_k2p5_lcx1_lin0_mlx0p2_mtz0_rec70_ruc2750_u0p4_ucx0", -61, FALSE),
  # Altiplano comp
  list("A1_G50_ch2_cs0p9_dmg1_ecx0p5_fr15_fs0p9_k2p5_lcx0p5_lin0_mlx0p5_mtz0_rec160_ruc2750_u0p8_ucx0p8", -35, TRUE)
)

comparisons_exx <- list(
  list("A1_G50_ch2_cs0p9_dmg1_ecx0p5_fr15_fs0p9_k2p5_lcx0_lin0_mlx0p7_mtz0_rec70_ruc2750_u0p4_ucx0p2", -84, TRUE)
)

p1 <- ggplot() +  
  scale_colour_identity() +
  geom_line(aes(x=-arizaro_tectonic_subs_t_Ma,y= -arizaro_tectonic_subs_depth_m/1000, colour="tomato"), linetype="41",linewidth=2.5,alpha=0.33) + 
  geom_line(aes(x=-altiplano_elev_t_Ma, y=altiplano_elev_km, colour="mediumorchid"), linetype="41", linewidth=2.5,alpha=0.33) + 
  geom_line(aes(x=-bidahochi_total_subs_t_Ma, y=-bidahochi_total_subs_depth_m/1000, colour="green3"), linetype="41",linewidth=2.5,alpha=0.33) +
  scale_alpha_identity()+
  new_scale_colour() + 
  scale_colour_manual(values=drip_type_colors) +
  scale_fill_manual(values=c("black","red")) +
  scale_y_continuous(breaks=seq(-2,6,0.5), limits=c(-2,4.5), expand=c(0,0)) +
  scale_x_continuous(breaks=seq(-100,0,5),limits=c(-30,0), expand=c(0,0)) +
  scale_linewidth_manual(values=c(0.75, 1.5)) +
  labs(x="Time (Myr)", y="Elev. (km)") +
  theme(legend.position = "none")
  
for(comp in comparisons){
  v = comp[[1]]
  tshift = comp[[2]]
  is_elev = comp[[3]]
  row = dat[dat$version==v, ]

  fname <- file.path("..","plots","deflections",paste0(v,".csv"))
  deflection = read_csv(fname,show_col_types = FALSE)
  
  extending_idx = row$depths[[1]]==row$m_extending
  time_extending = row$times[[1]][extending_idx]
  if(length(time_extending)==0) {
    time_extending = 1000 #Myr
  } else {
    time_extending = time_extending+0
  }
  
  in_time = deflection$t <= row$t_detach+2.5 & deflection$t >= row$t0-25
  time = deflection$t[in_time]
  elev = deflection$middle[in_time]
  depth = deflection$middle[in_time] - deflection$sides[in_time]
  topo = if (is_elev) elev else depth
  print(time_extending)
  p1 <- p1 + 
    geom_line(aes(x=!!time+!!tshift, y=!!topo, colour=!!row$drip_type),linewidth=0.75, linetype="solid") +
    geom_line(aes(x=!!time[!!time >= time_extending]+!!tshift, y=!!topo[!!time >= time_extending], colour=!!row$drip_type),linewidth=1.5, linetype="solid") +
    geom_point(aes(x=!!row$t_detach+!!tshift, y=!!topo[time==row$t_detach], colour=!!row$drip_type),size=3,shape=21,fill="white")
  
}

p2 <- ggplot() +  
  scale_x_continuous(breaks=seq(-30,0,5),limits=c(-30,0),expand=c(0,0)) +
  scale_y_continuous(breaks=seq(-3,3,1),limits=c(-2.5,2.5),expand=c(0,0)) +
  scale_linewidth_manual(values=c(0.75, 1.5)) +
  scale_colour_manual(values=drip_type_colors) +
  theme(legend.position = "none") +
  theme(axis.title.x=element_blank(), 
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.label.x.=element_blank()) +
  xlab(NULL) + ylab("E_xx")

for(comp in comparisons_exx) {
  
  v = comp[[1]]
  tshift = comp[[2]]
  row = dat[dat$version==v, ]
  print(row$exx_middle)
  exx0=1e-15
  fname <- file.path("..","plots","deflections",paste0(v,".csv"))

  deflection = read_csv(fname,show_col_types = FALSE)
  in_time = deflection$t <= row$t_detach+5 & deflection$t >= row$t0-2.5
  time = deflection$t[in_time] 
  exx = deflection$exx[in_time]
  p2 <- p2 + 
    geom_line(aes(x=!!time+!!tshift, y=!!exx/exx0, colour=!!row$drip_type), linewidth=0.75, linetype="solid") +
    geom_line(aes(x=!!time[exx>5e-17]+!!tshift, y=!!exx[exx>5e-17]/exx0, colour=!!row$drip_type), linewidth=1.5, linetype="solid") +
    geom_point(aes(x=!!row$t_detach+!!tshift, y=!!exx[time==row$t_detach]/exx0, colour=!!row$drip_type),size=3,shape=21,fill="white")
  
}

# Layout: p1 and p2 side-by-side, p3 spans the entire bottom row
p <- p2 / p1 +
  plot_layout(heights = c(1, 3))
p
ggsave(p1, file="real_world_comparison.png", width=6.5, height=4)
ggsave(p, file="real_world_comparison2.png", width=6.5, height=4)

