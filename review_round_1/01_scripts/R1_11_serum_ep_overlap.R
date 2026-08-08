# =====================================================================
# R1_11_serum_ep_overlap.R — serum/EP-overlapp med RIKTIG nevner
# Adresserer R5 punkt 2 («8 overlapping proteins (1.7%) has no
# interpretable denominator»), R2 concern 2
#
# Feilen i 201_Olink_vs_ECV.R: to signifikanslister fra to plattformer
# krysses direkte, uten å begrense til analytter som faktisk er MÅLBARE
# på begge. Nevneren blir da meningsløs.
#
# Her: overlapp beregnes innenfor de proteinene som er kvantifisert på
# begge plattformer, og suppleres med en terskelfri sammenligning av
# effektestimater — som ikke avhenger av signifikansgrenser i det hele tatt.
# Krever: R1_01_longterm.R
# =====================================================================
suppressPackageStartupMessages({library(data.table)})
source(file.path(Sys.getenv("ECV_ROOT", "/Users/sindrle/Research/Prosjekter/ECV"),
                 "review_round_1/01_scripts/R1_config.R"))
ROOT <- ECV_ROOT
OUT  <- file.path(ROOT, "review_round_1/02_output"); setwd(ROOT)
source(file.path(ROOT, "review_round_1/01_scripts/R1_theme.R"))

EP <- unique(readRDS(file.path(OUT,"R1_01_longterm.rds"))$PRIM[
  !is.na(Genes) & Genes != "", .(Genes, est_ep = Estimate, p_ep = p, q_ep = q)], by = "Genes")
SE <- fread(file.path(OLINK_DIR,"Exercise/MyoGlu_Olink_A1B1_All_Annotated.csv"), sep = ";")
SE <- unique(SE[!is.na(GeneName) & GeneName != "", .(Genes = GeneName, est_se = est, p_se = p, q_se = q)], by = "Genes")

cat(sprintf("EP-MS kvantifisert (etter filter): %d gener\n", nrow(EP)))
cat(sprintf("Serum-Olink målt              : %d gener\n", nrow(SE)))
both <- merge(EP, SE, by = "Genes")
cat(sprintf("MÅLBARE PÅ BEGGE              : %d gener  <- riktig nevner\n\n", nrow(both)))

sig_ep <- EP[q_ep < .05, Genes]; sig_se <- SE[q_se < .05, Genes]
ov <- intersect(sig_ep, sig_se)
cat("=== SLIK MANUSKRIPTET REGNER (hele lister, ulike nevnere) ===\n")
cat(sprintf("  EP sig %d | serum sig %d | overlapp %d = %.1f%% av EP-listen\n\n",
            length(sig_ep), length(sig_se), length(ov), 100*length(ov)/max(length(sig_ep),1)))

b <- both[, .(Genes, est_ep, q_ep, est_se, q_se)]
e_in <- b[q_ep < .05, Genes]; s_in <- b[q_se < .05, Genes]; o_in <- intersect(e_in, s_in)
exp_ov <- length(e_in) * length(s_in) / nrow(b)
cat("=== RIKTIG NEVNER: kun gener målbare på begge ===\n")
cat(sprintf("  av %d felles: EP-sig %d, serum-sig %d, overlapp %d\n",
            nrow(b), length(e_in), length(s_in), length(o_in)))
cat(sprintf("  forventet overlapp ved uavhengighet: %.1f  -> %.1fx anrikning (Fisher p = %.3g)\n",
            exp_ov, length(o_in)/max(exp_ov,1e-9),
            fisher.test(matrix(c(length(o_in), length(e_in)-length(o_in),
                                 length(s_in)-length(o_in),
                                 nrow(b)-length(e_in)-length(s_in)+length(o_in)), 2))$p.value))

cat("\n=== TERSKELFRITT (uavhengig av signifikansgrenser) ===\n")
ct <- cor.test(b$est_ep, b$est_se, method = "spearman")
cat(sprintf("  effektestimater EP vs serum: Spearman rho = %+.3f (p = %.3g), samme fortegn %.0f%%\n",
            ct$estimate, ct$p.value, 100*mean(b$est_ep*b$est_se > 0)))
cat("  Tolkning: rho nær null betyr at de to plattformene ikke måler samme\n  endring, og at et lavt listeoverlapp da er forventet uten biologi.\n")
fwrite(b, file.path(OUT, "R1_11_serum_ep_shared.csv"))

# ---- figur: effektestimater, EP vs serum ----
b[, klasse := fifelse(q_ep < .05 & q_se < .05, "Begge",
              fifelse(q_ep < .05, "Kun EP", fifelse(q_se < .05, "Kun serum", "Ingen")))]
p1 <- ggplot(b, aes(est_se, est_ep)) +
  geom_hline(yintercept = 0, linetype = "dotted") +
  geom_vline(xintercept = 0, linetype = "dotted") +
  geom_point(aes(colour = klasse), size = 0.9, alpha = 0.85) +
  scale_colour_manual(values = c("Ingen" = "grey75", "Kun EP" = cols[1],
                                 "Kun serum" = cols[2], "Begge" = "black")) +
  labs(x = expression(paste(Delta, "serum protein (NPX)")),
       y = expression(paste(Delta, "EP protein, log"[2], "(int)"))) +
  theme_ecv() + theme(legend.position = "right")
save_ecv("R1_11_EP_vs_serum_effects.pdf", p1, width = 12, height = 8)
