#!/usr/bin/env Rscript

# setwd("~/Projects/sasa_epi/methylUtil")

for (p in c("ggplot2", "ggbeeswarm", "cowplot", "nlme", "tidyverse")) {
    if (!suppressMessages(require(p, character.only = T))) {
        message(paste("Installing:", p))
        install.packages(p, repos = "https://mirror.its.dal.ca/cran", dependencies = T)
        suppressMessages(require(p, character.only = T))}
    rm(p)
}

theme_adjustments <- theme_linedraw() + theme(axis.text = element_text(size = 10, colour = "black"), 
                                              axis.title = element_text(size = 12, colour = "black"),
                                              panel.grid = element_blank(), panel.border = element_rect(colour = "black"))

x <- read.table("lab-epigenetics_phenos.txt", header = T)

x <- x %>% mutate(
    CF = 10^5 * WT / FL ^3,
    Crosstype = case_when(
        Crosstype == "SAS" ~ "SAS",
        Crosstype == "Wild" ~ "Wild",
        TRUE ~ NA_character_),
    FamilyID = str_c(Crosstype, str_remove(FatherID, ".*_"))
    )
# x <- x %>% filter(Crosstype != "SASxWild")
# x <- x %>% filter(FamilyID != "HOR5" & FamilyID != "NOR7")

cols <- c("goldenrod", paste0("goldenrod", 1:4), "dodgerblue", paste0("dodgerblue", 1:4))
names(cols) <- c(paste0("SAS", 1:5), paste0("Wild", c(1:2, 5 , 7:8)))

ggplot(x, aes(x=FL)) +
    geom_histogram(mapping = aes(fill=Crosstype), show.legend = FALSE) +
    geom_density(mapping = aes(FL, stat(count))) +
    facet_wrap(facets = ~FatherID, nrow = 2)

fl <- ggplot(x, aes(y=FL, x = Crosstype)) +
    geom_boxplot(outlier.shape = NA) +
    geom_beeswarm(aes(colour = FamilyID), cex = 1.2, show.legend = FALSE, size = 0.7, dodge.width = 0.5) +
    #geom_quasirandom(aes(colour = FamilyID)) +
    scale_color_manual(values = cols) +
    ylab("Fork length (mm)") +
    xlab("") +
    theme_adjustments
# summary(lm(FL~Crosstype,data=x))
summary(lme(FL ~ Crosstype, random = ~1|FamilyID, data = x))
anova.lme(lme(FL ~ Crosstype, random = ~1|FamilyID, data = x))
# car::Anova(lm(FL ~ Crosstype + FamilyID, data = x), type = "II")
# summary(aov(FL ~ Crosstype + Error(FamilyID), data = x))
x %>% group_by(Crosstype) %>% summarise(WTm = mean(WT), 
                                        WTsd = sd(WT),
                                        FLm = mean(FL),
                                        FLsd = sd(FL),
                                        CFm = mean(CF),
                                        CFsd = sd(CF))

cf <- ggplot(x, aes(y=CF, x = Crosstype)) +
    geom_boxplot(outlier.shape = NA) +
    geom_beeswarm(aes(colour = FamilyID), cex = 1.2, show.legend = FALSE, size = 0.7, dodge.width = 0.5) +
    #geom_quasirandom(aes(colour = FamilyID)) +
    scale_color_manual(values = cols) +
    ylab(expression("Condition Factor (10"^5 *"*W/L"^3*")")) +
    xlab("") +
    theme_adjustments
# summary(lm(CF~Crosstype,data=x))
summary(lme(CF ~ Crosstype, random = ~1|FamilyID, data = x))
anova.lme(lme(CF ~ Crosstype, random = ~1|FamilyID, data = x))

wt <- ggplot(x, aes(y=WT, x = Crosstype)) +
    geom_boxplot(outlier.shape = NA) +
    geom_beeswarm(aes(colour = FamilyID), cex = 1.2, show.legend = FALSE, size = 0.7, dodge.width = 0.5) +
    #geom_quasirandom(aes(colour = FamilyID)) +
    scale_color_manual(values = cols) +
    ylab("Weight (g)") +
    xlab("") +
    theme_adjustments
# summary(lm(WT~Crosstype,data=x))
summary(lme(WT ~ Crosstype, random = ~1|FamilyID, data = x))
anova.lme(lme(WT ~ Crosstype, random = ~1|FamilyID, data = x))


pg <- plot_grid(fl, wt, cf, labels = c("A", "B", "C"), ncol = 1)

save_plot("Figure3_juvenile-fl-wt.pdf", pg, base_width = 3.25, base_height = 9)

# ggplot(x, aes(x=FatherID, y = CF, fill = Crosstype)) +
#     geom_boxplot() +
#     geom_violin()
# anova(lm(CF~Crosstype,data=x))
# 
# x %>% summarise(WT_m = mean(WT), FL_m = mean(FL), CF_m = mean(CF))
# 
# y  <- x %>% mutate(CF_dev = abs(CF - 1.205))
# 
# y %>% group_by(FatherID) %>% sample_n(size = 8, weight = 1/CF_dev)
# 
