

# Calculate the line between the studies
d <- preplot %>%
  filter(Study != "C", !grepl("non-remiters", remission)) %>%
  group_by(Target, remission, Location) %>%
  mutate(
    d = max(meanAU) - min(meanAU),
    distU = (max(ymin) - min(ymax))/2,
    orig = max(ymin) - distU/10, # Make them not overlap witht the error bars.
    final = min(ymax) + distU/10,
    center = distU + final)

# Keep just one every two lines
empty <- seq(from = 1, to = nrow(d), by = 2)
d <- d[empty, ]


### W0 ####

# Compare studies by gene at w0 indepdendently of remission

preplot <- dff %>%
  group_by(Target, Study, remission, Time, General_location) %>%
  summarise(meanAU = mean(AU), sem = sd(AU)/sqrt(n())) %>%
  mutate(ymax = meanAU + sem, ymin = meanAU - sem,
         label = paste(Study, `remission`))

loc <- "colon"

pd <- position_dodge2(width = 0.5)
preplot %>%
  filter(Study != "C", Time == "w0", General_location == loc) %>%
  ggplot(aes(Study, meanAU)) +
  geom_point(aes(col = Study, group = Study, shape = remission),
             position = pd) +
  geom_errorbar(aes(col = Study, ymin = ymin, ymax = ymax), width = 0.2,
                position = pd)  +
  facet_wrap(~Target, scales = "free_y") +
  geom_hline(data = filter(preplot, Study == "C", General_location == loc),
             aes(yintercept = ymin), linetype = "dotted") +
  geom_hline(data = filter(preplot, Study == "C", General_location == loc),
             aes(yintercept = ymax), linetype = "dotted") +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black")) +
  labs(x = element_blank(), title = paste("Basal expression:", loc),
       subtitle = "Comparison between UPA and TNF")

preplot <- dff %>%
  group_by(Target, Study, General_location, remission, Time) %>%
  summarise(meanAU = mean(AU), sem = sd(AU)/sqrt(n())) %>%
  mutate(ymax = meanAU + sem, ymin = meanAU - sem)

st <- "UPA"

preplot %>%
  filter(Study != "C", General_location == loc, Study == st) %>%
  mutate(label = paste(Study, remission)) %>%
  ggplot(aes(Time, meanAU)) +
  # geom_point(aes(col = Study, shape = remission), position = pd) +
  geom_pointrange(aes(col = Study, ymin = ymin, ymax = ymax, shape = remission), position = pd) +
  geom_line(aes(group = label, col = Study), position = pd) +
  facet_wrap(~Target, scales = "free_y", ncol = 4) +
  geom_hline(data = filter(preplot, Study == "C", General_location == loc),
             aes(yintercept = ymin), linetype = "dotted") +
  geom_hline(data = filter(preplot, Study == "C", General_location == loc),
             aes(yintercept = ymax), linetype = "dotted") +
  theme_classic() +
  theme(strip.background = element_blank(),
        panel.border = element_rect(fill = NA, colour = "black")) +
  labs(x = element_blank(), title = "Expression of marker genes",
       subtitle = paste(loc, "at baseline and at week 14"))


preplot %>%
  filter(Study != "C", Target %in% c("AQP7", "DERL3", "OSM", "S100A8")) %>%
  ggplot(aes(remission, meanAU)) +
  geom_point(aes(col = Study, group = Study)) +
  expand_limits(y = 0) +
  geom_errorbar(aes(col = Study, group = Study, ymin = ymin, ymax = ymax), width = 0.2) +
  geom_line(aes(col = Study, group = Study)) +
  geom_hline(data = filter(preplot, Study == "C",
                           Target %in% c("AQP7", "DERL3", "OSM", "S100A8")),
             aes(yintercept = ymin), linetype = "dotted") +
  geom_hline(data = filter(preplot, Study == "C", Target %in% c("AQP7", "DERL3", "OSM", "S100A8")),
             aes(yintercept = ymax), linetype = "dotted") +
  geom_segment(data = filter(dw,  Target %in% c("AQP7", "DERL3", "OSM", "S100A8")),
               aes(x = remission, y = orig, xend = remission, yend = final)) +
  geom_text(data = filter(dw,  Target %in% c("AQP7", "DERL3", "OSM", "S100A8")),
            aes(x = remission, y = center), label = "\t*", size = 5) + # Dirty trick to dodge the symbol
  geom_text(data = filter(db,  Target %in% c("AQP7", "DERL3", "OSM", "S100A8")),
            aes(x = remission, y = meanAU), label = "\t+", size = 3) + # Dirty trick to dodge the symbol
  facet_wrap(Target ~ Location, scales = "free_y", ncol = 4) +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  geom_segment(aes(xend = remission, yend = meanAU, group = remission)) +
  ylab("AU (mean\u00B1SEM)")
