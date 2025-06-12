library(criteval)
seed <- 1
min.node.size <- 10
.n.skip <- 2
.grid.size <- 250
.tol.integrate <- 1e-3

n <- 2000
rho.W <- 0.75
kappa.ratio <- 1
q.W <- 1
p.X <- 3
sigma.eps <- 1

s <- 50
theta_FUN_list <- list(
  theta1 = function(x) plogis((x[,1] - 0.23), scale = 1/s),
  theta2 = function(x) plogis(-(x[,1] - 0.27), scale = 1/s),
  theta3 = function(x) 0.5 * plogis(x[,2] - 0.75, scale = 1/s),
  theta3 = function(x) 0.25 * plogis((x[,3] - 0.73), scale = 1/s),
  theta4 = function(x) 0.25 * plogis((x[,3] - 0.77), scale = 1/s)
)
nu_FUN <- function(x) {0 * x[,1]}

#--------------------------------------------------
#----- Start
#--------------------------------------------------
set.seed(1)

out <- run_criteval_sim(
  n = n,
  rho.W = rho.W,
  #R.W = R,
  kappa.ratio = kappa.ratio,
  q.W = q.W,
  p.X = p.X,
  sigma.eps = sigma.eps,
  theta_FUN_list = theta_FUN_list,
  nu_FUN = nu_FUN,
  min.node.size = min.node.size,
  .n.skip = .n.skip,
  .grid.size = .grid.size,
  .tol.integrate = .tol.integrate,
  seed = seed
)

#--------------------------------------------------
#----- Plot tests
#--------------------------------------------------
library(ggplot2)
df.res <- out$results

MY.ALPHA <- c(
  1, # Delta^*
  1, # Delta_V^*
  1, # Delta
  1, # Delta_V
  0, # GRF-grad
  0  # GRF-fpt
)
MY.COL <- c("#F8766D", "#00BFC4", "#F8766D", "#00BFC4", "#F8766D", "#00BFC4")
MY.LTY <- c("21", "21", "solid", "solid", "43", "43")
MY.LWD <- c(0.5, 0.5, 0.5, 0.5, 0.5, 0.5)

DIM.LABS <- setNames(
  paste0("Candidate~feature~italic(X)[", sort(unique(df.res$dim)), "]"),
  nm = sort(unique(df.res$dim))
)
CRIT.LABS <- list(
  "Delta.pop"  = expression(Delta^{"*"}),
  "Delta.V.pop" = expression({Delta[italic(V)]}^{"*"}),
  "Delta"     = expression(Delta),
  "Delta.V"    = expression(Delta[italic(V)]),
  "Delta.grad" = expression(widetilde(Delta)^{grad}),
  "Delta.fpt" = expression(widetilde(Delta)^{FPT})
)

df.max <- df.res %>%
  filter(normalized != "Normalized") %>%
  group_by(criterion) %>%
  dplyr::slice_max(value, n = 1, with_ties = FALSE) %>%
  ungroup()
df.max2 <- df.max %>%
  mutate(normalized = "Normalized") %>%
  bind_rows(df.max)


df.res %>%
  ggplot(aes(x = threshold, y = value, color = criterion, alpha = criterion,
             linewidth = criterion, linetype = criterion)) +
  labs(
    x = "Candidate threshold",
    y = "Criterion value",
    title = "Candidate split criterion values"
  ) +
  geom_hline(yintercept = 0, lwd = 1, col = "gray85") +
  geom_vline(
    data       = df.max2, #%>% filter(normalized == "Raw"),
    mapping    = aes(
      xintercept = threshold,
      color      = criterion,
      linetype   = criterion,
      alpha      = criterion
    ),
    linewidth = 0.5,
    show.legend = F
  ) +
  geom_line() +
  facet_grid(
    rows = vars(normalized),
    cols = vars(dim),
    scales = "free_y",
    labeller = labeller(
      normalized = label_value,
      dim        = as_labeller(DIM.LABS, label_parsed)
    )
  ) +
  scale_x_continuous(limits = c(0, 1)) +
  scale_color_manual(values = MY.COL, breaks = names(CRIT.LABS), labels = CRIT.LABS) +
  scale_linetype_manual(values = MY.LTY, breaks = names(CRIT.LABS), labels = CRIT.LABS) +
  scale_linewidth_manual(values = MY.LWD, breaks = names(CRIT.LABS), labels = CRIT.LABS) +
  scale_alpha_manual(values = MY.ALPHA, breaks = names(CRIT.LABS), labels = CRIT.LABS) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8, color = "gray50"),
    axis.ticks = element_blank(),
    axis.title = element_text(size =  11, color = "gray30"),
    axis.title.y = element_blank(),
    legend.key.width = unit(0.75, "cm"),
    legend.margin = margin(t = 0, r = 0, b = 0, l = -20),
    legend.text = element_text(size = 14),
    legend.title = element_blank(),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 5),
    plot.title = element_text(size = 11),
    strip.background = element_blank(),
    strip.text = element_text(size = 11, color = "gray30")
  )
