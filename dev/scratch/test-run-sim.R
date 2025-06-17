library(criteval)
seed <- 1
min.node.size <- 10
.n.skip <- 5
.grid.size <- 250
.tol.integrate <- 1e-3

n <- 1000
rho.W <- 0.7
kappa.ratio <- 1
q.W <- 1 # irrelevant if kappa.ratio = 1
p.X <- 1
sigma.eps <- 0

s <- 30
# theta_FUN_list <- list(
#   theta1 = function(x) plogis((x[,1] - 0.23), scale = 1/s),
#   theta2 = function(x) plogis(-(x[,1] - 0.27), scale = 1/s),
#   theta3 = function(x) 0.5 * plogis(x[,2] - 0.75, scale = 1/s),
#   theta4 = function(x) 0.25 * plogis((x[,3] - 0.73), scale = 1/s),
#   theta5 = function(x) 0.25 * plogis((x[,3] - 0.77), scale = 1/s)
# )
theta_FUN_list <- list(
  theta1 = function(x) plogis((x[,1] - 0.23), scale = 1/s),
  theta2 = function(x) plogis(-(x[,1] - 0.27), scale = 1/s),
  theta3 = function(x) 0.5 * plogis((x[,1] - 0.75), scale = 1/s)
)
# theta_FUN_list <- list(
#   theta1 = function(x) plogis((x[,1] - 0.48), scale = 1/s),
#   theta2 = function(x) plogis((x[,1] - 0.52), scale = 1/s)
# )
nu_FUN <- function(x) {0 * x[,1]}

#--------------------------------------------------
#----- Start
#--------------------------------------------------
set.seed(1)

out <- run_criteval_sim(
  n = n,
  rho.W = rho.W, #R.W = R,
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
#----- Plots: Criteria
#--------------------------------------------------
library(ggplot2)
library(dplyr)
library(tidyr)
df.res <- out$results

CRIT.ALPHA <- c(
  "Delta.pop"   = 1, # Delta^*
  "Delta.V.pop" = 1, # Delta_V^*
  "Delta"       = 1, # Delta
  "Delta.V"     = 1, # Delta_V
  "Delta.grad"  = 0, # GRF-grad
  "Delta.fpt"   = 0  # GRF-fpt
)

CRIT.COL <- c("#F8766D", "#00BFC4", "#F8766D", "#00BFC4", "#F8766D", "#00BFC4")
CRIT.LTY <- c("21", "21", "solid", "solid", "43", "43")
CRIT.LWD <- c(0.7, 0.7, 0.7, 0.7, 0.7, 0.7)

CRIT.LABS <- list(
  "Delta.pop"   = expression(Delta^{"*"}),
  "Delta.V.pop" = expression({Delta[italic(V)]^{"*"}}),
  "Delta"       = expression(Delta),
  "Delta.V"     = expression(Delta[italic(V)]),
  "Delta.grad"  = expression(widetilde(Delta)^{grad}),
  "Delta.fpt"   = expression(widetilde(Delta)^{FPT})
)
DIM.LABS <- setNames(
  paste0("Feature~italic(X)[", 1:out$dgp$params$p.X, "]"),
  nm = sort(unique(df.res$dim))
)

df.res.visible <- df.res %>%
  filter(criterion %in% names(which(CRIT.ALPHA > 0)))

df.max <- df.res.visible %>%
  filter(normalized != "Normalized") %>%
  group_by(criterion) %>%
  dplyr::slice_max(value, n = 1, with_ties = FALSE) %>%
  ungroup()
df.max2 <- df.max %>%
  mutate(normalized = "Normalized") %>%
  bind_rows(df.max)

plt.res <- df.res.visible %>%
  ggplot(aes(
    x = threshold,
    y = value,
    color = criterion,
    alpha = criterion,
    linewidth = criterion,
    linetype = criterion
  )) +
  labs(
    x = "Threshold",
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
  scale_color_manual(values = CRIT.COL, breaks = names(CRIT.LABS), labels = CRIT.LABS) +
  scale_linetype_manual(values = CRIT.LTY, breaks = names(CRIT.LABS), labels = CRIT.LABS) +
  scale_linewidth_manual(values = CRIT.LWD, breaks = names(CRIT.LABS), labels = CRIT.LABS) +
  scale_alpha_manual(values = CRIT.ALPHA, breaks = names(CRIT.LABS), labels = CRIT.LABS) +
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

#--------------------------------------------------
#----- Plots: Example theta functions
#--------------------------------------------------

x0.plot <- matrix(seq(0, 1, length.out = 250), ncol = 1)
X.plot <- kronecker(diag(1, out$dgp$params$p.X, out$dgp$params$p.X), x0.plot)
colnames(X.plot) <- colnames(out$data$X)
theta.plot <- out$dgp$theta_FUN(X.plot)

df.theta <- data.frame(X.plot, theta.plot) %>%
  pivot_longer(
    cols = colnames(theta.plot),
    names_to = "theta.feature",
    values_to = "theta"
  ) %>%
  pivot_longer(
    cols = colnames(X.plot),
    names_to = "X.feature",
    values_to = "X"
  ) %>%
  mutate(
    theta.feature = factor(theta.feature),
    X.feature = factor(X.feature)
  )

THETA.COLORS <- scales::hue_pal()(ncol(theta.plot))
THETA.LABS <- setNames(
  lapply(seq_len(ncol(theta.plot)), function(i) {
    call <- substitute({theta[ii]^"*"}({italic(x)}), list(ii = i))
    as.expression(call)
  }),
  nm = colnames(theta.plot)
)
X.LABS <- setNames(DIM.LABS, nm = colnames(X.plot))

#PLT.THETA.AXIS.TITLE.X <- expression(italic(x)[italic(ℓ)])
#PLT.THETA.AXIS.TITLE.Y <- expression(Effect~{theta[italic(k)]^"*"}(italic(x))~of~italic(W)[italic(k)]~on~italic(Y))
PLT.THETA.AXIS.TITLE.X <- expression(italic(x))
PLT.THETA.AXIS.TITLE.Y <- expression(Effect~of~italic(W)[italic(k)]~on~italic(Y))

plt.theta <- df.theta %>%
  #filter(X.feature == "X1") %>%
  ggplot(aes(x = X, y = theta, color = theta.feature)) +
  labs(
    x = PLT.THETA.AXIS.TITLE.X,
    y = PLT.THETA.AXIS.TITLE.Y,
    title = "Effect function curves (other dimensions of x held at zero)"
  ) +
  facet_wrap(
    vars(X.feature),
    labeller = as_labeller(X.LABS, label_parsed)
  ) +
  geom_hline(yintercept = 0, lwd = 1, col = "gray85") +
  geom_line(linewidth = 1) +
  scale_x_continuous(limits = c(0, 1)) +
  #scale_y_continuous(position = "right") +
  # scale_y_continuous(
  #   sec.axis = sec_axis(
  #     transform = ~ ., # Use the same transformation as the primary axis
  #     name = PLT.THETA.AXIS.TITLE.Y,
  #     breaks = NULL, # Hide tick marks
  #     labels = NULL  # Hide text labels
  #   )
  # ) +
  scale_color_manual(values = THETA.COLORS, breaks = names(THETA.LABS), labels = THETA.LABS) +
  theme_bw() +
  theme(
    axis.text = element_text(size = 8, color = "gray50"),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 10, color = "gray30"),
    #axis.title.y.left = element_blank(),
    legend.title = element_blank(),
    legend.margin = margin(t = -10, r = 0, b = 0, l = 0),
    legend.position = "bottom",
    legend.text = element_text(size = 10, color = "gray25"),
    plot.margin = margin(t = 5, r = 5, b = 5, l = 0),
    plot.title = element_text(size = 11),
    strip.background = element_blank(),
    strip.text = element_text(size = 11, color = "gray30"),
    strip.text.y = element_text(angle = 0)
  )




#--------------------------------------------------
#----- Plots: Empirical regressor distribution
#--------------------------------------------------
library(GGally)

W.COLORS <- setNames(
  scales::hue_pal()(ncol(out$data$W)),
  nm = colnames(out$data$W)
)
W.LABS <- setNames(
  paste0("italic(W)[", seq_along(names(W.COLORS)), "]"),
  nm = names(W.COLORS)
)
W.LIMS <- range(out$data$W)

# Grab the wrapped GGally helpers under their own names
gg_cor <- wrap(
  funcVal = "cor", # See GGally::ggcorr
  align_percent = 0.5,
  size = 4,
  digits = 3,
  stars = FALSE
)

#--- Custom GGally pairs plot templates
# Custom lower triangular plots: points + zero lines
custom_lower <- function(data, mapping, ...) {
  ggplot(data, mapping) +
    geom_hline(yintercept = 0, color = "gray85", linewidth = 1) +
    geom_vline(xintercept = 0, color = "gray85", linewidth = 1) +
    geom_point(shape = 16, stroke = 0, alpha = 0.35, size = 0.95) +
    scale_x_continuous(limits = W.LIMS, expand = expansion(mult = c(0.02, 0.02))) +
    scale_y_continuous(limits = W.LIMS, expand = expansion(mult = c(0.02, 0.02))) +
    #coord_fixed() +
    theme(
      panel.grid.major = element_line(color = "gray90", linewidth = 0.5)
    )
}

# Custom diagonal plots: modified "barDiag" + no grid
custom_diag <- function(data, mapping, ...) {
  varname   <- as_label(mapping$x)
  fill.col  <- ggplot2::alpha(W.COLORS[[varname]], 0.5)
  ggplot(data, mapping) +
    geom_histogram(
      bins = 20,
      fill = fill.col,
      color = "black"
    ) +
    scale_x_continuous(expand = expansion(mult = c(0.02, 0.02)), limits = W.LIMS) +
    scale_y_continuous(expand = expansion(mult = c(0.00, 0.05))) +
    theme(
      panel.grid = element_blank(),
      axis.text.y = element_blank()
    )
}

# Custom upper triangular plots: correlation text + no grid
custom_upper <- function(data, mapping, ...) {
  gg_cor(data = data, mapping = mapping, ...) +
    theme(panel.grid = element_blank())
}

set.seed(1)
W.plt.idx <- sample(
  NROW(out$data$W),
  min(NROW(out$data$W), 1000, na.rm = T),
  replace = F
)

#--- Call GGally::ggpairs
plt.W <- out$data$W %>%
  data.frame() %>%
  slice(W.plt.idx) %>%
  GGally::ggpairs(
    title        = "Regressor distribution",
    lower = list(continuous = custom_lower),
    diag  = list(continuous = custom_diag),
    upper = list(continuous = custom_upper),
    columnLabels = W.LABS,
    labeller = label_parsed
  ) +
  theme(
    axis.text = element_text(size = 8, color = "gray50"),
    axis.ticks = element_blank(),
    axis.title = element_text(size = 11, color = "gray30"),
    panel.background = element_rect(fill = NA),
    panel.border = element_rect(fill = NA, color = "black"),
    plot.title = element_text(size = 11),
    strip.background = element_blank(),
    strip.text = element_text(size = 11, color = "gray30"),
    strip.text.y = element_text(angle = 0)
  )

#--------------------------------------------------
#----- Arrange and display plots
#--------------------------------------------------
library(gridExtra)
# plt.W
# plt.theta
# plt.res

plt.W.grob <- GGally::ggmatrix_gtable(plt.W)
grid.arrange(
  plt.res,
  arrangeGrob(plt.theta, plt.W.grob, ncol = 1),
  ncol = 2
)

