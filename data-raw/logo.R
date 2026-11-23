# Hex logo for ComPoM. Run from the package root: Rscript data-raw/logo.R
# Fish (Poisson) as noteheads on a wavy staff (composition) over a spatial mesh (space).
suppressPackageStartupMessages({library(ggplot2); library(sf)})

pal <- list(bg = "#14213d", bg2 = "#1d3557", mesh = "skyblue", staff = "skyblue",
            fish = c("orange", "seagreen", "#b794f4"), text = "white", note = "#9fb0c8")

hex_pts <- function(r = 1, rot = 30) {
  a <- (seq(0, 300, 60) + rot) * pi / 180
  cbind(r * cos(a), r * sin(a))
}
hex <- hex_pts(1)
hex_poly <- st_polygon(list(rbind(hex, hex[1, ])))

# --- spatial mesh: Delaunay triangulation of jittered points, clipped to the hexagon
set.seed(7)
g <- expand.grid(x = seq(-1, 1, length.out = 9), y = seq(-1, 1, length.out = 10))
g$x <- g$x + rnorm(nrow(g), 0, 0.07); g$y <- g$y + rnorm(nrow(g), 0, 0.07)
tri <- st_collection_extract(st_triangulate(st_multipoint(as.matrix(g)), dTolerance = 0), "POLYGON")
tri <- st_intersection(st_sfc(tri), st_sfc(hex_poly))
tri_edges <- st_cast(st_boundary(tri), "LINESTRING")

# --- staff: five slightly wavy lines
xs <- seq(-0.78, 0.78, length.out = 200)
staff_y <- c(-0.02, 0.12, 0.26, 0.40, 0.54) - 0.05
staff <- do.call(rbind, lapply(seq_along(staff_y), function(i)
  data.frame(x = xs, y = staff_y[i] + 0.018 * sin(5 * xs + i * 0.8), line = i)))

# --- fish noteheads (body ellipse + triangular tail), stems and a beam
fish <- function(cx, cy, angle = 18, a = 0.125, b = 0.072, col) {
  th <- seq(0, 2 * pi, length.out = 60)
  body <- cbind(a * cos(th), b * sin(th))
  tail <- rbind(c(-a * 0.8, 0), c(-a * 1.55, b * 1.15), c(-a * 1.55, -b * 1.15))
  rot <- function(m) {
    r <- angle * pi / 180
    cbind(m[, 1] * cos(r) - m[, 2] * sin(r) + cx, m[, 1] * sin(r) + m[, 2] * cos(r) + cy)
  }
  list(body = data.frame(rot(body), col = col), tail = data.frame(rot(tail), col = col),
       eye = data.frame(x = cx + a * 0.55 * cos(angle * pi / 180) - 0.012,
                        y = cy + a * 0.55 * sin(angle * pi / 180) + 0.018))
}
notes <- data.frame(x = c(-0.50, -0.06, 0.36), y = c(0.03, 0.24, 0.12))
fs <- lapply(seq_len(3), function(i) fish(notes$x[i], notes$y[i], col = pal$fish[i]))
bodies <- do.call(rbind, lapply(seq_along(fs), function(i) cbind(fs[[i]]$body, id = i)))
tails  <- do.call(rbind, lapply(seq_along(fs), function(i) cbind(fs[[i]]$tail, id = i)))
eyes   <- do.call(rbind, lapply(fs, `[[`, "eye"))
names(bodies)[1:2] <- names(tails)[1:2] <- c("x", "y")
stem_x <- notes$x + 0.125 * cos(18 * pi / 180) * 0.92
# the stems end on a sloped beam
beam_y <- function(x) {
  y1 <- notes$y[1] + 0.50; y3 <- notes$y[3] + 0.50
  y1 + (y3 - y1) * (x - stem_x[1]) / (stem_x[3] - stem_x[1])
}
stems <- data.frame(x = stem_x, xend = stem_x, y = notes$y + 0.04, yend = beam_y(stem_x))
beam <- data.frame(x = stems$x[1], xend = stems$x[3], y = stems$yend[1], yend = stems$yend[3])
beam2 <- transform(beam, y = y - 0.06, yend = yend - 0.06)

# pom-poms: a fluffy cluster of small balls on the tail tip and on the back of each fish
lighten <- function(col, f = 0.45) {
  m <- grDevices::col2rgb(col) / 255
  grDevices::rgb(m[1] + (1 - m[1]) * f, m[2] + (1 - m[2]) * f, m[3] + (1 - m[3]) * f)
}
set.seed(3)
pompom <- function(cx, cy, r, col, id) {
  n  <- 16
  th <- runif(n, 0, 2 * pi); rr <- r * sqrt(runif(n))
  data.frame(x = cx + rr * cos(th), y = cy + rr * sin(th), col = col, id = id,
             fill = rep(c(col, lighten(col)), length.out = n))
}
pp <- do.call(rbind, lapply(seq_along(fs), function(i) {
  tail_tip <- fs[[i]]$tail[2:3, 1:2]
  rbind(pompom(mean(tail_tip[, 1]) - 0.015, mean(tail_tip[, 2]), 0.055, pal$fish[i], i),
        pompom(notes$x[i] - 0.01, notes$y[i] + 0.115, 0.045, pal$fish[i], i))
}))
hexdf <- data.frame(rbind(hex, hex[1, ])); names(hexdf) <- c("x", "y")

p <- ggplot() +
  geom_polygon(data = hexdf, aes(x, y), fill = pal$bg) +
  geom_sf(data = tri_edges, colour = pal$mesh, linewidth = 0.25, alpha = 0.28) +
  geom_path(data = staff, aes(x, y, group = line), colour = pal$staff, linewidth = 0.55, alpha = 0.9) +
  geom_polygon(data = tails, aes(x, y, group = id, fill = pal$fish[id])) +
  geom_polygon(data = bodies, aes(x, y, group = id, fill = pal$fish[id])) +
  geom_point(data = pp, aes(x, y, fill = fill), shape = 21, colour = pal$bg, stroke = 0.25,
             size = 2.9) +
  geom_point(data = eyes, aes(x, y), colour = pal$bg, size = 1.4) +
  geom_segment(data = stems, aes(x = x, xend = xend, y = y, yend = yend), colour = pal$note, linewidth = 0.9) +
  geom_segment(data = beam, aes(x = x, xend = xend, y = y, yend = yend), colour = pal$note, linewidth = 2.6, lineend = "butt") +
  geom_segment(data = beam2, aes(x = x, xend = xend, y = y, yend = yend), colour = pal$note, linewidth = 2.6, lineend = "butt") +
  annotate("text", x = 0, y = -0.50, label = "ComPoM", colour = pal$text, size = 17,
           fontface = "bold") +
  geom_polygon(data = hexdf, aes(x, y), fill = NA, colour = "skyblue", linewidth = 2.2,
               linejoin = "mitre") +
  scale_fill_identity() +
  coord_sf(xlim = c(-0.9, 0.9), ylim = c(-1.05, 1.05), expand = FALSE, datum = NA) +
  theme_void() +
  theme(plot.background = element_rect(fill = "transparent", colour = NA),
        panel.background = element_rect(fill = "transparent", colour = NA))

ggsave("man/figures/logo.png", p, width = 4.3, height = 4.98, dpi = 300, bg = "transparent")
cat("written man/figures/logo.png\n")
