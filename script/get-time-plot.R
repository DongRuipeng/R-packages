setwd("D:/Git/project/missing-data/result/time-result")
require(ggplot2)
require(ggpubr)
load("./time.RData")

p.p <- ggplot(T.p) +
  aes(x = p, y = time, colour = method) +
  geom_line(aes(linetype = method)) +
  geom_point(aes(shape = method)) +
  xlab("p") +
  ylab("Time (s)")

p.q <- ggplot(T.q) +
  aes(x = q, y = time, colour = method) +
  geom_line(aes(linetype = method)) +
  geom_point(aes(shape = method)) +
  xlab("q") +
  ylab("Time (s)")

# p.n <- ggplot(T.n) +
#   aes(x = n, y = time, colour = method) +
#   geom_line(aes(linetype = method)) +
#   geom_point(aes(shape = method)) +
#   xlab("n") +
#   ylab("Time (s)")

p <- ggarrange(p.p, p.q, ncol = 2, nrow = 1, common.legend = TRUE, legend = "right")

plot(p)
