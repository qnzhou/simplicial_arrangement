library(ggplot2)

data <- read.csv("arrangement_scaling.csv")
summary(data)

p <- ggplot(data)
p <- p + geom_errorbar(aes(x=num_planes, y=mean, ymax = mean + std_dev, ymin = mean - std_dev), color="darkgray")
p <- p + geom_point(aes(num_planes, mean))
p <- p + theme_minimal()
p <- p + scale_x_continuous(breaks = round(seq(10, 100, by = 10),1))
p <- p + ylab("Time (s)") + xlab("Number of functions") + ggtitle("Arrangement Scaling")
ggsave("arrangement_scaling.pdf", width=3, height=2, unit="in")

data <- read.csv("material_interface_scaling.csv")
summary(data)

p <- ggplot(data)
p <- p + geom_errorbar(aes(x=num_planes, y=mean, ymax = mean + std_dev, ymin = mean - std_dev), color="darkgray")
p <- p + geom_point(aes(num_planes, mean))
p <- p + theme_minimal()
p <- p + scale_x_continuous(breaks = round(seq(10, 100, by = 10),1))
p <- p + ylab("Time (s)") + xlab("Number of functions") + ggtitle("Material Interface Scaling")
ggsave("material_interface_scaling.pdf", width=3, height=2, unit="in")
