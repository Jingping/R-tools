par(mfrow=c(1,1))

x <- seq(-4, 4, length=100)
hx <- dnorm(x)

degf <- c(1, 3, 8, 30)
colors <- c("red", "blue", "darkgreen", "gold", "black")
labels <- c("df=1", "df=3", "df=8", "df=30", "normal")

plot(x, hx, type="l", lty=2, xlab="x value", ylab="Density", main="Comparison of t Distributions")

## shade
xmin<--1
xmax<-1
ymin<--.1
ymax<-.5
color<-"lightgray"
rect(xmin,ymin,xmax,ymax,col=color, density=100)

lines(x, hx, type="l", lty=2)

for (i in 1:4)
{	lines(x, dt(x,degf[i]), lwd=2, col=colors[i])
}

legend("topright", inset=.05, title="Distributions", labels, lwd=2, lty=c(1, 1, 1, 1, 2), col=colors)



