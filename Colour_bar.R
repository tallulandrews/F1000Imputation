color.bar <- function(colours, min=0, max=1, ticks.at=seq(min, max,len=min(10,length(colours)+1)), ticks.lab=ticks.at, title='', horiz=FALSE, add=FALSE) {
    scale = (length(colours))/(max-min)

    if (horiz) {
	if (!add) {
	        dev.new(width=5, height=1.75)
	}
        par(mar=c(4,0,1,0))
        plot(c(min,max),c(0,10), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main="")
        axis(1, at=ticks.at, labels=ticks.lab, las=1)
        for (i in 1:(length(colours))) {
             x = (i-1)/scale + min
             rect(x,0,x+1/scale,10, col=colours[i], border=NA)
        }
        mtext(side=1, line=2.5, title, font=2, cex=1.1)

    } else {
	if (!add) {
	        dev.new(width=1.75, height=5)
	}
	par(mar=c(0,4,1,0))
        plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
        axis(2, at=ticks.at, labels=ticks.lab, las=1)
        for (i in 1:(length(colours))) {
             y = (i-1)/scale + min
             rect(0,y,10,y+1/scale, col=colours[i], border=NA)
        }
    }
}

