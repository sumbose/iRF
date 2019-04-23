binplot.3d <- function(x, y, z, alpha=1, topcol="#ff0000",
                       sidecol="#aaaaaa", linecol="#000000")
{
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))
  x1 <- c(rep(c(x[1], x[2], x[2], x[1]), 3), rep(x[1], 4), rep(x[2], 4))
  z1 <- c(rep(0, 4), rep(c(0, 0, z, z), 4))
  y1 <- c(y[1], y[1], y[2], y[2], rep(y[1], 4), rep(y[2], 4), rep(c(y[1], y[2], y[2], y[1]), 2))
  x2 <- c(rep(c(x[1], x[1], x[2], x[2]), 2), rep(c(x[1], x[2], rep(x[1], 3), rep(x[2], 3)), 2))
  z2 <- c(rep(c(0, z), 4), rep(0, 8), rep(z, 8))
  y2 <- c(rep(y[1], 4), rep(y[2], 4), rep(c(rep(y[1], 3), rep(y[2], 3), y[1], y[2]), 2))

  rgl.quads(x1, z1, y1, col=rep(sidecol, each=4), alpha=alpha)
  rgl.quads(c(x[1], x[2], x[2], x[1]), rep(z, 4), c(y[1], y[1], y[2], y[2]), col=rep(topcol, each=4), alpha=1)
  rgl.lines(x2, z2, y2, col=linecol)
}

barplot3d <- function(z, alpha=1, col="#ff0000", scale=1, ints=NULL, main=NULL)
{
  save <- par3d(skipRedraw=TRUE)
  on.exit(par3d(save))
  z <- as.matrix(z)
  xy <- dim(z)
  x <- seq(xy[1])
  y <- seq(xy[2])
  mxz <- max(z, na.rm=TRUE)
  z <- z / max(z, na.rm=TRUE) * max(x, y) * scale
  for (i in x)
  {
    for (j in y)
    {
      binplot.3d(c(i, i+1), c(j, j+1), z[i,j], alpha=alpha, topcol=col)
    }
  }

  axis3d('xz', at=c(1.5, 2.5, 3.5), labels=c('AA', 'AD', 'DD'))
  if (!is.null(ints)) mtext3d(text=ints[1], edge='xz', line=2.5, size=5)
  axis3d('xz+', at=c(1.5, 2.5, 3.5), labels=c('AA', 'AD', 'DD'))
  if (!is.null(ints)) mtext3d(text=ints[1], edge='xz+', line=2.5, size=5)


  axis3d('z', at=c(1.5, 2.5, 3.5), labels=c('AA', 'AD', 'DD'))
  if (!is.null(ints)) mtext3d(text=ints[2], edge='z', line=2.5, size=5)
  axis3d('z+', at=c(1.5, 2.5, 3.5), labels=c('AA', 'AD', 'DD'))
  if (!is.null(ints)) mtext3d(text=ints[2], edge='z+', line=2.5, size=5)

  axis3d('y+', at=c(0, (1 / mxz) * max(z)), labels=c('0', '1'))
  mtext3d(text='P(Y=1)', edge='y+', line=1.5, size=5)
  title3d(main=main, size=8)
}

