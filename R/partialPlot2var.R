require(rgl)

partialPlot2var <- function(x1
                    , x2
                    , y
                    , gridlength=NULL
                    , x1_grid=NULL
                    , x2_grid=NULL
                    , x1lab='v1'
                    , x2lab='v2'
                    , ylab=NA
                    , range.color = NULL # range of values to be colored
                    , col.palette = c('blue', 'yellow') # a character vector of colors, used as input to colorRampPalette
                    , plot_quantile_scale = TRUE
                    , plot.colorbar = TRUE
                    , ... # additional arguments to pass to rgl::persp3d
                     ){
if ((is.null(x1_grid)|is.null(x2_grid)) & (is.null(gridlength))){
    stop('Either gridlength or both grids should be specified')
}
if (is.null(x1_grid))
  x1_grid=quantile(x1, prob=seq(0, 1, length.out=gridlength+1))
if (is.null(x2_grid))
  x2_grid=quantile(x2, prob=seq(0, 1, length.out=gridlength+1))

if (max(is.na(x1))|max(is.na(x2))|max(is.na(y)))
    stop('missing values in x1, x2 or y !!')

z=array(0,c(length(x1_grid)-1, length(x2_grid)-1))
for (i in 1:nrow(z))
   for (j in 1:ncol(z)){
      z[i,j]=100*mean(y[x1> x1_grid[i] & x1<= x1_grid[i+1] & 
                    x2> x2_grid[j] & x2<= x2_grid[j+1]])  # /sum(sum(y[x1 > x1_grid[i] & x1 <= x1_grid[i+1]]))
      if (is.na(z[i,j])){
          stop(paste('Not enough observation in cell', i, ':', j))
      }
   }
      
nbcol=100
zfacet = z

if (is.null(range.color)){
   range.color = range(zfacet)
}
else{
  range.color = range(range.color)
}
if ((range.color[1] > min(zfacet)) | (range.color[2] < max(zfacet)))
    warning('data out of range of range.color -- adjust range.color !!')
facetcol <- cut(c(range.color, as.vector(zfacet)), nbcol)[-seq(2)]


jet.colors <- colorRampPalette(col.palette)
color<- jet.colors(nbcol)

if (plot_quantile_scale == TRUE){
  x1_grid = seq(0, 1, length.out = gridlength+1)
  x2_grid = seq(0, 1, length.out = gridlength+1)
}

xtick = (x1_grid[-length(x1_grid)]+x1_grid[-1])/2
ytick = (x2_grid[-length(x2_grid)]+x2_grid[-1])/2

if (plot.colorbar)
    color.bar(colorRampPalette(col.palette)(nbcol), min = range.color[1], max=range.color[2])

persp3d(x=xtick
      , y=ytick
      , z
      , xlab=x1lab
      , ylab=x2lab
      , zlab=ylab
      , col=color[facetcol]
      , lit=TRUE
      , smooth = TRUE
      , ...
       )


return(z)
}

# Function to plot colorbar
# Courtesy:  http://www.colbyimaging.com/wiki/statistics/color-bars
color.bar <- function(lut, min, max=-min, nticks=11, ticks=round(seq(min, max, len=nticks), 2), title='') {
        scale = (length(lut)-1)/(max-min)

            dev.new(width=1.75, height=5)
            plot(c(0,10), c(min,max), type='n', bty='n', xaxt='n', xlab='', yaxt='n', ylab='', main=title)
            axis(2, ticks, las=1)
            for (i in 1:(length(lut)-1)) {
                     y = (i-1)/scale + min
                          rect(0,y,10,y+1/scale, col=lut[i], border=NA)
                         }
        }

