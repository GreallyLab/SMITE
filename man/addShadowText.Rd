\name{addShadowText}
\alias{addShadowText}
\alias{addShadowText,ANY-method}
\title{
Add shadow text (a second color bordering the text) to a plot
}
\description{
This is a usefule function to help text stand out on busy backgrounds
like gene networks
}
\usage{
addShadowText(x, y = NULL, labels, col = "white", bg = "black", 
theta = seq(pi/4, 2 * pi, length.out = 8), r = 0.1, ...)
}
\arguments{
    \item{x}{
A numeric vector of x coordinates
}
    \item{y}{
A numeric vector of y coordinates
}
    \item{labels}{
A character vector to be plotted at the specified coordinates
}
    \item{col}{
The text color
}
    \item{bg}{
The color of the outline
}
    \item{theta}{
The number of shadows to plot
}
    \item{r}{
The radius for the shadows
}
    \item{...}{
 Additional plotting arguments
}}
\details{
The function creates its effect by plotting theta shadows at r radius around the
text to create the illusion of a text shadow
}
\value{
Adds shadow text to plot
}
\references{
http://article.gmane.org/gmane.comp.lang.r.general/147787
}
\author{
Greg.Snow <at> imail.org
}
\note{
This function was adapted by N. Ari Wijetunga for SMITE.
}

\seealso{
text, mtext
}
\examples{
plot.new()

addShadowText(x = .5,y = .5,"TEST",col="white",bg="gray")


}
\keyword{ aplot }