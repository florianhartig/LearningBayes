

window.dream <-
    function(x,
             start = 1+(end(x$Sequences)-1)*(1-fraction),
             fraction = 0.5,
             ...)
{
    window(x$Sequences, start = start, ...)
}
