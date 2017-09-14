# test contplot (old cont)


# try with akima dataset
data(akima, package="akima")
akima
str(akima)

dat <- data.frame(cbind(a=akima$x, b=akima$y, c=akima$z))
#dat <- data.frame(cbind(x=akima$x, y=akima$y, z=akima$z))
dat
attributes(dat)

cont(x=a, y=b, z=c, data=dat, method="linear")
#cont(x=x, y=y, z=z, data=dat, method="linear")
# does not work. 

# this is what it should look like using the base plot devices:
surf<- akima::interp(x=dat$a, y=dat$b, z=dat$c,
                     xo=seq(min(dat$a),max(dat$a),length=length(dat$a)),
                       yo=seq(min(dat$b),max(dat$b),length=length(dat$b)))
with(surf, plot3D::image2D(x=x,y=y,z=z))
with(surf, points(akima$x,akima$y))
with(surf, plot3D::image3D(x=x,y=y,z=z))
