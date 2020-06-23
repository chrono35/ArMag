#  Licence ----
#
#  Copyright or © or Copr. CNRS	2019
#
#
# This package is under
# The “Creative Commons Attribution-ShareAlike International License” version 4.0

#  A copy of the License is available at
#  http://www.r-project.org/Licenses/
#  _________________________________________________________________________________

# Version 2020-01-15
#

#' @author "Philippe DUFRESNE"

# Equation du 3 degrées
# résolution du 3 eme degres pour calcul mcFadden importé de ARMAG
# modifier suivant livre photocopié
# Fonction interne
EQUATION_DEGRE_3 <- function(A1, A2, A3)
{

  Q <- (3*A2-A1*A1)/9
  R <- (9*A1*A2-27*A3-2*A1*A1*A1)/54
  D <- Q*Q*Q+R*R

 #  Discriminant équation du 3eme degre
  if (D>0) {
    S <- exp( log( R + sqrt(D)  )/3)
    if ( R>sqrt(D))
      T1 <- exp(log( R - sqrt(D)  )/3)
    else
      T1 <- -exp(log( -R + sqrt(D)  )/3)

    PZ1 <- S+ T1 -(A1/3)
     # Les deux autres solutions sont complexes
    PZ2 <- 0
    PZ3 <- 0
  }
  else if (D==0) {
    S <- exp( log(-Q/2)/3)
    PZ1 <- 2*S-(A1/3)
    PZ2 <- -S-(A1/3)
    PZ3 <- -S-(A1/3)
  }
  else if (D<0) {
    TETA <- acos( R/sqrt(-Q*Q*Q) )/3
    U <- 2*sqrt(-Q)
    PZ1 <- U*cos(TETA)-(A1/3);
    PZ2 <- U*cos(TETA + 2/3*pi)-(A1/3)
    PZ3 <- U*cos(TETA + 4/3*pi)-(A1/3)
  }

  rslt <- NULL
  rslt$PZ1 <- as.numeric(PZ1)
  rslt$PZ2 <- as.numeric(PZ2)
  rslt$PZ3 <- as.numeric(PZ3)

  return(rslt)
}

# Equation du 4 degrées
# résolution du 4 ème degres pour calcul mcFadden importé de ARMAG
# Fonction interne
EQUATION_DEGRE_4 <- function(B1, B2, B3, B4)
{
  E3 <- EQUATION_DEGRE_3(-B2, B1*B3-4*B4, 4*B2*B4-B3*B3-B1*B1*B4)
  TZ <- E3$PZ1
  if (E3$PZ2>TZ)
    TZ <- E3$PZ2
  if (E3$PZ3>TZ)
    TZ <- E3$PZ3

  SQ <- sqrt(B1*B1/4-B2+TZ)
  C1 <- B1/2-SQ
  C2 <- TZ/2-(B1*TZ/2-B3)/(2*SQ)
  if ((C1*C1-4*C2)>=0) {
    PX1 <- (-C1+ sqrt(C1*C1-4*C2))/2
    PX2 <- (-C1- sqrt(C1*C1-4*C2))/2
  } else {
    PX1 <- 0
    PX2 <- 0
  }

  D1 <- B1/2+SQ
  D2 <- TZ/2+(B1*TZ/2-B3)/(2*SQ)
  if ((D1*D1-4*D2)>=0)  {   # une solution double
    PX3 <- (-D1+ sqrt(D1*D1-4*D2))/2
    PX4 <- (-D1- sqrt(D1*D1-4*D2))/2
  }
  else {
    PX3 <- 0
    PX4 <- 0
  }
   rslt <- NULL
   rslt$PX1 <- as.numeric(PX1)
   rslt$PX2 <- as.numeric(PX2)
   rslt$PX3 <- as.numeric(PX3)
   rslt$PX4 <- as.numeric(PX4)

   return(rslt)

}

#' Calcul l'arc sinus d'un réel renvoie un angle entre -PI/2 et +PI/2
#' @return angle en radian
#' Fonction interne
n_arcsin <- function(x)
{
  n <- NA
  if ( abs(x) < 1) {
    n <- atan(x/sqrt(1-x*x))
  }
  else {
    if (x>=1)
    {n <- pi/2}
    else
    {n <- -pi/2}
  }
  return(n)
}


#' Calcul l'arc cosinus d'un réel renvoie un angle entre 0 et +PI
#' @return angle en radian
#' Fonction interne
n_arccos <- function(x)
{
  n <- NA
  if (abs(x)<1)
    n <- (pi/2-atan(x/sqrt(1-x*x)) )

  else
    if (x>=1)
      n <- 0

    else
      n <- pi

    return(n)
}

#' Calcul de l'angle D  entre  -pi <=  angD  < +3 pi  retourne en radian
#' @return angle en radian
#' Fonction interne
angleD <- function(X, Y)
{
  res <- NULL
  for (i in 1:length(X)) {
    if (X[i] == 0)
      res <- c(res, sign(Y[i])*pi/2)
    else {
      if (X[i]<0)
        res <- c(res, atan(Y[i]/X[i]) + pi)
      else
        res <- c(res, atan(Y[i]/X[i]) )
    }
  }

  return(res)
}

#' Statistique de mcFadden sur l'inclinaison seule à partir des coordonées XYZ
#' @seealso \code{\link{stat.mcFadden}}, \code{\link{stat.fisher}}
#' @export
stat.mcFadden.XYZ <- function(TabX, TabY, TabZ)
{
  N <- length(TabX)
  VP <- to.polar(TabX, TabY, TabZ)

  stat.mcFadden(inc, dec)

}

#' Statistique de mcFadden sur l'inclinaison seule à partir des coordonées I et D
#' @param Data liste des inclinaisons en degré ou une data.frame avec les variables $I et $D
#' @param dec liste des déclinaisons en degré
#' @param inc.absolue calcul avec la valeur absolue des inclinaisons
#' @return  en degré, un data.frame "n", "imoy.McFadden", "imoy.McElhinny", "a95.mcFad", "a95.eqFish", "Kb", "Kssb", "imin", "imax", "dmin", "dmax"
#' @seealso \code{\link{stat.mcFadden.XYZ}}, \code{\link{stat.fisher}}
#' @references https://doi.org/10.1111/j.1365-246X.1990.tb05683.x
#' @export
stat.mcFadden <- function(Data, dec = NULL, inc.absolue = TRUE)
{

  if (is.null(dec)) {
    inc <- Data$I
    dec <- Data$D
    Dmin <- min(dec)
    Dmax <- max(dec)
  } else {
    inc <- Data
    Dmin <- min(dec)
    Dmax <- max(dec)
  }

  if (inc.absolue)
    inc <- abs(inc)

  # Calcul des sommes
  N <- length(inc)

  if (N<=2)
  {
      warning("length(inc) < 3")
      return()
  }

  Imin <- min(inc)
  Imax <- max(inc)

  Imoy <- 0

  # Passage en radian pour les calculs
  i.rad <- as.numeric(inc*pi/180)

  A <- sum(sin(i.rad))
  B <- sum(cos(i.rad))


  P <- 2*(N+2*B)/A
  Q <- -6
  R <- 2*(N-2*B)/A
  S <- 1

  E4 <- EQUATION_DEGRE_4(P, Q, R, S)

  KJ <- N/(2*(N-A*sin(2*atan(E4$PX1))-B*cos(2*atan(E4$PX1))))

  Imoy <- 2*atan(E4$PX1)
  K <- KJ


  KJ <- N/(2*(N-A*sin(2*atan(E4$PX2))-B*cos(2*atan(E4$PX2))))
  if (KJ>K) {
    Imoy <- 2*atan(E4$PX2)
    K <- KJ
  }

  KJ <- N/(2*(N-A*sin(2*atan(E4$PX3))-B*cos(2*atan(E4$PX3))))
  if (KJ>K) {
    Imoy <- 2*atan(E4$PX3)
    K <- KJ
  }

  KJ <- N/(2*(N-A*sin(2*atan(E4$PX4))-B*cos(2*atan(E4$PX4))))
  if (KJ>K) {
    Imoy <- 2*atan(E4$PX4)
    K <- KJ
  }

  CFad <- sum(cos(i.rad-Imoy))
  SFad <- sum(sin(i.rad-Imoy))

  I.McElhinny <- Imoy /pi*180     # modif effectuée le 2017/01/25
  imoy.McFadden <- (Imoy + (SFad/CFad) ) /pi*180

  a95mcFadden <- 1 - ((SFad/CFad)^2)/2 - (qf(.95, 1, N-1)*(N-CFad)/(CFad*(N-1)))
  a95mcFadden <- acos(a95mcFadden) /pi*180
  Kssb <-  (N-1)/2/(N-CFad)

  # retour en degré #
  return(data.frame( n= N, imoy.McFadden = imoy.McFadden, imoy.McElhinny = I.McElhinny,
                     a95.mcFad = a95mcFadden, a95.eqFish = 2.4477/sqrt(N*K)/pi*180,
                         Kb = K, Kssb = Kssb,
                         imin = Imin, imax = Imax, Dmin = Dmin, Dmax = Dmax))

  }

#' Statistique de Fisher
#' modifié, pas de pondération calcul de A95 Vrai sans simplification tel que fisher 1953
#' @param Data liste des inclinaisons en degré ou une data.frame avec les variables $I et $D
#' @param dec liste des déclinaisons en degrés
#' @param aim liste des aimantations, facultatif
#' @param pfish pourcentage de confiance
#' @param inc.absolue calcul avec la valeur absolue des inclinaisons
#' @return  en degrés
#' @seealso \code{\link{stat.mcFadden}}
#' @keywords fisher
#' @export
stat.fisher <- function (Data, dec = NULL, aim = NA, pfish = 0.95, inc.absolue = TRUE)
{
  if (is.null(dec)) {
    inc <- Data$I
    dec <- Data$D
  } else {
    inc <- Data
    dec <- dec
  }

  n <- length(inc)
  if (length(dec) != n) {
    return("length (dec) diff length(inc)")
  }
  if (is.na(aim) || length(aim) != n) {
   # message("length(aim) diff length(inc)")
    aim <-  rep(1, n)
  }

  if (inc.absolue == TRUE)
    inc <- abs(inc)

  imin <- min(inc)
  imax <- max(inc)
  dmin <- min(dec)
  dmax <- max(dec)


  # Passage en radian pour les calculs
  i.rad <- inc/180*pi
  d.rad <- dec/180*pi

  sx <- sum(aim*cos(i.rad)*cos(d.rad))
  sy <- sum(aim*cos(i.rad)*sin(d.rad))
  sz <- sum(aim*sin(i.rad))
  sn <- sum(aim)

  r <- sqrt(sx*sx+sy*sy+sz*sz)

  imoy <- n_arcsin(sz/r)
  dmoy <- angleD(sx,sy)
  KF <- sn/(sn-r)
  # Calcul de A95 Vrai sans simplification tel que fisher 1953
  a95 <- exp( (1/(n-1))* log(1/(1-pfish)) )
  a95 <- (a95 - 1 )*(n-r)/r
  a95 <- acos(1-a95)

  # correction du biais
  KF <- ((n-1)/n) * KF

  if (KF<10) {
    delta <- log(1+ (1-pfish) * (exp(2*KF) -1) ) /KF
    delta <- acos(delta-1)
  } else {
    delta <- log(1-pfish)/KF
    delta <- acos(1 + delta)
  }


  # retour en degrés
  return(data.frame(n = n, imoy =imoy/pi*180, dmoy = dmoy/pi*180, alpha=a95/pi*180, pfish = pfish, delta = delta/pi*180, KF = KF,
                    imin = imin, imax = imax, dmin = dmin, dmax = dmax))
}


#' Valeurs maximun d'un tracer lambert
#' Permet de trouver les valeurs minimales et maximales pour tracer un lambert
#' Fonction interne
find.extremum <- function(i.min = 0, i.max = 90, d.min = -90, d.max = 270)
{
  d.seq <- seq( d.min, d.max, length.out = 300)
  x1.seq <- sin(d.seq*pi/180 ) * (1 - (i.min-i.min)/(90-i.min))
  x2.seq <- sin(d.seq*pi/180 ) * (1 - (i.max-i.min)/(90-i.min))

  y1.seq <- cos(d.seq*pi/180 ) * (1 - (i.min-i.min)/(90-i.min))
  y2.seq <- cos(d.seq*pi/180 ) * (1 - (i.max-i.min)/(90-i.min))

  x.range <- range(x1.seq, x2.seq)
  y.range <- range(y1.seq, y2.seq)

  return(c(x.range = x.range, y.range = y.range))
}

# Transformation des mesures en coordonnées graphiques ----------------------
#' Valeur de X graphique pour I et D
#' Fonction interne
#' @param i inclinaison en degrés
#' @param d déclinaison en degrés
X <- function (i, d, ray=1, i.min = 0, box.range)
{
  x <- sin(d*pi/180 ) * (1 - (abs(i)-i.min)/(90-i.min))  - (box.range[2] + box.range[1])/2
  c <- NA
  if ((box.range[2] - box.range[1]) > (box.range[4] - box.range[3]))
    c <- box.range[2] - box.range[1]
  else
    c <- box.range[4] - box.range[3]

  x /c *ray *2

}

#' Valeur de Y graphique pour I et D
#' Fonction interne
#' @param i inclinaison en degrés
#' @param d déclinaison en degrés
Y <- function(i, d, ray=1, i.min = 0, box.range)
{
  y <- cos(d*pi/180 ) *  (1 - (abs(i)-i.min)/(90-i.min))  - (box.range[4] + box.range[3])/2

  c<-NA
  if ((box.range[2] - box.range[1]) > (box.range[4] - box.range[3]))
    c <- box.range[2] - box.range[1]
  else
    c <- box.range[4] - box.range[3]

  y /c *ray *2

}

# Transformation de type de repère -----------------
## Conversion

#' Conversion de Degrés Minute Second en Degré Décimal
#' @param degre degrés entier
#' @param minute minute entière
#' @param second seconde en décimal
#' @seealso \code{\link{DD.to.DMS}}
#' @export
DMS.to.DD <- function(degre, minute, second = 0)
{
  if (degre>0) {
    return(as.numeric(degre + minute/60 + second/3600) )
  } else {
    return(as.numeric(degre - minute/60 - second/3600) )
  }

}

#' Conversion de Degré Décimal en Degrés Minute Second
#' @param degre degrés décimal
#' @seealso \code{\link{DMS.to.DD}}
#' @export
DD.to.DMS <- function(degre)
{
  absDegre <- abs(degre)
  res <- NULL
  res$Degre <- floor(absDegre)
  res$Minute <- floor((absDegre - res$Degre)*60)
  res$Second <- ((absDegre - res$Degre)*60 - res$Minute)*60
  if (degre<0)
    res$Degre <- res$Degre*(-1)
  return(res)
}

#' Convertie la déclinaison dans la convention AM
#' @param d déclinaison en degré
#' @return angle en degré en -90° et 270°
#' @export
D.AM <- function(d)
{
  res <- d
  if (d< (-90))
    res <- 360 + d

  if (d>270)
    res <- d-360

  return(res)
}

#' Convertie la déclinaison dans la convention paleomag
#' @param d déclinaison en degré
#' @return déclinaison entre 0 et 360°
#' @export
D.pal <- function(d)
{
  if ( d>270)
    return(d-360)
}


# Transformation de type de coordonnée ----

## to cartesian vector  ----
#' Valeur de Y pour I et D
#' @param inc liste des inclinaisons en degré
#' @param des liste des déclinaisons en degré
#' @seealso \code{\link{to.cartesian}} , \code{\link{to.polar}}
#' @export
to.cartesian.Y <- function(inc, dec, aim =1)
{
  as.numeric(aim*cos(inc/180*pi)*sin(dec/180*pi) )
}

#' Valeur de X pour I et D
#' @param inc liste des inclinaisons en degré
#' @param des liste des déclinaisons en degré
#' @seealso \code{\link{to.cartesian}} , \code{\link{to.polar}}
#' @export
to.cartesian.X <- function(inc, dec, aim=1)
{
  as.numeric(aim*cos(inc/180*pi)*cos(dec/180*pi) )
}

#' Valeur de Z pour I et D
#' @param inc liste des inclinaisons en degré
#' @param des liste des déclinaisons en degré
#' @seealso \code{\link{to.cartesian}} , \code{\link{to.polar}}
#' @export
to.cartesian.Z <- function(inc, dec, aim=1)
{
  as.numeric(aim*sin(inc/180*pi))
}


#' Valeur de X, Y et Z pour I et D
#' @param inc liste des inclinaisons en degrés
#' @param des liste des déclinaisons en degrés
#' @seealso  \code{\link{to.polar}}
#' @export
to.cartesian <- function(inc, dec, aim=1)
{
  res <- NULL
  resT <- NULL
  for (i in 1:length(inc)) {
    res$X <- as.numeric( aim*cos(inc[i]/180*pi)*cos(dec[i]/180*pi) )
    res$Y <- as.numeric( aim*cos(inc[i]/180*pi)*sin(dec[i]/180*pi) )
    res$Z <- as.numeric( aim*sin(inc[i]/180*pi) )
    res$I <- inc[i]
    res$D <- dec[i]
    res$F <- aim[i]
    resT <- c(resT, res)
  }
  return( resT )
}




#' Coordonnées polaires pour X, Y et Z
#' @return I, D, F en degré
#' @seealso \code{\link{to.polar}} , \code{\link{cartesien}}
#' @export
to.polar <- function(X, Y, Z)
{ #en A/m
  res <- NULL
  resT <- NULL
  for (i in 1:length(X)) {
    res$X <- X[i]
    res$Y <- Y[i]
    res$Z <- Z[i]
    aim <-  as.numeric( sqrt(X[i]*X[i] + Y[i]*Y[i] + Z[i]*Z[i]) )
    if ( aim== 0) {
      res$I <- as.numeric( 0 )
      res$D <- as.numeric( 0 )
    } else {
      res$I <- as.numeric( n_arcsin(Z[i]/aim) /pi*180 )
      res$D <- as.numeric( angleD(X[i],Y[i]) /pi*180 )
    }

    res$F <- aim

    resT <- c(resT, res)
  }
  return( resT )
}

# Calcul vecteur polaire ----
#' inclinaison pour X, Y et Z
#' @return angle en degré
#' @export
to.polar.I <- function(X, Y, Z)
{ #en A/m
  res <- NULL
  for (i in 1:length(X)) {
    A <- as.numeric( sqrt(X[i]*X[i]+Y[i]*Y[i]+Z[i]*Z[i]) )
    res <- c(res, n_arcsin(Z[i]/A)/pi*180)
  }
  return(res)
}

#' Déclinaison pour X, Y et Z
#' #' @return angle en degré
#' @export
to.polar.D <- function(X, Y, Z)
{ #en A/ms
  res <- NULL
  for (i in 1:length(X)) {
     res <- c(res, angleD(X[i],Y[i])/pi*180)
  }
  return(res)
}

#' champs pour X, Y et Z
#' @export
to.polar.F <- function(X, Y, Z)
{ #en A/m
  res <- NULL
  for (i in 1:length(X)) {
    res <- c(res, sqrt(X[i]*X[i]+Y[i]*Y[i]+Z[i]*Z[i]))
  }
  return(res)
}

# Trace des graphes ----

## Plot Lambert directionnel ID ----

#' lambert.ID.grid
#' Trace une grille dans le repère Lambert, fonctionne avec la fonction lambert()
#'
#' # Pour choisir les graduations
#' et pour paleomag : lab.pos$D = c(seq(270, 350, by=10), seq(0, 90, by=10))
#' Label.pos : doit être dans l'étendu des dex.min, dec.max, inc.min, inc.max
#'
#' @param  radlab écrit les labels sous forme d'étoile
#' @param  label.pos séquence de valeur à afficher en I et D, attention format particulier !!
#' @examples
#' label.pos = NULL
#' label.pos$I = seq(0, 90, by=20)
#' label.pos$D = seq(0, 90, by=10)
#' @export
lambert.ID.grid <- function (main = "", xlab = "", ylab = "", labels = NA, label.pos = NULL, radlab = FALSE,
                        start = 0, clockwise = FALSE, label.prop = 1.1,
                        grid.col = "gray", grid.bg = "transparent", show.radial.grid = TRUE, labels.precision = 0,
                        dec.min = -90, dec.max = 270, inc.min = 0, inc.max = 90, new = TRUE, ...)
{
  # setting up coord. system
  if (new == TRUE) {
    par( pty = "s")

    maxlength <- 100 # la valeur n'a pas d'influence, la fonction plot() calcul le reste
    plot(c(-maxlength, maxlength), c(-maxlength, maxlength), type = "n", axes = FALSE,main = main, xlab = xlab, ylab = ylab, new = new)
  }

  labelsD <- NULL
  maxlength <-  100
  par(xpd = TRUE)
  box.range <- find.extremum(inc.min, inc.max, dec.min, dec.max)

  anglesD <- seq(dec.min, dec.max, by = 1) # angle de deviation en degrée

  if (is.null(label.pos)) {
    labelI.pos <- seq(inc.min, inc.max, by = 10)
    labelD.pos <- seq(dec.min, dec.max, by = 10)
  }
  else {
    labelI.pos <- label.pos$I
    labelD.pos <- label.pos$D
  }
  # Supprime la superposition des labels pour D égale à -90 et 270
  if ((labelD.pos[1] == -90) && (labelD.pos[length(labelD.pos)] == 270))
    labelD.pos <- labelD.pos[-length(labelD.pos)]

  if (show.radial.grid) {
    # Trace un cercle  autour
    xpos <- X(inc.min, anglesD, maxlength, i.min = inc.min, box.range)
    ypos <- Y(inc.min, anglesD, maxlength, i.min = inc.min, box.range)
    lines(xpos, ypos, col = adjustcolor( grid.col, alpha.f = 0.5))

    # Trace les cercles radiaux concentriques
    if (length(labelI.pos)>0)
      for (i in seq(length(labelI.pos), 1, by = -1)) {
        xpos <- X(labelI.pos[i], anglesD, maxlength, i.min = inc.min, box.range)
        ypos <- Y(labelI.pos[i], anglesD, maxlength, i.min = inc.min, box.range)
        lines(xpos, ypos, col = adjustcolor( grid.col, alpha.f = 0.5)) #, border = grid.col)
      }


    if (!is.null(labels)) {
      if (is.na(labels[1]))
        labelsI <- as.character(round(labelI.pos, labels.precision))

      labelsD <- as.character(round(labelD.pos, labels.precision))
    }


    if (clockwise == FALSE)
      labelD.pos <- -labelD.pos
     if (start)
       labelD.pos <- labelD.pos + start
  # Trace les rayons

  #if (show.radial.grid) {
    for (i in 1: length(labelD.pos)) {
      xposA <- X(inc.min, labelD.pos[i], ray = maxlength, inc.min, box.range)
      yposA <- Y(inc.min, labelD.pos[i], ray = maxlength, inc.min, box.range)

      xposB <- X(inc.max, labelD.pos[i], ray = maxlength, inc.min, box.range)
      yposB <- Y(inc.max, labelD.pos[i], ray = maxlength, inc.min, box.range)
      segments( x0=xposA, y0=yposA, x1=xposB, y1=yposB, col = grid.col)

    }


    for (label in 2:length(labelI.pos)) {
      xpos <- X(labelI.pos[label], dec.min, ray = maxlength, inc.min, box.range)
      ypos <- Y(labelI.pos[label], dec.min, ray = maxlength, inc.min, box.range)
      #labelsrt <- dec.min + 90 # labelI.pos[label] + 90    #* label.prop
      corect<- 6
      if (dec.min >= -90 && dec.min < 0) {
        labelsrt <- ( 270 - dec.min)
        xpos <- xpos + label.prop*cos(dec.min/180*pi) * corect
        ypos <- ypos - label.prop*(2+sin(dec.min/180*pi)) * corect
      }
      else if (dec.min >= 0 && dec.min < 90) {
        labelsrt <- (90 - dec.min)
        xpos <- xpos - label.prop*cos(dec.min/180*pi) * corect
        ypos <- ypos + label.prop*sin(dec.min/180*pi) * corect
      }
      else if (dec.min >= 90 && dec.min < 180) {
        labelsrt <- (90 - dec.min)
        xpos <- xpos + label.prop*cos(dec.min/180*pi) * corect
        ypos <- ypos + label.prop*sin(dec.min/180*pi) * corect
      }
      else {
        labelsrt <- (270 - dec.min)
        xpos <- xpos - label.prop*cos(dec.min/180*pi) * corect
        ypos <- ypos - label.prop*sin(dec.min/180*pi) * corect
      }

      # if (labelD.pos[label] > 0 && labelD.pos[label] < 180)
      #   labelsrt <- (90 - labelD.pos[label])
      # else
      #   labelsrt <- (270 - labelD.pos[label]) labelsI[label],

      text(xpos, ypos,labelsI[label] , srt = labelsrt, cex = par("cex.axis"))
    #  boxed.labels(xpos, ypos, labelsI[label], ypad = par("cex.axis"), border = FALSE, cex =  par("cex.axis"))
    }

    # Ecrit les textes des angles
  }

  else {
    # Trace un cercle  autour
    #for (i in seq(dec.min, dec.max, by = 1)) {
      xpos <- X(inc.min, anglesD, maxlength, i.min = inc.min, box.range)
      ypos <- Y(inc.min, anglesD, maxlength, i.min = inc.min, box.range)
      lines(xpos, ypos, col = adjustcolor( grid.col, alpha.f = 0.5)) #, border = grid.col)
    #}
  }
    if (!is.null(labelsD)) {
      xpos <- X(inc.min, labelD.pos, ray = maxlength, inc.min, box.range) * label.prop
      ypos <- Y(inc.min, labelD.pos, ray = maxlength, inc.min, box.range) * label.prop



      if (radlab) { # radlab écrit les labels sous forme d'étoile
        for (label in 1:length(labelD.pos)) {
          if (labelD.pos[label] > 0 && labelD.pos[label] < 180)
            labelsrt <- (90 - labelD.pos[label])
          else
            labelsrt <- (270 - labelD.pos[label])

          text(xpos[label], ypos[label], labelsD[label], cex = par("cex.axis"), srt = labelsrt)
        }
      }
      else {  # boxed.labels(xpos, ypos, labelsD, ypad = 0.7, border = FALSE,  cex =  par("cex.axis"))
        for (label in 1:length(labelD.pos)) {
          text(xpos[label], ypos[label], labelsD[label], cex = par("cex.axis"), srt = 0)
        }
      }
    }


     par(xpd = FALSE)
}

#' Place des points I et D dans un repère Lambert avec un data.frame
#' @param data data.frame avec les variables $I d'inclinaison et $D de déclinaison
#' @param pt.names  Correspond à la liste des noms des points. Laissée vide n'affiche rien. Si on met pt.names = "", cela affiche les noms
#' @param label.pos  Séquence de valeur à afficher en I et D, voir fonction lambert.ID.grid
#' @param point.symbols défini la forme de points, correspond exactement au pch de la fonction points()
#' @param pch permet de changer la forme du symbole
#' @param show.grid permet d'afficher une grille en toile d'araigné sur le fond. Mettre à FALSE, si on superpose des diagrammes
#' @param show.grid.labels permet de changer échelle des graduations
#' @param inc.lim permet de restreindre l'affichage sur une étendue d'inclinaison. Ex: inc.lim = c(45, 90). Laissée à NULL, le diagramme s'addapte aux données
#' @param dec.min permet de restreindre l'étendue en déclinaison, borne minimale
#' @param dec.max permet de restreindre l'étendue en déclinaison, borne maximale
#' @param new permet d'initialiser la sortie graphique. Mettre à FALSE, si on superpose des diagrammes.
#' @export
lambert <- function (data , pt.names = NULL, labels = NA, label.pos = NULL,
                             radlab = FALSE, start = 0, clockwise = TRUE,
                             label.prop = 1.1, main = "", xlab = "", ylab = "", line.col = par("fg"),
                             lty = par("lty"), lwd = par("lwd"), mar = c(2, 2, 3, 2),
                             show.grid = TRUE, show.grid.labels = 10, show.radial.grid = TRUE,
                             grid.col = "gray", grid.bg = "transparent",
                             grid.left = FALSE, grid.unit = NULL, point.symbols = 1, pt.col = par("fg"), bg = pt.col,
                             inc.lim = NULL, radial.labels = NULL,
                             boxed.radial = TRUE, poly.col = NA,
                             dec.min = -90, dec.max = 270, new = TRUE, pch = 21, ...)
{
  if (is.null(pt.names))
    name <- NULL
  else
    name <- data$name

  lambert.ID (data$I, data$D , pt.names = name, labels = labels, label.pos = label.pos,
                               radlab = radlab, start = start, clockwise = clockwise,
                               label.prop = label.prop, main = main, xlab = xlab, ylab = ylab, line.col = line.col,
                               lty = lty, lwd = lwd, mar = mar,
                               show.grid = show.grid, show.grid.labels = show.grid.labels, show.radial.grid = show.radial.grid,
                               grid.col = grid.col, grid.bg = grid.bg,
                               grid.left = grid.left, grid.unit = grid.unit, point.symbols = point.symbols, pt.col = pt.col, bg = pt.col,
                               inc.lim = inc.lim, radial.labels = radial.labels,
                               boxed.radial = boxed.radial, poly.col = poly.col,
                               dec.min = dec.min, dec.max = dec.max, new = new, pch = pch, ...)
}

#' Place des points I et D dans un repère Lambert
#' @export
lambert.ID <- function (inc, dec , pt.names = NA, labels = NA, label.pos = NULL,
                                 radlab = FALSE, start = 0, clockwise = TRUE,
                                 label.prop = 1.1, main = "", xlab = "", ylab = "", line.col = par("fg"),
                                 lty = par("lty"), lwd = par("lwd"), mar = c(2, 2, 3, 2),
                                 show.grid = TRUE, show.grid.labels = 10, show.radial.grid = TRUE,
                                 grid.col = "gray", grid.bg = "transparent",
                                 grid.left = FALSE, grid.unit = NULL, point.symbols = 1, pt.col = par("fg"), bg = pt.col,
                                 inc.lim = NULL, radial.labels = NULL,
                                 boxed.radial = TRUE, poly.col = NA,
                                 dec.min = -90, dec.max = 270, new = TRUE, pch = 21, ...)
{

  if (is.null(inc.lim))
    inc.lim <- range(abs(inc))



  #lambert.ID(inclinaisons = inclinaisons, declinaisons = declinaisons, pt.names = pt.names, inc.lim = inc.lim,
  #           dec.min = dmin, dec.max = dmax, main = main, label.pos = label.pos, pt.col = pt.col, bg = bg, show.grid = show.grid, new = new)
  if (show.grid == TRUE) {

    if (is.null(label.pos)) {
      label.pos$I = seq(inc.lim[1], inc.lim[2], by=10)
      label.pos$D = seq(dec.min, dec.max, by=10)
    }

    lambert.ID.grid(main = main, label.pos = label.pos, labels = labels,
                         radlab = radlab,  start = start,
                         clockwise = clockwise, show.radial.grid = show.radial.grid,
                         dec.min = dec.min, dec.max = dec.max, inc.min = inc.lim[1] , inc.max = inc.lim[2], new  = new )
    if (new == TRUE)
      new <- FALSE
  }
  lambert.ID.point(inc, dec, dec.min = dec.min, dec.max = dec.max, inc.lim = inc.lim,
                        pt.names = pt.names, pt.col = pt.col, bg = bg, pch = pch, new = new)


}

#' Place des points I et D dans un repère Lambert et dessine le symbole
#' @export
lambert.ID.position <- function (data, declinaisons = NULL, position = "P",  pt.names = NA, labels = NA, label.pos = NULL,
                        radlab = FALSE, start = 0, clockwise = TRUE,
                        label.prop = 1.1, main = "Position auto", xlab = "", ylab = "", line.col = par("fg"),
                        lty = par("lty"), lwd = par("lwd"), mar = c(2, 2, 3, 2),
                        show.grid = TRUE, show.grid.labels = 10, show.radial.grid = TRUE,
                        grid.col = "gray", grid.bg = "transparent",
                        grid.left = FALSE, grid.unit = NULL, point.symbols = 1, pt.col = par("fg"), bg = pt.col,
                        inc.lim = c(0, 90), radial.labels = NULL,
                        boxed.radial = TRUE, poly.col = NA, add = FALSE,
                        dec.min = -90, dec.max = 90, new = TRUE, ...)
{

  if (is.null(declinaisons) ) {
    inclinaisons <- data$I
    declinaisons <- data$D
  }
  else {
    inclinaisons <- data
    declinaisons <- declinaisons
  }



  if (length(position) < length(inclinaisons))
    position <- rep(position, length(inclinaisons))

  if (is.null(inc.lim))
    inc.lim <- range(abs(inclinaisons))

  if (is.null(label.pos)) {
    lab.pos <- NULL
    lab.pos$I = seq(inc.lim[1], inc.lim[2], by = 20)
    lab.pos$D = seq(dec.min, dec.max, by = 10)
  } else {
    lab.pos <- label.pos
  }

  ne <- new
  if (length(which(position=="P")) > 0) {
    lambert.ID(inclinaisons[which(position=="P")], declinaisons[which(position=="P")],
             pt.names = pt.names[which(position=="P")], label.pos = lab.pos, main = main, show.grid = show.grid,
             inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, pt.col = pt.col,  pch = 23, new = ne)
    ne <- FALSE
    show.grid <- FALSE
    main <- ""
  }

  if (length(which(position=="C")) > 0) {
    lambert.ID(inclinaisons[which(position=="C")], declinaisons[which(position=="C")],
             pt.names = pt.names[which(position=="C")], label.pos = lab.pos, main = main, show.grid = show.grid,
             inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, pt.col = pt.col,  pch = 22, new = ne)
    ne <- FALSE
    show.grid <- FALSE
    main <- ""
  }


  if (length(which(position=="D")) > 0) {
    lambert.ID(inclinaisons[which(position=="D")], declinaisons[which(position=="D")],
             pt.names = pt.names[which(position=="D")], label.pos = lab.pos, main = main, show.grid = show.grid,
             inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, pt.col = pt.col, pch = 24, new = ne)
  }
}

#' Place des points I et D dans un repère Lambert
#' @export
lambert.ID.point <- function (inc, dec , pt.names = NA, labels = NA, label.pos = NULL,
                             start = 0, clockwise = TRUE,
                             lty = par("lty"),  mar = c(2, 2, 3, 2),
                             pt.col = par("fg"), bg = pt.col,
                             inc.lim = NULL,
                             dec.min = -90, dec.max = 270, new = TRUE, pch = 21, ...)
{

  maxlength <- 100 # la valeur n'a pas d'influence, la fonction plot() calcul le reste
  # radlab écrit les labels sous forme d'étoile
  if (is.null(inc.lim))
    inc.lim <- range(abs(inc))

  # Sélection des échantillons visibles dans le range
  index.supprim <- NULL
  for (i in 1:length(inc)) {
    if (abs(inc[i])< inc.lim[1] || abs(inc[i])> inc.lim[2] || dec[i]<dec.min || dec[i]>dec.max )
      index.supprim <- cbind(index.supprim, c(i))

  }
  if (!is.null(index.supprim)) {
    inc <- inc[-index.supprim]
    dec <- dec[-index.supprim]
  }


  inc.min = inc.lim[1]
  inc.max = inc.lim[2]

  nbpoints <- length(inc)

  if (clockwise == FALSE)
    dec <- -dec
  if (start)
    dec <- dec + start

  box.range <- find.extremum(inc.min, inc.max, dec.min, dec.max)

  oldpar <- par("xpd", "mar", "pty")

    # setting up coord. system
    if (new == TRUE) {
      par(mar = mar, pty = "s")

      maxlength <- 100 # la valeur n'a pas d'influence, la fonction plot() calcul le reste
      plot(c(-maxlength, maxlength), c(-maxlength, maxlength), type = "n", axes = FALSE, new = new)
    }

 # par(xpd = TRUE)

  if (length(pch) < nbpoints)
    pch <- rep(pch, length.out = nbpoints)

  if (length(pt.col) < nbpoints)
    pt.col <- rep(pt.col, length.out = nbpoints)


  xpos <- X(inc, dec, ray = maxlength, i.min = inc.min, box.range)
  ypos <- Y(inc, dec, ray = maxlength, i.min = inc.min, box.range)
  # print points names
  if (length(pt.names)>0 )
    for (i in 1: nbpoints )
      text(xpos[i], ypos[i], pt.names[i] , srt = 0, cex = par("cex"))



  # pch = 0, cercle
  # pch = 1, rond
  # pch = 2, triangle
  # pch = 3, plus
  # pch = 4, croix
  # pch = 5, losange
  # pch = 6, triangle vers le bas
  # pch = 7, carré avec croix
  # pch = 8, étoile
  # pch = 9, losange avec plus
  # pch = 10, cercle avec plus
  # pch = 11, triangles hauts et bas
  # pch = 12, carré avec plus
  # pch = 13, cercle avec croix
  # pch = 14, carré et triangle vers le bas
  # pch = 15, carré plein
  # pch = 16, cercle plein
  # pch = 17, triangle plein vers le haut
  # pch = 18, losange plein
  # pch = 19, cercle solide
  # pch = 20, petit rond plein
  # pch = 21, cercle plein bleu, accepte l'argument bg
  # pch = 22, carré plein bleu
  # pch = 23, losange plein bleu
  # pch = 24, triangle plein, pointe vers le haut bleu
  # pch = 25, triangle plein, pointe vers la bas bleu

  if (!is.null(pt.col))
    for (i in 1: nbpoints ) {
      if(pt.col[i] != "transparent") {
        if (inc[i]<0) {
          bg.col <- gray(0.95)
          points(xpos[i], ypos[i], pch = pch[i],  col = pt.col[i], bg = bg.col,...)
        }
        else
          points(xpos[i], ypos[i], pch = pch[i],  col = pt.col[i], bg = pt.col[i],  ...)
      }
    }


  invisible(oldpar)
}
## Plot lambert circle ----

## lambert.ID.circle.points
#' Calcul les points pour tracer un cercle
#' @param i.mean inclinaison du centre du cercle
#' @param d.mean déclinaiosn du centre du point
#' @param delta angle d'ouverture du cercle
#' @return une liste de I et D en degrés
#' @export
lambert.ID.circle.points <- function (i.mean, d.mean, delta)
{
  Imoy <- i.mean*pi/180
  Dmoy <- d.mean*pi/180
  Delta <- delta*pi/180
  pt <- NULL
  pt$I <- NULL
  pt$D <- NULL
  for (k in seq(0, pi, length.out = 180)) {
    I.tmp <- n_arcsin(sin(Imoy)*cos(Delta)+cos(k)*cos(Imoy)*sin(Delta))
    if (abs(Imoy)==(pi/2)) {
      D.tmp <- k
    } else {
      D.tmp <- Dmoy + n_arccos( (cos(Delta)-sin(Imoy)*sin(I.tmp))/(cos(Imoy)*cos(I.tmp)) )
    }

    pt$I<- c(pt$I , I.tmp*180/pi)
    pt$D <- c(pt$D, D.AM(D.tmp*180/pi))
  }
  for (k in seq(pi, 2*pi, length.out = 180)) {
    I.tmp <- n_arcsin(sin(Imoy)*cos(Delta)+cos(k)*cos(Imoy)*sin(Delta))
    if (abs(Imoy)==(pi/2)) {
      D.tmp <- k
    } else {
      D.tmp <- Dmoy - n_arccos( (cos(Delta)-sin(Imoy)*sin(I.tmp))/(cos(Imoy)*cos(I.tmp)) )
    }
    pt$I <- c(pt$I, I.tmp*180/pi)
    pt$D <- c(pt$D, D.AM(D.tmp*180/pi))
  }

  return(pt)
}

#' Trace un cercle sur le diagramme lambert
#' @param  col  définie la couleur de la ligne
#' @param absolue permet de tracer la partie inclinaison négative des cercles
#' @export
lambert.ID.circle <- function (i.mean, d.mean, delta, inc.lim = NULL, dec.min = -90, dec.max = 270,
                                    col = par("fg"), clockwise = TRUE, absolue = TRUE, ...)
{

  pt <- lambert.ID.circle.points(i.mean, d.mean, delta)

  maxlength <- 100 # la valeur n'a pas d'influence, la fonction plot() calcul le reste
  # radlab écrit les labels sous forme d'étoile
  if (is.null(inc.lim))
    inc.lim <- range(abs(pt$I))


  inc.min = inc.lim[1]
  inc.max = inc.lim[2]


  if (clockwise == FALSE)
    pt$D <- -pt$D


  box.range <- find.extremum(inc.min, inc.max, dec.min, dec.max)

  xpos <- X(pt$I, pt$D, ray = maxlength, i.min = inc.min, box.range)
  ypos <- Y(pt$I, pt$D, ray = maxlength, i.min = inc.min, box.range)

  before.bad.condition <- FALSE
  x0 <- xpos[1]
  y0 <- ypos[1]


  for (i in 1: length(pt$I) ) {
    if (absolue)
      bad.condition <- (abs(pt$I[i])< inc.lim[1] || abs(pt$I[i])> inc.lim[2] || pt$D[i]<dec.min || pt$D[i]>dec.max )
    else
      bad.condition <- (pt$I[i]< inc.lim[1] || pt$I[i]> inc.lim[2] || pt$D[i]<dec.min || pt$D[i]>dec.max )

    if (before.bad.condition == TRUE && bad.condition == FALSE) {
      x0 <- xpos[i]
      y0 <- ypos[i]
    }

    if (before.bad.condition == FALSE && bad.condition == FALSE) {
      x1 <- xpos[i]
      y1 <- ypos[i]
      if (pt$I[i]<0)
        segments(x0, y0, x1, y1, lty = 3, col = col, ...)
      else
        segments(x0, y0, x1, y1, lty = 1, col = col, ...)

      x0 <- xpos[i]
      y0 <- ypos[i]
    }

    before.bad.condition <- bad.condition

  }

}

#' Trace un cercle avec les champs de l'inclinaison possible
#' @param  field définit la zone possible du champ
#' @export
lambert.ID.field <- function (data, dec = NULL , pt.names = NA, field = c(50, 75), inc.lim = c(0, 90), dec.min = -90, dec.max = 270, col = par("fg"),
                              main = "", pt.col = par("fg"), bg = pt.col, label.pos= NULL, show.grid = TRUE, new = TRUE)
{
  #old.par <- par(no.readonly = TRUE) # all par settings which


  if (is.null(label.pos)) {
    label.pos$I = NA
    label.pos$D = seq(-90, 180, by=90)
  }


  if (is.null(dec) ) {
    inc <- data$I
    dec <- data$D
    if (is.null(pt.names))
      pt.names <- data$name
  }
  else {
    inc <- data
    dec <- dec
  }


  if (is.null(inc.lim))
    inc.lim <- range(abs(inc))

  if (is.null(label.pos)) {
    lab.pos$I = seq(inc.lim[1], inc.lim[2], by = 20)
    lab.pos$D = seq(dec.min, dec.max, by = 10)
  }



  if (show.grid == TRUE)
    lambert.ID.grid(main = main, label.pos = label.pos,
                        start = 0, dec.min = dec.min, dec.max = dec.max, inc.min = inc.lim[1], inc.max = inc.lim[2], new = new )

  lambert.ID.point(inc, dec, inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, pt.names = pt.names, pt.col = pt.col, bg = bg, new = FALSE)



  hinf <- 90 - field[1]
  hsup <- 90 - field[2]

  lambert.ID.circle(0, 180, hinf, col =  adjustcolor( col, alpha.f = 0.5) , inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, absolue = FALSE)
  lambert.ID.circle(0, 0, hinf, col = adjustcolor( col, alpha.f = 0.5) , inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, absolue = FALSE)
  lambert.ID.circle(0, 90, hinf, col = adjustcolor( col, alpha.f = 0.5) , inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, absolue = FALSE)
  lambert.ID.circle(0, -90, hinf, col = adjustcolor( col, alpha.f = 0.5) , inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, absolue = FALSE)
  lambert.ID.circle(90, -90, hinf, col = adjustcolor( col, alpha.f = 0.5) , inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, absolue = FALSE)


  lambert.ID.circle(0, 180, hsup, col = adjustcolor( col, alpha.f = 0.5) , inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, absolue = FALSE)
  lambert.ID.circle(0, 0, hsup, col = adjustcolor( col, alpha.f = 0.5) , inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, absolue = FALSE)
  lambert.ID.circle(0, 90, hsup, col = adjustcolor( col, alpha.f = 0.5) , inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, absolue = FALSE)
  lambert.ID.circle(0, -90, hsup, col = adjustcolor( col, alpha.f = 0.5) , inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, absolue = FALSE)
  lambert.ID.circle(90, -90, hsup, col = adjustcolor( col, alpha.f = 0.5) , inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, absolue = FALSE)

  #on.exit(par(old.par))
}

#' Place des points X, Y et Z dans un repère Lambert
#' @export
lambert.XYZ <- function( X, Y , Z, pt.names = NA, labels = NA, label.pos = NULL,
                                   radlab = FALSE, start = 0, clockwise = TRUE,
                                   label.prop = 1.1, main = "", xlab = "", ylab = "", line.col = par("fg"),
                                   lty = par("lty"), lwd = par("lwd"), mar = c(2, 2, 3, 2),
                                   show.grid = TRUE, show.grid.labels = 10, show.radial.grid = TRUE,
                                   grid.col = "gray", grid.bg = "transparent",
                                   grid.left = FALSE, grid.unit = NULL, point.symbols = 1, pt.col = "blue3", bg = pt.col,
                                   inc.lim = NULL, radial.labels = NULL,
                                   boxed.radial = TRUE, poly.col = NA, add = FALSE,
                                   dec.min = -90, dec.max = 270, new = TRUE, pch = 21, ...)
{
  I <- to.polar.I(X, Y, Z)
  D <- to.polar.D(X, Y, Z)
  lambert.ID(I, D, pt.names = pt.names, label.pos = label.pos, main = main, show.grid = show.grid,
             inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, pt.col = pt.col,  pch = 23, new = new)

}

#' Place des points X, Y et Z dans un repère Lambert avec un data.frame
#' @export
lambert.XYZ.specimen <- function( Data, Y , Z, pt.names = "", labels = NA, label.pos = NULL,
                         radlab = FALSE, start = 0, clockwise = TRUE,
                         label.prop = 1.1, main = NULL, xlab = "", ylab = "", line.col = "blue3",
                         lty = par("lty"), lwd = par("lwd"), mar = c(2, 2, 3, 2),
                         show.grid = TRUE, show.grid.labels = 10, show.radial.grid = TRUE,
                         grid.col = "lightgray", grid.bg = "transparent",
                         grid.left = FALSE, grid.unit = NULL, point.symbols = 1, pt.col = "blue", bg = pt.col,
                         inc.lim = c(0, 90), radial.labels = NULL,
                         boxed.radial = TRUE, poly.col = NA, add = FALSE,
                         dec.min = -90, dec.max = 270, new = TRUE, pch = 21, ...)
{
  eta <- NULL
  if(is.data.frame(Data)) {
    X <- Data$X
    Y <- Data$Y
    Z <- Data$Z
    if (is.null(pt.names))
      eta <- Data$step
    if (is.null(main))
      main <- as.character(Data$name[1])
  } else {
    X <- Data
    if (!is.null(pt.names))
      eta <- pt.names
  }
  I <- to.polar.I(X, Y, Z)
  D <- to.polar.D(X, Y, Z)

  if (is.null(label.pos)) {
    label.pos$I = c(90)
    label.pos$D = seq(-90, 270, by=90)
  }
  lambert.ID(I, D, pt.names = pt.names, label.pos = label.pos, main = main, show.grid = show.grid, type = "l",
                  grid.col = grid.col,
             inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, line.col = line.col, pt.col = pt.col, bg= par("fg"), new = new)
  lambert.ID(I, D, pt.names = "", label.pos = label.pos, main = "", show.grid = FALSE,
                  inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, line.col = line.col, pt.col = pt.col,  pch = 21, new = FALSE)
}

#' Place des points I et D dans un repère Lambert avec un data.frame
#' @export
lambert.ID.specimen <- function( Data, D , pt.names = "", labels = NA, label.pos = NULL,
                                       radlab = FALSE, start = 0, clockwise = TRUE,
                                       label.prop = 1.1, main = NULL, xlab = "", ylab = "", line.col = "blue3",
                                       lty = par("lty"), lwd = par("lwd"), mar = c(2, 2, 3, 2),
                                       show.grid = TRUE, show.grid.labels = 10, show.radial.grid = TRUE,
                                       grid.col = "lightgray", grid.bg = "transparent",
                                       grid.left = FALSE, grid.unit = NULL, point.symbols = 1, pt.col = "blue3", bg = pt.col,
                                       inc.lim = c(0, 90), radial.labels = NULL,
                                       boxed.radial = TRUE, poly.col = NA, add = FALSE,
                                       dec.min = -90, dec.max = 270, new = TRUE, pch = 21, ...)
{
  eta <- NULL
  if(is.data.frame(Data)) {
    I <- Data$I
    D <- Data$D
    if (is.null(pt.names))
      eta <- Data$step
    if (is.null(main))
      main <- as.character(Data$name[1])
  } else {
    I <- Data
    if (!is.null(pt.names))
      eta <- pt.names
  }

  if (is.null(label.pos)) {
    label.pos$I = c(90)
    label.pos$D = seq(-90, 270, by=90)
  }
  lambert.ID(I, D, pt.names = pt.names, label.pos = label.pos, main = main, show.grid = show.grid, type = "l",
                  grid.col = grid.col,
                  inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, line.col = line.col, pt.col = pt.col, bg= par("fg"), new = new)
  lambert.ID(I, D, pt.names = "", label.pos = label.pos, main = "", show.grid = FALSE,
                  inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, line.col = line.col, pt.col = pt.col,  pch = 21, new = FALSE)
}

#' Trace un diagramme de Zijderveld type 1
#' @param Data soit une data.frame avec les mesures X, Y et Z, soit que les valeurs de X
#' @param Y les valeurs de Y, si Data n''est pas une data.frame
#' @param Z les valeurs de Z, si Data n''est pas une data.frame
#' @param panel.first = grid() : affiche une grille
#' @param pt.name = "": n'affiche rien. pt.name = NULL: affiche les étapes si Data est une data.frame
#' @param legend.pos = NULL : n'affiche pas la legende. legend.pos = "topleft" affiche en haut à gauche
#' @export
zijderveld1<- function(Data, Y = NULL, Z = NULL, pt.names = "", main = NULL, panel.first = NULL, pt.col = c("forestgreen", "blue3"),
                            isometric = TRUE, ylim = NULL, legend.pos = NULL, legend.txt = c("(Y, X)", "(Y, Z)"), new = TRUE, ...)
{
  eta <- NULL
  if(is.data.frame(Data)) {
    X <- Data$X
    Y <- Data$Y
    Z <- Data$Z
    if (is.null(pt.names))
      eta <- Data$step
    else
      eta <- pt.names

    if (is.null(main))
      main <- as.character(Data$name[1])
  } else {
    X <- Data
    if (!is.null(pt.names))
      eta <- pt.names
  }

  if (isometric == TRUE) {
     asp <- 1
  } else {
     asp <- NA
  }

  if (is.null(ylim)) {
    Y.r <- range(X)
    Z.r <- range(Z)
    ylim <- c(min(Y.r[1], -Z.r[2]), max(Y.r[2], -Z.r[1]))
    if (isometric == TRUE) {
      if (ylim[1]>0)
        ylim[1]<-0
      if (ylim[2]<0)
        ylim[2]<-0

    }

  }

  if (new == TRUE) {
    plot(Y, X, type = "o", pch = 21, main = main, col = pt.col[1],  bg = adjustcolor( pt.col[1], alpha.f = 0.8), axes = FALSE,
         panel.first = panel.first, xlab = "", ylab = "", ylim = ylim, asp = asp, xaxt="n", yaxt="n", new = new, ...)



    if (!is.null(legend.pos))
      legend(legend.pos, legend.txt, pch = c(19, 21), col = pt.col, bg = c(par("bg"), adjustcolor( pt.col[1], alpha.f = 0.8), adjustcolor( pt.col, alpha.f = 0.05)),
             box.col = par("bg"), title = "")

    ax1 <- axis(1, pos = 0,  col = "darkgray")
    ax2 <- axis(2, pos = 0,  col = "darkgray") # Ordonnées

    text(0, ax2[length(ax2)], "+X", col = "gray5", adj = c(-.5, 1), cex = par("cex.lab"))
    text(0, ax2[1], "+Z", col = "gray5", adj = c(-.5, 0), cex = par("cex.lab"))
    text( ax1[length(ax1)], 0, "+Y", col = "gray5", adj = c(1, -.5), cex = par("cex.lab"))

  } else {
    lines(Y, X, type = "o", pch = 21, col = pt.col[1], bg = adjustcolor( pt.col[1], alpha.f = 0.8), ...)
  }

  lines(Y, -Z, type = "o", pch = 21, col = pt.col[2], bg = adjustcolor( pt.col[2], alpha.f = 0.05), ...)

  text(jitter(Y, 5, amount = 0), jitter(X, 5, amount = 0), eta, cex = par("cex.lab"))

}

#' Trace un diagramme de Zijderveld type 2
#' @export
zijderveld2<- function(Data, Y = NULL, Z = NULL, pt.names = "", main = NULL, panel.first = NULL, pt.col = c("forestgreen", "blue3"), isometric = TRUE,
                            ylim = NULL, legend.pos = NULL, new = TRUE, ...)
{
  eta <- NULL
  if(is.data.frame(Data)) {
    X <- Data$X
    Y <- Data$Y
    Z <- Data$Z
    if (is.null(main))
      main <- as.character(Data$name[1])
    if (is.null(pt.names))
      eta <- Data$step
  } else {
    X <- Data
    if (!is.null(pt.names))
      eta <- pt.names
  }

  if (isometric == TRUE) {
    asp <- 1
  } else {
    asp <- NA
  }

  if (is.null(ylim)) {
    Y.r <- range(Y)
    Z.r <- range(Z)
    ylim <- c(min(Y.r[1], -Z.r[2]), max(Y.r[2], -Z.r[1]))
    if (isometric == TRUE) {
      if (ylim[1]>0)
        ylim[1]<-0
      if (ylim[2]<0)
        ylim[2]<-0
    }
  }

  if (new == TRUE) {
    plot(-X, Y, main = main, type = "o", pch = 21, col = pt.col[1], bg = adjustcolor( pt.col[1], alpha.f = 0.8), axes = FALSE,
         panel.first = panel.first, xlab = "", ylab = "", ylim = ylim, asp = asp, new = new)

    text(jitter(-X, 5, amount = 0), jitter(Y, 5, amount = 0), eta, cex = par("cex.lab"))

    if (!is.null(legend.pos))
      legend(legend.pos, c("(-X, Y)", "(-X, Z)"), pch = c(19, 21), col = pt.col, bg = c(par("bg"), adjustcolor( pt.col[1], alpha.f = 0.8), adjustcolor( pt.col, alpha.f = 0.05)),
             box.col = par("bg"), title = "")

    ax1 <- axis(1, pos = 0, col = "darkgray")
    ax2 <- axis(2, pos = 0, col = "darkgray") # Ordonnées

    text(0, ax2[length(ax2)], "+Y", col = "gray5", adj = c(-.5, 1), cex = par("cex.lab"))
    text(0, ax2[1], "+Z", col = "gray5", adj = c(-.5, 0), cex = par("cex.lab"))
    text( ax1[length(ax1)] , 0, "-X", col = "gray5", adj = c(1, -.5), cex = par("cex.lab"))


  } else {
    lines(-X, Y, type = "o", pch = 21, col = pt.col[1], bg = adjustcolor( pt.col[1], alpha.f = 0.8), ...)
  }

  lines(-X, -Z, type = "o", pch = 21, col = pt.col[2], bg = adjustcolor( pt.col[2], alpha.f = 0.05))

  text(jitter(-X, 5, amount = 0), jitter(Y, 5, amount = 0), eta, cex = par("cex.lab"))

}

# Repliement ----

#' repliement
#' Calcul le repliement pour les trois position à plat "P", de chant "C" et debout "D"
#' pour les matériaux déplacés, la référence est le nord donc -angD
#' @param data valeur de l'inclinaison, ou une data.frame en degrés
#' @param dec en degrés
#' @param position trois valeurs posible à plat "P", de chant "C" et debout "D"
#' @return en degrés
#' @export
repliement <- function (data, dec = NULL, aim = 1, name = NULL, number = NULL,  position = "P")
{

  if (is.null(dec) ) {
    inc <- data$I
    dec <- data$D
    aim <- data$F
    nom <- data$name
    num <- data$number
  }
  else {
    inc <- data
    dec <- dec
    if (is.null(name))
      nom <- rep("", length(inc))
    else
      nom <- as.character(name)

    if (is.null(number))
      num <- c(1:length(inc))
    else
      num <- number
  }
  if (length(aim) < length(inc))
    aim <- rep(aim, length(inc))

  if (length(position) < length(inc))
    position <- rep(position, length(inc))

  X <- to.cartesian.X(inc, dec, aim)
  Y <- to.cartesian.Y(inc, dec, aim)
  Z <- to.cartesian.Z(inc, dec, aim)

  res <- NULL

 # Carottage à plat, repère conventionnel : Calcul des angles I-D pour les trois positions
  for (i in 1: length(inc)) {

    res$name <- c(res$name, as.character(nom[i]))
    res$number <- c(res$number, num[i])
    res$F <- c(res$F, aim[i])

    if (position[i] == "O") { # origine
      res$D <- c(res$D, angleD( X[i], Y[i] ) /pi * 180 )
      res$I <- c(res$I, n_arcsin( Z[i]/aim[i] ) /pi * 180)

      res$P <- c(res$P, "O")
    }

    if (position[i] == "P") { # à plat
      t.d <- angleD(abs(X[i]), -Y[i]) /pi * 180
      t.i <- n_arcsin(Z[i]/aim[i]) /pi * 180

      if (X[i]>= 0 && Z[i]>=0 ) {
        t.d <- -angleD(X[i], Y[i]) /pi * 180
      }
      if (X[i]> 0 && Z[i]<0 ) {
        t.d <- -angleD(X[i], -Y[i]) /pi * 180
      }

      if (X[i]< 0 && Z[i]>0 ) {
        t.d <- angleD(-X[i], Y[i]) /pi * 180
      }
      if (X[i]< 0 && Z[i]<0 ) {
        t.d <- angleD(-X[i], -Y[i]) /pi * 180
      }

      res$I <- c(res$I, t.i)
      res$D <- c(res$D, t.d )

      res$P <- c(res$P, "P")
    }

    if (position[i] == "C") { #  de chant -->> OK

      t.i <- sign(Z[i]) *abs(n_arcsin( -Y[i]/aim[i])) /pi * 180
      t.d <- angleD( X[i], Z[i]) /pi * 180

      if (X[i]> 0 && Y[i]<0 ) {
        t.d <- -t.d
      }

      if (X[i]< 0 && Y[i]>0 ) {
        t.d <- t.d - 180

      }


      if (X[i]<= 0 && Y[i]<=0 ) {
         t.d <- 180 - t.d
      }

      # normalisation de la déclinaison
      if (t.d> 270)
        t.d <- t.d - 360
      if (t.d < -90)
        t.d <- t.d + 360

      res$I <- c(res$I, t.i )
      res$D <- c(res$D, t.d )
      res$P <- c(res$P, "C")
    }

    if (position[i] == "D") { # debout

      t.i <- sign(Z[i]) * n_arcsin(abs(X[i])/aim[i]) /pi * 180
      t.d <- angleD(Y[i], -Z[i]) /pi * 180

      if (X[i]> 0 && Y[i]<=0 ) {
        t.d <- t.d - 180
      }

      if (X[i]< 0 && Y[i]>0 ) {
        t.d <- - t.d
      }
      if (X[i]< 0 && Y[i]<=0 ) {
        t.d <- 180 - t.d
      }

      res$I <- c(res$I, t.i)
      res$D <- c(res$D, t.d)
      res$P <- c(res$P, "D")
     }

    tmp <- to.cartesian(t.i, t.d, aim = aim[i])
    res$X <- c(res$X, tmp$X)
    res$Y <- c(res$Y, tmp$Y)
    res$Z <- c(res$Z, tmp$Z)
  }

  res.frame <- data.frame( number =  as.numeric(res$number), name = as.character(res$name),
                           I = as.numeric(res$I), D = as.numeric(res$D), F = as.numeric(res$F),
                           X = as.numeric(res$X), Y =as.numeric(res$Y), Z =as.numeric(res$Z), position = res$P, stringsAsFactors = FALSE)

  return(res.frame)
}

#' repliement automatique
#' Recherche la position qui permet d'avoir une position suivant les trois positions possible
#' donnant l'inclinaison la plus proche de la valeur inc.critique (en degrés)
#' @param  data valeur de l'inclinaison, ou une data.frame en degrés
#' @param  dec en degrés
#' @return I et D en degrés
#' @export
repliement.auto <- function (data, dec = NULL, aim = 1, name = NULL, number = NULL, inc.critique = 90)
{

  if (is.null(dec) ) {
    inc <- data$I
    dec <- data$D
    aim <- data$F
    nom <- data$name
    num <- data$number
  }
  else {
    inc <- data
    dec <- dec
    if (is.null(name))
      nom <- rep("", length(inc))
    else
      nom <- as.character(name)

    if (is.null(number))
      num <- c(1:length(inc))
    else
      num <- number
  }

  res <- NULL

  res.P <- repliement(inc, dec, aim, name = nom, number = num, position = "P")
  res.C <- repliement(inc, dec, aim, name = nom, number = num, position = "C")
  res.D <- repliement(inc, dec, aim, name = nom, number = num, position = "D")

  for (i in 1: length(inc)) {
    # inc.sup <- 90
    res.tmp <- NULL
    # Initialisation avec la position à Plat
    res.tmp$I <- res.P$I[i]
    res.tmp$D <- res.P$D[i]
    res.tmp$F <- res.P$F[i]
    res.tmp$X <- res.P$X[i]
    res.tmp$Y <- res.P$Y[i]
    res.tmp$Z <- res.P$Z[i]
    res.tmp$position <- res.P$position[i]

    inc.sup <- abs(inc.critique - abs(res.P$I[i]))

    # Comparaison avec la position de Chant
    if (abs(inc.critique - abs(res.C$I[i])) <= inc.sup) {
      res.tmp$I <- res.C$I[i]
      res.tmp$D <- res.C$D[i]
      res.tmp$F <- res.C$F[i]
      res.tmp$position <- res.C$position[i]
      res.tmp$X <- res.C$X[i]
      res.tmp$Y <- res.C$Y[i]
      res.tmp$Z <- res.C$Z[i]
      inc.sup <- abs(inc.critique - abs(res.C$I[i]))
    }
    # Comparaison avec la position Debout
    if (abs(inc.critique - abs(res.D$I[i])) <= inc.sup) {
      res.tmp$I <- res.D$I[i]
      res.tmp$D <- res.D$D[i]
      res.tmp$F <- res.D$F[i]
      res.tmp$position <- res.D$position[i]
      res.tmp$X <- res.D$X[i]
      res.tmp$Y <- res.D$Y[i]
      res.tmp$Z <- res.D$Z[i]

      #inc.sup <- abs(inc.critique-abs(res.D$I))
    }

    res$I <- c(res$I, res.tmp$I)
    res$D <- c(res$D, res.tmp$D)
    res$F <- c(res$F, res.tmp$F)

    res$name <- c(res$name, nom[i])
    res$number <- c(res$number, num[i])

    res$X <- c(res$X, res.tmp$X)
    res$Y <- c(res$Y, res.tmp$Y)
    res$Z <- c(res$Z, res.tmp$Z)

    res$position <- c(res$position, res.tmp$position)
  }

  res <- data.frame(number = res$number, name = res$name, I = as.numeric(res$I), D = as.numeric(res$D), F = as.numeric(res$F),
                    X = as.numeric(res$X), Y = as.numeric(res$Y), Z = as.numeric(res$Z), position = res$position, stringsAsFactors = FALSE)

  return(res)
}

#' repliement.tranche
#' Recherche la position qui permet d'avoir une position suivant la position debout ou dechant
#' donnant l'inclinaison la plus proche de la valeur inc.critique (en degrés)
#' @param  data valeur de l'inclinaison, ou une data.frame en degrés
#' @param  dec en degrés
#' @return I et D en degrés
#' @export
repliement.tranche <- function (data, dec = NULL, aim = 1,  name = NULL, number = NULL, inc.critique = 90)
{

  if (is.null(dec) ) {
    inc <- data$I
    dec <- data$D
    aim <- data$F
    nom <- data$name
    num <- data$number
  }
  else {
    inc <- data
    dec <- dec
    if (is.null(name))
      nom <- rep("", length(inc))
    else
      nom <- as.character(name)

    if (is.null(number))
      num <- c(1:length(inc))
    else
      num <- number
  }

  res <- NULL

  res.C <- repliement(inc, dec, aim, name = nom, number = num, position = "C")
  res.D <- repliement(inc, dec, aim, name = nom, number = num, position = "D")

  for (i in 1: length(inc)) {

    res.tmp <- NULL
    # Initialisation avec la position de Chant

    res.tmp$I <- res.C$I[i]
    res.tmp$D <- res.C$D[i]
    res.tmp$F <- res.C$F[i]
    res.tmp$position <- res.C$position[i]
    res.tmp$X <- res.D$X[i]
    res.tmp$Y <- res.D$Y[i]
    res.tmp$Z <- res.D$Z[i]
    inc.sup <- abs(inc.critique - abs(res.C$I[i]))

    # Comparaison avec la position Debout
    if (abs(inc.critique - abs(res.D$I[i])) <= inc.sup) {
      res.tmp$I <- res.D$I[i]
      res.tmp$D <- res.D$D[i]
      res.tmp$F <- res.D$F[i]
      res.tmp$position <- res.D$position[i]
      res.tmp$X <- res.D$X[i]
      res.tmp$Y <- res.D$Y[i]
      res.tmp$Z <- res.D$Z[i]
     }

    res$number <- c(res$number, num[i])
    res$name <- c(res$name, nom[i])
    res$I <- c(res$I, res.tmp$I)
    res$D <- c(res$D, res.tmp$D)
    res$F <- c(res$F, res.tmp$F)
    res$X <- c(res$X, res.tmp$X)
    res$Y <- c(res$Y, res.tmp$Y)
    res$Z <- c(res$Z, res.tmp$Z)
    res$position <- c(res$position, res.tmp$position)
  }
  res.Final <- data.frame( res , stringsAsFactors = FALSE)
  return(res.Final)
}

# Fonction sur fichiers ----

#' read.AM.mesures
#' Lecture des mesures d'un fichier AM
#' @param encoding  Pour les fichiers du magnétomètre, il faut "macroman" -> difficle à connaitre, peut être "latin1" ou "utf8".
#' @return une data.frame
#' @export
read.AM.mesures <- function(file.AM, encoding = "macroman")
{
  # Lecture et Copy du fichier
  lin<- NULL
  fil <- file(file.AM, "r", encoding = encoding)
  lin <- readLines(fil)
  close(fil)
  # Recherche position-ligne des noms
  g <- NULL
  for (i in 1:length(lin))
    if (grepl("Id:", lin[i])==TRUE)
      g<-c(g, i)

  # Lecture mesure par nom
  lname <- trimws(substr(lin[g],4, 15))
  list.mesure <- NULL
  lmes <- NULL

  first.mesure <- NULL
  last.mesure <- NULL
  for (i in 1:length(g)) {
    first.mesure <- g[i]+4
    if (i<length(g))
      last.mesure <- g[i+1]-2
    else
      last.mesure<- length(lin)

    tEtap <- NULL
    tEtap.val <- NULL
    tX <- NULL
    tY <- NULL
    tZ <- NULL
    tI <- NULL
    tD <- NULL
    tF <- NULL
    tQual <- NULL
    tApp <- NULL
    tSuscep <- NULL
    for (j in first.mesure:last.mesure) {
      lEtap <- substr(lin[j], 3, 7)
      lEtap.val <- substr(lin[j], 3, 5)
      lX <- substr(lin[j], 9, 18)
      lY <- substr(lin[j], 20, 29)
      lZ <- substr(lin[j], 31, 40)
      lQual <- substr(lin[j], 43, 44)
      lApp <- substr(lin[j], 46, 46)
      lSuscep <- substr(lin[j], 48, 52)
      if ( is.na(as.numeric(lX)) || is.na(as.numeric(lY)) || is.na(as.numeric(lY)) ) {
        lX <- 0
        lY <- 0
        lZ <- 0
        lI <- 0
        lD <- 0
        lF <- 0
      }
      else {
        lI <- to.polar.I(as.numeric(lX), as.numeric(lY), as.numeric(lZ))
        lD <- to.polar.D(as.numeric(lX), as.numeric(lY), as.numeric(lZ))
        lF <- to.polar.F(as.numeric(lX), as.numeric(lY), as.numeric(lZ))
      }

      tEtap <- c(tEtap, lEtap)
      tEtap.val <- c(tEtap.val, lEtap.val)
      tX <- c(tX, lX)
      tY <- c(tY, lY)
      tZ <- c(tZ, lZ)
      tI <- c(tI, lI)
      tD <- c(tD, lD)
      tF <- c(tF, lF)
      tQual <- c(tQual, lQual)
      tApp <- c(tApp, lApp)
      tSuscep <- c(tSuscep, lSuscep)
      lmes <- data.frame(number = i, name = lname[i], step = tEtap, step.value = as.numeric(tEtap.val),
                         X = as.numeric(tX), Y = as.numeric(tY), Z = as.numeric(tZ),
                         I = as.numeric(tI), D = as.numeric(tD), F = as.numeric(tF),
                         Quality = as.numeric(tQual), App = tApp, Suscep = as.numeric(tSuscep), stringsAsFactors = FALSE)
    }
    list.mesure <- rbind(list.mesure, lmes)

  }

  return(list.mesure)
}

#' Lecture des infos sur mesures d'un fichier AM
#' @return une data.frame avec les infos sur les spécimens
#' @export
read.AM.info <- function (file.AM, encoding="macroman")
{

  # Lecture et Copy du fichier
  lin<- NULL
  fil <- file(file.AM, "r", encoding=encoding)
  lin <- readLines(fil)
  close(fil)
  # Recherche position-ligne des noms
  g <- NULL
  for (i in 1:length(lin))
    if (grepl("Id:", lin[i])==TRUE)
      g<-c(g,i )

  # Lecture mesure par nom
  lname <- trimws(substr(lin[g], 4, 15))
  linc <- trimws(substr(lin[g], 20, 24))
  laz <- trimws(substr(lin[g], 29, 33))
  ltet <- trimws(substr(lin[g], 39, 43))
  lpsy <- trimws(substr(lin[g], 49, 53))
  lv <- trimws(substr(lin[g], 57, 61))
  lTH <- trimws(substr(lin[g], 69, 72))
  lshape <- trimws(substr(lin[g], 76, 78))

  lT1 <- trimws(substr(lin[g+2], 14, 17))
  lT2 <- trimws(substr(lin[g+2], 25, 28))
  lT3 <- trimws(substr(lin[g+2], 36, 39))
  lT4 <- trimws(substr(lin[g+2], 47, 50))

  list.mesure <- NULL
  list.mesure <- data.frame(number = c(1: length(g)), name = lname, inc = as.numeric(linc), az = as.numeric(laz), tet = as.numeric(ltet), psy = as.numeric(lpsy), TH =as.numeric(lTH),
                         shape = lshape, vol =as.numeric(lv),  T1 = as.numeric(lT1), T2 = as.numeric(lT2), T3 = as.numeric(lT3), T4 = as.numeric(lT4), stringsAsFactors = FALSE)


  return(list.mesure)
}
#' Lecture des infos d'orientation sur mesures d'un fichier AM
#' @return une data.frame avec les infos sur les spécimens
#' @export
read.AM.orient <- function (file.AM, encoding="macroman")
{

  # Lecture et Copy du fichier
  lin<- NULL
  fil <- file(file.AM, "r", encoding=encoding)
  lin <- readLines(fil)
  close(fil)
  # Recherche position-ligne des noms
  g <- NULL
  for (i in 1:length(lin))
    if (grepl("Id:", lin[i])==TRUE)
      g<-c(g,i+1 )

  # Lecture mesure par nom
  lname <-  lname <- trimws(substr(lin[g-1], 4, 15))
  lorient <- trimws(substr(lin[g], 8, 10))
  lj <- trimws(substr(lin[g], 14, 16))
  lm <- trimws(substr(lin[g], 20, 22))
  la <- trimws(substr(lin[g], 26, 30))
  lh <- trimws(substr(lin[g], 34, 37))
  lmin <- trimws(substr(lin[g], 41, 43))
  lsec <- trimws(substr(lin[g], 47, 49))
  lSM <- trimws(substr(lin[g], 54, 59))

  list.mesure <- NULL
  list.mesure <- data.frame(number = c(1: length(g)), name = lname, orient = as.character(lorient), jour = as.numeric(lj), mois = as.numeric(lm), annee = as.numeric(la), h = as.numeric(lh), min =as.numeric(lmin),
                            sec = as.numeric(lsec), SM =as.numeric(lSM), stringsAsFactors = FALSE)


  return(list.mesure)
}
#' Fonction de création de fichier pour les magnétomètres, pour des matériaux déplacés ----
#' @param encoding mettre "macroman" pour les mesures au molspin
#' @export
genere.AMD <- function(file.AMD = "fichier.AMD", list.ech, shape = "Cyl" , encoding = "macroman")
{
  entete<- c( "Spinner_Molspin 2008" ,
              "Commune : à définir",
              "Site : à définir",
              "Latitude  :   0°  0'  0\" ",
              "Longitude :   0°  0'  0\" IGRF:+00.0",
              "Prélèvements sur matériaux déplacés",
              "Type de carottage : à plat",
              paste0("Date de création : " , format(Sys.time(), "%d/%m/%Y")),
              "","")

  txt.mesures <- entete
  for (i in 1:length(list.ech)) {
    txt.mesures <- c( txt.mesures,
                      paste("Id:", format(list.ech[i], width = 13), "in:000.0 az:000.0 Tet:000.0 Psy:000.0 v:12.27 com:TH50.0µT ", shape, sep = ""),
                      "Repère:",
                      "CompDes:  T1:0000T+  T2:0000T+  T3:0000T-  T4:0000T-",
                      "")
  }


  # Ecriture du fichier

  filCon <- file(file.AMD, encoding = encoding)
  writeLines(txt.mesures, filCon)
  close(filCon)
}

#' Fonction de création de fichier pour les magnétomètres, pour des structures en place ----
#' @param encoding mettre "macroman" pour les mesures au molspin
#' @export
genere.AMP <- function(file.AMP = "fichier.AMP", list.ech, shape = "Cyl" , encoding = "macroman")
{
  entete<- c( "Spinner_Molspin 2008" ,
              "Commune : à définir",
              "Site : à Definir",
              "Latitude  :   0°  0'  0\" ",
              "Longitude :   0°  0'  0\" IGRF:+00.0",
              "Prélèvements sur matériaux déplacés",
              "Type de carottage : à plat",
              paste0("Date de création : " , format(Sys.time(), "%d/%m/%Y")),
              "","")

  txt.mesures <- entete
  for (i in 1:length(list.ech)) {
    txt.mesures <- c( txt.mesures,
                      paste("Id:", format(list.ech[i], width = 13), "in:000.0 az:000.0 Tet:000.0 Psy:000.0 v:12.27 com:TH50.0µT ", shape, sep = ""),
                      "Orient:NM  J: 1  M: 1  A:2000  H: 00  M: 0  S: 0  SM:000.0",
                      "CompDes:  T1:0000T+  T2:0000T+  T3:0000T-  T4:0000T-",
                      "")
  }


  # Ecriture du fichier

  filCon <- file(file.AMP, encoding = encoding)
  writeLines(txt.mesures, filCon)
  close(filCon)
}

#' Lecture des mesures d'un fichier Pal (*.txt)
#' @param encoding  Pour les fichiers du magnétomètre, il faut "macroman" -> difficle à connaitre, peut être "latin1" ou "utf8".
#' @return une data.frame
#' @export
read.Pal.mesures <- function(file.Pal, encoding = "macroman")
{
  # Lecture et Copy du fichier
  lin<- NULL
  fil <- file(file.Pal, "r", encoding = encoding)
  lin <- readLines(fil)
  close(fil)
  # Recherche position-ligne des noms
  g <- NULL
  for (i in 1:length(lin))
    if (grepl("Id:", lin[i])==TRUE)
      g<-c(g,i )

  # Lecture mesure par nom
  lname <- trimws(substr(lin[g],4, 15))
  list.mesure <- NULL
  lmes <- NULL

  first.mesure <- NULL
  last.mesure <- NULL
  for (i in 1:length(g)) {
    first.mesure <- g[i] + 2
    if (i<length(g))
      last.mesure <- g[i+1] - 1
    else
      last.mesure<- length(lin)

    tEtap <- NULL
    tEtap.val <- NULL
    tX <- NULL
    tY <- NULL
    tZ <- NULL
    tI <- NULL
    tD <- NULL
    tF <- NULL
    tQual <- NULL
    tApp <- NULL
    tSuscep <- NULL
    for (j in first.mesure:last.mesure) {
      lEtap <- substr(lin[j], 3, 7)
      lEtap.val <- substr(lin[j], 3, 5)
      lX <- substr(lin[j], 9, 18)
      lY <- substr(lin[j], 20, 29)
      lZ <- substr(lin[j], 31, 40)
      lQual <- substr(lin[j], 43, 44)
      lApp <- substr(lin[j], 46, 46)
      lSuscep <- substr(lin[j], 48, 52)
      if ( is.na(as.numeric(lX)) || is.na(as.numeric(lY)) || is.na(as.numeric(lY)) ) {
        lX <- 0
        lY <- 0
        lZ <- 0
        lI <- 0
        lD <- 0
        lF <- 0
      }
      else {
        lI <- to.polar.I(as.numeric(lX), as.numeric(lY), as.numeric(lZ))
        lD <- to.polar.D(as.numeric(lX), as.numeric(lY), as.numeric(lZ))
        lF <- to.polar.F(as.numeric(lX), as.numeric(lY), as.numeric(lZ))
      }

      tEtap <- c(tEtap, lEtap)
      tEtap.val <- c(tEtap.val, lEtap.val)
      tX <- c(tX, lX)
      tY <- c(tY, lY)
      tZ <- c(tZ, lZ)
      tI <- c(tI, lI)
      tD <- c(tD, lD)
      tF <- c(tF, lF)
      tQual <- c(tQual, lQual)
      tApp <- c(tApp, lApp)
      tSuscep <- c(tSuscep, lSuscep)
      lmes <- data.frame(number = i, name = lname[i], step = tEtap, step.value = as.numeric(tEtap.val),
                         X = as.numeric(tX), Y = as.numeric(tY), Z = as.numeric(tZ),
                         I = as.numeric(tI), D = as.numeric(tD), F = as.numeric(tF),
                         Quality = as.numeric(tQual), App = tApp, Suscep = as.numeric(tSuscep), stringsAsFactors = FALSE)
    }
    list.mesure <- rbind(list.mesure, lmes)

  }

  return(list.mesure)
}

#' Lecture des infos sur mesures d'un fichier AM
#' @return une data.frame avec les infos sur les spécimens
#' @export
read.Pal.info <- function (file.Pal, encoding="macroman")
{

  # Lecture et Copy du fichier
  lin<- NULL
  fil <- file(file.Pal, "r", encoding=encoding)
  lin <- readLines(fil)
  close(fil)
  # Recherche position-ligne des noms
  g <- NULL
  for (i in 1:length(lin))
    if (grepl("Id:", lin[i])==TRUE)
      g<-c(g,i )

  # Lecture mesure par nom
  lname <- trimws(substr(lin[g], 4, 15))
  linc <- trimws(substr(lin[g], 20, 24))
  laz <- trimws(substr(lin[g], 29, 33))
  ldip <- trimws(substr(lin[g], 39, 43))
  lstr <- trimws(substr(lin[g], 49, 53))
  lv <- trimws(substr(lin[g], 57, 61))
  lcom <- trimws(substr(lin[g], 69, 100))

  lL <- trimws(substr(lin[g+1], 3, 13))
  lG <- trimws(substr(lin[g+1], 16, 26))
  lH <- trimws(substr(lin[g+1], 33, 35))
  lT <- trimws(substr(lin[g+1], 38, 56))

  lazm <- trimws(substr(lin[g+1], 61, 66))
  lazs <- trimws(substr(lin[g+1], 71, 77))
  lOr <- trimws(substr(lin[g+1], 81, 92))
  lFm <- trimws(substr(lin[g+1], 96, 106))
  lLoc <- trimws(substr(lin[g+1], 110, 120))

  list.mesure <- NULL
  list.mesure <- data.frame(number = c(1: length(g)), name = lname, inc = as.numeric(linc), az = as.numeric(laz), dip = as.numeric(ldip), str = as.numeric(lstr),  v = as.numeric(lv), com =as.character(lcom),
                            L = as.numeric(lL), G = as.numeric(lG), H = as.numeric(lH), T = as.character(lT),
                            azm = as.numeric(lazm), azs = as.numeric(lazs), Or = as.character(lOr), Fm = as.character(lFm), Loc = as.character(lLoc),
                            stringsAsFactors = FALSE)


  return(list.mesure)
}

#' Lecture des mesures d'un fichier AM de l'ancien format
#' La valeur de la température d'anisotropie n'était pas informée, il est donc possible de le renseigner avec step.value
#' @param ani.step.value valeur de la température des manips d'anisotropie
#' @return une data.frame avec les infos sur les spécimens
#' @export
read.AM.oldType.mesures <- function(file.AM, encoding = "macroman", ani.step.value = 700)
{
  # Lecture et Copy du fichier
  lin<- NULL
  fil <- file(file.AM, "r", encoding = encoding)
  lin <- readLines(fil)
  close(fil)
  # Recherche position-ligne des noms
  g <- NULL

  # comptage et séparation des données, premier échantillon à ligne 12
  current.name <- "" #trimws(substr(lin[12], 1, 8))
  list.name <- NULL
  for (i in 11:length(lin)) {
    if ( trimws(substr(lin[i], 1, 8)) != current.name & trimws(substr(lin[i], 1, 8)) != "") {
      current.name <- trimws(substr(lin[i], 1, 8))
      list.name <- c(list.name, current.name)
      g <- c(g, i)
    }
  }


  # Lecture mesure par nom
  lname <- trimws(substr(lin[g],1, 8))
  list.mesure <- NULL
  lmes <- NULL

  first.mesure <- NULL
  last.mesure <- NULL
  for (i in 1:length(g)) {
    first.mesure <- g[i]+1
    if (i<length(g))
      last.mesure <- g[i+1]-2
    else
      last.mesure<- length(lin)

    tEtap <- NULL
    tEtap.val <- NULL
    tX <- NULL
    tY <- NULL
    tZ <- NULL
    tI <- NULL
    tD <- NULL
    tF <- NULL

    for (j in first.mesure:last.mesure) {
      lMesure <- trimws(substr(lin[j], 12, 14) ) # correspond au champ Mesure
      lEtap <- trimws(substr(lin[j], 19, 24) ) # correspond au champ Etape
      lEtap.val <- trimws(substr(lin[j], 19, 21) )
      lEtap.code <- trimws(substr(lin[j], 20, 23) )
      ordre <- trimws(substr(lin[j], 26, 30) ) # correspond au champ Ordre
      ordre.sens <- substr(ordre, 2, 2 )
      # Interpretation des codes
      if (lMesure == "ANI")
      {
        if (ordre.sens == "T")
        {
          lEtap.val <- ani.step.value
          code <- substr(lEtap.code, nchar(lEtap.code)-1, nchar(lEtap.code))
          lEtap.code <- code
        }
        else if (ordre.sens == "D" | ordre.sens == "I")
        {
          lEtap.val <- ani.step.value
          code <- substr(lEtap.code, nchar(lEtap.code), nchar(lEtap.code))
          if ( ordre.sens == "D") {
            lEtap.code <- paste(code, "+", sep = "")
          }
          else {
            lEtap.code <- paste(code, "-", sep = "")
          }

        }
        lEtap <- paste(lEtap.val, lEtap.code, sep = "")

      }

      # lecture des mesures
      lI <- trimws(substr(lin[j], 33, 43) )
      lD <- trimws(substr(lin[j], 47, 57) )
      lF <- trimws(substr(lin[j], 61, 71) )
      if ( is.na(as.numeric(lI)) || is.na(as.numeric(lD)) || is.na(as.numeric(lF)) ) {
        lX <- 0
        lY <- 0
        lZ <- 0
        lI <- 0
        lD <- 0
        lF <- 0
      }
      else {
        lX <- to.cartesian.X(as.numeric(lI), as.numeric(lD), as.numeric(lF))
        lY <- to.cartesian.Y(as.numeric(lI), as.numeric(lD), as.numeric(lF))
        lZ <- to.cartesian.Z(as.numeric(lI), as.numeric(lD), as.numeric(lF))
      }

      tEtap <- c(tEtap, lEtap)
      tEtap.val <- c(tEtap.val, lEtap.val)
      tX <- c(tX, lX)
      tY <- c(tY, lY)
      tZ <- c(tZ, lZ)
      tI <- c(tI, lI)
      tD <- c(tD, lD)
      tF <- c(tF, lF)

      lmes <- data.frame(number = i, name = lname[i], step = tEtap, step.value = as.numeric(tEtap.val),
                         X = as.numeric(tX), Y = as.numeric(tY), Z = as.numeric(tZ),
                         I = as.numeric(tI), D = as.numeric(tD), F = as.numeric(tF), stringsAsFactors = FALSE)
    }
    list.mesure <- rbind(list.mesure, lmes)

  }

  return(list.mesure)
}

#' Lecture des mesures d'un fichier AM de l'ancien format
#' @return une data.frame avec les infos sur les spécimens
#' @export
read.AM.oldType.info <- function (file.AM, encoding = "macroman")
{
  # Lecture et Copy du fichier
  lin<- NULL
  fil <- file(file.AM, "r", encoding=encoding)
  lin <- readLines(fil)
  close(fil)

  # Recherche position-ligne des noms
  g <- NULL

  # comptage et séparation des données, premier échantillon à ligne 12
  current.name <- "" #trimws(substr(lin[12], 1, 8))
  list.name <- NULL
  for (i in 11:length(lin)) {
    if ( trimws(substr(lin[i], 1, 8)) != current.name & trimws(substr(lin[i], 1, 8)) != "") {
      current.name <- trimws(substr(lin[i], 1, 8))
      list.name <- c(list.name, current.name)
      g <- c(g, i)
    }
  }


  # Lecture mesure par nom
  lname <- trimws(substr(lin[g], 1, 8))
  lop <- trimws(substr(lin[g], 10, 11))
  lLongueur <- trimws(substr(lin[g], 12, 17))
  lOrientation <- trimws(substr(lin[g], 22, 40))
  lH <- trimws(substr(lin[g], 45, 48))

  lT1 <- trimws(substr(lin[g], 53, 57))
  lT2 <- trimws(substr(lin[g], 61, 65))
  lT3 <- trimws(substr(lin[g], 69, 73))
  lT4 <- trimws(substr(lin[g], 77, 81))
  lPhi <- trimws(substr(lin[g], 86, 91))
  lTet <- trimws(substr(lin[g], 95, 100))

  list.mesure <- NULL
  list.mesure <- data.frame( number = c(1: length(g)), name = lname, Opt = lop, Long = as.numeric(lLongueur), Orient = lOrientation ,TH =as.numeric(lH),
                             T1 = as.numeric(lT1), T2 = as.numeric(lT2), T3 = as.numeric(lT3), T4 = as.numeric(lT4),
                             Phi = as.numeric(lPhi), Tet = as.numeric(lTet), stringsAsFactors = FALSE)



  return(list.mesure)
}


#' Extraction des mesures correspondant à un nom
#' @export
extract.mesures.specimen.name <- function( specimen.name, list.mesure)
{
  selec <- which(list.mesure$name == trimws(specimen.name))

  if (length(selec) == 0)
    warning("Pas de mesure avec ce nom !")

  res.list <- list.mesure[selec,]
  return(res.list)
}

#' Extraction des mesures correspondant à un numéro
#' @export
extract.mesures.specimen.number <- function( specimen.number, list.mesure)
{
  selec <- which(list.mesure$number == specimen.number)
  if (length(selec) == 0)
    warning("Pas de mesure avec ce numéro !")

  res.list <- list.mesure[selec,]
  return(res.list)
}

# Anisotropy Function ----

#' Calcul la matrice d'anisotropie symetrisée et normalisé pour un spécimen
#' qui sert à la correction
#' @param mesures data.frame contenant les mesures
#' @param step.value typiquement la température des mesures d'anisotropie
#' @param step.code code indiquant les étapes. les noms peuvent changer, mais pas l'ordre
#' @param volume la valeur du volume du spécimen
#' @param TH la valeur du champ appliqué
#' @return un data.frame avec les colonnes "L1", "L1.Inc", "L1.Dec", "L2", "L2.Inc", "L2.Dec", "L3", "L3.Inc", "L3.Dec", "F13", "F12", "F23"
#' @export
anisotropy.matrix.symetric <- function(mesures, step.value, step.code = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"), ...)
{

  ani.step <- trimws(paste(as.character(step.value), step.code, sep = "") )

  selec <- NULL
  for (i in 1:length(ani.step)) {
    selec <- c( selec, which(trimws(mesures$step) == trimws(ani.step[i])) )
  }

  res.list <- NULL
  res.list <- mesures[selec,]
  res.list <- res.list

  mat.plus <- matrix( c( res.list$X[3], res.list$X[5], res.list$X[1],
                         res.list$Y[3], res.list$Y[5], res.list$Y[1],
                         res.list$Z[3], res.list$Z[5], res.list$Z[1]) , 3, 3)


  mat.moins <- matrix( c( res.list$X[4], res.list$X[6], res.list$X[2],
                          res.list$Y[4], res.list$Y[6], res.list$Y[2],
                          res.list$Z[4], res.list$Z[6], res.list$Z[2]) , 3, 3)
  mat.reel <- (mat.plus - mat.moins)/2 #/ (volume * 1E-6)# / coef.norm

  # Symetrisation
  coef.norm <- 1 #TH* 10/ (4*pi)
  kxx <- mat.reel[1,1] / coef.norm
  kyy <- mat.reel[2,2] / coef.norm
  kzz <- mat.reel[3,3] / coef.norm
  kxy <- (mat.reel[1,2] + mat.reel[2,1]) / 2 / coef.norm
  kxz <- (mat.reel[1,3] + mat.reel[3,1]) / 2 / coef.norm
  kyz <- (mat.reel[2,3] + mat.reel[3,2]) / 2 / coef.norm

  # Normalisation
  suscept <- (kxx + kyy + kzz)/3
  mat.sym.norm <- matrix( c( kxx / suscept, kxy / suscept, kxz / suscept,
                             kxy / suscept, kyy / suscept, kyz / suscept,
                             kxz / suscept, kyz / suscept, kzz / suscept) , 3, 3)

  return(mat.sym.norm)

}


#' Calcul du vecteur propre et de la matrice d'anisotropie pour un spécimen
#' @param mesures data.frame contenant les mesures
#' @param step.value typiquement la température des mesures d'anisotropie
#' @param step.code code indiquant les étapes. les noms peuvent changer, mais pas l'ordre
#' @param volume la valeur du volume du spécimen
#' @param TH la valeur du champ appliqué
#' @return un data.frame avec les colonnes "L1", "L1.Inc", "L1.Dec", "L2", "L2.Inc", "L2.Dec", "L3", "L3.Inc", "L3.Dec", "F13", "F12", "F23"
#' @export
anisotropy.eigen <- function(mesures, step.value, step.code = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"), volume = 1, TH = 1,...)
{

  #step.code <- c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB")
  ani.step <- trimws(paste(as.character(step.value), step.code, sep = "") )

  selec <- NULL
  for (i in 1:length(ani.step)) {
    selec <- c( selec, which(trimws(mesures$step) == trimws(ani.step[i])) )
  }
  res.list <- NULL
  res.list <- mesures[selec,]
  res.list <- res.list

  mat.plus <- matrix( c( res.list$X[3], res.list$X[5], res.list$X[1],
                         res.list$Y[3], res.list$Y[5], res.list$Y[1],
                         res.list$Z[3], res.list$Z[5], res.list$Z[1]) , 3, 3)


  mat.moins <- matrix( c( res.list$X[4], res.list$X[6], res.list$X[2],
                          res.list$Y[4], res.list$Y[6], res.list$Y[2],
                          res.list$Z[4], res.list$Z[6], res.list$Z[2]) , 3, 3)
  mat.reel <- (mat.plus - mat.moins)/2 / (volume * 1E-6)# / coef.norm

  # Symetrisation
  coef.norm <- TH* 10/ (4*pi)
  kxx <- mat.reel[1,1] / coef.norm
  kyy <- mat.reel[2,2] / coef.norm
  kzz <- mat.reel[3,3] / coef.norm
  kxy <- (mat.reel[1,2] + mat.reel[2,1]) / 2 / coef.norm
  kxz <- (mat.reel[1,3] + mat.reel[3,1]) / 2 / coef.norm
  kyz <- (mat.reel[2,3] + mat.reel[3,2]) / 2 / coef.norm

  # normalisation
  suscept <- (kxx + kyy + kzz)/3
  mat.sym.norm <- matrix( c( kxx / suscept, kxy / suscept, kxz / suscept,
                             0,  kyy / suscept, kyz / suscept,
                             0, 0, kzz / suscept) , 3, 3)

  # vecteurs et valeurs propres
  v <- eigen(mat.sym.norm, symmetric = TRUE)

  # calcul des angles I,D des vecteurs propres}
  v1 <- NULL
  v1$I<-to.polar.I(v$vectors[1, 1], v$vectors[2, 1], v$vectors[3 ,1])
  v1$D<-to.polar.D(v$vectors[1, 1], v$vectors[2, 1], v$vectors[3 ,1])


  v2 <- NULL
  v2$I<-to.polar.I(v$vectors[1, 2], v$vectors[2, 2], v$vectors[3 ,2])
  v2$D<-to.polar.D(v$vectors[1, 2], v$vectors[2, 2], v$vectors[3 ,2])


  v3 <- NULL
  v3$I<-to.polar.I(v$vectors[1, 3], v$vectors[2, 3],v$vectors[3 ,3])
  v3$D<-to.polar.D(v$vectors[1, 3], v$vectors[2, 3],v$vectors[3 ,3])


  if (v3$I<0) {
    v3$I <- -v3$I
    v3$D <- D.AM(v3$D +180)
    v2$I <- -v2$I
    v2$D <- D.AM(v2$D +180)
  }
  F13 <- v$values[1]/v$values[3]
  F12 <- v$values[1]/v$values[2]
  F23 <- v$values[2]/v$values[3]


  ani <- c(eigen = v, L1 = v1, L2 = v2, L3 = v3, F13 = F13, F12 = F12, F23 = F23)
  return(ani)

}

#' Calcul les coordonnées des tenseurs propres matrice d'anisotropie
#' @param mesures data.frame contenant les mesures
#' @param step.value typiquement la température des mesures d'anisotropie
#' @param step.code code indiquant les étapes. les noms peuvent changer, mais pas l'ordre
#' @param volume la valeur du volume du spécimen
#' @param TH la valeur du champ appliqué
#' @return un data.frame avec les colonnes "L1", "L1.Inc", "L1.Dec", "L2", "L2.Inc", "L2.Dec", "L3", "L3.Inc", "L3.Dec", "F13", "F12", "F23"
#' @seealso \code{\link{anisotropy.eigen.tensor}}
#' @export
anisotropy.eigen.matrix <- function(mesures, step.value, step.code = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"), volume = 1, TH = 1,...)
{

  ani.step <- trimws(paste(as.character(step.value), step.code, sep = "") )

  selec <- NULL
  for (i in 1:length(ani.step)) {
    selec <- c( selec, which(trimws(mesures$step) == trimws(ani.step[i])) )
  }
  res.list <- NULL
  res.list <- mesures[selec,]
  res.list <- res.list

  mat.plus <- matrix( c( res.list$X[3], res.list$X[5], res.list$X[1],
                         res.list$Y[3], res.list$Y[5], res.list$Y[1],
                         res.list$Z[3], res.list$Z[5], res.list$Z[1]) , 3, 3)


  mat.moins <- matrix( c( res.list$X[4], res.list$X[6], res.list$X[2],
                          res.list$Y[4], res.list$Y[6], res.list$Y[2],
                          res.list$Z[4], res.list$Z[6], res.list$Z[2]) , 3, 3)
  mat.reel <- (mat.plus - mat.moins)/2 / (volume * 1E-6)# / coef.norm

  # Symetrisation
  coef.norm <- TH* 10/ (4*pi)
  kxx <- mat.reel[1,1] / coef.norm
  kyy <- mat.reel[2,2] / coef.norm
  kzz <- mat.reel[3,3] / coef.norm
  kxy <- (mat.reel[1,2] + mat.reel[2,1]) / 2 / coef.norm
  kxz <- (mat.reel[1,3] + mat.reel[3,1]) / 2 / coef.norm
  kyz <- (mat.reel[2,3] + mat.reel[3,2]) / 2 / coef.norm

  # normalisation
  suscept <- (kxx + kyy + kzz)/3
  mat.sym.norm <- matrix( c( kxx / suscept, kxy / suscept, kxz / suscept,
                             0,  kyy / suscept, kyz / suscept,
                             0, 0, kzz / suscept) , 3, 3)

  # vecteurs et valeurs propres
  v <- eigen(mat.sym.norm, symmetric = TRUE)

  return(as.matrix(v$vectors))

}

#' Calcul des directions et valeurs des vecteurs propres d'anisotropie partielle pour un spécimen à une certaine température
#' @param mesures data.frame contenant les mesures
#' @param step.value typiquement la température des mesures d'anisotropie
#' @param step.code code indiquant les étapes. Les noms peuvent changer, mais pas l'ordre
#' @seealso \code{\link{anisotropy.eigen.matrix}}
#' @return un data.frame avec les colonnes "L1", "L1.Inc", "L1.Dec", "L2", "L2.Inc", "L2.Dec", "L3", "L3.Inc", "L3.Dec", "F13", "F12", "F23"
#' @export
anisotropy.eigen.tensor <- function (mesures, step.value, step.code = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB") )
{
  ani.step <- trimws(paste(as.character(step.value), step.code, sep = "") )

  selec <- NULL
  for (i in 1:length(ani.step)) {
    selec <- c( selec, which(trimws(mesures$step) == trimws(ani.step[i])) )
  }

  res.list <- NULL
  res.list <- mesures[selec,]
  res.list <- res.list

  col.names <- c("L1", "L1.Inc", "L1.Dec", "L2", "L2.Inc", "L2.Dec", "L3", "L3.Inc", "L3.Dec", "F13", "F12", "F23")

  # vérification de l'existence de tous les termes
  verif <- 1
  for (i in 1 : 6)
  verif <- verif * res.list$Z[i]* res.list$Z[i]* res.list$Z[i]

  if ( is.null(verif) || is.na(verif) ) {
    Data <- c( L1 = 0, L2 = 0, L3 = 0,
               L1.Inc = 0, L1.Dec = 0,
               L2.Inc = 0, L2.Dec = 0,
               L3.Inc = 0, L3.Dec = 0,
               F13 = 0, F12 = 0, F23 = 0)

    return(as.data.frame(t(Data), col.names = col.names))

  }



  mat.plus <- matrix( c( res.list$X[3], res.list$X[5], res.list$X[1],
                         res.list$Y[3], res.list$Y[5], res.list$Y[1],
                         res.list$Z[3], res.list$Z[5], res.list$Z[1]) , 3, 3)


  mat.moins <- matrix( c( res.list$X[4], res.list$X[6], res.list$X[2],
                          res.list$Y[4], res.list$Y[6], res.list$Y[2],
                          res.list$Z[4], res.list$Z[6], res.list$Z[2]) , 3, 3)
  mat.reel <- (mat.plus - mat.moins)/2 #/ (volume * 1E-6)# / coef.norm

  # Symetrisation
  coef.norm <- 1 #TH* 10/ (4*pi)
  kxx <- mat.reel[1,1] / coef.norm
  kyy <- mat.reel[2,2] / coef.norm
  kzz <- mat.reel[3,3] / coef.norm
  kxy <- (mat.reel[1,2] + mat.reel[2,1]) / 2 / coef.norm
  kxz <- (mat.reel[1,3] + mat.reel[3,1]) / 2 / coef.norm
  kyz <- (mat.reel[2,3] + mat.reel[3,2]) / 2 / coef.norm

  # normalisation
  suscept <- (kxx + kyy + kzz)/3
  mat.sym.norm <- matrix( c( kxx / suscept, kxy / suscept, kxz / suscept,
                             0,  kyy / suscept, kyz / suscept,
                             0, 0, kzz / suscept) , 3, 3)

  # vecteurs et valeurs propres
  v <- eigen(mat.sym.norm, symmetric = TRUE)

  # calcul des angles I,D des vecteurs propres}
  v1 <- to.polar(v$vectors[1, 1], v$vectors[2, 1], v$vectors[3 ,1])
  v2 <- to.polar(v$vectors[1, 2], v$vectors[2, 2], v$vectors[3 ,2])
  v3 <- to.polar(v$vectors[1, 3], v$vectors[2, 3], v$vectors[3 ,3])


  if (v3$I<0) {
    v3$I <- -v3$I
    v3$D <- D.AM(v3$D +180)
    v2$I <- -v2$I
    v2$D <- D.AM(v2$D +180)
  }

  if (is.null(v$values[3]) || is.null(v$values[2])  ) {
    Data <- c( L1 = 0, L2 = 0, L3 = 0,
               L1.Inc = 0, L1.Dec = 0,
               L2.Inc = 0, L2.Dec = 0,
               L3.Inc = 0, L3.Dec = 0,
               F13 = 0, F12 = 0, F23 = 0)
  }
  else {

    F13 <- v$values[1]/v$values[3]
    F12 <- v$values[1]/v$values[2]
    F23 <- v$values[2]/v$values[3]


    Data <- c( L1 = v$values[1], L2 = v$values[2], L3 = v$values[3],
               L1.Inc = v1$I, L1.Dec = v1$D,
               L2.Inc = v2$I, L2.Dec = v2$D,
               L3.Inc = v3$I, L3.Dec = v3$D,
               F13 = F13, F12 = F12, F23 = F23)
  }


  return(as.data.frame(t(Data), col.names = col.names))


}

#' Calcul des tenseurs d'anisotropie pour une liste de numéro de spécimen
#' @param names liste de numéro de spécimen
#' @param Data.mesures data.frame contenant les mesures
#' @return un data.frame avec les colonnes "L1", "L2", "L3", "L1.Inc", "L1.Dec", "L2.Inc", "L2.Dec", "L3.Inc", "L3.Dec", "F13", "F12", "F23"
#' @export
anisotropy.eigen.tensors.numbers <- function(numbers, Data.mesures, step.value, step.code = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"))
{

  Data <- NULL
  for (i in 1:length(numbers) ) {
    mesures <-  extract.mesures.specimen.number(numbers[i], Data.mesures)
    ani <- anisotropy.eigen.tensor(mesures, step.value = step.value, step.code = step.code)
    Data <- rbind(Data, ani)
  }
  col.names <- c("L1", "L2", "L3", "L1.Inc", "L1.Dec", "L2.Inc", "L2.Dec", "L3.Inc", "L3.Dec", "F13", "F12", "F23")
  return(as.data.frame(Data, col.names = col.names))

}

#' Calcul des tenseurs d'anisotropie pour une liste de nom de spécimen
#' @param names liste de noms de spécimens
#' @param Data.mesures data.frame contenant les mesures
#' @return un data.frame avec les variables "L1", "L2", "L3", "L1.Inc", "L1.Dec", "L2.Inc", "L2.Dec", "L3.Inc", "L3.Dec", "F13", "F12", "F23"
#' @export
anisotropy.eigen.tensors.names <- function(names, Data.mesures, step.value, step.code = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"))
{

  Data <- NULL
  for (i in 1:length(names) ) {
    mesures <-  extract.mesures.specimen.name(names[i], Data.mesures)
    ani <- anisotropy.eigen.tensor(mesures, step.value = step.value, step.code = step.code)
    Data <- rbind(Data, ani)
  }
  col.names <- c("L1", "L2", "L3", "L1.Inc", "L1.Dec", "L2.Inc", "L2.Dec", "L3.Inc", "L3.Dec", "F13", "F12", "F23")
  return(as.data.frame(Data, col.names = col.names))
}

#' Calcul des tenseurs d'anisotropie pour tous les spécimens d'une liste de mesures
#' @return un data.frame avec les colonnes "L1", "L2", "L3", "L1.Inc", "L1.Dec", "L2.Inc", "L2.Dec", "L3.Inc", "L3.Dec", "F13", "F12", "F23"
#' @export
anisotropy.eigen.tensors.all <- function(Data.mesures, Data.number, step.value, step.code = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"))
{
  numbers <- NULL
  for (i in 1: length(Data.number) ) {
    if ( length(which( numbers == Data.number[i] )) == 0 )
      numbers <- c(numbers, Data.number[i])
  }

  Data <- NULL
  for (i in 1:length(numbers) ) {
    mesures <-  extract.mesures.specimen.number(numbers[i], Data.mesures)
    ani <- anisotropy.eigen.tensor(mesures, step.value = step.value, step.code = step.code )
    Data <- rbind(Data, ani)
  }

  col.names <- col.names <- c("L1", "L2", "L3", "L1.Inc", "L1.Dec", "L2.Inc", "L2.Dec", "L3.Inc", "L3.Dec", "F13", "F12", "F23")
  return(as.data.frame(Data, col.names = col.names))


}


#' Calcul de la matrice moyenne symétrisée d'anisotropie moyen
#' @param Data.mesures data.frame contenant les mesures
#' @param step.value typiquement la température des mesures d'anisotropie
#' @param Data.number numéro des échantillons
#' @seealso \code{\link{anisotropy.mean.eigen.tensor}}
#' @return une matrice carrée 3 x3
#' @export
anisotropy.mean.matrix <- function(Data.mesures, Data.number, step.value, step.code = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"))
{
  ani.moyen <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 3, ncol = 3)
  for (spe in 1: length(Data.number))
  {
    spe.mes <- extract.mesures.specimen.number(Data.number[spe], Data.mesures)
    ani.moyen <- ani.moyen + anisotropy.matrix.symetric(mesures = spe.mes, step.value = step.value, step.code = step.code)
  }

  ani.moyen <- ani.moyen/ length(Data.number)

  # Symetrisation
  kxx <- ani.moyen[1,1]
  kyy <- ani.moyen[2,2]
  kzz <- ani.moyen[3,3]
  kxy <- (ani.moyen[1,2] + ani.moyen[2,1]) / 2
  kxz <- (ani.moyen[1,3] + ani.moyen[3,1]) / 2
  kyz <- (ani.moyen[2,3] + ani.moyen[3,2]) / 2

  suscept <- (kxx + kyy + kzz)/3
  mat.sym <- matrix( c( kxx/suscept , kxy/suscept , kxz/suscept,
                        kxy/suscept , kyy/suscept , kyz/suscept,
                        kxz/suscept , kyz/suscept , kzz/suscept ) , 3, 3)

  return(mat.sym)

}

#' Calcul de la direction et des valeurs du vecteur propre d'anisotropie moyen
#' @param Data.mesures data.frame contenant les mesures
#' @param step.value typiquement la température des mesures d'anisotropie
#' @param Data.number numéro des échantillons
#' @seealso \code{\link{anisotropy.eigen.matrix}}, \code{\link{anisotropy.mean.matrix}}
#' @return un data.frame avec les colonnes "L1", "L2", "L3", "L1.Inc", "L1.Dec", "L2.Inc", "L2.Dec", "L3.Inc", "L3.Dec", "F13", "F12", "F23"
#' @export
anisotropy.mean.eigen.tensor <- function(Data.mesures, Data.number, step.value, step.code = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"))
{
  mat.sym <- anisotropy.mean.matrix(Data.mesures, Data.number, step.value = step.value, step.code = step.code)
  v <- eigen(mat.sym, symmetric = TRUE)

  # calcul des angles I,D des vecteurs propres}
  v1 <- to.polar(v$vectors[1, 1], v$vectors[2, 1], v$vectors[3 ,1])
  v2 <- to.polar(v$vectors[1, 2], v$vectors[2, 2], v$vectors[3 ,2])
  v3 <- to.polar(v$vectors[1, 3], v$vectors[2, 3], v$vectors[3 ,3])

  if (v3$I<0) {
    v3$I <- -v3$I
    v3$D <- D.AM(v3$D +180)
    v2$I <- -v2$I
    v2$D <- D.AM(v2$D +180)
  }
  F13 <- v$values[1]/v$values[3]
  F12 <- v$values[1]/v$values[2]
  F23 <- v$values[2]/v$values[3]

  Data <- c( L1 = v$values[1], L2 = v$values[2], L3 = v$values[3],
             L1.Inc = v1$I, L1.Dec = v1$D,
             L2.Inc = v2$I, L2.Dec = v2$D,
             L3.Inc = v3$I, L3.Dec = v3$D,
             F13 = F13, F12 = F12, F23 = F23)

  col.names <- c( "L1", "L2", "L3", "L1.Inc", "L1.Dec", "L2.Inc", "L2.Dec", "L3.Inc", "L3.Dec", "F13", "F12", "F23")
  return(as.data.frame(t(Data), col.names = col.names))

}

#' Corrige les mesures avec la matrice d'anisotropie donnée
#' Il faut une matrice 3 x3
#' @param mesures.frame data.frame avec les mesures et variables X, Y, Z, i, D, F
#' @param ani.matric matrice 3 x 3 correspondant à la matrice (symétrique et normalisé) à utiliser pour la correction
#' @return un data.frame du même type que mesures.frame
#' @export
correction.anisotropy <- function(mesures.frame, ani.matrix)
{

  if (!is.data.frame(mesures.frame))
    warning("mesures must be a data.frame")
  if (!is.matrix(ani.matrix))
    warning("ani.matrix must be a matrix")

  nbMesures <- length(mesures.frame[,1])

  # inversion de la matrice
  ani.inv <- solve(ani.matrix)

  r_cal_echt <- NULL
  # Calcul pour une mesure
  correct <- function(mesure, aniso.inv) {
    # Calcul de la correction
    re.cal.ani <- aniso.inv %*% c( mesure$X, mesure$Y, mesure$Z)
    # recopie des variables
    re.cal <- mesure
    re.cal$X <- re.cal.ani[1]
    re.cal$Y <- re.cal.ani[2]
    re.cal$Z <- re.cal.ani[3]
    re.cal.ID <- to.polar(re.cal$X, re.cal$Y, re.cal$Z)
    re.cal$I <- re.cal.ID$I
    re.cal$D <- re.cal.ID$D
    re.cal$F <- re.cal.ID$F

    return(re.cal)
  }

  resT <- NULL
  for (i in 1 : nbMesures)
    resT <- rbind(resT, correct(mesures.frame[i,], ani.inv))

  return(resT)
}

#' Tracer dans le diagramme Lambert/ Stéréo d'un tenseur
#' @param Data un data.frame avec les variables "L1.Inc", "L1.Dec", "L2.Inc", "L2.Dec", "L3.Inc", "L3.Dec"
#' @export
lambert.ID.tensors <- function(Data, pt.col = "blue3", new = TRUE, ...)
{
  # Restructuration des données
  L1.Inc <- as.numeric(Data$L1.Inc)
  L2.Inc <- as.numeric(Data$L2.Inc)
  L3.Inc <- as.numeric(Data$L3.Inc)

  L1.Dec <- as.numeric(Data$L1.Dec)
  L2.Dec <- as.numeric(Data$L2.Dec)
  L3.Dec <- as.numeric(Data$L3.Dec)

  if (length(L1.Inc) == 0)
    return(print("no DATA"))

  if (length(pt.col) < length(L1.Inc))
    pt.col <- rep(pt.col, length.out = length(L1.Inc))

  for (i in 1:length(L1.Inc) ) {
    Da.I <- c(L1.Inc[i], L2.Inc[i], L3.Inc[i])
    Da.D <- c(L1.Dec[i], L2.Dec[i], L3.Dec[i])
    lambert.ID(Da.I, Da.D,
                    inc.lim = c(0, 90), pch = c(22, 24, 21), pt.col = pt.col[i], new = new  )
    new <- FALSE
  }
}

#' Tracer d'un diagramme de Flinn
#' diagramme des paramètres d'anisotropie F12 en fonction de F23
#' @param Data soit correspond à un data.frame avec les variables $F23 et $F12, soit seulement des valeurs de F23 dans ce cas, renseignez Data.F12
#' @param Data.F12 des valeurs pour F12
#' @param absolue prend la valeur absolue des données
#' @export
flinn <- function( Data, Data.F12 = NULL, pt.names = NULL, pt.col = "blue3", pch = 21, type = "p",
                       xlab = "F23", ylab = "F12", main = "Flinn diagram", absolue = TRUE, new = TRUE)
{
   par(pty="s", "xaxp")

  if(is.data.frame(Data)) {
    if (absolue == TRUE) {
      X <- abs(Data$F23)
      Y <- abs(Data$F12)
    } else {
      X <- Data$F23
      Y <- Data$F12
    }

  } else {
    if (absolue == TRUE) {
      X <- abs(Data)
      Y <- abs(Data.F12)
    } else {
      X <- Data
      Y <- Data.F12
    }

  }
  X.lim <- range(X)
  if (X.lim[1]>1)
    X.lim[1] <- 1

  Y.lim <- range(Y)
  if (Y.lim[1]>1)
    Y.lim[1] <- 1

  XY.max <- max(X.lim[2], Y.lim[2]) * 1.05
  X.lim[2] <- XY.max
  Y.lim[2] <- XY.max

  if (new == TRUE) {
    plot(x = X, y = Y, xlab = xlab, ylab = ylab, xlim = X.lim, ylim = Y.lim, type = type, col = "gray50", bg = pt.col, pch = pch,
         xaxt = "n", yaxt = "n", asp = 1, bty ="n", main = main, new = TRUE )
    ax1 <- axis(1, pos = X.lim[1], col = "gray10")
    ax2 <- axis(2, pos = Y.lim[1], col = "gray10") # Ordonnées
    cc<-array(c(0,1), c(1,2))

    abline(  coef = cc, col = "gray90")

    main <- ""
  }  else {
    points(x = X, y = Y, xlim = X.lim, ylim = Y.lim,  type = type, col = "gray50", bg = pt.col, pch = pch,
         xaxt = "n", yaxt = "n", asp = 1, bty ="n", new = FALSE)
  }

  text(jitter(X, 5, amount = 0), jitter(Y, 5, amount = 0), pt.names)

}

#' Tracer d'un diagramme de désaimantation
#' @param normalize permet de comparer l'évolution de l'aimanation quelque soit l'amplitude en visualisant le résultat comme un pourcentage du maximum
#' @export
demag <- function( Data, F = NULL,  pt.col = "blue3", pch = 21, type = "b",
                             xlab = "°C", ylab = "", main = NULL,
                             names = NA, normalize = TRUE, step.J0 = NULL, new = TRUE, ...)
{

  tmp.frame <- NULL
  if (is.null(main))
    main <- " F vs step.value"

  if(is.data.frame(Data)) {
    tmp.frame$name <- Data$name[!is.na(Data$step.value)]
    tmp.frame$step.value <- Data$step.value[!is.na(Data$step.value)]
    tmp.frame$F <- Data$F[!is.na(Data$step.value)]


  } else {
    # Création du frame interne
    if (is.na(names))
      tmp.frame$name <- as.character(c(1: length(Data)))
    else
      tmp.frame$name <- names

    tmp.frame$step.value <- Data
    tmp.frame$F <- F
  }

  if (length(pt.col) < length(tmp.frame$step.value) )
    pt.col <- rep(pt.col, length(tmp.frame$step.value) )

  # comptage et séparation des données
  current <- tmp.frame$name[1]
  list.name <- current

  for (i in 1: length(tmp.frame$name)) {
    if ( tmp.frame$name[i] != current) {
      current <- tmp.frame$name[i]
      list.name<- c(list.name, current)
    }
  }

  if (is.null(step.J0) == FALSE) {
    if (length(step.J0) < length(tmp.frame$step.value) )
      step.J0 <- rep(step.J0, length(list.name) )
  } else {
    step.J0 <- rep(NA, length(list.name) )
  }

  xlim <- range(tmp.frame$step.value)

  # recherche du Ymax pour définir la taille du graphique
  Ymax <- 0
  J0 <- 0
  list.mesure.i <-  NULL
  if (normalize == TRUE) {
    for (i in 1 : length(list.name)) {
      list.mesure.i$F <- tmp.frame$F[tmp.frame$name == list.name[i]]
      if (is.na(step.J0[i]) == TRUE) {
        J0[i] <- list.mesure.i$F[1]
      } else {
        tmp.etp <- list.mesure.i$F[tmp.frame$step.value == step.J0[i]] # si il y a plusieurs valeurs pour la même étape, comme ani
        J0[i] <- tmp.etp[1]
       # J0[i] <- list.mesure.i$F[tmp.frame$step.value == step.J0[i]]
      }
      Ymax <- max(Ymax, tmp.frame$F[tmp.frame$name == list.name[i]]/J0[i] *100)
    }
  } else
    Ymax <- max(tmp.frame$F)

  list.mesure.i <-  NULL
  for (i in 1 : length(list.name)) {
    list.mesure.i$step.value <- tmp.frame$step.value[which(tmp.frame$name == list.name[i])]
    list.mesure.i$F <- tmp.frame$F[tmp.frame$name == list.name[i]]
    if (normalize == TRUE) {
      coefY <- 100 / J0[i]

    } else
      coefY <- 1

    Xi <- list.mesure.i$step.value
    Yi <- list.mesure.i$F  * coefY

    ylim <- c(0, Ymax)

    if (i == 1)
      plot(x = Xi, y = Yi, xlab = xlab, ylab = ylab, ylim = ylim, xlim = xlim, type = type,
           col = adjustcolor( pt.col[i], alpha.f = 0.7), bg = pt.col[i], pch = pch,
           yaxt = "n", bty = "n", main = main, new = new)
    else
      lines(x = Xi, y = Yi, type = type, col = adjustcolor( pt.col[i], alpha.f = 0.7), bg = pt.col[i], pch = pch, yaxt = "n", bty ="n")

  }

  axis(2, pos = 0, col = "darkgray") #cex.axis = 0.8) # Ordonnées

  if (normalize == TRUE) {
    mtext( paste(format(Ymax, digits = 3), "%"), side = 3, col = "gray5", adj = 0, cex = par("cex.lab"))
  } else {
    mtext( format(Ymax, digits = 3, scientific = TRUE), side = 3, col = "gray5", adj = 0, cex = par("cex.lab"))
  }


}


#' Correction de carottage sur tranche
#' Tourne les mesures de manière à retrouver le résultat suivant la convention de carottage à plat
#' @export
correction.carottage.tranche <- function(mesures.brutes)
{
  mesures.convent <- mesures.brutes
  mesures.convent$X <- mesures.brutes$X
  mesures.convent$Y <-  - mesures.brutes$Z
  mesures.convent$Z <- mesures.brutes$Y
  for (i in 1 : length(mesures.convent$X)) {
    mesures.convent$I[i] <- to.polar.I(mesures.convent$X[i], mesures.convent$Y[i], mesures.convent$Z[i])
    mesures.convent$D[i] <- to.polar.D(mesures.convent$X[i], mesures.convent$Y[i], mesures.convent$Z[i])
  }
  return(mesures.convent)
}

#' Correction de carottage en bout
#' Tourne les mesures de manière à retrouver le résultat suivant la convention de carottage à plat
#' @export
correction.carottage.enbout <- function(mesures.brutes)
{
  mesures.convent <- mesures.brutes

  for (i in 1 : length(mesures.brutes$X)) {
    mesures.convent$X[i] <- mesures.brutes$Z[i]
    mesures.convent$Y[i] <- mesures.brutes$X[i]
    mesures.convent$Z[i] <- mesures.brutes$Y[i]
    mesures.convent$I[i] <- to.polar.I(mesures.convent$X[i], mesures.convent$Y[i], mesures.convent$Z[i])
    mesures.convent$D[i] <- to.polar.D(mesures.convent$X[i], mesures.convent$Y[i], mesures.convent$Z[i])
  }
  return(mesures.convent)
}


# Partial component  ---

#' Calcul les directions du vecteur partiel pour une série d'étape
#' @param en0 permet de calculer la composante qui passe par l'origine (0, 0)
#' @return une data.frame "X", "Y", "Z", "I", "D", "F", "Sl", "MAD"
#' @export
partial.vector <- function(TabX, TabY, TabZ, en0 = TRUE)
{
  col.names <- c("X", "Y", "Z", "I", "D", "F", "Sl", "MAD")
  ntab <- length(TabX);

  if ( (ntab<1) || (ntab!=length(TabY)) || (ntab!=length(TabZ))) {
    warning("No mesure, to calculate the partial vector")
    Data <- c( X = 0, Y = 0, Z = 0,
               I = 0, D = 0, F = 0,
               Sl = NA, MAD = NA )

    return(as.data.frame(t(Data), col.names = col.names))  #il faut au moins 2 étapes
  }

  if ( ntab == 1) {
    message("Only one mesure, to calculate the partial vector")
    vp1 <- to.polar(TabX[1], TabY[1], TabZ[1])
    Data <- c( X = TabX[1], Y = TabY[1], Z = TabZ[1],
               I = vp1$I, D = vp1$D, F = vp1$F,
               Sl = 1, MAD = 0 )

    return(as.data.frame(t(Data), col.names = col.names))  #il faut au moins 2 étapes
  }


  # calcul du barycentre du nuage de points si 'en0' false
  somx <- 0; somy <- 0; somz <- 0

  if (en0 == FALSE) {
    somx <- sum(TabX)
    somy <- sum(TabY)
    somz <- sum(TabZ)
  }

  xm <- somx/ ntab
  ym <- somy/ ntab
  zm <- somz/ ntab

  Mom <-0; # variable pour Kirschink
  # calcul du moment suivant kischvink
  #  c.à.d. longueur du chemin total entre les points
  for ( j in 2 : ntab) {
    Mom <- Mom + sqrt( (TabX[j]-TabX[j-1])^2
                       + (TabY[j]-TabY[j-1])^2
                       + (TabZ[j]-TabZ[j-1])^2 )
  }

  # Longueur du plus court chemin
  Sl <- sqrt( (TabX[ntab]-TabX[1])^2
              + (TabY[ntab]-TabY[1])^2
              + (TabZ[ntab]-TabZ[1])^2);
  # Rapport des distances
  Sl <- Sl/Mom;


  sxx <- sum((TabX-xm)^2); sxy <- sum((TabX-xm)*(TabY-ym)); sxz <- sum((TabX-xm)*(TabZ-zm))
  syy <- sum((TabY-ym)^2);         syz <- sum((TabY-ym)*(TabZ-zm))
  szz <- sum((TabZ-zm)^2);
  # Matrice d inertie des points x, y, z pour le calcul de la droite
  Ixx = syy+szz;  Ixy = -sxy;     Ixz = -sxz;
  Iyy = sxx+szz;  Iyz = -syz;
  Izz = sxx+syy;


  #  MAD

  mat.sym.norm <- matrix( c( sxx, sxy, sxz,
                             0,  syy, syz,
                             0, 0, szz) , 3, 3)

  # vecteurs et valeurs propres pour MAD
  v.mad <- eigen(mat.sym.norm, symmetric = TRUE)


  if (( v.mad$values[1] == 0) || ((v.mad$values[2] + v.mad$values[3])/ v.mad$values[1] < 0)) {
    warning('MAD uncalculable')
    MAD <- NA
  } else {
    MAD <- atan(sqrt( ((v.mad$values[2] + v.mad$values[3])/ v.mad$values[1]) ));
  }

  #   Calcul composante partielle
  # Matrice d'inertie des points x,y,z

  # La direction correspond au vecteur propre minimum
  mat.sym.norm <- matrix( c( Ixx , Ixy , Ixz ,
                             0,  Iyy , Iyz ,
                             0, 0, Izz ) , 3, 3)

  # vecteurs et valeurs propres pour MAD
  v.cp <- eigen(mat.sym.norm, symmetric = TRUE)
  v3 <- to.polar(v.cp$vectors[1, 3],  v.cp$vectors[2, 3],  v.cp$vectors[3 ,3])

  vcorrect<- NULL
  vcorrect$I <- sign(TabZ[1] - TabZ[ntab]) * abs(v3$I)

  vcorrect$D <- v3$D

  if (sign(vcorrect$I) !=  sign(v3$I)) {
    vcorrect$D <- D.AM(v3$D +180)
  }


  # mise en forme du résultat
  MAD <- MAD * 180 /pi

  Data <- c( X = v.cp$vectors[1, 3], Y = v.cp$vectors[2, 3], Z = v.cp$vectors[3, 3],
             I = vcorrect$I, D = vcorrect$D, F = v3$F,
             Sl = Sl, MAD = MAD )

  return(as.data.frame(t(Data), col.names = col.names))
}

#' Calcul le vecteur de la composante partielle
#' et retourne aussi le MAD
#' @param en0 permet de calculer la composante qui passe par l'origine (0, 0)
#' @return une data.frame "X", "Y", "Z", "I", "D", "F", "Sl", "MAD", "DANG"
#' @seealso \code{\link{partial.component.T1T2}}
#' @export
partial.component <- function(TabX, TabY, TabZ, en0 = FALSE)
{
  vp <- partial.vector(TabX, TabY, TabZ, en0 = en0)
  ntab <- length(TabX)

  # Calcul de l'angle entre le vecteur passant par zéro et le vecteur ne passant pas par zéro
  # d'aprés Lisa Tauxe
  v.true <- partial.vector(TabX, TabY, TabZ, en0 = TRUE)
  v.false <- partial.vector(TabX, TabY, TabZ, en0 = FALSE)

  # produit scalaire x*xp+y*yp+z*zp=cos(teta)*(x2+y2+z2)*(xp2+yp2+zp2);
  DANG <- v.true$X * v.false$X + v.true$Y * v.false$Y + v.true$Z * v.false$Z
  DANG <- DANG / (sqrt( v.true$X^2 + v.true$Y^2 + v.true$Z^2)* sqrt( v.false$X^2+ v.false$Y^2+ v.false$Z^2) )
  DANG <- n_arccos(DANG) *180/pi
  return( cbind.data.frame(vp, DANG))
}

#' Calcul le vecteur de la composante partielle en calculant la composante entre les étapes T1 et T2
#' la fonction ne tient pas compte des étapes d'anisotropie, si step.code est bien rensiegné
#' @param en0 permet de calculer la composante qui passe par l'origine (0, 0)
#' @param corr.ani corrige de l'anisotropie
#' @return une data.frame "X", "Y", "Z", "I", "D", "F", "Sl", "MAD", "DANG"
#' @seealso \code{\link{partial.component, zijderveld1.T1T2}}
#' @export
partial.component.T1T2 <- function(Data, T1 = NULL, T2 = NULL, corr.ani = FALSE, ani.step.value = NULL,
                                      step.code = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"),
                                      en0 = FALSE )
{

  if (!is.data.frame(Data))
    warning("Data n'est pas une data.frame")

  res.list <- Data
  if (corr.ani == TRUE) {

    if (is.null(ani.step.value)) { # recherche des étapes d'anisotropie
      for (i in 1:length(Data$step)) {
        if (substring(Data$step[i], 4) == step.code[3])
          ani.step.value <- Data$step.value[i]
      }
    }
    ani <- anisotropy.matrix.symetric(Data, step.value = ani.step.value, step.code = step.code)
  }
  if (!is.null(ani.step.value))
    Data <- remove.step(Data, step.value = ani.step.value, step.code = step.code)

  # recherche étape en dessous de T1

  if (is.null(T1) | T1 == 0 | is.na(T1))
    T1 <- Data$step.value[1]


  if (is.null(T2) | is.na(T2))
    T2 <- Data$step.value[length(Data$step.value)]

  if (T2 < T1) {
    warning(" T2 < T1 ")
  }


  iT1 <- 1
  while(iT1 <= length(Data$step.value) &  Data$step.value[iT1] < T1 ) {
    iT1 <- iT1 + 1
  }

  # recherche étape au dessus de T2
  iT2 <- length(Data$step.value)
  while(iT2 > 0 & Data$step.value[iT2] > T2) {
    iT2 <- iT2 -1
  }

  # Selection de la composante partielle

  TabX <- Data$X[iT1:iT2]
  TabY <- Data$Y[iT1:iT2]
  TabZ <- Data$Z[iT1:iT2]

  # Correction des directions par l'anisotropie
  if (corr.ani == TRUE) {
    # inversion de la matrice
    aniso.inv <- solve(ani)
    for(i in 1: length(TabX)) {
      vp.cor <- aniso.inv %*% c( TabX[i], TabY[i], TabZ[i])
      TabX[i]  <- vp.cor[1]
      TabY[i] <- vp.cor[2]
      TabZ[i] <- vp.cor[3]
    }

  }

  ntab <- length(TabX)

  # Calcul de l'angle entre le vecteur passant par zéro et le vecteur ne passant pas par zéro
  # aprés correction d'anisotropie
  # d'aprés Lisa Tauxe
  v.true <- partial.vector(TabX, TabY, TabZ, en0 = TRUE)
  v.false <- partial.vector(TabX, TabY, TabZ, en0 = FALSE)


  # produit scalaire x*xp+y*yp+z*zp=cos(teta)*(x2+y2+z2)*(xp2+yp2+zp2);
  DANG <- (v.true$X * v.false$X) + (v.true$Y * v.false$Y) + (v.true$Z * v.false$Z)
  DANG <- DANG / (sqrt( v.true$X^2 + v.true$Y^2 + v.true$Z^2)* sqrt( v.false$X^2+ v.false$Y^2+ v.false$Z^2) )
  DANG <- n_arccos(DANG) *180/pi

  vp <- NULL
  if (en0 == TRUE) {
    vp <- v.true
  } else {
    vp <- v.false
  }

 return( cbind.data.frame(vp, DANG))
}

#' Trace un diagramme de zijderveld1 en calculant la composante entre les étapes T1 et T2
#' par défaut la fonction n'affiche pas les étapes d'anisotropie et ne corrige pas les directions de l'anisotropie
#' @param T1 correspond à step.value ou température la plus basse
#' @param T2 correspond à step.value ou température la plus haute
#' @param en0 booléen permettant de forcer la composante partielle à passer par l'origine
#' @param show.step booléen permettant l'affichage des étapes
#' @param withAni permet de voir les étapes d'anisotropie
#' @param ani.step.value correspond à l'step.value ou la température de la détermination de l'anisotropie
#' @param main titre de la figure
#' @param step.code chaîne de caractère représentant les étapes de l'anisotropie "Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB". Cette ordre est obligatoire
#' @seealso \code{\link{partial.component.T1T2, zijderveld2.T1T2}}
#' @export
zijderveld1.T1T2 <- function(Data, T1 = NULL, T2 = NULL, show.step = FALSE, ignore.ani = TRUE, ani.step.value = NULL, main = "",
                             step.code = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"),
                             legend.pos = NULL, legend.txt = c("(Y, X)", "(Y, Z)"),
                             en0 = FALSE )
{
  if (is.null(T1))
    T1 <- Data$step.value[1]
  if (is.null(T2))
    T2 <- Data$step.value[length(Data$step.value)]


  res.list <- Data

  if (ignore.ani == TRUE) { # suppression des étape d anisotropie
    Data <- remove.step(Data, step.value = ani.step.value, step.code = step.code)
  }
  # recherche étape en dessous de T1
  iT1 <- 1
  while(Data$step.value[iT1] < T1) {
    iT1 <- iT1 + 1
  }

  # recherche étape au dessus de T2
  iT2 <- length(Data$step.value)
  while(Data$step.value[iT2] > T2) {
    iT2 <- iT2 -1
  }

  # Calcul de la composante partielle

  TabX <- Data$X[iT1:iT2]
  TabY <- Data$Y[iT1:iT2]
  TabZ <- Data$Z[iT1:iT2]

  vp <- partial.component(TabX, TabY, TabZ, en0 = en0)
  ntab <- length(TabX)

  somx <- 0; somy <- 0; somz <- 0

  if (en0 == FALSE) {
    somx <- sum(TabX)
    somy <- sum(TabY)
    somz <- sum(TabZ)
  }

  xm <- somx/ ntab
  ym <- somy/ ntab
  zm <- somz/ ntab

  # a, b : Valeurs indiquant le point d’interception sur l’axe des y et la pente de la droite -> y = a + bx
  # droite sur les déviations axe (Y, X)
  bd <- vp$X/vp$Y
  ad <- xm - bd*ym
  # droite sur les inclinaisons axe (Y, Z), +Z est négatif
  bi <- -vp$Z/vp$Y
  ai <- -zm - bi*ym

  Y.r <- range(Data$X)
  Z.r <- range(Data$Z)
  ylim <- c(min(Y.r[1], -Z.r[2]), max(Y.r[2], -Z.r[1]))

  if (ylim[1]>0)
    ylim[1]<-0

  if (ylim[2]<0)
    ylim[2]<-0

  if (show.step == TRUE) {
    pt.names <- Data$step
  } else {
    pt.names <- rep("", length(Data$X) )
  }

  zijderveld1(Data$X[1:iT1], Data$Y[1:iT1], Data$Z[1:iT1],  ylim = ylim, pt.names = pt.names[1:iT1], legend.pos = legend.pos, legend.txt = legend.txt, main = main )
  if (iT2 != length(Data$X))
    zijderveld1(Data$X[iT2:length(Data$X)], Data$Y[iT2:length(Data$X)], Data$Z[iT2:length(Data$X)], pt.names = pt.names[iT2:length(Data$X)], new = FALSE)

  zijderveld1(Data$X[iT1:iT2], Data$Y[iT1:iT2], Data$Z[iT1:iT2], pt.col = c("red", "red") , pt.names = pt.names[iT1:iT2], new = FALSE)

  abline(ad, bd)
  abline(ai, bi)

}

#' Trace un diagramme de zijderveld en calculant la composante entre les step.value T1 et T2
#' @param T1 correspond à step.value ou température la plus basse
#' @param T2 correspond à step.value ou température la plus haute
#' @param en0 booléen permettant de forcer la composante partielle de passer par l'origine
#' @param show.step booléen permettant l'affichage des étapes
#' @param TRUE permet de voir les étapes d'anisotropie
#' @param ani.step.value correspond à l'step.value ou la température de la détermination de l'anisotropie
#' @param step.code chaîne de caractère représentant les étapes de l'anisotropie "Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB". Cette ordre est obligatoire
#' @seealso \code{\link{partial.component.T1T2, zijderveld2.T1T2}}
#' @export
zijderveld2.T1T2 <- function(Data, T1 = NULL, T2 = NULL, show.step = FALSE, ignore.ani = TRUE, ani.step.value = NULL, main = "",
                             step.code = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"),
                             legend.pos = NULL, legend.txt = c("(Y, X)", "(Y, Z)"),
                             en0 = FALSE )
{
  if (is.null(T1))
    T1 <- Data$step.value[1]
  if (is.null(T2))
    T2 <- Data$step.value[length(Data$step.value)]


  res.list <- Data

  if (ignore.ani == TRUE) { # suppression des étape d anisotropie
    Data <- remove.step(Data, step.value = ani.step.value, step.code = step.code)
  }
  # recherche étape en dessous de T1
  iT1 <- 1
  while(Data$step.value[iT1] < T1) {
    iT1 <- iT1 + 1
  }

  # recherche étape au dessus de T2
  iT2 <- length(Data$step.value)
  while(Data$step.value[iT2] > T2) {
    iT2 <- iT2 -1
  }

  # Calcul de la composante partielle

  TabX <- Data$X[iT1:iT2]
  TabY <- Data$Y[iT1:iT2]
  TabZ <- Data$Z[iT1:iT2]

  vp <- partial.component(TabX, TabY, TabZ, en0 = en0)
  ntab <- length(TabX)

  somx <- 0; somy <- 0; somz <- 0

  if (en0 == FALSE) {
    somx <- sum(TabX)
    somy <- sum(TabY)
    somz <- sum(TabZ)
  }

  xm <- somx/ ntab
  ym <- somy/ ntab
  zm <- somz/ ntab

  # a, b : Valeurs indiquant le point d’interception sur l’axe des y et la pente de la droite -> y = a + bx
  # droite sur les déviations axe (Y, X)
  bd <- vp$Y/(-vp$X)
  ad <- ym + bd*xm
  # droite sur les inclinaisons axe (Y, Z), +Z est négatif
  bi <- -vp$Z/(-vp$X)
  ai <- -zm + bi*xm

  if (show.step == TRUE) {
    pt.names <- Data$step
  } else {
    pt.names <- rep("", length(Data$X) )
  }

  zijderveld2(Data$X[1:iT1], Data$Y[1:iT1], Data$Z[1:iT1], pt.names = pt.names[1:iT1], legend.pos = legend.pos, legend.txt = legend.txt, main = main )
  if (iT2 != length(Data$X))
    zijderveld2(Data$X[iT2:length(Data$X)], Data$Y[iT2:length(Data$X)], Data$Z[iT2:length(Data$X)], pt.names = pt.names[iT2:length(Data$X)], new = FALSE)

  zijderveld2(Data$X[iT1:iT2], Data$Y[iT1:iT2], Data$Z[iT1:iT2], pt.col = c("red", "red") , pt.names = pt.names[iT1:iT2], new = FALSE)

  abline(ad, bd)
  abline(ai, bi)
}

#' Supprime les étapes à la valeur step.value et avec les codes step.code
#' Utiliser par défaut pour enlever les étapes d'anisotropie
#' @param Data un data.frame possédant la variable $step de type chr
#' @param step.value valeur de la température . Par defaut NUll, alors la tempértature est retrouvée automatiquement avec les step.code définis. Si step.value = '', alors suppression de toutes les étapes avec le step.code
#' @param step.code chaîne de caractère représentant par exemple les étapes de l'anisotropie "Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB", ou des erreurs "??"
#' @param verbose affiche des commentaires et avertissements
#' @export
remove.step <- function(Data, step.value = NULL, step.code = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"), verbose =TRUE )
{
  selec <- NULL
  if (is.null(step.value)) {
    for (i in 1:length(Data$step)) {
      if (substring(Data$step[i], 4) == step.code[3])
        step.value <- Data$step.value[i]
    }
  } else  if (step.value == "") {
    for (i in 1:length(mes$step))
    {
      for (j in 1: length(step.code)) {
        selec <- c( selec, which(substring(Data$step, 4) == step.code[j]) )
      }
    }
  }

  remove.step <- trimws(paste(step.value, step.code, sep = ""))

  for (i in 1:length(remove.step)) {
    selec <- c( selec, which(trimws(Data$step) == trimws(remove.step[i])) )
  }

  if (verbose == TRUE) {
    if (length(selec)== 0) {
      warning("No step to remove")
      return(Data)
    } else {
      print(paste(Data[selec,]$step))
    }
  }

  if (length(selec) > 0) {
    res.list <- Data[-selec,]
  }
  else {
    res.list <- Data
  }

  return(res.list)
}

# Calcul rotation ----

#' Fait tourner les mesures selon deviation
#' @param Data un data.frame possédant la variable $step de type chr
#' @param deviation valeur de l'angle de rotation (en degrés).
#' @export
rotation.mesure <- function(Data, deviation)
{
  res <- Data

  for (i in 1 : length(res$D)) {
    res$D[i] <- as.numeric(res$D[i] + deviation)
    res$X[i] <- as.numeric(to.cartesian.X(res$I[i], res$D[i], res$F[i]) )
    res$Y[i] <- as.numeric(to.cartesian.Y(res$I[i], res$D[i], res$F[i]) )
    res$Z[i] <- as.numeric(to.cartesian.Z(res$I[i], res$D[i], res$F[i]) )
  }

  return(res)
}

# Calcul astronomique ----

# Le jour julien 0 commence le 24 novembre -4713 (4712 BC) à 12h
#' The number of Julian days for astronomical calculations
#'  Julian Day 0 starts November 24th -4713 (4712 BC) at 12:00 pm
#' @seealso \code{\link{https://codes-sources.commentcamarche.net/source/31774-calcul-de-la-position-du-soleil-declinaison-angle-horaire-altitude-et-azimut-altaz-solaire}}
#' @export
julian.day <- function( day, month, year, hour, minute, seconde)
{
  day.hour <- day + hour/24.0 + minute/1440.0 + seconde/86400.0

  if (month == 1 || month == 2) {
    year <- year-1.0
    month <- month+12.0
  }

  a <- trunc(year/100.0)
  b <- 2 - a + trunc(a/4.0)

  julian <- trunc(365.25*(year+4716.0)) + trunc(30.6001*(month+1.0)) + day.hour + b - 1524.5

  return (as.numeric(julian))
}


#' Calculate the azimuth of the sun in a place at a given date at a given time "UTC".
#' @seealso \code{\link{https://codes-sources.commentcamarche.net/source/31774-calcul-de-la-position-du-soleil-declinaison-angle-horaire-altitude-et-azimut-altaz-solaire}}
#' @seealso \code{\link{https://fr.planetcalc.com/320/}}
#' @export
sun.azimuth <- function(day, month, year, hour, minute, seconde=0, longdeg, longmin=0, longsec=0, latdeg, latmin=0, latsec=0)
{
  longitude <- DMS.to.DD(longdeg, longmin, longsec)
  latitude <- DMS.to.DD(latdeg, latmin, latsec)
  #     Heure d'hiver ou d'été
  correction_heure <- 0

  jj <- julian.day(day, month, year, hour, minute, seconde) - correction_heure/24.0 - 2451545.0

  #     Calculs ascension droite et déclinaison
  g <- 357.529 + 0.98560028*jj
  q <- 280.459 + 0.98564736*jj
  l <- q + 1.915 * sin(g*pi/180.0) + 0.020*sin(2*g*pi/180.0) # Ellipticité
  e <- 23.439 - 0.00000036*jj

  ascension_droite <- atan(cos(e*pi/180.0)*sin(l*pi/180.0)/cos(l*pi/180.0))*(180.0/pi)/15.0
  if (cos(l*pi/180.0) < 0) {
    ascension_droite <- 12.0 + ascension_droite
  }
  if (cos(l*pi/180.0)>0 && sin(l*pi/180.0) < 0) {
    ascension_droite <- ascension_droite + 24.0
  }
  declinaison <- asin(sin(e*pi/180.0)*sin(l*pi/180.0))*180.0/pi


  nb_siecle <- jj/36525.0
  heure_siderale1 <- (24110.54841 + (8640184.812866*nb_siecle) + (0.093104*(nb_siecle*nb_siecle)) - (0.0000062*(nb_siecle*nb_siecle*nb_siecle)))/3600.0
  heure_siderale2 <- ((heure_siderale1/24.0) - trunc(heure_siderale1/24.0)) * 24.0

  angleH <- 360.0*heure_siderale2/23.9344
  angleT <- (hour - correction_heure - 12.0 + minute/60.0 + seconde/3600.0)*360.0/23.9344 # jour sidéraux = 23h56min 4,0989s -> 23,93447h -> 0,9972696 jour solaire
  angle <- angleT + angleH

  angle_horaire <- angle - ascension_droite*15.0 + longitude

  #       calculs altitude et azimut

  altitude <- asin( sin(declinaison*pi/180.0)*sin(latitude*pi/180.0) - cos(declinaison*pi/180.0)*cos(latitude*pi/180.0)*cos(angle_horaire*pi/180.0) )*180.0/pi

  azimut <- acos( (sin(declinaison*pi/180.0) - sin(latitude*pi/180.0)*sin(altitude*pi/180.0)) / (cos(latitude*pi/180.0)*cos(altitude*pi/180.0)) )*180.0/pi
  sinazimut <- (cos(declinaison*pi/180.0)*sin(angle_horaire*pi/180.0)) / cos(altitude*pi/180.0)
  if (sinazimut < 0) {
    azimut <- 360 - azimut;
  }

  return(as.numeric(azimut))
}

# Reduction ----

#' Geocentric Axial Dipole
#' Reduction suivant le Dipôle Axiale Centré
#' @param lat.reloc Latitude de réduction. On peut aussi utiliser "Paris" pour 48.85
#' "Madrid" pour 40.4; "Meriden" pour 52.43; "Athenes" pour 37.96
#' @return I.reloc en degré
#' @examples
#' relocate.GAD(65, 47.12)
#' @seealso \code{\link{relocate.VADM}}, \code{\link{relocate.VDM}}, \code{\link{relocate.VGP} }
#' @export
relocate.GAD <- function( I.site, lat.site,  lat.reloc = "Paris")
{
  if (lat.reloc == "Paris") {
    lat.reloc <- 48.85
  }
  if (lat.reloc == "Madrid") {
    lat.reloc <- 40.4
  }
  if (lat.reloc == "Meriden") {
    lat.reloc <- 52.43
  }
  if (lat.reloc == "Athenes") {
    lat.reloc <- 37.96
  }

  I.reloc <- I.site + 0.5 * (3 * cos(I.site * pi /180)^2 + 1) * (lat.reloc - lat.site)

  return(data.frame(I.reloc = I.reloc) )
}

#' Virtual Axial Dipole Moment VADM
#' @examples
#' Reduction par defaut à Paris
#'  relocate.VADM(55, 47.12)
#'  relocate.VADM(55, 47.12, "Paris")
#'  relocate.VADM(55, 47.12, lat.reloc = "48.85")
#' @seealso \code{\link{relocate.GAD}}, \code{\link{relocate.VDM}}, \code{\link{relocate.VGP} }
#' @export
relocate.VADM <- function( F.site,  lat.site,  lat.reloc = "Paris")
{
  if (lat.reloc == "Paris") {
    lat.reloc <- 48.85
  }
  if (lat.reloc == "Madrid") {
    lat.reloc <- 40.4
  }
  if (lat.reloc == "Meriden") {
    lat.reloc <- 52.43
  }
  if (lat.reloc == "Athenes") {
    lat.reloc <- 37.96
  }
  F.reloc <-  as.numeric(F.site * sqrt((3 * sin(lat.reloc * pi /180)^2 + 1) / (3 * sin(lat.site * pi /180)^2 + 1)) )
  return(data.frame(F.reloc = F.reloc))
}

#' Virtual Dipôle Moment
#' @examples
#' relocate.VDM(55, 65,5, 47.12, 2.175,"Paris")
#' @seealso \code{\link{relocate.GAD}}, \code{\link{relocate.VADM}}, \code{\link{relocate.VGP} }
#' @export
relocate.VDM <- function( F.site, I.site, D.site, lat.site, lon.site, lat.reloc = "Paris", lon.reloc = 0)
{
  if (lat.reloc == "Paris") {
    lon.reloc <- 2.3
    lat.reloc <- 48.85
  }
  if (lat.reloc == "Madrid") {
    lon.reloc <- -3.68
    lat.reloc <- 40.4
  }
  if (lat.reloc == "Meriden") {
    lon.reloc <- -1.62
    lat.reloc <- 52.43
  }
  if (lat.reloc == "Athenes") {
    lon.reloc <- 23.72
    lat.reloc <- 37.96
  }

  rad <- pi /180
  LAR <- lat.reloc * rad
  LGR <- lon.reloc * rad
  LA <- lat.site * rad
  LG <- lon.site * rad

  # Calcul de la distance p du pôle
  pol <- (pi / 2) - atan(tan(I.site * rad) / 2)
  lap <- asin(sin(LA) * cos(pol) + cos(LA) * sin(pol) * cos(D.site * rad))
  beta <- asin(sin(pol) * sin(D.site * rad) / cos(lap))
  if (cos(pol) >= (sin(LA) * sin(lap)))
  {
    lgp <- LG + beta
  }
  else
  {
    lgp <- LG + pi - beta
  }

  # Calcul de Iij et Dij corrigés
  pol <- acos(sin(LAR) * sin(lap) + cos(LAR) * cos(lap) * cos(lgp - LGR))
  I.reloc <-  as.numeric(atan(2 * tan((pi / 2) - pol)) )
  F.reloc <-  as.numeric( F.site * sqrt((3 * cos(I.site * rad)^2 + 1) / (3 *  cos(I.reloc)^2 + 1)) )

  return(data.frame(I.reloc = I.reloc /rad , F.reloc = F.reloc, lat.pole = lap / rad , lon.pole = lgp / rad))
}

#' Correction Virtual Geomagnetic Pole
#' Paris   lat:48.85 lon:2.3
#' Madrid  lat:40.4 lon:-3.68
#' Meriden lat:52.43 lon:-1.62
#' Athenes lat:37.96 lon:23.72
#' @examples
#' relocate.VGP(65,5,DMS.to.DD(47,7,15),DMS.to.DD(2,10,12),"Paris")
#' @seealso \code{\link{relocate.GAD, relocate.VADM, relocate.VDM }}
#' @references Westphal, p.117 et Merill McElhinny, p.80
#' @export
relocate.VGP <- function(I.site, D.site, lat.site, lon.site,  lat.reloc = "Paris", lon.reloc = 0)
{
  rad <- as.numeric(pi /180)
  degre <- as.numeric(180 /pi)

  if (lat.reloc == "Paris") {
    lon.reloc <- 2.3
    lat.reloc <- 48.85
  }
  if (lat.reloc == "Madrid") {
    lon.reloc <- -3.68
    lat.reloc <- 40.4
  }
  if (lat.reloc == "Meriden") {
    lon.reloc <- -1.62
    lat.reloc <- 52.43
  }
  if (lat.reloc == "Athenes") {
    lon.reloc <- 23.72
    lat.reloc <- 37.96
  }
  LAR <- lat.reloc * rad
  LGR <- lon.reloc * rad
  LA <- lat.site * rad
  LG <- lon.site * rad
  # Calcul de la distance p du pôle
  pol <- (pi / 2) - atan(tan(I.site * rad) / 2)
  lap <- asin(sin(LA) * cos(pol) + cos(LA) * sin(pol) * cos(D.site * rad))
  beta <- asin(sin(pol) * sin(D.site * rad) / cos(lap))
  if (cos(pol) >= (sin(LA) * sin(lap)))
  {
    lgp <- LG + beta
  }
  else
  {
    lgp <- LG + pi - beta
  }
  # Calcul de Iij et Dij corrigés
  pol <- acos(sin(LAR) * sin(lap) + cos(LAR) * cos(lap) * cos(lgp - LGR))
  I.reloc <- atan(2 * tan((pi / 2) - pol)) * degre
  D.reloc <- asin(sin(lgp - LGR) * cos(lap) / sin(pol)) * degre

  return(data.frame(I.reloc = as.numeric(I.reloc),  D.reloc = as.numeric(D.reloc), lat.pole = as.numeric(lap*degre) , lon.pole = as.numeric(lgp*degre) ) )
}

# igrf13syn
# les coefs
g0 <- c(-31543,-2298, 5922, -677, 2905,-1061,  924, 1121,1022,-1469, -330, 1256,  3,  572,  523,  876, 628,  195,  660,  -69, -361, -210,  134,  -75, -184,  328, -210,  264, 53,  5,  -33,  -86, -124,  -16,  3, 63, 61, -9,  -11, 83, -217,  2,  -58,  -35, 59, 36,  -90,  -69, 70,  -55,  -45,  0,  -13, 34,  -10,  -41, -1,  -21, 28, 18,  -12,  6,  -22, 11,  8,  8, -4,  -14, -9,  7,  1,  -13,  2,  5, -9, 16,  5, -5,  8,  -18,  8, 10,  -20,  1, 14,  -11,  5, 12, -3,  1, -2, -2,  8,  2, 10, -1, -2, -1,  2, -3, -4,  2,  2,  1, -5,  2, -2,  6,  6, -4,  4,  0,  0, -2,  2,  4,  2,  0,  0, -6.)

g1 <- c(-31464,-2298, 5909, -728, 2928,-1086, 1041, 1065, 1037,-1494, -357, 1239, 34,  635,  480,  880,  643,  203,  653,  -77, -380, -201,  146,  -65, -192,  328, -193,  259, 56, -1,  -32,  -93, -125,  -26, 11, 62, 60, -7,  -11, 86, -221,  4,  -57,  -32, 57, 32,  -92,  -67, 70,  -54,  -46,  0,  -14, 33,  -11,  -41,  0,  -20, 28, 18,  -12,  6,  -22, 11,  8,  8, -4,  -15, -9,  7,  1,  -13,  2,  5, -8, 16,  5, -5,  8,  -18,  8, 10,  -20,  1, 14,  -11,  5, 12, -3,  1, -2, -2,  8,  2, 10,  0, -2, -1,  2, -3, -4,  2,  2,  1, -5,  2, -2,  6,  6, -4,  4,  0,  0, -2,  2,  4,  2,  0,  0, -6.)

g2 <- c(-31354,-2297, 5898, -769, 2948,-1128, 1176, 1000, 1058,-1524, -389, 1223, 62,  705,  425,  884,  660,  211,  644,  -90, -400, -189,  160,  -55, -201,  327, -172,  253, 57, -9,  -33, -102, -126,  -38, 21, 62, 58, -5,  -11, 89, -224,  5,  -54,  -29, 54, 28,  -95,  -65, 71,  -54,  -47,  1,  -14, 32,  -12,  -40,  1,  -19, 28, 18,  -13,  6,  -22, 11,  8,  8, -4,  -15, -9,  6,  1,  -13,  2,  5, -8, 16,  5, -5,  8,  -18,  8, 10,  -20,  1, 14,  -11,  5, 12, -3,  1, -2, -2,  8,  2, 10,  0, -2, -1,  2, -3, -4,  2,  2,  1, -5,  2, -2,  6,  6, -4,  4,  0,  0, -2,  2,  4,  2,  0,  0, -6.)

g3 <- c(-31212,-2306, 5875, -802, 2956,-1191, 1309,  917, 1084,-1559, -421, 1212, 84,  778,  360,  887,  678,  218,  631, -109, -416, -173,  178,  -51, -211,  327, -148,  245, 58,  -16,  -34, -111, -126,  -51, 32, 61, 57, -2,  -10, 93, -228,  8,  -51,  -26, 49, 23,  -98,  -62, 72,  -54,  -48,  2,  -14, 31,  -12,  -38,  2,  -18, 28, 19,  -15,  6,  -22, 11,  8,  8, -4,  -15, -9,  6,  2,  -13,  3,  5, -8, 16,  6, -5,  8,  -18,  8, 10,  -20,  1, 14,  -11,  5, 12, -3,  1, -2, -2,  8,  2, 10,  0, -2, -1,  2, -3, -4,  2,  2,  1, -5,  2, -2,  6,  6, -4,  4,  0,  0, -2,  1,  4,  2,  0,  0, -6.)

g4 <- c(-31060,-2317, 5845, -839, 2959,-1259, 1407,  823, 1111,-1600, -445, 1205,  103,  839,  293,  889,  695,  220,  616, -134, -424, -153,  199,  -57, -221,  326, -122,  236, 58,  -23,  -38, -119, -125,  -62, 43, 61, 55,  0,  -10, 96, -233, 11,  -46,  -22, 44, 18, -101,  -57, 73,  -54,  -49,  2,  -14, 29,  -13,  -37,  4,  -16, 28, 19,  -16,  6,  -22, 11,  7,  8, -3,  -15, -9,  6,  2,  -14,  4,  5, -7, 17,  6, -5,  8,  -19,  8, 10,  -20,  1, 14,  -11,  5, 12, -3,  1, -2, -2,  9,  2, 10,  0, -2, -1,  2, -3, -4,  2,  2,  1, -5,  2, -2,  6,  6, -4,  4,  0,  0, -2,  1,  4,  3,  0,  0, -6.)

g5 <- c(-30926,-2318, 5817, -893, 2969,-1334, 1471,  728, 1140,-1645, -462, 1202,  119,  881,  229,  891,  711,  216,  601, -163, -426, -130,  217,  -70, -230,  326,  -96,  226, 58,  -28,  -44, -125, -122,  -69, 51, 61, 54,  3, -9, 99, -238, 14,  -40,  -18, 39, 13, -103,  -52, 73,  -54,  -50,  3,  -14, 27,  -14,  -35,  5,  -14, 29, 19,  -17,  6,  -21, 11,  7,  8, -3,  -15, -9,  6,  2,  -14,  4,  5, -7, 17,  7, -5,  8,  -19,  8, 10,  -20,  1, 14,  -11,  5, 12, -3,  1, -2, -2,  9,  2, 10,  0, -2, -1,  2, -3, -4,  2,  2,  1, -5,  2, -2,  6,  6, -4,  4,  0,  0, -2,  1,  4,  3,  0,  0, -6.)

g6 <- c(-30805,-2316, 5808, -951, 2980,-1424, 1517,  644, 1172,-1692, -480, 1205,  133,  907,  166,  896,  727,  205,  584, -195, -422, -109,  234,  -90, -237,  327,  -72,  218, 60,  -32,  -53, -131, -118,  -74, 58, 60, 53,  4, -9,  102, -242, 19,  -32,  -16, 32,  8, -104,  -46, 74,  -54,  -51,  4,  -15, 25,  -14,  -34,  6,  -12, 29, 18,  -18,  6,  -20, 11,  7,  8, -3,  -15, -9,  5,  2,  -14,  5,  5, -6, 18,  8, -5,  8,  -19,  8, 10,  -20,  1, 14,  -12,  5, 12, -3,  1, -2, -2,  9,  3, 10,  0, -2, -2,  2, -3, -4,  2,  2,  1, -5,  2, -2,  6,  6, -4,  4,  0,  0, -2,  1,  4,  3,  0,  0, -6.)

g7 <- c(-30715,-2306, 5812,-1018, 2984,-1520, 1550,  586, 1206,-1740, -494, 1215,  146,  918,  101,  903,  744,  188,  565, -226, -415,  -90,  249, -114, -241,  329,  -51,  211, 64,  -33,  -64, -136, -115,  -76, 64, 59, 53,  4, -8,  104, -246, 25,  -25,  -15, 25,  4, -106,  -40, 74,  -53,  -52,  4,  -17, 23,  -14,  -33,  7,  -11, 29, 18,  -19,  6,  -19, 11,  7,  8, -3,  -15, -9,  5,  1,  -15,  6,  5, -6, 18,  8, -5,  7,  -19,  8, 10,  -20,  1, 15,  -12,  5, 11, -3,  1, -3, -2,  9,  3, 11,  0, -2, -2,  2, -3, -4,  2,  2,  1, -5,  2, -2,  6,  6, -4,  4,  0,  0, -1,  2,  4,  3,  0,  0, -6.)
g8 <- c(-30654,-2292, 5821,-1106, 2981,-1614, 1566,  528, 1240,-1790, -499, 1232,  163,  916, 43,  914,  762,  169,  550, -252, -405,  -72,  265, -141, -241,  334,  -33,  208, 71,  -33,  -75, -141, -113,  -76, 69, 57, 54,  4, -7,  105, -249, 33,  -18,  -15, 18,  0, -107,  -33, 74,  -53,  -52,  4,  -18, 20,  -14,  -31,  7, -9, 29, 17,  -20,  5,  -19, 11,  7,  8, -3,  -14,  -10,  5,  1,  -15,  6,  5, -5, 19,  9, -5,  7,  -19,  8, 10,  -21,  1, 15,  -12,  5, 11, -3,  1, -3, -2,  9,  3, 11,  1, -2, -2,  2, -3, -4,  2,  2,  1, -5,  2, -2,  6,  6, -4,  4,  0,  0, -1,  2,  4,  3,  0,  0, -6.)
g9 <- c(-30594,-2285, 5810,-1244, 2990,-1702, 1578,  477, 1282,-1834, -499, 1255,  186,  913,  -11,  944,  776,  144,  544, -276, -421,  -55,  304, -178, -253,  346,  -12,  194, 95,  -20,  -67, -142, -119,  -82, 82, 59, 57,  6,  6,  100, -246, 16,  -25, -9, 21,  -16, -104,  -39, 70,  -40,  -45,  0,  -18,  0,  2,  -29,  6,  -10, 28, 15,  -17, 29,  -22, 13,  7, 12, -8,  -21, -5,  -12,  9, -7,  7,  2,  -10, 18,  7,  3,  2,  -11,  5,  -21,  -27,  1, 17,  -11, 29,  3, -9, 16,  4, -3,  9, -4,  6, -3,  1, -4,  8, -3, 11,  5,  1,  1,  2,  -20, -5, -1, -1, -6,  8,  6, -1, -4, -3, -2,  5,  0, -2, -2.)
ga <- c(-30554,-2250, 5815,-1341, 2998,-1810, 1576,  381, 1297,-1889, -476, 1274,  206,  896,  -46,  954,  792,  136,  528, -278, -408,  -37,  303, -210, -240,  349,  3,  211,  103,  -20,  -87, -147, -122,  -76, 80, 54, 57, -1,  4, 99, -247, 33,  -16,  -12, 12,  -12, -105,  -30, 65,  -55,  -35,  2,  -17,  1,  0,  -40, 10, -7, 36,  5,  -18, 19,  -16, 22, 15,  5, -4,  -22, -1,  0, 11,  -21, 15, -8,  -13, 17,  5, -4, -1,  -17,  3, -7,  -24, -1, 19,  -25, 12, 10,  2,  5,  2, -5,  8, -2,  8,  3,  -11,  8, -7, -8,  4, 13, -1, -2, 13,  -10, -4,  2,  4, -3, 12,  6,  3, -3,  2,  6, 10, 11,  3,  8)

gb <- c(-30500,-2215, 5820,-1440, 3003,-1898, 1581,  291, 1302,-1944, -462, 1288,  216,  882,  -83,  958,  796,  133,  510, -274, -397,  -23,  290, -230, -229,  360, 15,  230,  110,  -23,  -98, -152, -121,  -69, 78, 47, 57, -9,  3, 96, -247, 48, -8,  -16,  7,  -12, -107,  -24, 65,  -56,  -50,  2,  -24, 10, -4,  -32,  8,  -11, 28,  9,  -20, 18,  -18, 11,  9, 10, -6,  -15,  -14,  5,  6,  -23, 10,  3, -7, 23,  6, -4,  9,  -13,  4,  9,  -11, -4, 12, -5,  7,  2,  6,  4, -2,  1, 10,  2,  7,  2, -6,  5,  5, -3, -5, -4, -1,  0,  2, -8, -3, -2,  7, -4,  4,  1, -2, -3,  6,  7, -2, -1,  0, -3.)
gc <- c(-30421,-2169, 5791,-1555, 3002,-1967, 1590,  206, 1302,-1992, -414, 1289,  224,  878, -130,  957,  800,  135,  504, -278, -394,  3,  269, -255, -222,  362, 16,  242,  125,  -26, -117, -156, -114,  -63, 81, 46, 58,  -10,  1, 99, -237, 60, -1,  -20, -2,  -11, -113,  -17, 67,  -56,  -55,  5,  -28, 15, -6,  -32,  7, -7, 23, 17,  -18,  8,  -17, 15,  6, 11, -4,  -14,  -11,  7,  2,  -18, 10,  4, -5, 23, 10,  1,  8,  -20,  4,  6,  -18,  0, 12, -9,  2,  1,  0,  4, -3, -1,  9, -2,  8,  3,  0, -1,  5,  1, -3,  4,  4,  1,  0,  0, -1,  2,  4, -5,  6,  1,  1, -1, -1,  6,  2,  0,  0, -7.)
gd <- c(-30334,-2119, 5776,-1662, 2997,-2016, 1594,  114, 1297,-2038, -404, 1292,  240,  856, -165,  957,  804,  148,  479, -269, -390, 13,  252, -269, -219,  358, 19,  254,  128,  -31, -126, -157,  -97,  -62, 81, 45, 61,  -11,  8,  100, -228, 68,  4,  -32,  1, -8, -111, -7, 75,  -57,  -61,  4,  -27, 13, -2,  -26,  6, -6, 26, 13,  -23,  1,  -12, 13,  5,  7, -4,  -12,  -14,  9,  0,  -16,  8,  4, -1, 24, 11, -3,  4,  -17,  8, 10,  -22,  2, 15,  -13,  7, 10, -4, -1, -5, -1, 10,  5, 10,  1, -4, -2,  1, -2, -3,  2,  2,  1, -5,  2, -2,  6,  4, -4,  4,  0,  0, -2,  2,  3,  2,  0,  0, -6.)
ge <- c(-30220,-2068, 5737,-1781, 3000,-2047, 1611, 25, 1287,-2091, -366, 1278,  251,  838, -196,  952,  800,  167,  461, -266, -395, 26,  234, -279, -216,  359, 26,  262,  139,  -42, -139, -160,  -91,  -56, 83, 43, 64,  -12, 15,  100, -212, 72,  2,  -37,  3, -6, -112,  1, 72,  -57,  -70,  1,  -27, 14, -4,  -22,  8, -2, 23, 13,  -23, -2,  -11, 14,  6,  7, -2,  -15,  -13,  6, -3,  -17,  5,  6,  0, 21, 11, -6,  3,  -16,  8, 10,  -21,  2, 16,  -12,  6, 10, -4, -1, -5,  0, 10,  3, 11,  1, -2, -1,  1, -3, -3,  1,  2,  1, -5,  3, -1,  4,  6, -4,  4,  0,  1, -1,  0,  3,  3,  1, -1, -4.)
gf<- c(-30100,-2013, 5675,-1902, 3010,-2067, 1632,  -68, 1276,-2144, -333, 1260,  262,  830, -223,  946,  791,  191,  438, -265, -405, 39,  216, -288, -218,  356, 31,  264,  148,  -59, -152, -159,  -83,  -49, 88, 45, 66,  -13, 28, 99, -198, 75,  1,  -41,  6, -4, -111, 11, 71,  -56,  -77,  1,  -26, 16, -5,  -14, 10,  0, 22, 12,  -23, -5,  -12, 14,  6,  6, -1,  -16,  -12,  4, -8,  -19,  4,  6,  0, 18, 10,  -10,  1,  -17,  7, 10,  -21,  2, 16,  -12,  7, 10, -4, -1, -5, -1, 10,  4, 11,  1, -3, -2,  1, -3, -3,  1,  2,  1, -5,  3, -2,  4,  5, -4,  4, -1,  1, -1,  0,  3,  3,  1, -1, -5.)

gg <- c(-29992,-1956, 5604,-1997, 3027,-2129, 1663, -200, 1281, -2180, -336, 1251,  271,  833, -252,  938, 782,  212,  398, -257, -419, 53,  199, -297, -218,  357, 46,  261,  150,  -74, -151, -162,  -78,  -48, 92, 48, 66,  -15, 42, 93, -192, 71,  4,  -43, 14, -2, -108, 17, 72,  -59,  -82,  2,  -27, 21, -5,  -12, 16,  1, 18, 11,  -23, -2,  -10, 18,  6,  7,  0,  -18,  -11,  4, -7,  -22,  4,  9,  3, 16,  6,  -13, -1,  -15,  5, 10,  -21,  1, 16,  -12,  9,  9, -5, -3, -6, -1,  9,  7, 10,  2, -6, -5,  2, -4, -4,  1,  2,  0, -5,  3, -2,  6,  5, -4,  3,  0,  1, -1,  2,  4,  3,  0,  0, -6.)
gi <- c(-29873,-1905, 5500,-2072, 3044,-2197, 1687, -306, 1296,-2208, -310, 1247,  284,  829, -297,  936,  780,  232,  361, -249, -424, 69,  170, -297, -214,  355, 47,  253,  150,  -93, -154, -164,  -75,  -46, 95, 53, 65,  -16, 51, 88, -185, 69,  4,  -48, 16, -1, -102, 21, 74,  -62,  -83,  3,  -27, 24, -2, -6, 20,  4, 17, 10,  -23,  0, -7, 21,  6,  8,  0,  -19,  -11,  5, -9,  -23,  4, 11,  4, 14,  4,  -15, -4,  -11,  5, 10,  -21,  1, 15,  -12,  9,  9, -6, -3, -6, -1,  9,  7,  9,  1, -7, -5,  2, -4, -4,  1,  3,  0, -5,  3, -2,  6,  5, -4,  3,  0,  1, -1,  2,  4,  3,  0,  0, -6.)
gj <- c(-29775,-1848, 5406,-2131, 3059,-2279, 1686, -373, 1314,-2239, -284, 1248,  293,  802, -352,  939,  780,  247,  325, -240, -423, 84,  141, -299, -214,  353, 46,  245,  154, -109, -153, -165,  -69,  -36, 97, 61, 65,  -16, 59, 82, -178, 69,  3,  -52, 18,  1,  -96, 24, 77,  -64,  -80,  2,  -26, 26,  0, -1, 21,  5, 17,  9,  -23,  0, -4, 23,  5, 10, -1,  -19,  -10,  6,  -12,  -22,  3, 12,  4, 12,  2,  -16, -6,  -10,  4,  9,  -20,  1, 15,  -12, 11,  9, -7, -4, -7, -2,  9,  7,  8,  1, -7, -6,  2, -3, -4,  2,  2,  1, -5,  3, -2,  6,  4, -4,  3,  0,  1, -2,  3,  3,  3, -1,  0, -6.)

gk <- c(-29692,-1784, 5306,-2200, 3070,-2366, 1681, -413, 1335,-2267, -262, 1249,  302,  759, -427,  940,  780,  262,  290, -236, -418, 97,  122, -306, -214,  352, 46,  235,  165, -118, -143, -166,  -55,  -17,  107, 68, 67,  -17, 68, 72, -170, 67, -1,  -58, 19,  1,  -93, 36, 77,  -72,  -69,  1,  -25, 28,  4,  5, 24,  4, 17,  8,  -24, -2, -6, 25,  6, 11, -6,  -21, -9,  8,  -14,  -23,  9, 15,  6, 11, -5,  -16, -7, -4,  4,  9,  -20,  3, 15,  -10, 12,  8, -6, -8, -8, -1,  8, 10,  5, -2, -8, -8,  3, -3, -6,  1,  2,  0, -4,  4, -1,  5,  4, -5,  2, -1,  2, -2,  5,  1,  1, -2,  0, -7.)
gk <- c(gk, rep(0, 75))

gl <- c(-29619.4,-1728.2, 5186.1,-2267.7, 3068.4,-2481.6, 1670.9, -458.0, 1339.6,-2288.0, -227.6, 1252.1,  293.4,  714.5, -491.1,  932.3,  786.8,  272.6,  250.0, -231.9, -403.0,  119.8,  111.3, -303.8, -218.8,  351.4, 43.8,  222.3,  171.9, -130.4, -133.1, -168.6,  -39.3,  -12.9,  106.3, 72.3, 68.2,  -17.4, 74.2, 63.7, -160.9, 65.1, -5.9,  -61.2, 16.9,  0.7,  -90.4, 43.8, 79.0,  -74.0,  -64.6,  0.0,  -24.2, 33.3,  6.2,  9.1, 24.0,  6.9, 14.8,  7.3,  -25.4, -1.2, -5.8, 24.4,  6.6, 11.9, -9.2,  -21.5, -7.9,  8.5,  -16.6,  -21.5,  9.1, 15.5,  7.0,  8.9, -7.9,  -14.9, -7.0, -2.1,  5.0,  9.4,  -19.7,  3.0, 13.4, -8.4, 12.5,  6.3, -6.2, -8.9, -8.4, -1.5,  8.4,  9.3,  3.8, -4.3, -8.2, -8.2,  4.8, -2.6, -6.0,  1.7,  1.7,  0.0, -3.1,  4.0, -0.5,  4.9,  3.7, -5.9,  1.0, -1.2,  2.0, -2.9,  4.2,  0.2,  0.3, -2.2, -1.1, -7.4,  2.7, -1.7,  0.1, -1.9,  1.3,  1.5, -0.9, -0.1, -2.6,  0.1,  0.9, -0.7, -0.7,  0.7, -2.8,  1.7, -0.9,  0.1, -1.2,  1.2, -1.9,  4.0, -0.9, -2.2, -0.3, -0.4,  0.2,  0.3,  0.9,  2.5, -0.2, -2.6,  0.9,  0.7, -0.5,  0.3,  0.3,  0.0, -0.3,  0.0, -0.4,  0.3, -0.1, -0.9, -0.2, -0.4, -0.4,  0.8, -0.2, -0.9, -0.9,  0.3,  0.2,  0.1,  1.8, -0.4, -0.4,  1.3, -1.0, -0.4, -0.1,  0.7,  0.7, -0.4,  0.3,  0.3,  0.6, -0.1,  0.3,  0.4, -0.2,  0.0, -0.5,  0.1, -0.9)

gm <- c(-29554.63,-1669.05, 5077.99,-2337.24, 3047.69,-2594.50, 1657.76, -515.43, 1336.30,-2305.83, -198.86, 1246.39,  269.72,  672.51, -524.72,  920.55,  797.96,  282.07,  210.65, -225.23, -379.86,  145.15,  100.00, -305.36, -227.00,  354.41, 42.72,  208.95,  180.25, -136.54, -123.45, -168.05,  -19.57,  -13.55,  103.85, 73.60, 69.56,  -20.33, 76.74, 54.75, -151.34, 63.63,  -14.58,  -63.53, 14.58,  0.24,  -86.36, 50.94, 79.88,  -74.46,  -61.14, -1.65,  -22.57, 38.73,  6.82, 12.30, 25.35,  9.37, 10.93,  5.42,  -26.32,  1.94, -4.64, 24.80,  7.62, 11.20,  -11.73,  -20.88, -6.88,  9.83,  -18.11,  -19.71, 10.17, 16.22,  9.36,  7.61,  -11.25,  -12.76, -4.87, -0.06,  5.58,  9.76,  -20.11,  3.58, 12.69, -6.94, 12.67,  5.01, -6.72,  -10.76, -8.16, -1.25,  8.10,  8.76,  2.92, -6.66, -7.73, -9.22,  6.01, -2.17, -6.12,  2.19,  1.42,  0.10, -2.35,  4.46, -0.15,  4.76,  3.06, -6.58,  0.29, -1.01,  2.06, -3.47,  3.77, -0.86, -0.21, -2.31, -2.09, -7.93,  2.95, -1.60,  0.26, -1.88,  1.44,  1.44, -0.77, -0.31, -2.27,  0.29,  0.90, -0.79, -0.58,  0.53, -2.69,  1.80, -1.08,  0.16, -1.58,  0.96, -1.90,  3.99, -1.39, -2.15, -0.29, -0.55,  0.21,  0.23,  0.89,  2.38, -0.38, -2.63,  0.96,  0.61, -0.30,  0.40,  0.46,  0.01, -0.35,  0.02, -0.36,  0.28,  0.08, -0.87, -0.49, -0.34, -0.08,  0.88, -0.16, -0.88, -0.76,  0.30,  0.33,  0.28,  1.72, -0.43, -0.54,  1.18, -1.07, -0.37, -0.04,  0.75,  0.63, -0.26,  0.21,  0.35,  0.53, -0.05,  0.38,  0.41, -0.22, -0.10, -0.57, -0.18, -0.82)

gp <- c(-29496.57,-1586.42, 4944.26,-2396.06, 3026.34,-2708.54, 1668.17, -575.73, 1339.85,-2326.54, -160.40, 1232.10,  251.75,  633.73, -537.03,  912.66,  808.97,  286.48,  166.58, -211.03, -356.83,  164.46, 89.40, -309.72, -230.87,  357.29, 44.58,  200.26,  189.01, -141.05, -118.06, -163.17, -0.01, -8.03,  101.04, 72.78, 68.69,  -20.90, 75.92, 44.18, -141.40, 61.54,  -22.83,  -66.26, 13.10,  3.02,  -78.09, 55.40, 80.44,  -75.00,  -57.80, -4.55,  -21.20, 45.24,  6.54, 14.00, 24.96, 10.46,  7.03,  1.64,  -27.61,  4.92, -3.28, 24.41,  8.21, 10.84,  -14.50,  -20.03, -5.59, 11.83,  -19.34,  -17.41, 11.61, 16.71, 10.85,  6.96,  -14.05,  -10.74, -3.54,  1.64,  5.50,  9.45,  -20.54,  3.45, 11.51, -5.27, 12.75,  3.13, -7.14,  -12.38, -7.42, -0.76,  7.97,  8.43,  2.14, -8.42, -6.08,  -10.08,  7.01, -1.94, -6.24,  2.73,  0.89, -0.10, -1.07,  4.71, -0.16,  4.44,  2.45, -7.22, -0.33, -0.96,  2.13, -3.95,  3.09, -1.99, -1.03, -1.97, -2.80, -8.31,  3.05, -1.48,  0.13, -2.03,  1.67,  1.65, -0.66, -0.51, -1.76,  0.54,  0.85, -0.79, -0.39,  0.37, -2.51,  1.79, -1.27,  0.12, -2.11,  0.75, -1.94,  3.75, -1.86, -2.12, -0.21, -0.87,  0.30,  0.27,  1.04,  2.13, -0.63, -2.49,  0.95,  0.49, -0.11,  0.59,  0.52,  0.00, -0.39,  0.13, -0.37,  0.27,  0.21, -0.86, -0.77, -0.23,  0.04,  0.87, -0.09, -0.89, -0.87,  0.31,  0.30,  0.42,  1.66, -0.45, -0.59,  1.08, -1.14, -0.31, -0.07,  0.78,  0.54, -0.18,  0.10,  0.38,  0.49,  0.02,  0.44,  0.42, -0.25, -0.26, -0.53, -0.26, -0.79)

gq <- c(-29441.46,-1501.77, 4795.99,-2445.88, 3012.20,-2845.41, 1676.35, -642.17, 1350.33,-2352.26, -115.29, 1225.85,  245.04,  581.69, -538.70,  907.42,  813.68,  283.54,  120.49, -188.43, -334.85,  180.95, 70.38, -329.23, -232.91,  360.14, 46.98,  192.35,  196.98, -140.94, -119.14, -157.40, 15.98,  4.30,  100.12, 69.55, 67.57,  -20.61, 72.79, 33.30, -129.85, 58.74,  -28.93,  -66.64, 13.14,  7.35,  -70.85, 62.41, 81.29,  -75.99,  -54.27, -6.79,  -19.53, 51.82,  5.59, 15.07, 24.45,  9.32,  3.27, -2.88,  -27.50,  6.61, -2.32, 23.98,  8.89, 10.04,  -16.78,  -18.26, -3.16, 13.18,  -20.56,  -14.60, 13.33, 16.16, 11.76,  5.69,  -15.98, -9.10, -2.02,  2.26,  5.33,  8.83,  -21.77,  3.02, 10.76, -3.22, 11.74,  0.67, -6.74,  -13.20, -6.88, -0.10,  7.79,  8.68,  1.04, -9.06, -3.89,  -10.54,  8.44, -2.01, -6.26,  3.28,  0.17, -0.40,  0.55,  4.55, -0.55,  4.40,  1.70, -7.92, -0.67, -0.61,  2.13, -4.16,  2.33, -2.85, -1.80, -1.12, -3.59, -8.72,  3.00, -1.40,  0.00, -2.30,  2.11,  2.08, -0.60, -0.79, -1.05,  0.58,  0.76, -0.70, -0.20,  0.14, -2.12,  1.70, -1.44, -0.22, -2.57,  0.44, -2.01,  3.49, -2.34, -2.09, -0.16, -1.08,  0.46,  0.37,  1.23,  1.75, -0.89, -2.19,  0.85,  0.27,  0.10,  0.72,  0.54, -0.09, -0.37,  0.29, -0.43,  0.23,  0.22, -0.89, -0.94, -0.16, -0.03,  0.72, -0.02, -0.92, -0.88,  0.42,  0.49,  0.63,  1.56, -0.42, -0.50,  0.96, -1.24, -0.19, -0.10,  0.81,  0.42, -0.13, -0.04,  0.38,  0.48,  0.08,  0.48,  0.46, -0.30, -0.35, -0.43, -0.36, -0.71)
gr <- c(-29404.8, -1450.9,  4652.5, -2499.6,  2982.0, -2991.6, 1677.0,  -734.6,  1363.2, -2381.2, -82.1,  1236.2,  241.9, 525.7,  -543.4, 903.0, 809.5, 281.9, 86.3,  -158.4,  -309.4, 199.7,  48.0,  -349.7, -234.3, 363.2,  47.7, 187.8, 208.3,  -140.7, -121.2,  -151.2,  32.3,  13.5,  98.9,  66.0, 65.5, -19.1,  72.9,  25.1,  -121.5,  52.8,  -36.2, -64.5,  13.5, 8.9, -64.7,  68.1, 80.6, -76.7, -51.5,  -8.2, -16.9,  56.5,  2.2,  15.8,  23.5, 6.4,  -2.2,  -7.2,  -27.2, 9.8,  -1.8,  23.7, 9.7, 8.4,  -17.6, -15.3,  -0.5,  12.8, -21.1, -11.7, 15.3,  14.9,  13.7, 3.6, -16.5,  -6.9, -0.3, 2.8, 5.0, 8.4, -23.4, 2.9, 11.0,  -1.5, 9.8,  -1.1,  -5.1, -13.2, -6.3, 1.1, 7.8, 8.8, 0.4,  -9.3, -1.4, -11.9, 9.6,  -1.9,  -6.2, 3.4, -0.1,  -0.2, 1.7, 3.6,  -0.9, 4.8,  0.7,  -8.6,  -0.9,  -0.1, 1.9,  -4.3,  1.4,  -3.4,  -2.4,  -0.1,  -3.8,  -8.8,  3.0,  -1.4, 0.0,  -2.5, 2.5, 2.3, -0.6,  -0.9,  -0.4, 0.3, 0.6,  -0.7, -0.2,  -0.1,  -1.7, 1.4,  -1.6,  -0.6, -3.0, 0.2,  -2.0, 3.1,  -2.6,  -2.0, -0.1,  -1.2, 0.5, 0.5, 1.3, 1.4, -1.2,  -1.8, 0.7, 0.1, 0.3, 0.8,  0.5,  -0.2,  -0.3, 0.6,  -0.5, 0.2,  0.1,  -0.9,  -1.1, 0.0,  -0.3, 0.5,  0.1,  -0.9,  -0.9, 0.5, 0.6, 0.7,  1.4,  -0.3,  -0.4, 0.8,  -1.3, 0.0, -0.1, 0.8, 0.3, 0.0,  -0.1, 0.4,  0.5, 0.1, 0.5, 0.5,  -0.4,  -0.5, -0.4,  -0.4,  -0.6)

gs <- c(  5.7, 7.4, -25.9, -11.0,  -7.0, -30.2,  -2.1, -22.4, 2.2,  -5.9, 6.0, 3.1,  -1.1, -12.0, 0.5,  -1.2,  -1.6,  -0.1,  -5.9, 6.5, 5.2, 3.6,  -5.1,  -5.0,  -0.3, 0.5, 0.0,  -0.6, 2.5, 0.2,  -0.6, 1.3, 3.0, 0.9, 0.3,  -0.5,  -0.3, 0.0, 0.4,  -1.6, 1.3,  -1.3,  -1.4, 0.8, 0.0, 0.0, 0.9, 1.0,  -0.1,  -0.2, 0.6, 0.0, 0.6, 0.7,  -0.8, 0.1,  -0.2,  -0.5,  -1.1,  -0.8, 0.1, 0.8, 0.3, 0.0, 0.1,  -0.2,  -0.1, 0.6, 0.4,  -0.2,  -0.1, 0.5, 0.4,  -0.3, 0.3,  -0.4,  -0.1, 0.5, 0.4, 0.0)
gs <- c(gs, rep(0, 115))

gh <- c( g0, g1, g2, g3, g4, g5, g6, g7, g8, g9, ga, gb, gc, gd, ge, gf, gg, gi, gj, gk, gl, gm, gp, gq, gr, gs)
p <- rep(0, 105)
q <- rep(0, 105)
cl <- rep(0, 13)
sl <- rep(0, 13)

#'     This is a synthesis routine for the 13th generation IGRF as agreed
#'    in December 2019 by IAGA Working Group V-MOD. It is valid 1900.0 to
#'     2025.0 inclusive. Values for dates from 1945.0 to 2015.0 inclusive are
#'     definitive, otherwise they are non-definitive.
#'   INPUT
#'@param isv   = 0 if main-field values are required,  = 1 if secular variation values are required
#'@param date  = year A.D. Must be greater than or equal to 1900.0 and less than or equal to 2030.0. Warning message is given
#'             for dates greater than 2025.0. Must be double precision.
#'@param itype = 1 if geodetic (spheroid);  = 2 if geocentric (sphere)
#'@param alt   = height in km above sea level if itype = 1;   = distance from centre of Earth in km if itype = 2 (>3485 km)
#'@param colat = colatitude (0-180)
#'@param elong = east-longitude (0-360)
#'     alt, colat and elong must be double precision.
#'   OUTPUT
#'     x     = north component (nT) if isv = 0, nT/year if isv = 1
#'     y     = east component (nT) if isv = 0, nT/year if isv = 1
#'     z     = vertical component (nT) if isv = 0, nT/year if isv = 1
#'     f     = total intensity (nT) if isv = 0, rubbish if isv = 1
#'
#'     To get the other geomagnetic elements (D, I, H and secular
#'     variations dD, dH, dI and dF) use routines ptoc and ptocsv.
#'
#'     Adapted from 8th generation version to include new maximum degree for
#'     main-field models for 2000.0 and onwards and use WGS84 spheroid instead
#'     of International Astronomical Union 1966 spheroid as recommended by IAGA
#'     in July 2003. Reference radius remains as 6371.2 km - it is NOT the mean
#'     radius (= 6371.0 km) but 6371.2 km is what is used in determining the
#'     coefficients. Adaptation by Susan Macmillan, August 2003 (for
#'     9th generation), December 2004, December 2009 & December 2014;
#'     by William Brown, December 2019, February 2020.
#'
#'     Coefficients at 1995.0 incorrectly rounded (rounded up instead of
#'     to even) included as these are the coefficients published in Excel
#'     spreadsheet July 2005.
#' igrf13syn(isv=isv, date=date, itype=itype , alt=alt, colat= colat, elong=elong)
#' @seealso \code{\link{https://www.ngdc.noaa.gov/IAGA/vmod/igrf.html}}
#' @references Thébault, E., Finlay, C.C., Beggan, C.D. et al. International Geomagnetic Reference Field: the 12th generation. Earth Planet Sp 67, 79 (2015).
#' https://doi.org/10.1186/s40623-015-0228-9 (\link{https://earth-planets-space.springeropen.com/articles/10.1186/s40623-015-0228-9}): the 12th generation, Erwan Thébault, Christopher C Finlay, Ciarán D Beggan, Patrick Alken, Julien Aubert, Olivier Barrois, Francois Bertrand, Tatiana Bondar, Axel Boness, Laura Brocco, Elisabeth Canet, Aude Chambodut, Arnaud Chulliat, Pierdavide Coïsson, Francois Civet, Aimin Du, Alexandre Fournier, Isabelle Fratter, Nicolas Gillet, Brian Hamilton, Mohamed Hamoudi, Gauthier Hulot, Thomas Jager, Monika Korte, Weijia Kuang, Xavier Lalanne, Benoit Langlais, Jean-Michel Léger, Vincent Lesur, Frank J Lowes et al. Earth, Planets and Space 2015, 67:79 (27 May 2015)
#'
#' @examples
#' isv <- 1
#' date <- 2000
#' itype <- 1
#' alt <- 50
#' colat <- 90 - 45
#' elong <- 10
#'
#' @export
igrf13syn <- function(isv, date, itype, alt, colat, elong) {
  # set initial values

  x     = 0.0
  y     = 0.0
  z     = 0.0

  if (date < 1900.0 || date > 2025.0) {
    message(' Date must be in the range [1900.0..2030.0]')
    break;
  }

  if (date > 2020.0) {
    t     = date - 2020.0
    tc    = 1.0
    if (isv == 1) {
      t = 1.0
      tc = 0.0
    }
    #     pointer for last coefficient in pen-ultimate set of MF coefficients...

    ll    = 3255
    nmx   = 13
    nc    = nmx*(nmx+2)
    kmx   = (nmx+1)*(nmx+2)/2


  } else {
    t    = 0.2*(date - 1900.0)
    ll    = t
    one   = ll
    t     = t - one

    #
    #SH models before 1995.0 are only to degree 10
    #
    if (date < 1995.0) {
      nmx   = 10
      nc    = nmx*(nmx+2)
      ll    = nc*ll
      kmx   = (nmx+1)*(nmx+2)/2
    }
    else {
      nmx   = 13
      nc    = nmx*(nmx+2)
      ll    = 0.2*(date - 1995.0)
      #
      #     19 is the number of SH models that extend to degree 10
      #
      ll    = 120*19 + nc*ll
      kmx   = (nmx+1)*(nmx+2)/2
    }

    tc    = 1.0 - t
    if (isv== 1) {
      tc = -0.2
      t = 0.2
    }
  }


  r     = alt
  one   = as.numeric(colat*0.017453292) # con version deg to rad
  ct    = as.numeric( cos(one) )
  st    = as.numeric( sin(one) )
  one   = elong*0.017453292
  cl[1] = as.numeric(cos(one) )
  sl[1] = as.numeric(sin(one) )
  cd    = 1.0
  sd    = 0.0
  l     = 1
  m     = 1
  n     = 0
  if (itype != 2) {

    #     conversion from geodetic to geocentric coordinates
    #     (using the WGS84 spheroid)

    a2    = 40680631.6
    b2    = 40408296.0
    one   = a2*st*st
    two   = b2*ct*ct
    three = one + two
    rho   = as.numeric( sqrt(three))
    r     = as.numeric( sqrt(alt*(alt + 2.0*rho) + (a2*one + b2*two)/three) )
    cd    = (alt + rho)/r
    sd    = (a2 - b2)/rho*ct*st/r
    one   = ct
    ct    = ct*cd -  st*sd
    st    = st*cd + one*sd

  }

  ratio = 6371.2/r
  rr    = ratio*ratio

  #     computation of Schmidt quasi-normal coefficients p and x(=q)

  p[1]  = 1.0
  p[3]  = st
  q[1]  = 0.0
  q[3]  =  ct
  for (k in seq(2, kmx) ) {

    if (n < m) {
      m     = 0
      n     = n + 1
      rr    = rr*ratio
      fn    = n
      gn    = n - 1
    }

    fm    = m
    if (m != n) {
      gmm    = m*m
      one   = as.numeric( sqrt(fn*fn - gmm))
      two   = as.numeric( sqrt(gn*gn - gmm)/one )
      three = (fn + gn)/one
      i     = k - n
      j     = i - n + 1
      p[k]  = three*ct*p[i] - two*p[j]
      q[k]  = three*(ct*q[i] - st*p[i]) - two*q[j]
    } else if (k != 3) {
      one   = as.numeric(sqrt(1.0 - 0.5/fm) )
      j     = k - n - 1

      p[k]  = one*st*p[j]
      q[k]  = one*(st*q[j] + ct*p[j])
      cl[m] = cl[m-1]*cl[1] - sl[m-1]*sl[1]
      sl[m] = sl[m-1]*cl[1] + cl[m-1]*sl[1]
    }

    #     synthesis of x, y and z in geocentric coordinates

    lm    = ll + l
    one   = (tc*gh[lm] + t*gh[lm+nc])*rr
    if (m == 0) {
      x     = x + one*q[k]
      z     = z - (fn + 1.0)*one*p[k]
      l     = l + 1
    } else {
      two   = (tc*gh[lm+1] + t*gh[lm+nc+1])*rr
      three = one*cl[m] + two*sl[m]
      x     = x + three*q[k]
      z     = z - (fn + 1.0)*three*p[k]
      if (st == 0.0) {
        y     = y + (one*sl[m] - two*cl[m])*q[k]*ct
      } else {
        y     = y + (one*sl[m] - two*cl[m])*fm*p[k]/st
      }

      l     = l + 2

    }


    m     = m + 1
  }
  #     conversion to coordinate system specified by itype

  one   = x
  x     = x*cd +   z*sd
  z     = z*cd - one*sd
  f     = sqrt(x*x + y*y + z*z)
  ret<-NULL
  ret$x <- x
  ret$y <- y
  ret$z <- z
  ret$f <-f

  if (isv == 0) {
    ret$type <- 'Main-Field'
  }  else {
    ret$type <- 'Secular Variation'
  }

  hsq = x*x + y*y
  hoz  = sqrt(hsq)
  dec = atan2(y,x)
  inc = atan2(z,hoz)
  ret$dec <- dec/pi*180
  ret$inc <- inc/pi*180
  ret$hoz <- hoz


  return(ret)
}


#'   Standart calcul en plot to study paleo and archeo magnetic magnétisation
#'   by convention we use the two last character of the variable step to indicate the type of manip
#'   example 100RA , 100 is the temperature (step.value) , R is the sens of the magnetisation and A is the name of the step (step.name)
#' @references  Coe 1978 : DOI: 10.1029/JB083iB04p01740
#' Prévost et Al. 1985 DOI: 10.1029/JB090iB12p10417
#'
#'@param mesures data.frame with the package convention format
#'@param relative plot with a relative value in percent
#'@param verbose show comment
#'@param show.plot  display the plot
#'@param TH lab field
#'@param aim.coef used to correct the measurement often 1E-10x1E6xvolume
#'@param show.step.value display on the plot the value of each step
#'@param R.mark = 'R', V.mark = 'V' conventional notation to mark the sens of the magnetisation with the step name
#'@param P.mark = 'P' mark the pTRM check loop
#'@param  L.mark = "L", Q.mark = "Q" mark to indiquate the slow cooling step and the quick (fast-normal) cooling step
#'@param step.J0 is the step of the natural magnetisation eg "0N0"
#'@param begin.step.value the first step.value (temperature) used to determinate the magnetic field
#'@param end.step.value the last step.value (temperature) used to determinate the magnetic field
#'@param loop.col the color of the line showing the the check loop process
#'@param pt.col the color of the line and the plot
#' @export
arai <- function(mesures, relative = TRUE, verbose = TRUE, show.plot = TRUE, TH = 60, aim.coef = 1E-10*1E6, step.J0 = "20N0", show.step.value = FALSE, R.mark = 'R', V.mark = 'V', P.mark = 'P', L.mark = "L", Q.mark = "Q", pt.col = "blue", loop.col = "forestgreen", begin.step.value = 0, end.step.value = 1000) {
  # __________________
  if (is.null(step.J0)) {
    step.J0 <- mes.sel$step[1]
  }

  ATRR <- mesures[which(substr(mesures$step, nchar(mesures$step)-1, nchar(mesures$step)-1 ) == R.mark),]
  ATRV <- mesures[which(substr(mesures$step, nchar(mesures$step)-1, nchar(mesures$step)-1 ) == V.mark),]

  J0 <- mesures [which(trimws(mesures$step) ==trimws(step.J0) ), ]
  if (J0$F < 0 ) {
    warning("error on step.J0, F must be positive")

  } else {
    ATRR <- rbind.data.frame(J0, ATRR)
    ATRV <- rbind.data.frame(J0, ATRV)
  }

  for (i in 1:length(ATRR$step) ) {
    ATRR$step.name[i] <- substr(ATRR$step[i], nchar(ATRR$step[i]), nchar(ATRR$step[i]) )
  }
  for (i in 1:length(ATRV$step) ) {
    ATRV$step.name[i] <- substr(ATRV$step[i], nchar(ATRV$step[i]), nchar(ATRV$step[i]) )
  }


  ARN <- NULL
  ATR <- NULL
  RN <- NULL
  TR <- NULL
  pt.col.res <- NULL
  # tableau ATR vs ARN
  for (i in ATRR$step.name) {
    iATRR <- which(ATRR$step.name == i)
    iATRV <- which(ATRV$step.name == i)

    atrXR <- ATRR[iATRR, ]$X *aim.coef
    atrXV <- ATRV[iATRV, ]$X *aim.coef

    atrYR <- ATRR[iATRR, ]$Y *aim.coef
    atrYV <- ATRV[iATRV, ]$Y *aim.coef

    atrZR <- ATRR[iATRR, ]$Z *aim.coef
    atrZV <- ATRV[iATRV, ]$Z *aim.coef

    if (ATRR[iATRR, ]$step.value != ATRV[iATRV, ]$step.value) {
      warning(paste0("Error in step.name within step.value = ", ATRR[iATRR, ]$step, ATRV[iATRV, ]$step) )
      next
    }

    RN$X <- (atrXR + atrXV)/2
    RN$Y <- (atrYR + atrYV)/2
    RN$Z <- (atrZR + atrZV)/2

    RN <- to.polar(RN$X, RN$Y, RN$Z)
    RN$step.value <- ATRR[iATRR, ]$step.value
    RN$step.name <- i

    ARN <-rbind.data.frame(ARN, RN, stringsAsFactors = FALSE)

    TR$X <- (atrXR - atrXV)/2
    TR$Y <- (atrYR - atrYV)/2
    TR$Z <- (atrZR - atrZV)/2
    TR <- to.polar(TR$X, TR$Y, TR$Z)
    TR$step.value <- ATRR[iATRR, ]$step.value
    TR$step.name <- i

    ATR <-rbind.data.frame(ATR, TR, stringsAsFactors = FALSE)
    if ((RN$step.value >= begin.step.value) && (RN$step.value <= end.step.value)) {
      pt.col.res <- c(pt.col.res, "red")
    } else {
      pt.col.res <- c(pt.col.res, pt.col)
    }
  }
  xlim <- range(ATR$F)
  ylim <- range(ARN$F)

  if (relative == TRUE) {
    ATRmax <- max(ATR$F) / 100
    ARNmax <- max(ARN$F) / 100
  } else {
    ATRmax <- 1
    ARNmax <- 1
  }

  # Plot Arai
  if (show.plot == TRUE) {

    xlim <- range(ATR$F/ATRmax)
    ylim <- range(ARN$F/ARNmax)
    if (relative == TRUE) {
      xlab <- "% ATR"
      ylab <- "% ARN"
    } else {
      xlab <- "ATR"
      ylab <- "ARN"
    }

    plot(y = ARN$F/ARNmax, x = ATR$F/ATRmax, type='o', col = pt.col.res, pch = 20, xlim = xlim, ylim = ylim, xlab = xlab, ylab = ylab)
    if (show.step.value == TRUE) {
      text(y= ARN$F, x=ATR$F, ARN$step.value, adj = 0)
    }
  }

  # stat with all points
  Xmoy <- mean(ATR$F)
  Ymoy <- mean(ARN$F)
  n <- length(ATR$F)

  # calcul with a population, not with a sample
  covYX <- cov(ATR$F, ARN$F) * (n-1)/n
  S2x <- var(ATR$F) * (n-1)/n
  S2y <- var(ARN$F) * (n-1)/n

  a_tot <- -sqrt(S2y/S2x)
  b_tot <- Ymoy - a_tot*Xmoy
  JTRM <- (-b_tot/a_tot)

  # loop pTRM check
  ATRP <- mesures[which(substr(mesures$step, nchar(mesures$step)-1, nchar(mesures$step)-1 ) == P.mark),]
  for (i in 1:length(ATRP$step) ) {
    ATRP$step.name[i] <- substr(ATRP$step[i], nchar(ATRP$step[i]), nchar(ATRP$step[i]) )
  }

  # Plot loop
  TP <- NULL
  ATP <- NULL
  for (i in ATRP$step.name) {
    iATRV <- which(ATRV$step.name == i)
    iATRP <- which(ATRP$step.name == i)

    atrXV <- ATRV[iATRV, ]$X *aim.coef
    atrXP <- ATRP[iATRP, ]$X *aim.coef

    atrYV <- ATRV[iATRV, ]$Y *aim.coef
    atrYP <- ATRP[iATRP, ]$Y *aim.coef

    atrZV <- ATRV[iATRV, ]$Z *aim.coef
    atrZP <- ATRP[iATRP, ]$Z *aim.coef


    TP$X <-(atrXP - atrXV)/2
    TP$Y <- (atrYP - atrYV)/2
    TP$Z <- (atrZP - atrZV)/2
    TP <- to.polar(TP$X, TP$Y, TP$Z)
    TP$step.value <- ATRP[iATRP, ]$step.value
    TP$step.name <- i

    ATP <-rbind.data.frame(ATP, TP, stringsAsFactors = FALSE)

    arnLoop.name <- ARN[which(ARN$step.name == i), ]
    arnLoop.value <- ARN[which(ARN$step.value == TP$step.value), ]
    atrLoop <- ATR[which(ATR$step.name == i), ]
    if (show.plot == TRUE) {
      points( x= TP$F/ATRmax, y= arnLoop.value$F/ARNmax, col = loop.col, pch = 4)
      lines( x= c(atrLoop$F/ATRmax, TP$F/ATRmax, TP$F/ATRmax), y = c(arnLoop.name$F/ARNmax, arnLoop.name$F/ARNmax, arnLoop.value$F/ARNmax), col = loop.col)
    }
  }

  res<- NULL
  res$ARN <- ARN
  res$ATR <- ATR
  res$ATP <- ATP

  # Statistic

  n <- 0
  Crm <- CrmMax <- 0
  tabX <- NULL
  tabY <- NULL
  for (i in 1:length(ARN$step.name))  {
    if ((ARN$step.value[i] >= begin.step.value) && (ARN$step.value[i] <= end.step.value)) {
      tabX <- c( tabX, ATR[which(ATR$step.value == ARN$step.value[i]),]$F  /ATRmax)
      tabY <- c( tabY, ARN$F[i]  /ARNmax)
      n <- n+1

      if (n==1) {
        incl_ET1.rad <- ARN$I[i]/180*pi
      } else {
        Crm <- ( sin( incl_ET1.rad -ARN$I[i]/180*pi) /sin(incl_ET1.rad -pi/2))*ARN$F[i]
        if (abs(Crm)>CrmMax) {
          CrmMax <- abs(Crm)
        }
      }


    }
  }
  CrmMax <- CrmMax/ ( tabX[length(tabX)] - tabX[1])*100

  if (n >= 2) {
    Xmoy <- mean(tabX)
    Ymoy <- mean(tabY)

    # calcul sur une population, pas sur un echantillon
    covYX <- cov(tabY, tabX) * (n-1)/n
    S2x <- var(tabX) * (n-1)/n
    S2y <- var(tabY) * (n-1)/n

    a_aff <- -sqrt(S2y/S2x)
    b_aff <- Ymoy - a_aff*Xmoy

    abline(b_aff, a_aff, col= "red")

    # Calcul du coef. de corrélation linéaire de la droite entre les 2 étapes
    R <- cor(tabY, tabX)

    # calcul de la variance eq n°3 : Prévost et Al. 1985 DOI: 10.1029/JB090iB12p10417

    s <-  covYX/sqrt(S2x*S2y)
    if (s < -1) {
      warning('Error on sigma Prévost =', s);
      s<- abs(s);

    } else {
      s <- sqrt( 2+(2*s));
      sigmaPrevost85 <- s/sqrt(n-2);
      # Calcul de sigma pour Coe 1978 : DOI: 10.1029/JB083iB04p01740
      sigmaCoe <- ((2*S2y)-2*a_aff*covYX) /((n-2)*S2x)
      sigmaCoe <- sqrt(sigmaCoe)

    }
    # fraction of pTRM used
    fs1s2 <- function(tx,ty, i1, i2) {
      (a_aff*tabX[i1] + b_aff  + tabY[i1] )/2 - (a_aff*tabX[i2] + b_aff + tabY[i2])/2
    }

    nStep <- length(tabX)
    # Coe et Al. 1978)
    fCoe78 <- fs1s2(tabX*ATRmax, tabY*ATRmax, 1, nStep) /b_aff

    g<-0
    for (i in 2:length(tabX)) {
      g <- g + fs1s2(tabX*ATRmax, tabY*ATRmax, i-1, i)^2
    }
    g <- 1 - g/(fCoe78*b_aff)^2

    qCoe78 <- abs(a_aff)*fCoe78*g/sigmaCoe;

    qPrevost85 <- fCoe78*g/sigmaCoe

    SigFe <- sigmaCoe*TH
    Fe <- -a_aff*TH

    if (verbose == TRUE) {
      Commentaire <- paste0('Coef. Corr. lin. R= ', format(R, digits = 3), ' pour ', n,' points')
      Commentaire <- c(Commentaire, paste0(' Sigmab ( Coe 1978) = ',format(sigmaCoe, digits = 3)) )
      Commentaire <- c(Commentaire, paste0('with lab field : ', TH, ' µT => Fe= ', format(Fe, digits = 3), ' ± ', format(SigFe, digits = 3), ' µT (Coe et Al. 1978)') )
      Commentaire <- c(Commentaire, paste0('f = ', format(fCoe78*100, digits = 4), '% (Coe et Al. 1978)') )
      Commentaire <- c(Commentaire, paste0('q = ', format(qCoe78, digits = 3)) )
      Commentaire <- c(Commentaire, paste0('g  = ', format(g, digits = 3)) )
      Commentaire <- c(Commentaire, paste0('Crm  = ', format(CrmMax, digits = 4), ' (Coe 1984)') )


      Commentaire <- c(Commentaire, paste0('q (Prévost 85) = ', format(qPrevost85, digits = 3)) )
      Commentaire <- c(Commentaire, paste0('sigma (Prévost 85) = ', format(sigmaPrevost85, digits = 3)) )
      print(Commentaire)
    }

    # Cooling rate
    # slow step
    ATRL <- mes.sel[which(substr(mes.sel$step, nchar(mes.sel$step)-1, nchar(mes.sel$step)-1 ) == L.mark),]
    for (i in 1:length(ATRL$step) ) {
      ATRL$step.name[i] <- substr(ATRL$step[i], nchar(ATRL$step[i]), nchar(ATRL$step[i]) )
    }
    ATRQ <- mes.sel[which(substr(mes.sel$step, nchar(mes.sel$step)-1, nchar(mes.sel$step)-1 ) == Q.mark),]
    for (i in 1:length(ATRQ$step) ) {
      ATRQ$step.name[i] <- substr(ATRQ$step[i], nchar(ATRQ$step[i]), nchar(ATRQ$step[i]) )
    }

    # quick step
    for (i in ATRL$step.name) {
      iARN <- which(ARN$step.name == i)
      iATR <- which(ATR$step.name == i)

      iATRL <- which(ATRL$step.name == i)
      atrXL <- ATRL[iATRL, ]$X*aim.coef  - ARN[iARN, ]$X
      atrYL <- ATRL[iATRL, ]$Y *aim.coef - ARN[iARN, ]$Y
      atrZL <- ATRL[iATRL, ]$Z*aim.coef  - ARN[iARN, ]$Z
      trL <- to.polar(atrXL, atrYL, atrZL)
      rateL <- (trL$F - ATR[iATR, ]$F) / (ATR[iATR, ]$F)
      FeL <- Fe*(1-rateL)
      SigFeL <- SigFe * (1-rateL)

      iATRQ <- which(ATRQ$step.name == i)
      atrXQ <- ATRQ[iATRQ, ]$X *aim.coef - ARN[iARN, ]$X
      atrYQ <- ATRQ[iATRQ, ]$Y *aim.coef - ARN[iARN, ]$Y
      atrZQ <- ATRQ[iATRQ, ]$Z *aim.coef - ARN[iARN, ]$Z
      trQ <- to.polar(atrXQ, atrYQ, atrZQ)
      rateQ <- (trQ$F - ATR[iATR, ]$F) / (ATR[iATR, ]$F)

      if (verbose == TRUE) {
        Commentaire <- paste0('Speed Rate with slow step: ', format(rateL*100, digits = 4), ' %')
        Commentaire <- c( Commentaire, paste0('Derive speed Rate with quick step: ', format(rateQ*100, digits = 4), ' %') )
        Commentaire <- c( Commentaire, paste0('Fe Lent= ', format(FeL, digits = 3), ' ± ', format(SigFeL, digits = 3), ' µT') )

        print(Commentaire)

      }
    }

    res$stat <- data.frame(Fe, SigFe, FeL, SigFeL, rateL, rateQ, sigmaCoe, JTRM, fCoe78, qCoe78, g, CrmMax, qPrevost85, sigmaPrevost85)
  }
  return(res)
}


