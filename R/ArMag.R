#  Licence ------------------------------------------------------------------
#
#  Copyright or © or Copr. CNRS	2014 - 2018
#
#
# This package is under
# The “Creative Commons Attribution-ShareAlike International License” version 4.0

#  A copy of the License is available at
#  http://www.r-project.org/Licenses/
#  ------------------------------------------------------------------

# Version 2019-09-27
#

#' @author "Philippe DUFRESNE"

# Equation du 3 degrées
# résolution du 3 eme degres pour calcul mcFadden importé de ARMAG
#  modifier suivant livre photocopié
#  Fonction interne
EQUATION_DEGRE_3 <- function(A1, A2, A3)
{

  Q <- (3*A2-A1*A1)/9
  R <- (9*A1*A2-27*A3-2*A1*A1*A1)/54
  D <- Q*Q*Q+R*R;

 #  Discriminant équation du 3eme degre
  if (D>0) {
    S <- exp( log( R + sqrt(D)  )/3)
    if ( R>sqrt(D))
      T1 <- exp(log( R - sqrt(D)  )/3)
    else
      T1 <- -exp(log( -R + sqrt(D)  )/3)

    PZ1 <- S+ T1 -(A1/3);
     # Les deux autres solutions sont complexes
    PZ2 <- 0;
    PZ3 <- 0;
  }
  else if (D==0) {
    S <- exp( log(-Q/2)/3)
    PZ1 <- 2*S-(A1/3)
    PZ2 <- -S-(A1/3)
    PZ3 <- -S-(A1/3)
  }
  else if (D<0) {
    TETA <- acos( R/sqrt(-Q*Q*Q) )/3;
    U <- 2*sqrt(-Q);
    PZ1 <- U*cos(TETA)-(A1/3);
    PZ2 <- U*cos(TETA + 2/3*pi)-(A1/3);
    PZ3 <- U*cos(TETA + 4/3*pi)-(A1/3);
  }

  rslt <- NULL
  rslt$PZ1 <- PZ1
  rslt$PZ2 <- PZ2
  rslt$PZ3 <- PZ3

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
   rslt$PX1 <- PX1
   rslt$PX2 <- PX2
   rslt$PX3 <- PX3
   rslt$PX4 <- PX4

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

#' Statistique de mcFadden sur l'inclinaison seul à partir des coordonées XYZ
#' @seealso stat.mcFadden, stat.fisher
#' @export
stat.mcFadden.XYZ <- function(TabX, TabY, TabZ)
{
  N <- length(TabX)
  VP <- Calcul_Vecteur_Polaire(TabX, TabY, TabZ)

  stat.mcFadden(inc, dec)

}

#' Statistique de mcFadden sur l'inclinaison seul à partir des coordonées I et D
#' @param Data liste des inclinaisons en degrésou une data.frame avec les variables $I et $D
#' @param des liste des déclinaisons en degrés
#' @param inc.absolue calcul avec la valeur absolue des inclinaisons
#' @return  en degrés, un data.frame "n", "imoy.McFadden", "imoy.McElhinny", "a95.mcFad", "a95.eqFish", "Kb", "Kssb", "imin", "imax", "dmin", "dmax"
#' @seealso stat.mcFadden.XYZ, stat.fisher
#' @references https://doi.org/10.1111/j.1365-246X.1990.tb05683.x
#' @export
stat.mcFadden <- function(Data, dec = NULL, inc.absolue = TRUE)
{

  if (is.null(dec)) {
    inc <- Data$I
    dec <- Data$D
  } else {
    inc <- Data
    dec <- dec
  }

  if (inc.absolue)
    inc <- abs(inc)

  # Calcul des sommes
  N <- length(inc)
  if ((N<=2) || (N != length(dec)) ) return()


  Imin <- min(inc)
  Imax <- max(inc)
  Dmin <- min(dec)
  Dmax <- max(dec)
  Imoy <- 0

  # Passage en radian pour les calculs
  i.rad <- inc/180*pi
  d.rad <- dec/180*pi

  A <- sum(sin(i.rad))
  B <- sum(cos(i.rad))


  P <- 2*(N+2*B)/A
  Q <- -6
  R <- 2*(N-2*B)/A
  S <- 1

  E4 <- EQUATION_DEGRE_4(P, Q, R, S)

  K <- 0
  KJ <- N/(2*(N-A*sin(2*atan(E4$PX1))-B*cos(2*atan(E4$PX1))));

  if (KJ>K) {
    Imoy <- 2*atan(E4$PX1)
    K <- KJ
  }

  KJ <- N/(2*(N-A*sin(2*atan(E4$PX2))-B*cos(2*atan(E4$PX2))));

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

  I.McElhinny <- Imoy      # modif effectuée le 2017/01/25
  Imoy <- Imoy + (SFad/CFad)

  a95mcFadden <- acos(1-(SFad/CFad/4)^2-(qf(.95, 1, N-1)*(N-CFad)/(CFad*(N-1))))
  Kssb <-  (N-1)/2/(N-CFad);

  # retour en degrés
  Data <- c( n= N, imoy.McFadden = Imoy/pi*180, imoy.McElhinny = I.McElhinny/pi*180,
            a95.mcFad = a95mcFadden/pi*180, a95.eqFish = 2.4477/sqrt(N*K)/pi*180,
            Kb = K, Kssb = Kssb,
            imin = Imin, imax = Imax, dmin = Dmin, dmax = Dmax
           )

  col.names <- c("n", "imoy.McFadden", "imoy.McElhinny", "a95.mcFad", "a95.eqFish", "Kb", "Kssb", "imin", "imax", "dmin", "dmax")
  return(as.data.frame(t(Data), col.names = col.names, stringsAsFactors = FALSE))


  }

#' Statistique de Fisher
#' modifié, pas de pondération calcul de A95 Vrai sans simplification tel que fisher 1953
#' @param inc liste des inclinaisons en degrés
#' @param des liste des déclinaisons en degrés
#' @param aim liste des aimantations, facultatif
#' @param pfish pourcentage de confiance
#' @param inc.absolue calcul avec la valeur absolue des inclinaisons
#' @return  en degrés
#' @seealso stat.mcFadden
#' @keywords fisher
#' @export
stat.fisher <- function (inc, dec, aim=NA, pfish = 0.95, inc.absolue = TRUE)
{
  n <- length(inc)
  if (length(dec) != n) {
    return("length (dec) diff length(inc)")
  }
  if (is.na(aim) || length(aim) != n) {
   # message("length(aim) diff length(inc)")
    aim <-  rep(1, n)
  }

  imin <- min(inc)
  imax <- max(inc)
  dmin <- min(dec)
  dmax <- max(dec)

  if (inc.absolue == TRUE)
    inc <- abs(inc)

  # Passage en radian pour les calculs
  i.rad <- inc/180*pi
  d.rad <- dec/180*pi

  sx <- sum(aim*cos(i.rad)*sin(d.rad))
  sy <- sum(aim*cos(i.rad)*cos(d.rad))
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

  delta <- log(1+(1-pfish) * (exp(2*KF)-1)) /KF
  delta <- acos(delta-1)

  return( )

  # retour en degrés
  Data <- c(n = n, imoy =imoy/pi*180, dmoy = dmoy/pi*180, alpha=a95/pi*180, pfish = pfish, delta = delta/pi*180, KF = KF,
            imin = imin, imax = imax, dmin = dmin, dmax = dmax)

  col.names <- c("n", "imoy", "dmoy", "a95", "pFish", "Delta", "KF", "imin", "imax", "dmin", "dmax")
  return(as.data.frame(t(Data), col.names = col.names, stringsAsFactors = FALSE))

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
{ # angle en degré

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

#' Conversion de Degres minute second en Degré Décimal
#' @param degre degrés entier
#' @param minute minute entière
#' @param second seconde en décimal
#' @seealso DD.to.DMS
#' @export
DMS.to.DD <- function(degre, minute, second)
{
  if (degre>0) {
    return(degre + minute/60 + second/3600)
  } else {
    return(degre - minute/60 - second/3600)
  }

}

#' Conversion de Degré Décimal en Degres minute second
#' @param degre degrés décimal
#' @seealso DMS.to.DD
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
#' @param d déclinaison en degrés
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
#' @param d déclinaison en degrés
#' @return déclinaison entre 0 et 360°
#' @export
D.pal <- function(d)
{
  if ( d>270)
    return(d-360)
}

## Calcul vecteur cartesien ----
#' Valeur de Y pour I et D
#' @param inc liste des inclinaisons en degré
#' @param des liste des déclinaisons en degré
#' @export
calcul.vecteur.cartesien.Y <- function(inc, dec, aim =1)
{
  aim*cos(inc/180*pi)*sin(dec/180*pi)
}

#' Valeur de X pour I et D
#' @param inc liste des inclinaisons en degré
#' @param des liste des déclinaisons en degré
#' @export
calcul.vecteur.cartesien.X <- function(inc, dec, aim=1)
{
  aim*cos(inc/180*pi)*cos(dec/180*pi)
}

#' Valeur de Z pour I et D
#' @param inc liste des inclinaisons en degré
#' @param des liste des déclinaisons en degré
#' @export
calcul.vecteur.cartesien.Z <- function(inc, dec, aim=1)
{
  aim*sin(inc/180*pi)
}

#' Valeur de X, Y et Z pour I et D
#' @param inc liste des inclinaisons en degré
#' @param des liste des déclinaisons en degré
#' @seealso cartesien , polaire
#' @export
calcul.vecteur.cartesien <- function(inc, dec, aim=1)
{
  res <- NULL
  resT <- NULL
  for (i in 1:length(inc)) {
    res$X <- aim*cos(inc[i]/180*pi)*cos(dec[i]/180*pi)
    res$Y <- aim*cos(inc[i]/180*pi)*sin(dec[i]/180*pi)
    res$Z <- aim*sin(inc[i]/180*pi)
    resT <- c(resT, res)
  }
  return( resT )
}

#' Valeur de X, Y et Z pour I et D
#' @param inc liste des inclinaisons en degrés
#' @param des liste des déclinaisons en degrés
#' @seealso calcul.vecteur.cartesien , polaire
#' @export
cartesien <- function(inc, dec, aim=1)
{
  res <- NULL
  resT <- NULL
  for (i in 1:length(inc)) {
    res$X <- aim*cos(inc[i]/180*pi)*cos(dec[i]/180*pi)
    res$Y <- aim*cos(inc[i]/180*pi)*sin(dec[i]/180*pi)
    res$Z <- aim*sin(inc[i]/180*pi)
    resT <- c(resT, res)
  }
  return( resT )
}


# Transformation de type de coordonnée ----
#' Coordonnées polaires pour X, Y et Z
#' @return angle en degré
#' @seealso polaire
#' @export
calcul.vecteur.polaire <- function(X, Y, Z)
{ #en A/m
  res <- NULL
  resT <- NULL
  for (i in 1:length(X)) {
    res$F <- sqrt(X[i]*X[i]+Y[i]*Y[i]+Z[i]*Z[i])
    res$I <- n_arcsin(Z[i]/res$F)/pi*180
    res$D <- angleD(X[i],Y[i])/pi*180
    resT <- c(resT, res)
  }
  return( resT )
}

# Transformation de type de coordonnée ----
#' Coordonnées polaires pour X, Y et Z
#' @return degrés
#' @seealso calcul.vecteur.polaire , cartesien
#' @export
polaire <- function(X, Y, Z)
{ #en A/m
  res <- NULL
  resT <- NULL
  for (i in 1:length(X)) {
    res$F <- sqrt(X[i]*X[i] + Y[i]*Y[i] + Z[i]*Z[i])
    res$I <- n_arcsin(Z[i]/res$F) /pi*180
    res$D <- angleD(X[i],Y[i]) /pi*180
    resT <- c(resT, res)
  }
  return( resT )
}

# Calcul vecteur polaire ----
#' inclinaison pour X, Y et Z
#' @return angle en degré
#' @export
calcul.vecteur.polaire.I <- function(X, Y, Z)
{ #en A/m
  res <- NULL
  for (i in 1:length(X)) {
    A <- sqrt(X[i]*X[i]+Y[i]*Y[i]+Z[i]*Z[i])
    res <- c(res, n_arcsin(Z[i]/A)/pi*180)
  }
  return(res)
}

#' Déclinaison pour X, Y et Z
#' #' @return angle en degré
#' @export
calcul.vecteur.polaire.D <- function(X, Y, Z)
{ #en A/ms
  res <- NULL
  for (i in 1:length(X)) {
     res <- c(res, angleD(X[i],Y[i])/pi*180)
  }
  return(res)
}

#' champs pour X, Y et Z
#' @export
calcul.vecteur.polaire.F = function(X, Y, Z)
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
#' @param  radlab écrit les labels sous forme d'étoile
#' @param  label.pos séquence de valeur à afficher en I et D, attention format particulier !!

#' # Pour choisir les graduations
#' label.pos = NULL
#' label.pos$I = seq(0, 90, by=20)
#' label.pos$D = seq(0, 90, by=10)
#' et pour paleomag : lab.pos$D = c(seq(270, 350, by=10), seq(0, 90, by=10))
#' Label.pos : doit être dans l'étendu des dex.min, dec.max, inc.min, inc.max
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
#' @param inclinaisons, declinaisons : les listes des données d'inclinaison et de déclinaison
#' @param pt.names : Correspond à la liste des noms des points. Laissée vide n'affiche rien. Si on met pt.names = "", cela affiche les noms
#' @param label.pos : Séquence de valeur à afficher en I et D, voir fonction lambert.ID.grid
#' @param point.symbols : défini la forme de points, correspond exactement au pch de la fonction points()
#' @param pch : permet de changer la forme du symbole
#' @param show.grid :  permet d'afficher une grille en toile d'araigné sur le fond. Mettre à FALSE, si on superpose des diagrammes
#' @param show.grid.labels = 10 : permet de changer échelle des graduations
#' @param new : permet d'initialiser la sortie graphique. Mettre à FALSE, si on superpose des diagrammes.
#' @param inc.lim : permet de restraire l'affichage sur une étendue d'inclinaison. Ex: inc.lim = c(45, 90). Laissée à NULL, le diagramme s'addapte aux données
#' @param dec.min = -90, dec.max = 180 : permet de restraire l'étendue en déclianison
#' @export
lambert <- function (data , pt.names = NULL, labels = NA, label.pos = NULL,
                             radlab = FALSE, start = 0, clockwise = TRUE,
                             label.prop = 1.1, main = "", xlab = "", ylab = "", line.col = par("fg"),
                             lty = par("lty"), lwd = par("lwd"), mar = c(2, 2, 3, 2),
                             show.grid = TRUE, show.grid.labels = 10, show.radial.grid = TRUE,
                             grid.col = "gray", grid.bg = "transparent",
                             grid.left = FALSE, grid.unit = NULL, point.symbols = 1, point.col = par("fg"), bg = point.col,
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
                               grid.left = grid.left, grid.unit = grid.unit, point.symbols = point.symbols, point.col = point.col, bg = point.col,
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
                                 grid.left = FALSE, grid.unit = NULL, point.symbols = 1, point.col = par("fg"), bg = point.col,
                                 inc.lim = NULL, radial.labels = NULL,
                                 boxed.radial = TRUE, poly.col = NA,
                                 dec.min = -90, dec.max = 270, new = TRUE, pch = 21, ...)
{

  if (is.null(inc.lim))
    inc.lim <- range(abs(inc))



  #lambert.ID(inclinaisons = inclinaisons, declinaisons = declinaisons, pt.names = pt.names, inc.lim = inc.lim,
  #           dec.min = dmin, dec.max = dmax, main = main, label.pos = label.pos, point.col = point.col, bg = bg, show.grid = show.grid, new = new)
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
                        pt.names = pt.names, point.col = point.col, bg = bg, pch = pch, new = new)


}

#' Place des points I et D dans un repère Lambert et dessine le symbole
#' @export
lambert.ID.position <- function (data, declinaisons = NULL, position = "P",  pt.names = NA, labels = NA, label.pos = NULL,
                        radlab = FALSE, start = 0, clockwise = TRUE,
                        label.prop = 1.1, main = "Position auto", xlab = "", ylab = "", line.col = par("fg"),
                        lty = par("lty"), lwd = par("lwd"), mar = c(2, 2, 3, 2),
                        show.grid = TRUE, show.grid.labels = 10, show.radial.grid = TRUE,
                        grid.col = "gray", grid.bg = "transparent",
                        grid.left = FALSE, grid.unit = NULL, point.symbols = 1, point.col = par("fg"), bg = point.col,
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
             inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, point.col = point.col,  pch = 23, new = ne)
    ne <- FALSE
    show.grid <- FALSE
    main <- ""
  }

  if (length(which(position=="C")) > 0) {
    lambert.ID(inclinaisons[which(position=="C")], declinaisons[which(position=="C")],
             pt.names = pt.names[which(position=="C")], label.pos = lab.pos, main = main, show.grid = show.grid,
             inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, point.col = point.col,  pch = 22, new = ne)
    ne <- FALSE
    show.grid <- FALSE
    main <- ""
  }


  if (length(which(position=="D")) > 0) {
    lambert.ID(inclinaisons[which(position=="D")], declinaisons[which(position=="D")],
             pt.names = pt.names[which(position=="D")], label.pos = lab.pos, main = main, show.grid = show.grid,
             inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, point.col = point.col, pch = 24, new = ne)
  }
}

#' Place des points I et D dans un repère Lambert
#' @export
lambert.ID.point <- function (inc, dec , pt.names = NA, labels = NA, label.pos = NULL,
                             start = 0, clockwise = TRUE,
                             lty = par("lty"),  mar = c(2, 2, 3, 2),
                             point.col = par("fg"), bg = point.col,
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

  if (length(point.col) < nbpoints)
    point.col <- rep(point.col, length.out = nbpoints)


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

  if (!is.null(point.col))
    for (i in 1: nbpoints ) {
      if(point.col[i] != "transparent") {
        if (inc[i]<0) {
          bg.col <- gray(0.95)
          points(xpos[i], ypos[i], pch = pch[i],  col = point.col[i], bg = bg.col,...)
        }
        else
          points(xpos[i], ypos[i], pch = pch[i],  col = point.col[i], bg = point.col[i],  ...)
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
                              main = "", point.col = par("fg"), bg = point.col, label.pos= NULL, show.grid = TRUE, new = TRUE)
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

  lambert.ID.point(inc, dec, inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, pt.names = pt.names, point.col = point.col, bg = bg, new = FALSE)



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
                                   grid.left = FALSE, grid.unit = NULL, point.symbols = 1, point.col = "blue3", bg = point.col,
                                   inc.lim = NULL, radial.labels = NULL,
                                   boxed.radial = TRUE, poly.col = NA, add = FALSE,
                                   dec.min = -90, dec.max = 270, new = TRUE, pch = 21, ...)
{
  I <- calcul.vecteur.polaire.I(X, Y, Z)
  D <- calcul.vecteur.polaire.D(X, Y, Z)
  lambert.ID(I, D, pt.names = pt.names, label.pos = label.pos, main = main, show.grid = show.grid,
             inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, point.col = point.col,  pch = 23, new = new)

}

#' Place des points X, Y et Z dans un repère Lambert avec un data.frame
#' @export
lambert.XYZ.specimen <- function( Data, Y , Z, pt.names = "", labels = NA, label.pos = NULL,
                         radlab = FALSE, start = 0, clockwise = TRUE,
                         label.prop = 1.1, main = NULL, xlab = "", ylab = "", line.col = "blue3",
                         lty = par("lty"), lwd = par("lwd"), mar = c(2, 2, 3, 2),
                         show.grid = TRUE, show.grid.labels = 10, show.radial.grid = TRUE,
                         grid.col = "lightgray", grid.bg = "transparent",
                         grid.left = FALSE, grid.unit = NULL, point.symbols = 1, point.col = "blue", bg = point.col,
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
      eta <- Data$etape
    if (is.null(main))
      main <- as.character(Data$name[1])
  } else {
    X <- Data
    if (!is.null(pt.names))
      eta <- pt.names
  }
  I <- calcul.vecteur.polaire.I(X, Y, Z)
  D <- calcul.vecteur.polaire.D(X, Y, Z)

  if (is.null(label.pos)) {
    label.pos$I = c(90)
    label.pos$D = seq(-90, 270, by=90)
  }
  lambert.ID(I, D, pt.names = pt.names, label.pos = label.pos, main = main, show.grid = show.grid, type = "l",
                  grid.col = grid.col,
             inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, line.col = line.col, point.col = point.col, bg= par("fg"), new = new)
  lambert.ID(I, D, pt.names = "", label.pos = label.pos, main = "", show.grid = FALSE,
                  inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, line.col = line.col, point.col = point.col,  pch = 21, new = FALSE)
}

#' Place des points I et D dans un repère Lambert avec un data.frame
#' @export
lambert.ID.specimen <- function( Data, D , pt.names = "", labels = NA, label.pos = NULL,
                                       radlab = FALSE, start = 0, clockwise = TRUE,
                                       label.prop = 1.1, main = NULL, xlab = "", ylab = "", line.col = "blue3",
                                       lty = par("lty"), lwd = par("lwd"), mar = c(2, 2, 3, 2),
                                       show.grid = TRUE, show.grid.labels = 10, show.radial.grid = TRUE,
                                       grid.col = "lightgray", grid.bg = "transparent",
                                       grid.left = FALSE, grid.unit = NULL, point.symbols = 1, point.col = "blue3", bg = point.col,
                                       inc.lim = c(0, 90), radial.labels = NULL,
                                       boxed.radial = TRUE, poly.col = NA, add = FALSE,
                                       dec.min = -90, dec.max = 270, new = TRUE, pch = 21, ...)
{
  eta <- NULL
  if(is.data.frame(Data)) {
    I <- Data$I
    D <- Data$D
    if (is.null(pt.names))
      eta <- Data$etape
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
                  inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, line.col = line.col, point.col = point.col, bg= par("fg"), new = new)
  lambert.ID(I, D, pt.names = "", label.pos = label.pos, main = "", show.grid = FALSE,
                  inc.lim = inc.lim, dec.min = dec.min, dec.max = dec.max, line.col = line.col, point.col = point.col,  pch = 21, new = FALSE)
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
      eta <- Data$etape
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

    ax1 <- axis(1, pos = 0, cex.axis = 0.8, col = "darkgray")
    ax2 <- axis(2, pos = 0, cex.axis = 0.8, col = "darkgray") # Ordonnées

    text(0, ax2[length(ax2)], "+X", col = "gray5", adj = c(-.5, 1), cex = par("cex.lab"))
    text(0, ax2[1], "+Z", col = "gray5", adj = c(-.5, 0), cex = par("cex.lab"))
    text( ax1[length(ax1)], 0, "+Y", col = "gray5", adj = c(1, -.5), cex = par("cex.lab"))

  } else {
    lines(Y, X, type = "o", pch = 21, col = pt.col[1], bg = adjustcolor( pt.col[1], alpha.f = 0.8), ...)
  }

  lines(Y, -Z, type = "o", pch = 21, col = pt.col[2], bg = adjustcolor( pt.col[2], alpha.f = 0.05), ...)

  text(jitter(Y, 5, amount = 0), jitter(X, 5, amount = 0), eta)

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
      eta <- Data$etape
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

    text(jitter(-X, 5, amount = 0), jitter(Y, 5, amount = 0), eta)

    if (!is.null(legend.pos))
      legend(legend.pos, c("(-X, Y)", "(-X, Z)"), pch = c(19, 21), col = pt.col, bg = c(par("bg"), adjustcolor( pt.col[1], alpha.f = 0.8), adjustcolor( pt.col, alpha.f = 0.05)),
             box.col = par("bg"), title = "")

    ax1 <- axis(1, pos = 0,  col = "darkgray") # cex.axis = 0.8,
    ax2 <- axis(2, pos = 0,  col = "darkgray") # Ordonnées

    text(0, ax2[length(ax2)], "+Y", col = "gray5", adj = c(-.5, 1), cex = par("cex.lab"))
    text(0, ax2[1], "+Z", col = "gray5", adj = c(-.5, 0), cex = par("cex.lab"))
    text( ax1[length(ax1)] , 0, "-X", col = "gray5", adj = c(1, -.5), cex = par("cex.lab"))


  } else {
    lines(-X, Y, type = "o", pch = 21, col = pt.col[1], bg = adjustcolor( pt.col[1], alpha.f = 0.8), ...)
  }

  lines(-X, -Z, type = "o", pch = 21, col = pt.col[2], bg = adjustcolor( pt.col[2], alpha.f = 0.05))


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

  X <- calcul.vecteur.cartesien.X(inc, dec, aim)
  Y <- calcul.vecteur.cartesien.Y(inc, dec, aim)
  Z <- calcul.vecteur.cartesien.Z(inc, dec, aim)

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

    tmp <- calcul.vecteur.cartesien(t.i, t.d, aim = aim[i])
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
repliement.tranche <- function (data, dec, aim = 1,  name = NULL, number = NULL, inc.critique = 90)
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
  res <- data.frame(number = as.numeric(res$number), name = as.character(res$name), I = as.numeric(res$I), D = as.numeric(res$D), F = as.numeric(res$F),
                    X = as.numeric(res$X), Y =as.numeric(res$Y), Z =as.numeric(res$Z), position = res$positon, stringsAsFactors = FALSE)

  return(res)
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
      g<-c(g,i )

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
        lI <- calcul.vecteur.polaire.I(as.numeric(lX), as.numeric(lY), as.numeric(lZ))
        lD <- calcul.vecteur.polaire.D(as.numeric(lX), as.numeric(lY), as.numeric(lZ))
        lF <- calcul.vecteur.polaire.F(as.numeric(lX), as.numeric(lY), as.numeric(lZ))
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
      lmes <- data.frame(number = i, name = lname[i], etape = tEtap, etape.value = as.numeric(tEtap.val),
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

#' Fonction de création de fichier pour les magnétomètres ----
#' @param encoding mettre "macroman" pour les mesures au molspin
#' @export
genere.AMD <- function(file.AMD = "fichier.AMD", list.ech, shape = "Cyl" , encoding = "macroman")
{
  entete<- c( "Spinner_Molspin 2008" ,
              "Commune : Laval",
              "Site : St Pierre-le-Potier",
              "Latitude  :   0°  0'  0\" ",
              "Longitude :   0°  0'  0\" IGRF:+00.0",
              "Prélèvements sur matériaux déplacés",
              "Type de carottage : à plat",
              "Date de création : 27/05/2019", "","")

  txt.mesures <- entete
  for (i in 1:length(list.ech)) {
    txt.mesures <- c( txt.mesures,
                      paste("Id:", format(list.ech[i], width = 13), "in:000.0 az:000.0 Tet:000.0 Psy:000.0 v:12.27 com:TH50.0µT ", shape, sep = ""),
                      "Repère:",
                      "CompDes:  T1:0000T+  T2:0000T+  T3:0000T-  T4:0000T-",
                      "")
  }


  # Ecriture du fichier

  filCon <- file(file.AM, encoding = encoding)
  writeLines(txt.mesures, filCon)
  close(filCon)
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

# Fonction Anisotropie ----

#' Calcul la matrice d'anisotropie symetrisée et normalisé pour un spécimen
#' qui sert à la correction
#' @param mesures data.frame contenant les mesures
#' @param etape.value typiquement la température des mesures d'anisotropie
#' @param etape.sigle sigle indiquant les étapes. les noms peuvent changer, mais pas l'ordre
#' @param volume la valeur du volume du spécimen
#' @param TH la valeur du champ appliqué
#' @return un data.frame avec les colonnes "L1", "L1.Inc", "L1.Dec", "L2", "L2.Inc", "L2.Dec", "L3", "L3.Inc", "L3.Dec", "F13", "F12", "F23"
#' @export
anisotropie.matrix.symetric <- function(mesures, etape.value, etape.sigle = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"), ...)
{

  ani.etape <- trimws(paste(as.character(etape.value), etape.sigle, sep = "") )

  selec <- NULL
  for (i in 1:length(ani.etape)) {
    selec <- c( selec, which(trimws(mesures$etape) == trimws(ani.etape[i])) )
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
#' @param etape.value typiquement la température des mesures d'anisotropie
#' @param etape.sigle sigle indiquant les étapes. les noms peuvent changer, mais pas l'ordre
#' @param volume la valeur du volume du spécimen
#' @param TH la valeur du champ appliqué
#' @return un data.frame avec les colonnes "L1", "L1.Inc", "L1.Dec", "L2", "L2.Inc", "L2.Dec", "L3", "L3.Inc", "L3.Dec", "F13", "F12", "F23"
#' @export
anisotropie.eigen <- function(mesures, etape.value, etape.sigle = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"), volume = 1, TH = 1,...)
{

  #etape.sigle <- c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB")
  ani.etape <- trimws(paste(as.character(etape.value), etape.sigle, sep = "") )

  selec <- NULL
  for (i in 1:length(ani.etape)) {
    selec <- c( selec, which(trimws(mesures$etape) == trimws(ani.etape[i])) )
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
  v1$I<-calcul.vecteur.polaire.I(v$vectors[1, 1], v$vectors[2, 1], v$vectors[3 ,1])
  v1$D<-calcul.vecteur.polaire.D(v$vectors[1, 1], v$vectors[2, 1], v$vectors[3 ,1])


  v2 <- NULL
  v2$I<-calcul.vecteur.polaire.I(v$vectors[1, 2], v$vectors[2, 2], v$vectors[3 ,2])
  v2$D<-calcul.vecteur.polaire.D(v$vectors[1, 2], v$vectors[2, 2], v$vectors[3 ,2])


  v3 <- NULL
  v3$I<-calcul.vecteur.polaire.I(v$vectors[1, 3], v$vectors[2, 3],v$vectors[3 ,3])
  v3$D<-calcul.vecteur.polaire.D(v$vectors[1, 3], v$vectors[2, 3],v$vectors[3 ,3])


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
#' @param etape.value typiquement la température des mesures d'anisotropie
#' @param etape.sigle sigle indiquant les étapes. les noms peuvent changer, mais pas l'ordre
#' @param volume la valeur du volume du spécimen
#' @param TH la valeur du champ appliqué
#' @return un data.frame avec les colonnes "L1", "L1.Inc", "L1.Dec", "L2", "L2.Inc", "L2.Dec", "L3", "L3.Inc", "L3.Dec", "F13", "F12", "F23"
#' @seealso anisotropie.eigen.tensor
#' @export
anisotropie.eigen.matrix <- function(mesures, etape.value, etape.sigle = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"), volume = 1, TH = 1,...)
{

  ani.etape <- trimws(paste(as.character(etape.value), etape.sigle, sep = "") )

  selec <- NULL
  for (i in 1:length(ani.etape)) {
    selec <- c( selec, which(trimws(mesures$etape) == trimws(ani.etape[i])) )
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
#' @param etape.value typiquement la température des mesures d'anisotropie
#' @param etape.sigle sigle indiquant les étapes. Les noms peuvent changer, mais pas l'ordre
#' @seealso anisotropie.eigen.matrix
#' @return un data.frame avec les colonnes "L1", "L1.Inc", "L1.Dec", "L2", "L2.Inc", "L2.Dec", "L3", "L3.Inc", "L3.Dec", "F13", "F12", "F23"
#' @export
anisotropie.eigen.tensor <- function (mesures, etape.value, etape.sigle = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB") )
{
  ani.etape <- trimws(paste(as.character(etape.value), etape.sigle, sep = "") )

  selec <- NULL
  for (i in 1:length(ani.etape)) {
    selec <- c( selec, which(trimws(mesures$etape) == trimws(ani.etape[i])) )
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
  v1 <- polaire(v$vectors[1, 1], v$vectors[2, 1], v$vectors[3 ,1])
  v2 <- polaire(v$vectors[1, 2], v$vectors[2, 2], v$vectors[3 ,2])
  v3 <- polaire(v$vectors[1, 3], v$vectors[2, 3], v$vectors[3 ,3])


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
anisotropie.eigen.tensors.numbers <- function(numbers, Data.mesures, etape.value, etape.sigle = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"))
{

  Data <- NULL
  for (i in 1:length(numbers) ) {
    mesures <-  extract.mesures.specimen.number(numbers[i], Data.mesures)
    ani <- anisotropie.eigen.tensor(mesures, etape.value = etape.value, etape.sigle = etape.sigle)
    Data <- rbind(Data, ani)
  }
  col.names <- c("L1", "L2", "L3", "L1.Inc", "L1.Dec", "L2.Inc", "L2.Dec", "L3.Inc", "L3.Dec", "F13", "F12", "F23")
  return(as.data.frame(Data, col.names = col.names))

}

#' Calcul des tenseurs d'anisotropie pour une liste de nom de spécimen
#' @param names liste de noms de spécimen
#' @param Data.mesures data.frame contenant les mesures
#' @return un data.frame avec les colonnes "L1", "L2", "L3", "L1.Inc", "L1.Dec", "L2.Inc", "L2.Dec", "L3.Inc", "L3.Dec", "F13", "F12", "F23"
#' @export
anisotropie.eigen.tensors.names <- function(names, Data.mesures, etape.value, etape.sigle = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"))
{

  Data <- NULL
  for (i in 1:length(names) ) {
    mesures <-  extract.mesures.specimen.names(names[i], Data.mesures)
    ani <- anisotropie.eigen.tensor(mesures, etape.value = etape.value, etape.sigle = etape.sigle)
    Data <- rbind(Data, ani)
  }
  col.names <- c("L1", "L2", "L3", "L1.Inc", "L1.Dec", "L2.Inc", "L2.Dec", "L3.Inc", "L3.Dec", "F13", "F12", "F23")
  return(as.data.frame(Data, col.names = col.names))
}

#' Calcul des tenseurs d'anisotropie pour tous les spécimens d'une liste de mesures
#' @return un data.frame avec les colonnes "L1", "L2", "L3", "L1.Inc", "L1.Dec", "L2.Inc", "L2.Dec", "L3.Inc", "L3.Dec", "F13", "F12", "F23"
#' @export
anisotropie.eigen.tensors.all <- function(Data.mesures, Data.number, etape.value, etape.sigle = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"))
{
  numbers <- NULL
  for (i in 1: length(Data.number) ) {
    if ( length(which( numbers == Data.number[i] )) == 0 )
      numbers <- c(numbers, Data.number[i])
  }

  Data <- NULL
  for (i in 1:length(numbers) ) {
    mesures <-  extract.mesures.specimen.number(numbers[i], Data.mesures)
    ani <- anisotropie.eigen.tensor(mesures, etape.value = etape.value, etape.sigle = etape.sigle )
    Data <- rbind(Data, ani)
  }

  col.names <- col.names <- c("L1", "L2", "L3", "L1.Inc", "L1.Dec", "L2.Inc", "L2.Dec", "L3.Inc", "L3.Dec", "F13", "F12", "F23")
  return(as.data.frame(Data, col.names = col.names))


}


#' Calcul de la matrice moyenne symétrisée d'anisotropie moyen
#' @param Data.mesures data.frame contenant les mesures
#' @param etape.value typiquement la température des mesures d'anisotropie
#' @param Data.number numéro des échantillons
#' @seealso anisotropie.mean.eigen.tensor
#' @return une matrice carrée 3 x3
#' @export
anisotropie.mean.matrix <- function(Data.mesures, Data.number, etape.value, etape.sigle = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"))
{
  ani.moyen <- matrix(c(0, 0, 0, 0, 0, 0, 0, 0, 0), nrow = 3, ncol = 3)
  for (spe in 1: length(Data.number))
  {
    spe.mes <- extract.mesures.specimen.number(Data.number[spe], Data.mesures)
    ani.moyen <- ani.moyen + anisotropie.matrix.symetric(mesures = spe.mes, etape.value = etape.value, etape.sigle = etape.sigle)
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
#' @param etape.value typiquement la température des mesures d'anisotropie
#' @param Data.number numéro des échantillons
#' @seealso anisotropie.eigen.matrix, anisotropie.mean.matrix
#' @return un data.frame avec les colonnes "L1", "L2", "L3", "L1.Inc", "L1.Dec", "L2.Inc", "L2.Dec", "L3.Inc", "L3.Dec", "F13", "F12", "F23"
#' @export
anisotropie.mean.eigen.tensor <- function(Data.mesures, Data.number, etape.value, etape.sigle = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"))
{
  mat.sym <- anisotropie.mean.matrix(Data.mesures, Data.number, etape.value = etape.value, etape.sigle = etape.sigle)
  v <- eigen(mat.sym, symmetric = TRUE)

  # calcul des angles I,D des vecteurs propres}
  v1 <- polaire(v$vectors[1, 1], v$vectors[2, 1], v$vectors[3 ,1])
  v2 <- polaire(v$vectors[1, 2], v$vectors[2, 2], v$vectors[3 ,2])
  v3 <- polaire(v$vectors[1, 3], v$vectors[2, 3], v$vectors[3 ,3])

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
correction.anisotropie <- function(mesures.frame, ani.matrix)
{

  if (!is.data.frame(mesures.frame))
    warning("mesures n'est pas une data.frame")
  if (!is.matrix(ani.matrix))
    warning("ani.matrix n'est pas une matrice")

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
    re.cal.ID <- polaire(re.cal$X, re.cal$Y, re.cal$Z)
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
lambert.ID.tensors <- function(Data, point.col = "blue3", new = TRUE, ...)
{

  # restructuration des données
  L1.Inc <- as.numeric(Data$L1.Inc)
  L2.Inc <- as.numeric(Data$L2.Inc)
  L3.Inc <- as.numeric(Data$L3.Inc)

  L1.Dec <- as.numeric(Data$L1.Dec)
  L2.Dec <- as.numeric(Data$L2.Dec)
  L3.Dec <- as.numeric(Data$L3.Dec)

  if (length(L1.Inc) == 0)
    return(print("no DATA"))

  if (length(point.col) < length(L1.Inc))
    point.col <- rep(point.col, length.out = length(L1.Inc))

  for (i in 1:length(L1.Inc) ) {
    Da.I <- c(L1.Inc[i], L2.Inc[i], L3.Inc[i])
    Da.D <- c(L1.Dec[i], L2.Dec[i], L3.Dec[i])
    lambert.ID(Da.I, Da.D,
                    inc.lim = c(0, 90), pch = c(22, 24, 21), point.col = point.col[i], new = new  )
    new <- FALSE
  }
}

#' Tracer d'un diagramme de Flin
#' diagramme des paramètres d'anisotropie F12 en fonction de F23
#' @param Data soit correspond à un data.frame avec les variables $F23 et $F12, soit seulement des valeurs de F23 dans ce cas renseignez Data.F12
#' @param Data.F12 des valeurs pour F12
#' @param absolue prend la valeur absolue des données
#' @export
flin <- function( Data, Data.F12=NULL, pt.names= NULL, point.col = "blue3", pch = 21, type= "p",
                       xlab = "F23", ylab = "F12", main = "Flin diagram", absolue = TRUE, new = TRUE)
{

   par(pty="s")

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

  if (new == TRUE) {
    plot(x = X, y = Y, xlab = xlab, ylab = ylab, type = type, col = "gray50", bg = point.col, pch = pch,
         xaxt = "n", yaxt = "n", asp = 1, bty ="n", main = main, new = TRUE)
    ax1 <- axis(1, pos = 1, cex.axis = 0.8, col = "gray10")
    ax2 <- axis(2, pos = 1, cex.axis = 0.8, col = "gray10") # Ordonnées
    cc<-array(c(0,1), c(1,2))

    abline(  coef = cc, col = "gray90")

    main <- ""
  }  else {
    points(x = X, y = Y, type = type, col = "gray50", bg = point.col, pch = pch,
         xaxt = "n", yaxt = "n", asp = 1, bty ="n", new = FALSE)
  }

  text(jitter(X, 5, amount = 0), jitter(Y, 5, amount = 0), pt.names)

}

#' Tracer d'un diagramme de désaimantation
#' @param normalize permet de comparer l'évolution de l'aimanation quelque soit l'amplitude en visualisant le résultat comme un pourcentage du maximum
#' @export
desaim <- function( Data, F = NULL,  point.col = "blue3", pch = 21, type = "b",
                             xlab = "°C", ylab = "", main = NULL,
                             names = NA, normalize = TRUE, etape.J0 = NULL, new = TRUE, ...)
{

  tmp.frame <- NULL
  if (is.null(main))
    main <- " F vs etape.value"

  if(is.data.frame(Data)) {
    tmp.frame$name <- Data$name[!is.na(Data$etape.value)]
    tmp.frame$etape.value <- Data$etape.value[!is.na(Data$etape.value)]
    tmp.frame$F <- Data$F[!is.na(Data$etape.value)]


  } else {
    # Création du frame interne
    if (is.na(names))
      tmp.frame$name <- as.character(c(1: length(Data)))
    else
      tmp.frame$name <- names

    tmp.frame$etape.value <- Data
    tmp.frame$F <- F
  }

  if (length(point.col) < length(tmp.frame$etape.value) )
    point.col <- rep(point.col, length(tmp.frame$etape.value) )

  # comptage et séparation des données
  current <- tmp.frame$name[1]
  list.name <- current

  for (i in 1: length(tmp.frame$name)) {
    if ( tmp.frame$name[i] != current) {
      current <- tmp.frame$name[i]
      list.name<- c(list.name, current)
    }
  }

  if (is.null(etape.J0) == FALSE) {
    if (length(etape.J0) < length(tmp.frame$etape.value) )
      etape.J0 <- rep(etape.J0, length(list.name) )
  } else {
    etape.J0 <- rep(NA, length(list.name) )
  }

  xlim <- range(tmp.frame$etape.value)

  # recherche du Ymax pour définir la taille du graphique
  Ymax <- 0
  J0 <- 0
  list.mesure.i <-  NULL
  if (normalize == TRUE) {
    for (i in 1 : length(list.name)) {
      list.mesure.i$F <- tmp.frame$F[tmp.frame$name == list.name[i]]
      if (is.na(etape.J0[i]) == TRUE) {
        J0[i] <- list.mesure.i$F[1]
      } else {
        tmp.etp <- list.mesure.i$F[tmp.frame$etape.value == etape.J0[i]] # si il y a plusieurs valeurs pour la même étape, comme ani
        J0[i] <- tmp.etp[1]
       # J0[i] <- list.mesure.i$F[tmp.frame$etape.value == etape.J0[i]]
      }
      Ymax <- max(Ymax, tmp.frame$F[tmp.frame$name == list.name[i]]/J0[i] *100)
    }
  } else
    Ymax <- max(tmp.frame$F)

  list.mesure.i <-  NULL
  for (i in 1 : length(list.name)) {
    list.mesure.i$etape.value <- tmp.frame$etape.value[which(tmp.frame$name == list.name[i])]
    list.mesure.i$F <- tmp.frame$F[tmp.frame$name == list.name[i]]
    if (normalize == TRUE) {
      coefY <- 100 / J0[i]

    } else
      coefY <- 1

    Xi <- list.mesure.i$etape.value
    Yi <- list.mesure.i$F  * coefY

    ylim <- c(0, Ymax)

    if (i == 1)
      plot(x = Xi, y = Yi, xlab = xlab, ylab = ylab, ylim = ylim, xlim = xlim, type = type,
           col = adjustcolor( point.col[i], alpha.f = 0.7), bg = point.col[i], pch = pch,
           yaxt = "n", bty = "n", main = main, new = new)
    else
      lines(x = Xi, y = Yi, type = type, col = adjustcolor( point.col[i], alpha.f = 0.7), bg = point.col[i], pch = pch, yaxt = "n", bty ="n")

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
    mesures.convent$I[i] <- calcul.vecteur.polaire.I(mesures.convent$X[i], mesures.convent$Y[i], mesures.convent$Z[i])
    mesures.convent$D[i] <- calcul.vecteur.polaire.D(mesures.convent$X[i], mesures.convent$Y[i], mesures.convent$Z[i])
  }
  return(mesures.convent)
}

# Composante partielle ----

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
    mesures.convent$I[i] <- calcul.vecteur.polaire.I(mesures.convent$X[i], mesures.convent$Y[i], mesures.convent$Z[i])
    mesures.convent$D[i] <- calcul.vecteur.polaire.D(mesures.convent$X[i], mesures.convent$Y[i], mesures.convent$Z[i])
  }
  return(mesures.convent)
}

#' Calcul les directions du vecteur partiel pour une série d'étape
#' @return une data.frame "X", "Y", "Z", "I", "D", "F", "Sl", "MAD"
#' @export
vecteur.partiel <- function(TabX, TabY, TabZ, en0 = TRUE)
{
  col.names <- c("X", "Y", "Z", "I", "D", "F", "Sl", "MAD")
  ntab <- length(TabX);

  if ( (ntab<1) || (ntab!=length(TabY)) || (ntab!=length(TabZ))) {
    Data <- c( X = 0, Y = 0, Z = 0,
               I = 0, D = 0, F = 0,
               Sl = NA, MAD = NA )

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
    warning('MAD non calculable')
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
  v3 <- polaire(v.cp$vectors[1, 3],  v.cp$vectors[2, 3],  v.cp$vectors[3 ,3])

  vcorrect<- NULL
  vcorrect$I <- sign(TabZ[1] - TabZ[ntab]) * abs(v3$I)

  vcorrect$D <- v3$D

  if (sign(vcorrect$I) !=  sign(v3$I)) {
    vcorrect$D <- D.AM(v3$D +180)
  }


  # mise en forme du résultat
  MAD <- MAD * 180 /pi

  col.names <- c("X", "Y", "Z", "I", "D", "F", "Sl", "MAD")
  Data <- c( X = v.cp$vectors[1, 3], Y = v.cp$vectors[2, 3], Z = v.cp$vectors[3, 3],
             I = vcorrect$I, D = vcorrect$D, F = v3$F,
             Sl = Sl, MAD = MAD )

  return(as.data.frame(t(Data), col.names = col.names))
}

#' Calcul le vecteur de la composante partielle
#' et retourne aussi le MAD
#' @param en0 permet de calculer la composante qui passe par l'origine (0, 0)
#' @return une data.frame "X", "Y", "Z", "I", "D", "F", "Sl", "MAD", "DANG"
#' @seealso composante.partielle.T1T2
#' @export
composante.partielle <- function(TabX, TabY, TabZ, en0 = FALSE)
{
  vp <- vecteur.partiel(TabX, TabY, TabZ, en0 = en0)
  ntab <- length(TabX)

  # Calcul de l'angle entre le vecteur passant par zéro et le vecteur ne passant pas par zéro
  # d'aprés Lisa Tauxe
  v.true <- vecteur.partiel(TabX, TabY, TabZ, en0 = TRUE)
  v.false <- vecteur.partiel(TabX, TabY, TabZ, en0 = FALSE)

  # produit scalaire x*xp+y*yp+z*zp=cos(teta)*(x2+y2+z2)*(xp2+yp2+zp2);
  DANG <- v.true$X * v.false$X + v.true$Y * v.false$Y + v.true$Z * v.false$Z
  DANG <- DANG / (sqrt( v.true$X^2 + v.true$Y^2 + v.true$Z^2)* sqrt( v.false$X^2+ v.false$Y^2+ v.false$Z^2) )
  DANG <- n_arccos(DANG) *180/pi
  return( cbind.data.frame(vp, DANG))
}

#' Calcul le vecteur de la composante partielle en calculant la composante entre les étapes T1 et T2
#' par défaut la fonction n'affiche pas les étapes d'anisotropie
#' @param en0 permet de calculer la composante qui passe par l'origine (0, 0)
#' @param corr.ani corrige de l'anisotropie
#' @return une data.frame "X", "Y", "Z", "I", "D", "F", "Sl", "MAD", "DANG"
#' @seealso composante.partielle, zijderveld1.T1T2
#' @export
composante.partielle.T1T2 <- function(Data, T1=NULL, T2=NULL, corr.ani = FALSE, ani.etape.value= NULL,
                                      etape.sigle = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"),
                                      en0 = FALSE )
{

  if (!is.data.frame(Data))
    warning("Data n'est pas une data.frame")

  if (is.null(T1))
    T1 <- Data$etape.value[1]
  if (is.null(T2))
    T2 <- Data$etape.value[length(Data$etape.value)]


  res.list <- Data
  if (corr.ani == TRUE) {

    if (is.null(ani.etape.value)) { # recherche des étapes d'anisotropie
      for (i in 1:length(Data$etape)) {
        if (substring(Data$etape[i], 4) == etape.sigle[3])
          ani.etape.value <- Data$etape.value[i]
      }
    }
    ani <- anisotropie.matrix.symetric(Data, etape.value = ani.etape.value, etape.sigle = etape.sigle)
  }

  Data <- supprime.etape(Data, ani.etape.value = ani.etape.value, etape.sigle = etape.sigle)

  # recherche étape en dessous de T1
  iT1 <- 1
  while(Data$etape.value[iT1] < T1) {
    iT1 <- iT1 + 1
  }

  # recherche étape au dessus de T2
  iT2 <- length(Data$etape.value)
  while(Data$etape.value[iT2] > T2) {
    iT2 <- iT2 -1
  }

  # Selection de la composante partielle

  TabX <- Data$X[iT1:iT2]
  TabY <- Data$Y[iT1:iT2]
  TabZ <- Data$Z[iT1:iT2]

  # Correction des direction par l'anisotropie
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
  v.true <- vecteur.partiel(TabX, TabY, TabZ, en0 = TRUE)
  v.false <- vecteur.partiel(TabX, TabY, TabZ, en0 = FALSE)


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
#' @param T1 correspond à etape.value ou température la plus basse
#' @param T2 correspond à etape.value ou température la plus haute
#' @param en0 booléen permettant de forcer la composante partielle de passer par l'origine
#' @param show.etape booléen permettant l'affichage des étapes
#' @param withAni permet de voir les étapes d'anisotropie
#' @param ani.etape.value correspond à l'etape.value ou la température de la détermination de l'anisotropie
#' @param etape.sigle chaîne de caractère représentant les étapes de l'anisotropie "Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB". Cette ordre est obligatoire
#' @seealso composante.partielle.T1T2, zijderveld2.T1T2
#' @export
zijderveld1.T1T2 <- function(Data, T1 = NULL, T2 = NULL, show.etape = FALSE, withAni = FALSE, ani.etape.value = NULL,
                             etape.sigle = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"),
                             legend.pos = NULL, legend.txt = c("(Y, X)", "(Y, Z)"),
                             en0 = FALSE )
{
  if (is.null(T1))
    T1 <- Data$etape.value[1]
  if (is.null(T2))
    T2 <- Data$etape.value[length(Data$etape.value)]


  res.list <- Data

  if (withAni == FALSE) { # suppression des étape d anisotropie
    Data <- supprime.etape(Data, ani.etape.value = ani.etape.value, etape.sigle = etape.sigle)
  }
  # recherche étape en dessous de T1
  iT1 <- 1
  while(Data$etape.value[iT1] < T1) {
    iT1 <- iT1 + 1
  }

  # recherche étape au dessus de T2
  iT2 <- length(Data$etape.value)
  while(Data$etape.value[iT2] > T2) {
    iT2 <- iT2 -1
  }

  # Calcul de la composante partielle

  TabX <- Data$X[iT1:iT2]
  TabY <- Data$Y[iT1:iT2]
  TabZ <- Data$Z[iT1:iT2]

  vp <- composante.partielle(TabX, TabY, TabZ, en0 = en0)
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

  if (show.etape == TRUE) {
    pt.names <- Data$etape
  } else {
    pt.names <- rep("", length(Data$X) )
  }

  zijderveld1(Data$X[1:iT1], Data$Y[1:iT1], Data$Z[1:iT1],  ylim = ylim, pt.names = pt.names[1:iT1], legend.pos = legend.pos, legend.txt = legend.txt )
  if (iT2 != length(Data$X))
    zijderveld1(Data$X[iT2:length(Data$X)], Data$Y[iT2:length(Data$X)], Data$Z[iT2:length(Data$X)], pt.names = pt.names[iT2:length(Data$X)], new = FALSE)

  zijderveld1(Data$X[iT1:iT2], Data$Y[iT1:iT2], Data$Z[iT1:iT2], pt.col = c("red", "red") , pt.names = pt.names[iT1:iT2], new = FALSE)

  abline(ad, bd)
  abline(ai, bi)

  Data$etape
}

#' trace un diagramme de zijderveld en calculant la composante entre les étapes T1 et T2
#' @param T1 correspond à etape.value ou température la plus basse
#' @param T2 correspond à etape.value ou température la plus haute
#' @param en0 booléen permettant de forcer la composante partielle de passer par l'origine
#' @param show.etape booléen permettant l'affichage des étapes
#' @param withAni permet de voir les étapes d'anisotropie
#' @param ani.etape.value correspond à l'etape.valeu ou la température de la détermination de l'anisotropie
#' @param etape.sigle chaîne de caractère représentant les étapes de l'anisotropie "Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB". Cette ordre est obligatoire
#' @seealso composante.partielle.T1T2, zijderveld2.T1T2
#' @export
zijderveld2.T1T2 <- function(Data, T1 = NULL, T2 = NULL, show.etape = FALSE, withAni = FALSE, ani.etape.value = NULL,
                             etape.sigle = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"),
                             legend.pos = NULL, legend.txt = c("(Y, X)", "(Y, Z)"),
                             en0 = FALSE )
{
  if (is.null(T1))
    T1 <- Data$etape.value[1]
  if (is.null(T2))
    T2 <- Data$etape.value[length(Data$etape.value)]


  res.list <- Data

  if (withAni == FALSE) { # suppression des étape d anisotropie
    Data <- supprime.etape(Data, ani.etape.value = ani.etape.value, etape.sigle = etape.sigle)
  }
  # recherche étape en dessous de T1
  iT1 <- 1
  while(Data$etape.value[iT1] < T1) {
    iT1 <- iT1 + 1
  }

  # recherche étape au dessus de T2
  iT2 <- length(Data$etape.value)
  while(Data$etape.value[iT2] > T2) {
    iT2 <- iT2 -1
  }

  # Calcul de la composante partielle

  TabX <- Data$X[iT1:iT2]
  TabY <- Data$Y[iT1:iT2]
  TabZ <- Data$Z[iT1:iT2]

  vp <- composante.partielle(TabX, TabY, TabZ, en0 = en0)
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

  Y.r <- range(Data$X)
  Z.r <- range(Data$Z)
  ylim <- c(min(Y.r[1], -Z.r[2]), max(Y.r[2], -Z.r[1]))

  if (ylim[1]>0)
    ylim[1]<-0

  if (ylim[2]<0)
    ylim[2]<-0

  if (show.etape == TRUE) {
    pt.names <- Data$etape
  } else {
    pt.names <- rep("", length(Data$X) )
  }

  zijderveld2(Data$X[1:iT1], Data$Y[1:iT1], Data$Z[1:iT1], pt.names = pt.names, ylim = ylim, legend.pos = legend.pos, legend.txt = legend.txt )
  if (iT2 != length(Data$X))
    zijderveld2(Data$X[iT2:length(Data$X)], Data$Y[iT2:length(Data$X)], Data$Z[iT2:length(Data$X)], pt.names = pt.names, new = FALSE)

  zijderveld2(Data$X[iT1:iT2], Data$Y[iT1:iT2], Data$Z[iT1:iT2], pt.col = c("red", "red") , pt.names = pt.names, new = FALSE)

  abline(ad, bd)
  abline(ai, bi)
}

#' Supprime les étapes à la valeur ani.etape.value et avec les sigles etape.sigle
#' Utiliser par défaut pour supprimer les étapes d'anisotropie
#' @param Data un data.frame possédant la variable $etape de type chr
#' @param ani.etape.value valeur de la température d'anisotropie. Par defaut NUll, alors la tempértature est retrouvée automatiquement avec les etape.sigle définis
#' @param etape.sigle chaîne de caractère représentant les étapes de l'anisotropie "Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB"
#' @export
supprime.etape <- function(Data, ani.etape.value= NULL, etape.sigle = c("Z+", "Z-", "X+", "X-", "Y+", "Y-", "ZB") )
{
  if (is.null(ani.etape.value)) {
    for (i in 1:length(Data$etape)) {
      if (substring(Data$etape[i], 4) == etape.sigle[3])
        ani.etape.value <- Data$etape.value[i]
    }
  }
  ani.etape <- trimws(paste(ani.etape.value, etape.sigle, sep = ""))


  selec <- NULL
  for (i in 1:length(ani.etape)) {
    selec <- c( selec, which(trimws(Data$etape) == trimws(ani.etape[i])) )
  }

  if (length(selec) > 0)
    res.list <- Data[-selec,]
  else
    res.list <- Data

  return(res.list)
}

