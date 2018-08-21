# R-SKRIPT ZUR ANPASSUNG VON GW-FLIESSRATEN AN T-PROFILE
#    NACH SCHMIDT ET AL. (2006)


# Lese Datei mit gemessenen T-Profilen
#   Messtiefen (in cm!!) in den Spaltenueberschriften
Tmess = read.table("T-Messungen.csv", header=T, 
                    na.strings="-9999", sep=";")

# Fehlermass (RMSE), das minimiert wird
rmse = function(x1,x2) {
  if(length(x1) != length(x2)) {
    stop("vectors not of equal length")
  }
  return( sqrt( sum((x1-x2)^2)/length(x2)  ) )
} 

# Modellparameter
rho_f = 1000 # Dichte von Wasser (kg/m3)
TL = 11      # GW-Temperatur (Grad C) in Tiefe L
c_f = 4190   # Spezifische Waermekapazitaet von Wasser (J/kg/K)
K_fs = 2     # Waermeleitfaehigkeit des gesaettigten Sediments (J/s/m/K)
L = 10       # Tiefe ab welcher GW-Temperatur als konstant angenommen wird (m)

# Moeglicher Wertebereich fuer q_z (fuer die Optimierung)
lower = -1.e-4
upper = -1.e-8

# Vektor der Messtiefen (m) gemaess Spaltennamen
tiefen = as.numeric(gsub("z", "", names(Tmess)[4:ncol(Tmess)])) / 100.
ix  = 4:(length(tiefen) + 3)

# Modellgleichung fuer stationaeres T-Profil (nach Schmidt et al. 2006, Gl. 2)
Tmod = function(p, z, T0, L) {
  if (p["q_z"]==0) {
      return( z*(TL-T0)+T0 ) 
    } else {
        return(
          ( exp(p["q_z"] * rho_f * c_f * z / K_fs) - 1 ) * (TL - T0) /
            ( exp(p["q_z"] * rho_f * c_f * L / K_fs) - 1 ) + T0
        )
      }
    }

# Nun wird fuer jedes gemessene T-Profil das optimale q_z gesucht
result = data.frame()
for (i in 1:nrow(Tmess)) {
  T0 = Tmess$z0[i]

  Temp = which(!is.na(Tmess[i,ix]))+3
  depth = tiefen[which(!is.na(Tmess[i,ix]))]
  minfun <- function(p) rmse(Tmess[i,Temp],
                             Tmod(p=c(q_z=as.numeric(p)), 
                                  z=depth, T0=T0, L))
  optimized = optim(par=param, fn=minfun, method="Brent",
                    control=list(maxit=1000),
                    lower=-1.e-5, upper=-1.e-7)
  result=rbind(result, 
               data.frame(TNR=Tmess[i,1], 
                          WEITE=Tmess[i,2],
                          RMSE=optimized$value, 
                          OPTIM=optimized$par))
}

# Ergebnistabelle
result

# Plotten der Ergebnisse
windows()

op <- par(mfrow=c(5,5), mar = c(2,2,2,2))

for (i in 1:nrow(Tmess)) {
  tns=result$TNR[i]
  shore_dist=result$WEITE[i]
  T0=Tmess$z0[i]
  plot(as.numeric(Tmess[i,ix]), tiefen, 
       main=paste(tns, shore_dist, sep=" "), xlab="Temperatur deg C", 
       ylim=c(0.6,0), xlim=c(13,25), pch=16)
  lines(Tmod(p=c(q_z=result$OPTIM[i]), z=depth, T0=T0, L=L), depth, col="black", lwd=2)
  if (i==25) break
}

par(op)

