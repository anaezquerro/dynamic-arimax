
library(dplyr)   

column_names <- c('ccaa', 'idcaso', 'anno', 'semana', 'edad', 'sexo', 'codigomc', 'sdgripal', 'sarscov2')


read_csv_filter <- function(filename, sep=',') {
    dat <- read.csv(paste0('iras/', filename, '.csv'), sep=sep)
    dat <- dat[, column_names]
    return(dat)
}


#' --------------------- LECTURA DE LOS DATOS INDIVIDUALES --------------------
#' En los archivos DatosIras...csv se encuentran los casos individuales de 
#' infecciones respiratorias (de interés, gripe y COVID). Debemos cargar todos 
#' los csv (separados en distintos archivos por semana).
#' Anotamos los archivos en los que se encuentren los datos
archivos <- c('CasosIRAS_S1', 'CasosIRAS_S13 v2', 'CasosIRAS_S40-43_2021',
              'CasosIRAS_S46', 'CasosIRAS_S47', 'CasosIRAS_S48',
              'CasosIRAS_S49', 'CasosIRAS_S50', 'CasosIRAS_S51',
              'CasosIRAS_S52')

#' Los juntamos en un mismo data.frame y eliminamoslos duplicados
individual <- read_csv_filter('CasosIRAS_temp2020-21', sep=';')
dats <- lapply(archivos, read_csv_filter)

for (i in 1:length(dats)) {   # concatenamos todos los data.frames
    individual <- rbind(individual, dats[[i]])
    individual <- unique(individual)
}


cat(paste0('Número total de columnas del data.frame: ', ncol(individual), '\n'))
cat(paste0('Número total de filas del data.frame: ', nrow(individual)))



#' En el dataframe completo añadimos una nueva columna para representar la 
#' el año y semana
individual$fecha <- paste0(individual$anno, '-', 
                           ifelse(individual$semana < 10, "0", ""), individual$semana)

head(individual[, c('anno', 'semana', 'fecha')])   # obsérvese el resultado


#' En el objeto `dat` almacenaremos los datos agregados a nivel semanal. Primero 
#' agregamos los casos de gripe y COVID y ordenamos por fecha.
#' 
dat <- aggregate(individual[, c('sdgripal', 'sarscov2')],
                 by=list(fecha=individual$fecha), FUN=sum, na.rm=T)
dat <- dat[order(dat$fecha),]

head(dat)



# ------------------ AGREGACIÓN DE DATOS POR EDADES ----------------------------
#' Una vez agregados los datos a nivel "completo", procedemos a separar por edades.
#' Los rangos de edad se escogen en base a la separación que hacen los archivos 
#' FitxerPoblacio...csv. La idea es tener, para cada semana, el número de 
#' individuos que figuran en FitxerPoblacio...csv con la edad de su respectivo rango.


rangos_edad <-  list("edad04" = c(0, 4), "edad514" = c(5, 14), "edad1544" = c(15, 44),
                     "edad4564" = c(45, 64), "edad65" = c(65, Inf))

#' Para aplicar la agregación por edades vamos a hacer uso de la siguiente función, 
#' que le asigna un valor discreto de las claves de `rangos_edad` según la edad de 
#' un individuo.
get_ranges <- function(n) {
  if (is.na(n)) {return(NA) }
  for (k in names(rangos_edad)) {
    if (between(n, rangos_edad[[k]][1], rangos_edad[[k]][2]))
      return(k)  
  }
}

get_ranges(10)      # un ejemplo


#' Aplicamos esto a los valores del data.frame de datos individuales
individual$rango_edad <- unlist(lapply(individual$edad, get_ranges))

head(individual[, c('edad', 'rango_edad')])   # obsérvese el resultado


#' Una vez clasificada cada muestra individual en un rango de edad, podemos 
#' proceder a agregar los datos por semana y por grupo de edad
for (rango in names(rangos_edad)) {
    individual_edad <- individual[individual$rango_edad == rango,]
    dat_edad <- aggregate(individual_edad[, c('sdgripal', 'sarscov2')], 
                          by=list(fecha=individual_edad$fecha), FUN=sum, na.rm=T)

    dat[[paste0(rango, '.sdgripal')]] <- dat_edad$sdgripal
    dat[[paste0(rango, '.sarscov2')]] <- dat_edad$sarscov2
}


#' Ahora podemos observar que se han añadido columnas que indican, para cada semana, 
#' el número de casos de gripe y COVID en cada grupo de edad
head(dat)


#' Una vez hecho esto, se añade la información de los datos de FitxerPoblacio...csv, 
#' que indican el número de individuos de cada rango de edad. Para leer estos datos 
#' repetimos el mismo proceso que realizamos para concatenar las filas de los csv 
#' de DatosIRAS


archivos <- c('FitxerPoblacio_S1', 'FitxerPoblacio_S46', 
              'FitxerPoblacio_S47', 'FitxerPoblacio_S48', 'FitxerPoblacio_S49', 
              'FitxerPoblacio_S50', 'FitxerPoblacio_S51', 'FitxerPoblacio_S52')


dats <- lapply(archivos, function(x) read.csv(paste0('iras/', x, '.csv')))

pob <- read.csv('iras/FitxerPoblacio_temp2020-21.csv', sep=';')


for (i in 1:length(dats)) {
    print(ncol(dats[[i]]))
    pob <- rbind(pob, dats[[i]])
    pob <- unique(pob)
}

#' Los datos de FitxerPoblacio_S13.csv tienen otra división de edades que se debe 
#' corregir "a mano".
pob13 <- read.csv('iras/FitxerPoblacio_S13.csv')
head(pob13)
pob13$pob04 <- pob13$pob0_4
pob13$pob514 <- pob13$pob5_9 + pob13$pob10_14
pob13$pob1544 <- apply(pob13[, c('pob15_19', 'pob20_24', 'pob25_29', 'pob30_34', 'pob35_39', 'pob40_44')], 1, sum)
pob13$pob4564 <- apply(pob13[, c('pob45_49', 'pob50_54', 'pob55_59', 'pob60_64')], 1, sum)
pob13$pob65 <- apply(pob13[, c('pob65_69', 'pob70_74', 'pob75_79', 'pob80_84', 'pob85_89', 'pob90_94', 'pob95_99', 'pob100mas')], 1, sum)


# lo añadimos al resto de datos
pob <- unique(rbind(pob, pob13[, colnames(pob)]))

# Añadimos la fecha de la misma forma que con `individual`
pob$fecha <- paste0(pob$anno, '-', ifelse(pob$semana < 10, "0", ""), pob$semana)
pob <- pob[order(pob$fecha), c('fecha', 'pob04', 'pob514', 'pob1544', 'pob4564', 'pob65')]
head(pob)

# Lo añadimos a los datos agregados
dat <- merge(dat, pob, by=c('fecha'))
head(dat)


# ---------------------- AGREGACIÓN DE DATOS POR REGIONES ----------------------------

conv <- read.csv('iras/conversió regió sanitària.csv')
names(conv) <- c("region", "codigomc")

# Cambiamos los nombres de las regiones 
conv$region[conv$region == "ALT PIRINEU - ARAN"] <- "PIRINEU"
conv$region[conv$region == "BARCELONA CIUTAT"] <- "BARCELONA"
conv$region[conv$region == "CATALUNYA CENTRAL"] <- "CANTALUNYA"
conv$region[conv$region == "CATALUNYA CENTRAL"] <- "CANTALUNYA"
conv$region[conv$region == "METROPOLITANA NORD"] <- "METR_NORD"
conv$region[conv$region == "METROPOLITANA SUD"] <- "METR_SUD"
conv$region[conv$region == "TERRES DE L'EBRE"] <- "TERRES_EBRE"
  

head(conv)

#' Ahora en el data.frame `individual` cambiamos el códigomc por la región

get_region <- function(x) {
    if (x %in% conv$codigomc) {
        return(conv$region[conv$codigomc == x])
    }
    return('UNK')
}

get_region(individual$codigomc[1])   # ejemplo

individual$region <- unlist(lapply(individual$codigomc, get_region))
head(individual[, c('codigomc', 'region')])


#' Ahora agregamos por región de la misma forma que hicimos por edades

for (region in unique(conv$region)) {
    individual_region <- individual[(!is.na(individual$region)) & (individual$region == region),]
    dat_region <- aggregate(individual_region[, c('sdgripal', 'sarscov2')],
                            by=list(fecha=individual_region$fecha), FUN=sum, na.rm=T)
    
    dat_region <- merge(dat, dat_region, all.x=T, by=c('fecha'))[, c('fecha', 'sdgripal.y', 'sarscov2.y')]
    dat[[paste0(region, '.sdgripal')]] <- dat_region$sdgripal
    dat[[paste0(region, '.sarscov2')]] <- dat_region$sarscov2
    
}

head(dat)


# --------------------------- VACUNACIÓN.-----------------------------
vac1 <- read.csv('iras/CobVac_2021-11-03.csv')
vac2 <- read.csv('iras/CobVac_2021-12-09.csv')
vac <- unique(rbind(vac1, vac2))
colnames(vac) <- c('fecha', 'vac1218', '% [12,18)', 'vac1845', '% [18,45)',
                  'vac4565', '% [45,65)', 'vac6580', '% [65,80)', 'vac80',
                  '% [80,Inf)', 'vactotal', '% total')

head(vac)

# Añadimos los datos al objeto final
vac <- vac[, c('fecha', colnames(vac)[startsWith(colnames(vac), 'vac')])]
head(vac)


dat <- merge(dat, vac, by=c('fecha'), all.x=T)
head(dat)


write.csv(dat, file="data/evolucion_gripe_covid.csv", row.names=F)



