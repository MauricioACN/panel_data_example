library(tidyverse)
library(sampleSelection)
library(plm)

##### Modelos de datos panel ####

### valuidacion de supuestos
## 1. Normalidad
## 2. Homocedasticidad
## 3. Independencia de los residuos

# Prueba de robustez

# data_cafe = read.csv('data/data_cafe.csv',sep = ';', dec = ',')
# data_modelo = read.csv('data/data_cania.csv',sep = ';', dec = ',')
data_modelo = read.csv('data/data_papa.csv',sep = ';', dec = ',')
# data_modelo_platano =

#### Analisis Exploratorio de Datos
##

for (variable in names(data_modelo)) {

  if (!is.numeric(data_modelo[,variable])) {
    data_modelo[,variable] = as.numeric(data_modelo[,variable])
  }

}

summary(data_modelo)

data_modelo = drop_na(data_modelo)

#eliminar variables con NA's

# columnas_na = c('granp_q','granp_vlr')
# data_clean = data_modelo %>% select(-all_of(columnas_na))
data_clean = data_modelo

cor_matrix = cor(data_clean)
correlacionadas <- abs(cor_matrix) > 0.7
variables_correlacionadas <- names(data_modelo)[colSums(correlacionadas) > 0]

## seleccion final de variables
## bajo 2 criterios, el primero correlacion absoluta superior al 70% y variables independientes
variables_modelo = c(
  "cod_mun",
  "año",
  "protvd",
  "qtotal_cred",
  # "vlrtotal_cred",
  # "granp_q",
  # "granp_vlr",
  # "medp_q",
  # "medp_vlr",
  # "peqp_q",
  "peqp_vlr",
  # "hom_q",
  # "hom_vlr",
  # "muj_q",
  # "muj_vlr",
  # "persjur_q",
  # "persjur_vlr",
  # "cagrop_q",
  "cagrop_vlr",
  # "credes_q",
  "credes_vlr",
  # "csustit_q",
  "csustit_vlr",
  # "lkdtrab_q",
  "lkdtrab_vlr",
  # "linver_q",
  "linver_vlr",
  # "lnorm_q",
  "lnorm_vlr",
  'pcc'
)

independientes = variables_modelo[-c(1,2,3)]
dependiente = 'protvd'
municipio = 'cod_mun'
anio = 'año'

input_modelo = data_clean %>% select(all_of(variables_modelo))
input_modelo$acceso_credito = ifelse(input_modelo$qtotal_cred>0,TRUE,FALSE)

# 1. Primer modelo de selección de Heckman
heckman_modelo <- heckit(acceso_credito~pcc,protvd ~ peqp_vlr+linver_vlr+lnorm_vlr,input_modelo)
resumen_heckman = summary(heckman_modelo)


# 2. Segundo modelo de sample selection
selection_modelo <- selection(acceso_credito~pcc,protvd ~ peqp_vlr+linver_vlr+lnorm_vlr,input_modelo)
resumen_selection = summary(selection_modelo)

# 3. Modelo clasico de datos panel de efectos aleatorios

independientes_ = independientes[-c(1,9)]
lista_modelos_cafe = do.call("c", lapply(seq_along(independientes_), function(i) combn(independientes, i, FUN = list)))

lista_formulas_modelos = lapply(lista_modelos_cafe, function(x){

  as.formula(paste(dependiente,'~',paste(x,collapse = "+")))
})

list_modelos_reg_aleatorios = lapply(lista_formulas_modelos, function(x){

  list_modelos_reg_aleatorios = plm(x, index=c(anio,municipio), model="random", data=input_modelo)
  list_modelos_reg_aleatorios = summary(list_modelos_reg_aleatorios)
  df_salida_modelos = as.data.frame(t(list_modelos_reg_aleatorios$r.squared))
  df_salida_modelos = cbind(df_salida_modelos,paste(names(list_modelos_reg_aleatorios$model),collapse = '-'))
  names(df_salida_modelos) = c('rsq','adjrsq','variables')
  df_salida_modelos
})


full_modelos_salida = do.call('rbind',list_modelos_reg_aleatorios)

# bajo la logica anterior el mejor modelo es:

modelo_clasico = plm(protvd~qtotal_cred+peqp_vlr+csustit_vlr+lkdtrab_vlr, index=c(anio,municipio), model="random", data=input_modelo)
resumen_modelo_clasico = summary(modelo_clasico)

## comparacion de modelos

# stargazer
stargazer::stargazer(heckman_modelo, selection_modelo,modelo_clasico,
                     se=list(NULL,NULL,NULL),
                     title="Resumen de modelos", type="text",
                     df=FALSE, digits=4)

################################## supuestos ######################################################

valor_significancia = 0.05


########### test de normalidad de los residuos ##########
test = shapiro.test(resid(heckman_modelo))
ifelse(test$p.value>valor_significancia,"pasa normalidad en los residuos","no pasa normalidad en los residuos")

test = shapiro.test(resid(selection_modelo))
ifelse(test$p.value>valor_significancia,"pasa normalidad en los residuos","no pasa normalidad en los residuos")

test = shapiro.test(resid(modelo_clasico))
ifelse(test$p.value>valor_significancia,"pasa normalidad en los residuos","no pasa normalidad en los residuos")

# Heckman pasa los test de normalidad en los residuos pero el modelo selection no

########## test de homocedasticidad ##########

#modelo de heckman
n = nrow(input_modelo)
ei <- resid(heckman_modelo)
fit = heckit(acceso_credito~pcc,ei^2 ~ peqp_vlr+cagrop_vlr+lkdtrab_vlr+linver_vlr+lnorm_vlr,input_modelo)
R2 <- summary(fit)$rSquared$R2
k <- 2
estadistico <- n * R2
valorP <- pchisq(q=estadistico, df=k, lower.tail=FALSE)
cbind(estadistico, valorP)

ifelse(valorP>valor_significancia,"pasa test de homocedasticidad","no pasa test de homocedasticidad")


#modelo de seleccion
ei <- resid(selection_modelo)
fit <- lm(ei^2 ~ peqp_vlr+cagrop_vlr+lkdtrab_vlr+linver_vlr+lnorm_vlr, data=input_modelo) # Modelando ei^2 ~ x1 + x2
R2 <- summary(fit)$r.squared
k <- 2
estadistico <- n * R2
valorP <- pchisq(q=estadistico, df=k, lower.tail=FALSE)
cbind(estadistico, valorP)

ifelse(valorP>valor_significancia,"pasa test de homocedasticidad","no pasa test de homocedasticidad")


#modelo clasico
test = lmtest::bptest(modelo_clasico)
ifelse(test$p.value>valor_significancia,"pasa test de homocedasticidad","no pasa test de homocedasticidad")

######### test de independencia de los residuos #########

#modelo heckman
test = Box.test(resid(heckman_modelo))
ifelse(test$p.value>valor_significancia,"pasa test de independencia de los residuos","no pasa test de independencia de los residuos")

#modelo selection
test = Box.test(resid(selection_modelo))
ifelse(test$p.value>valor_significancia,"pasa test de independencia de los residuos","no pasa test de independencia de los residuos")

#modelo clasico
test = Box.test(resid(modelo_clasico))
ifelse(test$p.value>valor_significancia,"pasa test de independencia de los residuos","no pasa test de independencia de los residuos")
