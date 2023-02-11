library(tidyverse)
library(sampleSelection)
library(plm)

##### Modelos de datos panel ####

### valuidacion de supuestos
## 1. Normalidad
## 2. Homocedasticidad
## 3. Independencia de los residuos

# cargue de la base de datos
data_cafe = read.csv('data/data_cafe.csv',sep = ';', dec = ',')

#### Analisis Exploratorio de Datos

#eliminar variables con NA's

#al contener nas, estas variables no se inlcuyen en el analisis
columnas_na = c('cagrop_q.1',
                'cagrop_vlr.1',
                'credes_q.1',
                'credes_vlr.1',
                'csustit_q.1',
                'csustit_vlr.1',
                'lkdtrab_q.1',
                'lkdtrab_vlr.1',
                'linver_q.1',
                'linver_vlr.1',
                'lnorm_q.1',
                'lnorm_vlr.1')
data_cafe_clean = data_cafe %>% select(-all_of(columnas_na))

# Matriz de correlaciones
cor_matrix = cor(data_cafe_clean)
correlacionadas <- abs(cor_matrix) > 0.7
variables_correlacionadas <- names(data_cafe_clean)[colSums(correlacionadas) > 0]

## seleccion final de variables
## bajo 2 criterios, el primero correlacion absoluta superior al 70% y variables independientes
variables_modelo_cafe = c(
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

independientes = variables_modelo_cafe[-c(1,2,3)]
dependiente = 'protvd'
municipio = 'cod_mun'
anio = 'año'

#dataset de insumo para el modelo
input_modelo_cafe = data_cafe_clean %>% select(all_of(variables_modelo_cafe))
# creacion de variable dicotomica para los modelos heckman y sample selection
input_modelo_cafe$acceso_credito = ifelse(input_modelo_cafe$qtotal_cred>0,TRUE,FALSE)

# 1. Primer modelo de selección de Heckman
heckman_cafe_modelo <- heckit(acceso_credito~pcc,protvd ~ cagrop_vlr+linver_vlr,input_modelo_cafe)
resumen_heckman_cafe = summary(heckman_cafe_modelo)


# 2. Segundo modelo de sample selection
#+lkdtrab_vlr+linver_vlr+lnorm_vlr
selection_cafe_modelo <- selection(acceso_credito~pcc,protvd ~ peqp_vlr+cagrop_vlr,input_modelo_cafe)
resumen_selection_cafe = summary(selection_cafe_modelo)

# 3. Modelo clasico de datos panel de efectos aleatorios

independientes_ = independientes[-c(1,9)]
lista_modelos_cafe = do.call("c", lapply(seq_along(independientes_), function(i) combn(independientes, i, FUN = list)))

lista_formulas_modelos = lapply(lista_modelos_cafe, function(x){
  as.formula(paste(dependiente,'~',paste(x,collapse = "+")))
})

list_modelos_reg_aleatorios = lapply(lista_formulas_modelos, function(x){
  list_modelos_reg_aleatorios = plm(x, index=c(anio,municipio), model="random", data=input_modelo_cafe)
  list_modelos_reg_aleatorios = summary(list_modelos_reg_aleatorios)
  df_salida_modelos = as.data.frame(t(list_modelos_reg_aleatorios$r.squared))
  df_salida_modelos = cbind(df_salida_modelos,paste(names(list_modelos_reg_aleatorios$model),collapse = '-'))
  names(df_salida_modelos) = c('rsq','adjrsq','variables')
  df_salida_modelos
})

# lista de todos los modelos posubles iterativos sobre las variables disponibles
# tener en cuenta que plm limita el uso de maximo 9 variables como explicativas
full_modelos_salida = do.call('rbind',list_modelos_reg_aleatorios)

# bajo la logica anterior el mejor modelo es:
modelo_clasico = plm(protvd~peqp_vlr+cagrop_vlr+credes_vlr+linver_vlr+pcc, index=c(anio,municipio), model="random", data=input_modelo_cafe)
resumen_modelo_clasico = summary(modelo_clasico)

## comparacion de modelos

# stargazer
stargazer::stargazer(heckman_cafe_modelo, selection_cafe_modelo,modelo_clasico,
          se=list(NULL,NULL,NULL),
          title="Comparacion de modelos", type="text",
          df=FALSE, digits=4)


################################## supuestos ######################################################

valor_significancia = 0.05
heckman_modelo = heckman_cafe_modelo
selection_modelo = selection_cafe_modelo
input_modelo = input_modelo_cafe

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

