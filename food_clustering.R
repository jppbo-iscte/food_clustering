### Importação de bibliotecas

install.packages("DataExplorer")
install.packages("corrplot")
install.packages("psych")

library(DataExplorer)
library(corrplot)
library(psych)

### Importar dataset 

df <- read.csv('/work/food.csv')
t(introduce(df))

df_full <- df

# Removendo variáveis de ID e descrição (não se repetem)
df <- drop_columns(df, c('ID', 'Descrip'))
# Removendo variáveis com correlação 1
df <- drop_columns(df, c('Calcium_USRDA', 'Copper_USRDA', 'Folate_USRDA',
                         'Magnesium_USRDA', 'Niacin_USRDA', 'Phosphorus_USRDA',
                         'Riboflavin_USRDA', 'Selenium_USRDA', 'Thiamin_USRDA',
                         'VitA_USRDA', 'VitB12_USRDA', 'VitB6_USRDA',
                         'VitC_USRDA', 'VitE_USRDA', 'Zinc_USRDA'))

t(introduce(df))

### Análise Exploratória

# Verificando a quantidade de zeros no data frame
zero_df <- data.frame()
row_count <- nrow(df)
for (col in names(df)) {
    zero_count <- nrow(df[df[, col]==0,])
    zero_df <- rbind(zero_df, cbind(col, zero_count, round(zero_count/row_count, 3)))
}
colnames(zero_df) <- c('col_name', 'zero_rows', 'percentage')
zero_df
# Apesar dehaver uma grande quantidade de zeros, é algo esperado deste tipo de dados
# Afinal, nem todos os alimentos possuem todos estes componentes em sua composição
# Dificilmente encontrará Açúcar em carnes e peixes
# Da mesma forma forma que não encontrará Gordura em vegetais como alface

# O que causou mais estranheza foram as calorias, porém, vendoa descrição dos item fazem sentido
df_full[df_full$Energy_kcal==0,]

## Análise de variáveis originais

# Profiler
plot_bar(df)

# Input
plot_histogram(df)

# Matriz de correlação
correlation <- cor(df[,-1])
par(oma = c(2, 2, 2, 2)) # space around for text
corrplot.mixed(correlation,
    order = "hclust", #order of variables
    tl.pos = "lt", #text left + top
    upper = "ellipse", 
    tl.cex = 0.5,
    number.cex = 0.5)

# Viabilidade do PCA
cortest.bartlett(correlation)
KMO(correlation)

## Análise de variáveis transformadas

# Aplicando log
df_log <- df
df_log[, -1] <- log(df_log[, -1] + 0.01)

# Input
plot_histogram(df_log)

# Matriz de Correlação
correlation_log <- cor(df_log[,-1])
par(oma = c(2, 2, 2, 2)) # space around for text
corrplot.mixed(correlation_log,
    order = "hclust", #order of variables
    tl.pos = "lt", #text left + top
    upper = "ellipse", 
    tl.cex = 0.5,
    number.cex = 0.4)

# Viabilidade do PCA
cortest.bartlett(correlation_log)
KMO(correlation_log)

### PCA

# O PCA pode ser feito com o dataset original e o dataset transformado com a função logaritmo
dfZ <- scale(df[,-1], center=TRUE, scale=TRUE)
df_logZ <- scale(df_log[,-1], center=TRUE, scale=TRUE)

# Obtenção do modelo sem rotação para determinar o número de componentes para o dataset estandartizado
pc <- principal(dfZ, nfactors=length(names(dfZ)), rotate="none", scores=TRUE)
# Print dos valores próprios - 7 componentes com valores maiores que 1
round(pc$values,3)
# Scree Plot - No método a risca, apenas 1 componente, porém podemos flexibilizar para 2 ou 5 dependendo da variância explicada. Método inconclusivo.
plot(pc$values, type="b", main="Scree Plot", xlab="Number of PC", ylab="Eigenvalue")
# Pela variância acumulada, seriam necessário pelo menos 6 componentes para reter mais de 60% da informação. 
pc$loadings
# Recalculando o PCA definindo o número de componentes para estimar quanto foi retido de cada variável 
pc_5 <- principal(dfZ, nfactors=5, rotate="none", scores=TRUE)
# Selenium_mcg, Manganese_mg, VitE_mg e VitC_mg obtiveram uma retenção menor de 40%
pc_5$communality

# Obtenção do modelo sem rotação para determinar o número de componentes para o dataset transformado com log e estandartizado
pcLog <- principal(df_logZ, nfactors=length(names(df_logZ)), rotate="none", scores=TRUE)
# Print dos valores próprios - 5 componentes com valores maiores que 1
round(pcLog$values,3)
# Scree Plot - No método a risca, apenas 1 componente, porém podemos flexibilizar para 3 dependendo da variância explicada. Método inconclusivo.
plot(pcLog$values, type="b", main="Scree Plot", xlab="Number of PC", ylab="Eigenvalue")
# Pela variância acumulada, seriam necessário pelo menos 3 componentes para reter mais de 60% da informação. 
pcLog$loadings
# Recalculando o PCA definindo o número de componentes para estimar quanto foi retido de cada variável 
pcLog_3 <- principal(df_logZ, nfactors=3, rotate="none", scores=TRUE)
# VitE_mg e VitA_mg obtiveram uma retenção menor de 40%
pcLog_3$communality
# Aumentando o número de componentes pois haviam atributos mal representados
pcLog_5 <- principal(df_logZ, nfactors=5, rotate="none", scores=TRUE)
# A variável mais mal representada foi Calcium_mg com 41%
pcLog_5$communality

# Identificação das componentes principais
pcLog_5r <- principal(df_logZ, nfactors=5, rotate="varimax", scores=TRUE)
pcLog_5r$loadings

pca_loadings <- data.frame(matrix(nrow=23, ncol=0))
pca_loadings <- cbind(pca_loadings, as.data.frame(pcLog_5r$loadings[, 1]))
pca_loadings <- cbind(pca_loadings, as.data.frame(pcLog_5r$loadings[, 2]))
pca_loadings <- cbind(pca_loadings, as.data.frame(pcLog_5r$loadings[, 3]))
pca_loadings <- cbind(pca_loadings, as.data.frame(pcLog_5r$loadings[, 4]))
pca_loadings <- cbind(pca_loadings, as.data.frame(pcLog_5r$loadings[, 5]))
colnames(pca_loadings) <- c('RC1', 'RC2', 'RC5', 'RC3', 'RC4')

round(pca_loadings[order(pca_loadings$RC1), ], 2) # Alterar a componente de ordenação para facilitar a identificação

# RC1 - Riboflavin (B2), Niacin (B3), Thiamin (B1), VitB6, Iron, Zinc, Protein, Phosphorous
# RC1: Complexo_B - Aumenta hemoglobina, glóbulos vermelhos e imunidade

# RC2 - Carb, Fiber, Sugar, Low_VitB12
# RC2: Carboidratos

# RC3 - Fat, Energy, Low_VitC
# RC3: Lipidos

# RC4 - VitE, VitA
# RC4: Vitaminas_Liposoluveis

# RC5 - Copper, Magnesium, Folate, Manganese
# RC5: Minerais

### Dataset final

df_final <- df_full[, 2:3]
df_final

df_final$Complexo_B <- pcLog_5$scores[,1]
df_final$Carboidratos <- pcLog_5$scores[,2]
df_final$Lipidos <- pcLog_5$scores[,3]
df_final$Vitaminas_Liposoluveis <- pcLog_5$scores[,4]
df_final$Minerais <- pcLog_5$scores[,5]

df_final

write.csv(df_final, '/work/pca.csv', row.names=F)

### Importação de bibliotecas

install.packages('cluster')
install.packages('mclust')

library(cluster)
library(mclust)

### Importar dataset ou continuar

# df <- df_final
df <- read.csv('/work/pca.csv')

### Escolhendo o melhor método de linkage para definir os clusteres

dist_food <- dist(df)

for (i in c('single', 'complete', 'average', 'ward.D2')) {
    hclust_food <- hclust(dist_food,method=i)
    plot(hclust_food, main=i, hang=-1)
}

# O melhor método foi o de Ward
hclust_ward <- hclust(dist_food, method='ward.D2')
# Neste dendrograma, o ponto de corte foi por volta de y=80, resultando em 5 clusteres
plot(hclust_ward, main='Ward - 5 Clusters', hang=-1)
rect.hclust(hclust_ward, k=5, border="red")
# Valores médios por cluster
groups_k5 <- cutree(hclust_ward, k=5)
aggregate(df[, -1:-2],list(groups_k5), mean)
# Silhueta dos clusteres
plot(silhouette(groups_k5, dist_food), border=NA)
# Plot results
pairs(df[, -1:-2], col=groups_k5)

# Método do cotovelo para definir o número de clusteres no Kmeans
wssplot <- function(xx, nc=15, seed=1234){
    wss <- (nrow(xx)-1)*sum(apply(xx,2,var))
    for (i in 2:nc){
        set.seed(seed)
        wss[i] <- sum(kmeans(xx, centers=i)$withinss)
    }
    plot(1:nc, wss, type="b", xlab="Number of Clusters",
    ylab="Within groups sum of squares")
}
wssplot(df[, -1:-2], nc=15)

# Por este método, também obteve-se 6 clusteres
set.seed(100)
kmeans_k6 <- kmeans(df[, -1:-2], 6, nstart=100)
# Valores médios por cluster
kmeans_k6$centers
# Silhueta dos clusteres
plot(silhouette(kmeans_k6$cluster,dist_food), border=NA)
# Plot results
pairs(df[, -1:-2], col=kmeans_k6$cluster)

# Método de PAM com 6 clusteres
set.seed(100)
pam_k6 <- pam(df[,-1:-2], 6)
# Valores médios por cluster
aggregate(df[, -1:-2], list(pam_k6$cluster), mean)
# Silhueta dos clusteres
plot(silhouette(pam_k6$cluster,dist_food), border=NA)
# Plot results
pairs(df[, -1:-2], col=pam_k6$cluster)

# Método de GMM com 6 clusteres
set.seed(100)
gmm_k6 <- Mclust(df[, -1:-2], G=6)
# Valores médios por cluster
summary(gmm_k6, parameters=T)
aggregate(df[, -1:-2], list(gmm_k6$classification), mean)
# Silhueta dos clusteres
plot(silhouette(gmm_k6$classification,dist_food), border=NA)
# Plot results
pairs(df[, -1:-2], col=gmm_k6$classification)
plot(gmm_k6, what="density") 
plot(gmm_k6, what="uncertainty")

# Validação do número de clusters
BIC <- mclustBIC(df[, -1:-2], G=seq(from=1,to=9, by=1), penalty=1)
plot(BIC)
BIC <- mclustBIC(df[, -1:-2], G=seq(from=1,to=15, by=1), penalty=1)
plot(BIC)
BIC <- mclustBIC(df[, -1:-2], G=seq(from=1,to=30, by=3), penalty=1)
plot(BIC)
BIC <- mclustBIC(df[, -1:-2], G=seq(from=1,to=15, by=1), penalty=100)
plot(BIC)

### Observações que alteraram de cluster
# ward vs. kmeans
table(groups_k5, kmeans_k6$cluster)
# ward vs. pam
table(groups_k5, pam_k6$clustering)
# ward vs. gmm
table(groups_k5, gmm_k6$classification)

df_clusters <- df
df_clusters$ward_group <- as.factor(groups_k5)
new_kmeans <- c()
for (i in kmeans_k6$cluster) {
    new_i <- ifelse(i==2, 1,
    ifelse(i==3, 2,
    ifelse(i==4, 3,
    ifelse(i==6, 4,
    ifelse(i==1, 5, 6)))))
    new_kmeans <- c(new_kmeans, new_i)
}
df_clusters$kmeans_group <- as.factor(new_kmeans)
df_clusters$pam_group <- as.factor(pam_k6$cluster)
new_gmm <- c()
for (i in gmm_k6$classification) {
    new_i <- ifelse(i==2, 1,
    ifelse(i==4, 2,
    ifelse(i==3, 3,
    ifelse(i==1, 4,
    ifelse(i==5, 5, 6)))))
    new_gmm <- c(new_gmm, new_i)
}
df_clusters$gmm_group <- as.factor(new_gmm)

for (i in c('ward_group', 'kmeans_group', 'pam_group', 'gmm_group')) {
    par(mai=c(1,2,1,1))
    barplot(
        apply(
            table(df_clusters$FoodGroup, df_clusters[, i]),
            1, function(x){x/sum(x,na.rm=T)}
        ),
        col=c('#6ECCAF', '#A8E890', '#F0E9D2', '#E5BA73', '#256D85', '#913678'),
        border="white", 
        space=0.04, 
        font.axis=1, 
        xlab="group",
        horiz=T,
        las=1,
        cex.names=.5,
        cex.axis=.5,
        main=i,
        legend=T,
        args.legend = list(x ='right', bty='n', inset=c(-0.15,0)),
        ylim=c(25,1)
    )
}

boxplot(Carboidratos ~ kmeans_group, df_clusters)
boxplot(Complexo_B ~ kmeans_group, df_clusters)
boxplot(Lipidos ~ kmeans_group, df_clusters)
boxplot(Minerais ~ kmeans_group, df_clusters)
boxplot(Vitaminas_Liposoluveis ~ kmeans_group, df_clusters)


# kmeans vs. pam
table(df_clusters$kmeans_group, df_clusters$pam_group)

# Grupo 6
# Ricos em carboidratos, poucos lipidos e mediano em vitaminas e minerais
# A maior parte dos vegetais e frutas pertencem a este grupo e alguma parte significativa de sopas, doces, bebidas pertencem também
# Sweets, Beverages e Dairy and Egg products parecem estranhos, então vamos verificá-los
# Sweets: Doces com base em fruta ou chocolate
# Beverages: Bebidas de baixo valor calórico ou baseado em frutas
# Dairy and Egg Products: Leites e yogurtes na versão desnatada e baixa gordura
df_clusters[df_clusters$kmeans_group==6 & df_clusters$FoodGroup=='Vegetables and Vegetable Products', ]
df_clusters[df_clusters$kmeans_group==6 & df_clusters$FoodGroup=='Beverages', ]
df_clusters[df_clusters$kmeans_group==6 & df_clusters$FoodGroup=='Fruits and Fruit Juices', ]

# Nome: Produtos Leves - Produtos naturais, fonte de energia e vitaminas com baixo teor de gordura

# Grupo 5
# Poucos carboidratos e vitaminas do complexo B, mediano em lipidos, vitaminas lipossolúveis e minerais
# A maior parte dos óleos e gorduras pertencem a este grupo e alguma parte significativa de sopas, doces, bebidas e legumes pertencem também
# Sweets, Legumes and Legume Products e Soups, Sauces, and Gravies parecem estranhos, então vamos verificá-los
# Sweets: Doces genéricos como balas, pirulitos, chicletes e etc.
# Legumes and Legume Products: Produtos de soja
# Soups, Sauces, and Gravies: Sopas genéricas
df_clusters[df_clusters$kmeans_group==5 & df_clusters$FoodGroup=='Sweets', ]
df_clusters[df_clusters$kmeans_group==5 & df_clusters$FoodGroup=='Soups, Sauces, and Gravies', ]
# Por ser difícil encontrar um padrão entre eles, veremos o grupo predominante também
df_clusters[df_clusters$kmeans_group==5 & df_clusters$FoodGroup=='Fats and Oils', ]
# Predominante por óleo de soja e de girassol
# Este cluster aparentemente não há um agrupamento muito óbivio e o principal motívo foi a ausência ou baixa quantidade de vitaminas do complexo B.

# Nome: Gorduras Pobres. Produtos com baixo teor de vitaminas e sem ganho nutricional.

#Grupo 4
#Alimentos ricos em vitaminas, tendo um valor mediano de carboidratos para determinado grupo de alimentos
df_clusters[df_clusters$kmeans_group==4 & df_clusters$FoodGroup=='Baby Foods', ]
df_clusters[df_clusters$kmeans_group==4 & df_clusters$FoodGroup=='Dairy and Egg Products', ]
df_clusters[df_clusters$kmeans_group==4 & df_clusters$FoodGroup=='Restaurant Foods', ]

#Nome: Ricos em vitaminas - Produtos mais relacionados com o equilibrio do organismo, evitar doenças, saudáveis ... 

#Grupo 3
#Alimentos ricos em lipidos e minerais
df_clusters[df_clusters$kmeans_group==3 & df_clusters$FoodGroup=='Baked Products', ]
df_clusters[df_clusters$kmeans_group==3 & df_clusters$FoodGroup=='Breakfast Cereals', ]
df_clusters[df_clusters$kmeans_group==3 & df_clusters$FoodGroup=='Legumes and Legume Products', ]
df_clusters[df_clusters$kmeans_group==3 & df_clusters$FoodGroup=='Meals, Entrees, and Side Dishes', ]

#Nome: Gorduras Ricas

# Grupo 2

df_clusters[df_clusters$kmeans_group==2 & df_clusters$FoodGroup=='Beef Products', ]
df_clusters[df_clusters$kmeans_group==2 & df_clusters$FoodGroup=='Finfish and Shellfish Products', ]
df_clusters[df_clusters$kmeans_group==2 & df_clusters$FoodGroup=='Lamb, Veal, and Game Products', ]
df_clusters[df_clusters$kmeans_group==2 & df_clusters$FoodGroup=='Pork Products', ]
df_clusters[df_clusters$kmeans_group==2 & df_clusters$FoodGroup=='Poultry Products', ]
df_clusters[df_clusters$kmeans_group==2 & df_clusters$FoodGroup=='Sausages and Luncheon Meats', ]

#Quantidades consideraveis de vitaminas do Complexo_B, medianas de Vitaminas_Liposoluveis e Minerais
#Quantidades bastante reduzidas de Carboidratros e de Lipidos
#Na distribuição do PCA, a proteina, que é predominante em produtos animais, ficou inserida no grupo das vitaminas do Complexo B, uma vez que são correlacionaveis
#Dai o numero bastante elevado de vitaminas do Complexo B

#A grande maioria de produtos de origem animal, sejam eles processados ou naturais pertence a este grupo, estando assim este grupo tambem presente em Fast Foods e nas comidas nativas Americanas/Indianas 


# Nome: Produtos Proteicos. Toda a carne de seres vivos esta presente neste cluster


# Grupo 1


df_clusters[df_clusters$kmeans_group==1 & df_clusters$FoodGroup=='Cereal Grains and Pasta', ]
df_clusters[df_clusters$kmeans_group==1 & df_clusters$FoodGroup=='Legumes and Legume Products', ]
df_clusters[df_clusters$kmeans_group==1 & df_clusters$FoodGroup=='Nut and Seed Products', ]
df_clusters[df_clusters$kmeans_group==1 & df_clusters$FoodGroup=='Snacks', ]
df_clusters[df_clusters$kmeans_group==1 & df_clusters$FoodGroup=='Spices and Herbs', ]


#Quantidade elevadas de Carboidratos, Lipidos e de vitaminas do complexo B
#Quantidades bastante reduzidas de Vitaminas Liposoluveis e de Minerais

#Alimentos como as massas e arroz, que são conhecidos como sendo carboidratos, estão presentes neste cluster
#Os frutos secos são ricos em carboidratos, mas tambem em proteina e em lipidos
#como observado no grupo 1, a proteina está correlacionada com as vitaminas do complexo B.
#Desta forma, as quantidades deste cluster fazem sentido, uma vez que os cereais, massas e frutos secos são os elementos predominantes do cluster
#Muitos legumes são ricos nestes 3 componentes, Carboidratos, Proteina e Lipidos, dai fazer sentido terem uma grande inserção neste cluster
#O grupo das ervas e especiarias é tambem ele rico em proteina, carboidratos e Lipidos.
#O grupo dos snacks tambem está predominante neste cluster, e tambem faz sentido, uma vez que por exemplo barras de cereais, batatas fritas entre outros snacks utilizam produtos para a sua produção que são pertencentes a este cluster


#verifica-se presenças de outros grupos neste cluster, sobretudo os Baked Products e Breakfast Cereal, uma vez que são grupos que tem base nos alimentos deste cluster
#por exemplo os cereais para o pequeno almoço ou o pão, são produtos que utilizam cereais, mas não são cereais puros
#É um cereal processado

# Nome: Produtos completos naturais
# Carboidratos, proteínas e lipídos (gorduras) são os principais substratos de energia para o nosso corpo
#Todos os componentes 

#https://www.ufrgs.br/laranjanacolher/2020/10/19/carboidratos-lipideos-e-proteinas-afinal-o-que-sao/


# Testes de implementação

kmeans_names <- c()
for (i in df_clusters$kmeans_group) {
    kmeans_names <- c(kmeans_names, switch(as.integer(i), 'Produtos Completos Naturais', 'Produtos Proteicos', 'Gorduras Ricas', 'Ricos em Vitaminas', 'Gorduras Pobres', 'Produtos Leves'))
}

df_results <- df_clusters[, 1:2]
df_results$cluster <- kmeans_names

set.seed(111)
df_1 <- df_results[df_results$cluster=='Produtos Completos Naturais',]
print('Produtos Completos Naturais')
df_1[sample(1:nrow(df_1), 10, replace=FALSE), ]

set.seed(222)
df_2 <- df_results[df_results$cluster=='Produtos Proteicos',]
print('Produtos Proteicos')
df_2[sample(1:nrow(df_2), 10, replace=FALSE), ]

set.seed(333)
df_3 <- df_results[df_results$cluster=='Gorduras Ricas',]
print('Gorduras Ricas')
df_3[sample(1:nrow(df_3), 10, replace=FALSE), ]

set.seed(444)
df_4 <- df_results[df_results$cluster=='Ricos em Vitaminas',]
print('Ricos em Vitaminas')
df_4[sample(1:nrow(df_4), 10, replace=FALSE), ]

set.seed(555)
df_5 <- df_results[df_results$cluster=='Gorduras Pobres',]
print('Gorduras Pobres')
df_5[sample(1:nrow(df_5), 10, replace=FALSE), ]

set.seed(666)
df_6 <- df_results[df_results$cluster=='Produtos Leves',]
print('Produtos Leves')
df_6[sample(1:nrow(df_6), 10, replace=FALSE), ]

dist_food1 <- as.matrix(dist(df, diag=T, upper=T))
set.seed(100)
for (i in df[sample(1:nrow(df), 4, replace=FALSE), 'Descrip']) {
    print(i)
    View(df[tail(head(order(dist_food1[df$Descrip==i, ]), 11), 10), ])
}
