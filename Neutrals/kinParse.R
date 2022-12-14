# Kinetic Parser ############################################################

nbReac = 0
reactants = products = params = type = orig =  list()

## Data from Dobrijevic/Loison Google docs
## https://docs.google.com/spreadsheets/d/1i48b1zjYQX85uhssyWkR1CI-D5OBz6wrEs1YsU8X53U/edit#gid=13
filename = 'Titan - Réactions bimoléculaires.csv'
scheme  = read.csv(
  file = paste0(sourceDir, filename),
  header = FALSE,
  stringsAsFactors = FALSE
)
comments = scheme[, ncol(scheme)]
scheme  = t(apply(scheme, 1, function(x) gsub(" ", "", x)))
for (i in 1:nrow(scheme)) {
  if(substr(scheme[i,1],1,1)=='#') next
  nbReac = nbReac + 1
  terms = scheme[i, 1:3]
  reactants[[nbReac]] = terms[!is.na(terms) & terms != ""]
  terms = scheme[i, 4:8]
  products[[nbReac]]  = terms[!is.na(terms) & terms != ""]
  terms = scheme[i, 9:13]
  params[[nbReac]]    = terms[!is.na(terms) & terms != ""]
  params[[nbReac]][6] = 'kooij'
  type[[nbReac]]      = 'kooij'
  orig[[nbReac]]      = filename
}

filename = 'Titan - Réactions trimoléculaires.csv'
scheme  = read.csv(
  file = paste0(sourceDir, filename),
  header = FALSE,
  stringsAsFactors = FALSE
)
comments = c(comments, scheme[, ncol(scheme)])
scheme  = t(apply(scheme, 1, function(x) gsub(" ", "", x)))
for (i in 1:nrow(scheme)) {
  if(substr(scheme[i,1],1,1)=='#') next
  nbReac = nbReac + 1
  terms = scheme[i, 1:3]
  reactants[[nbReac]] = terms[!is.na(terms) & terms != ""]
  terms = scheme[i, 4:8]
  products[[nbReac]]  = terms[!is.na(terms) & terms != ""]
  terms = scheme[i, 9:24]
  params[[nbReac]]    = terms[!is.na(terms) & terms != ""]
  params[[nbReac]][17]= 'assocMD'
  type[[nbReac]]      = 'assocMD'
  orig[[nbReac]]      = filename
}

## Additional bimolecular data from misc. sources

filename = 'bimol_supp.csv'
scheme  = read.csv(
  file = paste0(sourceDir, filename),
  header = FALSE,
  stringsAsFactors = FALSE
)
comments = c(comments, scheme[, ncol(scheme)])
scheme  = t(apply(scheme, 1, function(x) gsub(" ", "", x)))
for (i in 1:nrow(scheme)) {
  if(substr(scheme[i,1],1,1)=='#') next
  nbReac = nbReac + 1
  terms = scheme[i, 1:3]
  reactants[[nbReac]] = terms[!is.na(terms) & terms != ""]
  terms = scheme[i, 4:8]
  products[[nbReac]]  = terms[!is.na(terms) & terms != ""]
  terms = scheme[i, 9:13]
  params[[nbReac]]    = terms[!is.na(terms) & terms != ""]
  params[[nbReac]][6] = 'kooij'
  type[[nbReac]]      = 'kooij'
  orig[[nbReac]]      = filename
}

## Additional trimolecular data from Vuitton2019
## (With specific parameterization)

filename = 'trimol_VV.csv'
scheme  = read.csv(
  file = paste0(sourceDir, filename),
  header = FALSE,
  stringsAsFactors = FALSE
)
comments = c(comments, scheme[, ncol(scheme)])
scheme  = t(apply(scheme, 1, function(x) gsub(" ", "", x)))
for (i in 1:nrow(scheme)) {
  if(substr(scheme[i,1],1,1)=='#') next
  nbReac = nbReac + 1
  terms = scheme[i, 1:3]
  reactants[[nbReac]] = terms[!is.na(terms) & terms != ""]
  terms = scheme[i, 4:8]
  products[[nbReac]]  = terms[!is.na(terms) & terms != ""]
  terms = scheme[i, 9:24]
  params[[nbReac]]    = terms[!is.na(terms) & terms != ""]
  params[[nbReac]][17]= 'assocVV'
  type[[nbReac]]      = 'assocVV'
  orig[[nbReac]]      = filename
}

# filename = 'Titan - Réactions bimol_trimol_association.csv'
# scheme  = read.csv(file=paste0(sourceDir,filename),header=FALSE)
# scheme  = t(apply(scheme,1,function(x) gsub(" ","",x)))
# for (i in 1:nrow(scheme)) {
# if(substr(scheme[i,1],1,1)=='#') next
#   nbReac = nbReac + 1
#   terms=scheme[i,1:3]
#   reactants[[nbReac]] = terms[!is.na(terms) & terms!="" & terms!="HV"]
#   terms=scheme[i,4:8]
#   products[[nbReac]]  = terms[!is.na(terms) & terms!="" & terms!="HV"]
#   terms=scheme[i,9:23]
#   params[[nbReac]]    = terms[!is.na(terms) & terms!=""]
#   params[[nbReac]][16] = 'assoc'
#   type[[nbReac]]      = 'assoc'
#   orig[[nbReac]]      = filename
# }
