# Startup message
.onAttach <- function(lib, pkg) {
  packageStartupMessage(
    paste0(
      '\n >===================================================================================<',
      '\n "DRrequiredAgeing" is developed by International Mouse Phenotyping Consortium (IMPC) ',
      '\n More details  https://www.mousephenotype.org/                                        ',
      '\n Contact us hamedhm@ebi.ac.uk                                                         ',
      '\n Source code and issues:  https://git.io/fj9WP                                        ',
      '\n >===================================================================================<'
    ),
    domain = NULL,
    appendLF = TRUE
  )
}

noSexes0 = function(obj){
  if(is.null((obj)))
    return(1)
  if('Sex' %in% colnames(obj@datasetPL)){
    r = nlevels(obj@datasetPL$Sex)
  }else{
    r = 1
  }
  return(r)
}

# Behaviour parameters
BehviourParamters = function(x, ...) {
  r =  c(
    "IMPC_ACS_006_001",
    "IMPC_ACS_033_001",
    "IMPC_ACS_034_001",
    "IMPC_ACS_035_001",
    "IMPC_ACS_036_001",
    "IMPC_ACS_037_001",
    "IMPC_CSD_026_001",
    "MGP_CSD_064_001",
    "IMPC_CSD_029_001",
    "IMPC_CSD_030_001",
    "MGP_CSD_007_001",
    "IMPC_CSD_032_001",
    "IMPC_CSD_033_001",
    "MGP_CSD_005_001",
    "IMPC_CSD_036_001",
    "MGP_CSD_008_001",
    "IMPC_CSD_037_001",
    "IMPC_CSD_038_001",
    "MGP_CSD_003_001",
    "IMPC_CSD_039_001",
    "MGP_CSD_002_001",
    "IMPC_CSD_077_001",
    "MGP_CSD_011_001",
    "IMPC_CSD_078_001",
    "IMPC_CSD_079_001",
    "IMPC_CSD_080_001",
    "MGP_CSD_012_001",
    "IMPC_FEA_003_001",
    "IMPC_FEA_009_001",
    "IMPC_FEA_015_001",
    "IMPC_FEA_021_001",
    "IMPC_GRS_008_001",
    "IMPC_GRS_009_001",
    "JAX_HBD_001_001",
    "IMPC_LDT_001_001",
    "JAX_LDT_001_001",
    "IMPC_LDT_006_001",
    "IMPC_LDT_007_001",
    "JAX_LDT_008_001",
    "JAX_LDT_006_001",
    "IMPC_OFD_007_001",
    "IMPC_OFD_009_001",
    "IMPC_OFD_014_001",
    "IMPC_OFD_020_001",
    "IMPC_OFD_022_001",
    "JAX_ROT_002_001",
    "ICS_ROT_002_001",
    "HMGU_ROT_002_001",
    "JAX_ROT_004_001",
    "JAX_TLS_001_001"
  )
  return(r)
}



# Objects that must be removed at each iterations
ObjectsThatMustBeRemovedInEachIteration = function(x = NULL, ...) {
  if (is.null(x)) {
    rm0(
      c(
        'note',
        'minSampRequired',
        'n3.4',
        'n3.5',
        'lvls',
        'counts',
        'n3.5.1',
        'note',
        'SubSubDirOrdFileName',
        'outDir',
        'outpfile',
        'AllLevels',
        'GenePageURL',
        'BodyWeightCurvURL',
        'depVariable',
        'depVar',
        'n3.5_summary',
        'n3.5.1_F_list',
        'n3.5.1_v_list',
        'isException',
        'ExtraCols',
        'a',
        'PhenListSpecIds',
        'OrgSpecIds',
        'a_summary_before_concurrent',
        'a_phenlist_concurrent_summary',
        'method',
        'args',
        'aTmp',
        'c.ww0',
        'c.ww.vec',
        'outP',
        'optFail'
      )
    )
  } else{
    rm0(x)
  }
}

getNotEmptyValue = function(x, head = Inf) {
  ux    = unique(x)
  neux  = which(ux != ''     &
                  !is.na(ux)   &
                  !is.null(ux))
  if (length(neux) > 0) {
    r = head(ux [neux], head)
  } else{
    r = 'NotExist'
  }
  return(as.character(r))
}

sortDataset = function(x = NULL, BatchCol = NULL) {
  if (is.null(x)        ||
      is.null(BatchCol) ||
      length(BatchCol) < 1 ||
      any(dim(x) < 1))
    return(x)
  if (BatchCol %in% names(x)) {
    message0('Sorting dataset based on ', BatchCol, '...')
    x = x[order(Date2Integer(x[, BatchCol])),]
  }
  return(x)
}

SimplesortDataset = function(x = NULL, cols = NULL) {
  if (is.null(x)        ||
      is.null(cols) ||
      length(cols) < 1 ||
      any(dim(x) < 1))
    return(x)
  for (col in cols) {
    if (col %in% names(x) && sum(!is.na(x[, col])) > 0) {
      message0('Sorting dataset based on `', col, '`...')
      x = x[order(x[, col]), ]
    } else{
      message0('Column `', col, ' does not exist or all NA')
    }
  }
  return(x)
}
# Only works for the data from Solr
# Automatic select the corresponding column for the type of data (only for category and data_point)
getResponseColumn = function(x, activate = TRUE) {
  if (activate) {
    if (length(na.omit(x)) > 0) {
      lbl      = unique(na.omit(x))
      column   = paste(lbl, collapse = '-', sep = '-')
      accepted = TRUE

      if (length(lbl) > 1)
        accepted = FALSE

      if (lbl %in% 'categorical')
        column = 'category'
      else if (lbl %in% 'unidimensional')
        column = 'data_point'
      else if (lbl %in% 'time_series')
        column = 'data_point'
      else if (lbl %in% 'text')
        column =  'text_value'
      else
        accepted = FALSE
    } else{
      accepted = FALSE
      column = lbl = 'Unknown'
    }
  } else{
    lbl    = NA
    column = lbl = NA
    accepted = FALSE
  }
  return(list(
    column   = column,
    accepted = accepted,
    lbl      = lbl
  ))
}

removeOffweightsFromPLobject = function (plObject = NULL,
                                         threshold = sqrt(.Machine$double.eps) * 10) {
  if (is.null(plObject) || is.null(threshold))
    return(plObject)

  if (!'AllModelWeights' %in% names(plObject@datasetPL))
    return(plObject)

  message0('Removing off threshold weights from PL object for windowing only. Threshold: ', threshold)
  plObject@datasetPL = plObject@datasetPL[plObject@datasetPL$AllModelWeights >
                                            threshold,]

  return(plObject)

}

MeanVarOverTime = function(mm, tt, data = phenlistObject@datasetPL) {
  if (length(tt) < 1 || length(mm) < 1)
    return(1)
  v = sapply(unique(tt[mm]), function(i) {
    ind      = 1:length(tt)
    CriTeria = (tt %in% i) & (ind %in% mm)
    if (sum(CriTeria) > 1) {
      sd0(data[CriTeria], na.rm = TRUE)
    } else{
      NA
    }
  })
  if (length(v[!is.na(v)]) > 0) {
    return(mean(v, na.rm = TRUE))
  } else{
    return(1)
  }
}

IsInBlackListCategories = function(x, len = 1, blackList = NULL) {
  note = NULL
  if (!is.null(blackList) && !is.numeric(x) && nlevels(x) == len) {
    r    = levels(x) %in% blackList
    note = 'The categorical variable has only one level that is found in the skip list'
  } else{
    r    = FALSE
  }
  return(list(result = r, note = note))
}


#
MergeLevels = function(x                         ,
                       listOfLevelMaps           ,
                       parameter_stable_id       ,
                       AllowedParametersList     ,
                       report         = FALSE    ,
                       reportFileName = 'MergedCategoriesParameterStableIds.txt',
                       fileColumnSep  = '\t') {
  note = NULL
  if (all(parameter_stable_id %in% AllowedParametersList) &&
      !is.null(x)    &&
      !is.numeric(x) &&
      length(x) > 0  &&
      nlevels(as.factor(x)) > 0) {
    x  = as.factor(x)
    nl = nlevels(x)
    l  = levels(x)
    for (i in 1:nl) {
      if (l[i] %in% names(listOfLevelMaps)) {
        newCat = as.character(listOfLevelMaps[l[i] == names(listOfLevelMaps)])
        levels(x)[which(levels(x) == l[i])] = ifelse(newCat != 'NA', newCat, NA)
        note = c(note, paste0(as.character(l[i]), '-->', newCat))
      }
    }
  }
  if (report && !is.null(note)) {
    write(
      paste(
        Sys.time()                   ,
        parameter_stable_id          ,
        paste(note, collapse = ', ') ,
        sep      = fileColumnSep     ,
        collapse = fileColumnSep
      ),
      file = reportFileName,
      append = ifelse(file.exists(reportFileName),
                      TRUE                       ,
                      FALSE)                     ,
      ncolumns = 3500
    )
  }
  return(list(x = x, note = if (is.null(note)) {
    note
  } else{
    paste(note, collapse = ', ')
  }))
}


# Gene page URL
GenePageURL = function(obj) {
  if (any(dim(obj) == 0) || is.null(obj)) {
    return ('Empty dataset has no URL!')
  }
  activation  = c(
    'gene_accession_id',
    'allele_accession_id',
    'zygosity',
    'parameter_stable_id',
    'pipeline_stable_id',
    'phenotyping_center',
    'metadata_group'    ,
    'biological_sample_group'
  ) %in% names(obj)
  ##########
  if (prod(activation) > 0) {
    obj2 = subset(obj, obj$biological_sample_group %in% 'experimental')
    r = paste(
      'https://www.mousephenotype.org/data/charts?'  ,
      'accession='                                   ,
      getNotEmptyValue(obj2$gene_accession_id, 1)    ,
      '&allele_accession_id='                        ,
      getNotEmptyValue(obj2$allele_accession_id, 1)  ,
      '&zygosity='                                   ,
      getNotEmptyValue(obj2$zygosity, 1)             ,
      '&parameter_stable_id='                        ,
      getNotEmptyValue(obj2$parameter_stable_id, 1)  ,
      '&pipeline_stable_id='                         ,
      getNotEmptyValue(obj2$pipeline_stable_id, 1)   ,
      '&phenotyping_center='                         ,
      getNotEmptyValue(obj2$phenotyping_center, 1)   ,
      '&metadata_group='                             ,
      getNotEmptyValue(obj2$metadata_group, 1)       ,
      sep = ''                                       ,
      collapse = ''
    )
    r = URLencode(r, repeated = TRUE)
  } else{
    r = NA
  }
  return(r)
}

# BWT URL
BodyWeightCurvURL = function(obj) {
  if (any(dim(obj) == 0) || is.null(obj)) {
    return ('Empty dataset has no BWT URL!')
  }
  # allele_accession_id to make sure that we are safe for external application of the R package
  activation = c('allele_accession_id', 'gene_accession_id','biological_sample_group') %in% names(obj)
  obj2 = subset(obj, obj$biological_sample_group %in% 'experimental')
  if (prod(activation) > 0) {
    r = paste(
      'https://www.mousephenotype.org/data/charts?'    ,
      'accession='                                     ,
      getNotEmptyValue(obj2$gene_accession_id, 1)       ,
      '&parameter_stable_id='                          ,
      'IMPC_BWT_008_001'                               ,
      '&chart_type=TIME_SERIES_LINE'                   ,
      sep = ''                                         ,
      collapse = ''
    )
    r = URLencode(r, repeated = TRUE)
  } else{
    r = NA
  }
  return(r)
}

#### create dir
dir.create0 = function(x, ...) {
  if (!dir.exists(x))
    dir.create(x, ...)
  return(x)
}
###
LogicalToNote = function(x) {
  if (x) {
    r = 'Enabled'
  } else{
    r = 'Disabled'
  }
  return(r)
}
#### Concurrent sample selection
concurrentContSelect = function(activate = TRUE,
                                PhenListObj,
                                minSmp  = 6,
                                depVar  = 'data_point',
                                control = 'control'   ,
                                sexColName      = 'Sex'     ,
                                GenotypeColName = 'Genotype',
                                BatchColName    = 'Batch'   ,
                                minSampRequired,

                                plot = FALSE) {
  note = NULL
  PhenListObjBck = PhenListObj
  # Old keys
  # cun_con_sel | cun_con_selection | concurrent_control_selection
  # Part 1 messages
  if (is.numeric(PhenListObj@datasetPL[, depVar])) {
    m  =  (PhenListObj@datasetPL[PhenListObj@datasetPL[, GenotypeColName] != control,])
    m1 =  (PhenListObj@datasetPL[PhenListObj@datasetPL[, BatchColName] %in% m[, BatchColName],])
    if (length(unique(m[, BatchColName])) == 1) {
      tb = table(m1[, sexColName], m1[, GenotypeColName])
      if (control %in% colnames(tb)) {
        if (min0(tb[, which(colnames(tb) == control)]) >= minSmp) {
          if (activate) {
            message0 ('Concurrent control selection in progress ...')
            PhenListObj@datasetPL = m1
            note$'Concurrent control selection' = paste0(LogicalToNote(activate),
                                                         ". Single Batch in use")
          } else{
            note$'Concurrent control selection' = paste0(LogicalToNote(activate),
                                                         '. Single Batch possible')
          }
        } else{
          note$'Concurrent control selection' = paste0(
            LogicalToNote(activate),
            '. There are not enough controls in the dataset at the mutant dates (min = ',
            minSmp,
            ' but there are ',
            min0(tb[, which(colnames(tb) == control)]),
            ')'
          )
        }
      } else{
        note$'Concurrent control selection' = paste0(
          LogicalToNote(activate),
          ". The is no control in the dataset at the same day as mutants DOE"
        )
      }
    } else{
      note$'Concurrent control selection' = paste0(
        LogicalToNote(activate),
        ". Concurrent control selection not possible. Mutants are scattered on different dates"
      )
    }
  } else {
    tb  = table(PhenListObj@datasetPL[, BatchColName])
    if (length(tb) > 1) {
      note$'Concurrent control selection' = paste0(LogicalToNote(activate),
                                                   '. Multi batch found in data (categorical data)')
    } else if (length(tb) == 1) {
      note$'Concurrent control selection' = paste0(LogicalToNote(activate),
                                                   '. Original mutants in single batch (categorical data)')
    } else{
      note$'Concurrent control selection' = paste0(LogicalToNote(activate),
                                                   '. But batch not clear. Please check data')
    }
  }
  # Part 2 extra checks
  if (activate && is.numeric(PhenListObj@datasetPL[, depVar])) {
    ### Check whether there is any variation in data
    # No worries it only applies to MM framework
    checkNoZeroVar = RemoveZerovarCategories(
      x = PhenListObj@datasetPL,
      depVar = depVar,
      sex = sexColName,
      genotype = GenotypeColName
    )
    #PhenListObj@datasetPL = checkNoZeroVar$x
    PhenListObj = OpenStats::OpenStatsList(
      dataset = checkNoZeroVar$x,
      testGenotype = PhenListObjBck@testGenotype,
      refGenotype = PhenListObjBck@refGenotype
    )
    note =  c(note, checkNoZeroVar$note)


    conditions =  (
      !is.null(PhenListObj@datasetPL) &&
        # data.frame is not zero
        min0(dim(PhenListObj@datasetPL)) > 0 &&
        # is it  really exist!
        length(unique(PhenListObj@datasetPL[, GenotypeColName])) > 1 &&
        # include mut and cont
        min0(table(PhenListObj@datasetPL[, GenotypeColName])) >= minSampRequired &&
        # include at least 4/2 of each genotype
        checkGenInCol(
          PhenListObj@datasetPL,
          sex = sexColName,
          genotype = GenotypeColName
        ) &&
        # each sex and genotype
        # there must be variation in data
        NonZeroVariation(PhenListObj@datasetPL[, depVar])
    )

    if (!conditions) {
      PhenListObj = PhenListObjBck
      note$'Concurrent control selection error' = paste0(
        LogicalToNote(activate),
        'Cuncurrent control selection would cause failure in model then DISABLED'
      )
    }
    ####################
  }
  if (plot)
    plot(PhenListObj, depVariable = depVar)

  return (list(obj = PhenListObj, note = note))
}



# Variation in data
NonZeroVariation = function(x) {
  x = na.omit(x)
  result = ifelse(length(x) > 0 &&
                    length(unique(x)) > 1, TRUE, FALSE)
  return(result)
}

# Is this an ABR screen
is.ABR = function(x) {
  length(grep("_ABR", x, fixed = TRUE)) > 0
}

# Is this a 3I screen
is.3i = function(x) {
  length(
    grep(
      "MGP_BCI|MGP_PBI|MGP_ANA|MGP_CTL|MGP_EEI|MGP_BMI|MGP_MLN|MGP_IMM|MGP_BHP|MGP_MIC",
      x,
      fixed = FALSE
    )
  ) > 0

}

# Exception map
getException =  function(parameter,
                         exceptionMap = NULL)
{
  matched001 = sapply(names(exceptionMap), grepl, parameter)
  if (sum(matched001) > 0) {
    result = FALSE
  } else{
    result = TRUE
  }
  return(result)
}

# Method map
getMethodi =  function(var,
                       type = 'numeric',
                       methodMap = NULL)
{
  matched001 = sapply(names(methodMap), grepl, var)
  if (sum(matched001) > 0) {
    if (sum(matched001) > 1) {
      methodMapReduced = methodMap[matched001]
      mappedPatternLengths = nchar(names(methodMapReduced))
      method =
        unlist(methodMapReduced[which(mappedPatternLengths ==  max(mappedPatternLengths))])
    } else {
      method = unlist(methodMap[matched001])
    }
  } else{
    if (type == 'numeric') {
      method = 'MM'
    } else{
      method = 'FE'
    }
  }
  return(method)
}

# Equation map
getEquation =  function(var,
                        equationMap = NULL)
{
  matched001 = sapply(names(equationMap), grepl, var)
  if (sum(matched001) > 0) {
    if (sum(matched001) > 1) {
      equationMapReduced = equationMap[matched001]
      mappedPatternLengths = nchar(names(equationMapReduced))
      equation =
        unlist(equationMapReduced[which(mappedPatternLengths ==  max(mappedPatternLengths))])
    } else {
      equation = unlist(equationMap[matched001])
    }
  } else{
    equation = 'withWeight'
  }
  return(equation)
}

local =  function(x=NULL){
  r= system.file("extdata", package = "DRrequiredAgeing")
  return(r)
}


getEarlyAdultsFromParameterStableIds = function(LA_parameter_stable_id,
                                                map         = read.csv(file = file.path(local(),
                                                                                        'EA2LA_parameter_mappings_2019-09-24.csv')),
                                                LA_data     = NULL,
                                                solrBaseURL = 'http://hx-noah-74-10:8090') {
  message0('EA-LA binding in progress ...')
  if (!any(c('LA_parameter', 'EA_parameter') %in% names(map))) {
    stop(
      'Check the LA/EA mapping column names! It must have at least two columns: EA_parameter and LA_parameter '
    )
  }
  plist = as.character(unique(map$EA_parameter[map$LA_parameter %in% LA_parameter_stable_id]))

  if (length(plist) < 1) {
    message0('Skip this parameter. No match for the EA corresponded to this parameter, ',
             LA_parameter_stable_id)
    return(NULL)
  }
  q = URLencode(
    paste0(
      solrBaseURL,
      '/solr/experiment/select?q=*:*',
      '&fl=',
      paste(sort(requiredDataColumns(0)), collapse = ','),
      '&fq=datasource_name:IMPC*&fq=observation_type:(',
      paste(
        '"',
        unique(LA_data$observation_type),
        '"',
        collapse = ' OR ',
        sep = ''
      ),
      ')&fq=parameter_stable_id:(',
      paste('"', plist, '"', collapse = ' OR ', sep = ''),
      ')&rows=500000000&wt=csv'
    )
  )
  message0(
    'EA data found for the parameter: ',
    LA_parameter_stable_id,
    '\n\t => ',
    paste(plist, sep = ', ', collapse = ', '),
    '\n\t => ',
    q
  )
  ######################################
  EA_data = read.csv(q)
  #EA_data = subset(EA_data,EA_data$colony_id %in% LA_data$colony_id )
  if (any(dim(LA_data) < 1) || any(dim(EA_data) < 1)) {
    message0('Skip this parameter. No EA data to match')
    return(NULL)
  } else{
    message0(nrow(EA_data), ' data point added to the dataset')
  }
  ######################################
  requireNamespace("plyr")
  message0('Binding the EA and LA data ...')
  df      = plyr::rbind.fill(EA_data, LA_data)
  ######################################
  df      = df[df$colony_id %in%
                 unique(c(LA_data$colony_id, EA_data$colony_id[EA_data$biological_sample_group %in% 'control'])), ]
  ######################################
  message0('Remove any possible duplicate from the combined EA & LA dataset.')
  df = df[!duplicated(df),]
  ######################################
  plist = unique(c(plist,LA_parameter_stable_id))
  ######################################
  f = function(str,
               cut = 1,
               split = '_',
               partIndex = NULL) {
    r = sapply(str, function(x) {
      s = unlist(strsplit(x, split = split))
      if (is.null(partIndex)) {
        r0 = paste(head(s, length(s) - cut),
                   sep = '_',
                   collapse = '_')
      } else{
        r0 = paste(s[partIndex],
                   sep = '_',
                   collapse = '_')
      }
      return(r0)
    })
    return(unique(r))
  }
  ######################################
  df$LifeStage = ifelse(grepl(pattern = '(IP_)|(LA_)', df$parameter_stable_id),
                        'Late',
                        'Early')
  df$LifeStage = as.factor(df$LifeStage)
  ######################################
  df$parameter_stable_id_renamed = df$parameter_stable_id
  df$parameter_stable_id = gsub(
    pattern     = paste0('(', paste(f(plist, cut = 0), collapse  = ')|('), ')'),
    replacement = f(plist[1], cut = 0),
    x           = df$parameter_stable_id
  )

  df$procedure_stable_id_renamed = df$procedure_stable_id
  df$procedure_stable_id         = gsub(
    pattern     = paste0('(', paste(f(plist, cut = 2), collapse  = ')|('), ')'),
    replacement =  df$procedure_stable_id[grepl(pattern = f(str = plist[1],
                                                            cut = 2),
                                                x = df$procedure_stable_id)][1],
    x            =  df$procedure_stable_id
  )

  df$procedure_group_renames = df$procedure_group
  df$procedure_group         = gsub(
    pattern = paste0('(', paste(f(plist, cut = 2), collapse  = ')|('), ')'),
    replacement = f(plist[1], cut = 2),
    df$procedure_group
  )
  if (sum(df$age_in_weeks <= 16 & df$LifeStage %in% 'Late')) {
    message0('LA must have age_in_weeks > 16')
    write(
      x = paste(Sys.Date(), LA_parameter_stable_id, sep = '\t'),
      file = 'LAWithAgeLessThan16Weeks.txt',
      append = TRUE
    )
  }
  # Return them into factors ...
  df$parameter_stable_id  = as.factor(df$parameter_stable_id)
  df$procedure_stable_id  = as.factor(df$procedure_stable_id)
  df$procedure_group      = as.factor(df$procedure_group )
  return(df)
}


# Read config files
readConf = function(file, path = NULL, ...) {
  if (is.null(path))
    path = system.file("extdata", package = "DRrequiredAgeing")
  message0('Reading the config file:\n\t ~> ',file)

  r   = read.dcf(file.path(path, file), all = TRUE,  ...)
  return(r)
}
# Read config files
readFile = function(file,
                    path = NULL,
                    FUN = function(x, ...) {
                      as.vector(read.table(file = x)[, 1])
                    },
                    ...) {
  if (is.null(path))
    path = system.file("extdata", package = "DRrequiredAgeing")
  message0('Reading the config file:\n\t ~> ',file)
  r   = FUN(file.path(path, file),  ...)
  return(r)
}

# Replace blank elements of a vector with user define character
RepBlank = function(x, match = '', replace = 'control') {
  x[x %in% match] = replace
  return(x)
}

# To prevent:
#                   female male
# control           4       0
# experimental      0       5
checkGenInCol = function(x, sex = 'sex', genotype = 'biological_sample_group') {
  x  = droplevels0(x)
  r = FALSE
  # not null
  if (min0(dim(x)) > 0) {
    tb = table(x[, genotype], x[, sex])
    if (min0(dim(tb)) > 0) {
      if (dim(tb)[2] > 1)
        r = (prod(tb[, 1]) + prod(tb[, 2])) > 0
      else
        r = (prod(tb[, 1])) > 0
    }
  }
  return(r)
}




### check for null or errorneous object
NullOrError = function(obj) {
  if (is.null(obj) ||
      class(obj) %in% c("simpleWarning", "warning", "condition")) {
    return(TRUE)
  } else{
    return(FALSE)
  }
}


# Generate Batch file only for EBI
BatchGenerator = function(file                       ,
                          dir = NULL                 ,
                          procedure = NULL           ,
                          parameter = NULL           ,
                          center = NULL              ,
                          cpu = 1                    ,
                          memory = "8G"              ,
                          time = "10:00:00"          ,
                          extraBatchParameters = NULL) {
  dirOut = file.path(dir, 'ClusterOut')
  dirErr = file.path(dir, 'ClusterErr')
  dir.create0(dirOut)
  dir.create0(dirErr)
  
  logfile_basename <- basename(file)
  oname = file.path(dirOut, logfile_basename)
  ename = file.path(dirErr, logfile_basename)

  ro = paste(' -o ', paste0('"', oname, '.ClusterOut', '"'), sep = '')
  re = paste(' -e ', paste0('"', ename, '.ClusterErr', '"'), sep = '')
  rf = paste(
    "sbatch --job-name=impc_stats_pipeline_job --mem=", memory,
    " --time=", time,
    extraBatchParameters  ,
    ' --cpus-per-task='                ,
    cpu                   ,
    ' '                   ,
    re                    ,
    ' '                   ,
    ro                    ,
    " --wrap=\'Rscript function.R ",
    paste('"', file, '"', sep = ''),
    "\'",
    sep = ''
  )
  return(rf)
}


# remove special characters from a string
RemoveSpecialChars = function(x,
                              what = '[^0-9A-Za-z]',
                              replaceBy = '_',
                              message = FALSE) {
  what = gsub(
    pattern = ']',
    x = what,
    replacement = paste0(replaceBy, ']')
  )
  if (message)
    message0('pattern: ',
             what ,
             '; replaced by ',
             replaceBy)
  r = gsub(what, replaceBy , x , ignore.case = TRUE)
  ###
  if (any(nchar(r) < 1)) {
    RN = RandomRegardSeed(1, stringOutput = TRUE, round = 8)
    r[nchar(r) < 1] = paste('no_name',
                            RN,
                            sep = '_',
                            collapse = '_')
  }
  return(r)
}


# Create sub directory from a url if does not exist
CreateSubDirIfNotExist = function(x) {
  base = dirname(x)
  dir.create0(base)
  return(x)
}

# Unique and not NULL
UniqueAndNNull = function(x,
                          removeSpecials = TRUE,
                          removeNewLine  = TRUE,
                          newlineRepChar = ' ' ,
                          collapse = '~') {
  if (length(na.omit(x)) > 0)
    x = na.omit(x)
  if (length(x[!is.null(x)]) > 0)
    x = x[!is.null(x)]
  if (length(x[!x %in% '']) > 0)
    x = x[!x %in% '']
  if (length(x) > 0) {
    x = as.character(unique(x))
    x = paste0(x, collapse = collapse)
    if(removeNewLine)
      x = gsub(pattern = '\n',
               x = x,
               replacement = newlineRepChar)
    if (removeSpecials)
      x = RemoveSpecialChars(x = x, replaceBy = ' ')
  } else{
    x = 'NotSpecified'
  }
  return(x)
}

# remove useless characters
removeNAPipe = function(x,
                        pattern = 'NA',
                        replacement = '') {
  x = na.omit(x)
  if (length(x) > 0)
    x = gsub(
      pattern     = pattern,
      replacement = replacement,
      x           = x,
      fixed       = FALSE
    )
  return(x)
}


# as.numeric with quote!
as.numeric0 = function(x, ...) {
  r = suppressWarnings(as.numeric(x, ...))
  return(r)
}

# is.numeric with quote!
is.numeric0 = function(x) {
  r0  = is.na(as.numeric0(x))
  r1  = is.na(x) | (x %in% c('na', 'NA', 'TRUE', 'FALSE'))
  r   = r0 | r1
  return(!r)
}

## solve "}"]  and ["
RemoveTwoSpecialFromOutput  =  function(x) {
  ind = grepl(pattern = '{',
              x = x,
              fixed = TRUE)
  for (i in 1:length(x)) {
    if (!ind[i]) {
      x[i] = paste(
        ifelse(length(x[i]) > 1, '[', ''),
        ifelse(is.numeric0(x[i]), as.numeric0(x[i]), paste0('"', x[i], '"')),
        ifelse(length(x[i]) > 1, ']', ''),
        sep = ''
      )
    } else{
      x[i] = paste('{', x[i], '}', sep = '')
    }
  }
  return(x)
}

### removeObject if exist
rm0 = function(x, silent = FALSE) {
  x      = unique(x)
  obList = ls(all.names = TRUE, pos = ".GlobalEnv")
  InObj  = obList %in% x
  if (sum(InObj) > 0) {
    rm(list = as.character(obList[InObj]), pos = ".GlobalEnv")
    if (!silent)
      message0(
        'Cleaning the environment. Removed objects from memory:\n\t ~~~~>>>> ',
        paste(obList[InObj], collapse = ', ')
      )
  } else{
    if (!silent)
      message0('Cleaning the environment. No object removed from the memory')
  }
}



# Generate random numbers even when seed is set (max 10^6 numbers)
RandomRegardSeed = function(n = 1,
                            max = 10 ^ 9,
                            decimal = TRUE,
                            round = 5,
                            stringOutput = FALSE,
                            what = '[^0-9A-Za-z]',
                            replaceBy = '_') {
  r = runif(n, 1, as.numeric(Sys.time()) * 10000) %% max
  if (decimal)
    r = r %% 1
  r = round(r, digits = round)
  if (stringOutput) {
    r = RemoveSpecialChars(as.character(r), what = what, replaceBy = replaceBy)
  }
  return(r)
}

complete.cases0 = function(x, ...) {
  if (!is.null(x) && length(x) > 0) {
    r = complete.cases(x)
    if (is.null(r))
      r = NULL
  } else{
    message0('Null input in the "complete.case". Null retured!')
    r = NULL
  }
  return(r)
}
#From PhenStat
columnLevelsVariationRadio = function (dataset,
                                       columnName,
                                       genotypeCol = 'biological_sample_group',
                                       sexCol = 'sex')
{
  message0(paste('depVar = ', columnName))
  if (is.numeric(dataset[, c(columnName)])) {
    columnOfInterest <- na.omit(dataset[, c(columnName)])
    values <- c(length(columnOfInterest))
    Genotype_levels <- levels(factor(dataset[, genotypeCol]))
    Sex_levels <- levels(factor(dataset[, sexCol]))
    values <- c(values, length(levels(factor(
      columnOfInterest
    ))))
    values <-
      c(values, length(Genotype_levels) * length(Sex_levels))
    for (i in 1:length(Genotype_levels)) {
      GenotypeSubset <- subset(dataset, dataset[, genotypeCol] ==
                                 Genotype_levels[i])
      for (j in 1:length(Sex_levels)) {
        GenotypeSexSubset <-
          subset(GenotypeSubset, GenotypeSubset[, sexCol] ==
                   Sex_levels[j])
        columnOfInterestSubset <- na.omit(GenotypeSexSubset[,
                                                            c(columnName)])
        values <- c(values, length(columnOfInterestSubset))
      }
    }
    if (length(values) > 1 && values[1] != 0) {
      ratio = values[2] / values[1]
      message0('Variability (PhenStat function) = ', ratio)
    } else{
      ratio = 1
    }
  } else{
    ratio = 1
  }
  return(as.vector(ratio))
}
# remove zero count categories (filter on sex and genotype)
RemoveZeroFrequencyCategories = function(x,
                                         minSampRequired,
                                         depVar = 'data_point',
                                         drop = TRUE,
                                         totalLevels = 4,
                                         sexCol = 'sex',
                                         genotypeCol = 'biological_sample_group',
                                         sep = '_',
                                         activated = TRUE) {
  note   = NULL
  if (all(dim(x) > 0) && !is.null(complete.cases0(x[, depVar]))) {
    x = x[complete.cases0(x[, depVar]), ,drop=FALSE]
    if (is.numeric(x[, depVar])) {
      lvls   = interaction(x[, sexCol], x[, genotypeCol], sep = sep, drop = drop)
    } else{
      lvls   = interaction(x[, sexCol], x[, genotypeCol], sep = sep, drop = drop)
      #lvls   = interaction(x[, sexCol], x[, genotypeCol], x[, depVar], sep = sep, drop = drop)
    }
    lvls   = lvls[!is.na(x[, depVar])]
    counts = table(lvls)
    newx   = (x[lvls %in% names(counts)[counts >= minSampRequired],])
    #####
    if ((length(names(counts[counts >= minSampRequired])) != totalLevels ||
         !identical(droplevels0(newx), droplevels0(x))) &&
        length(lvls) > 0) {
      note$'Sex-Genotype data included in the analysis' = c(
        paste(
          'Not all Sex*Genotype interactions exist in data. Threshold: ',
          paste0(minSampRequired),
          ' sample per group. INCLUDED categories: ',
          paste(names(counts)[counts >= minSampRequired],
                sep = ',',
                collapse = ', '),
          '',
          sep = '',
          collapse = ''
        )
      )
    }
  } else{
    newx = NULL
  }
  if (any(dim(newx) == 0) || any(dim(x) == 0))
    newx = NULL

  return(list(x = newx, note = note))
}

# Droplevels does not work with null/empty objects
droplevels0 = function(x, ...) {
  if (!is.null(x) && all(dim(x) > 0)) {
    r = droplevels(x = x, ...)
    message0('droplevel applied.')
  } else{
    r = x
    message0('droplevel did not applied.')
  }
  return(r)
}

normalisePhenList =   function(phenlist, colnames = NULL) {
  message0('Normalising the PhenList object in progress ...')
  if (!is.null(colnames)) {
    colnames = colnames[colnames %in% names(phenlist@datasetPL)]

    if (length(colnames) < 1) {
      message0('No variable name found in the data. Please check "colnames" ...')
      return(phenlist)
    }
  } else{
    message0('No variable selected for normalisation. All numerical variables will be normalised.')
    colnames = names(phenlist@datasetPL)
  }
  ######
  phenlist@datasetPL  [, colnames]      = as.data.frame(lapply(
    phenlist@datasetPL[, colnames, drop = FALSE],
    FUN = function(x) {
      if (is.numeric(x) && length(unique(x)) > 1) {
        sdx = sd0(x, na.rm = TRUE)
        r   =    (x - mean(x, na.rm = TRUE)) / ifelse(sdx > 0, sdx, 1)
      } else{
        r = x
      }
      return(r)
    }
  ))
  return(phenlist)
}



# remove zero count categories (filter on sex and genotype)
RemoveZerovarCategories = function(x,
                                   depVar,
                                   minvar = 0,
                                   sex = 'sex',
                                   genotype = 'biological_sample_group',
                                   sep = '_',
                                   drop = TRUE,
                                   method = 'MM') {
  note   = NULL
  # do not move me
  if (any(dim(x) == 0) || is.null(complete.cases0(x[, depVar])))
    return(list(x = NULL, note = note))

  #x      = droplevels0(x)
  x = x[complete.cases(x[, depVar]), , drop = FALSE]
  newx   =  x
  if (is.numeric(x[, depVar]) && (method %in% 'MM')) {
    lvls   = interaction(x[, sex], x[, genotype], sep = sep, drop = drop)
    vars   = tapply(x[, depVar], INDEX = lvls, function(xx) {
      if (is.numeric(xx) && length(na.omit(xx)) > 1) {
        r = var(xx, na.rm = TRUE)
      } else{
        r = 0
      }
      return(r)
    })

    if (!is.null(vars)) {
      Zvars  = names(vars)[vars <= minvar]
      newx   = x[!lvls %in% Zvars, ]
      if (!identical(droplevels0(newx), droplevels0(x))) {
        note$'Sex-Genotype interactions zero variation' = paste0(
          'Some categories have zero variation and then removed from the raw data. List of zero variation categories: ',
          paste0(Zvars, collapse = ', '),
          ''
        )
      }
    }
  }else{
    message0('Checking variation only applied to the MM method. The input mmethod: ',
             method)
  }
  if (is.null(x) ||
      is.null(newx) || any(dim(newx) == 0) || any(dim(x) == 0))
    newx = NULL

  return(list(x = newx, note = note))
}


sd0 = function(x, ...) {
  if (!is.numeric(x))
    return(NA)
  r = if (length(na.omit(x)) > 1)
    sd(x, ...)
  else
    0
  return(r)
}

cleanNULLkeys = function(list){
  requireNamespace('rlist')
  if (!is.null(list)) {
    message0 ('Removing the NULL keys ...')
    list = rlist::list.clean(
      list,
      fun = function(x) {
        length(x) == 0L || is.null(x)
      },
      recursive = TRUE
    )
  }
  return(list)
}

NA2LabelInFactor = function(x, label = 'control') {
  if (is.null(x) || !is.factor(x) || length(x)<1)
    return(x)

  x = factor(
    x,
    levels  = levels(addNA(x)),
    labels  = c(levels(x), label),
    exclude = NULL
  )
  return(x)
}


varsInColsOrReturn = function(x, vars, retValue = NULL) {
  b.vars = vars
  r =  if (!all(vars %in% names(x)) || nrow(x) < 1) {
    retValue
  } else{
    vars[vars %in% names(x)]
  }
  if (!identical(b.vars, vars)) {
    message0('Not all variables in the data: ',
             paste(
               setdiff(b.vars, vars),
               sep      = ', ',
               collapse = ','
             ))
  }
  return(r)
}

FlatteningTheSummary = function(x, name = 'Raw data summary statistics') {
  rs = lapply(x[[name]], function(y) {
    y[1]
  })
  rs1  = unlist(rs)
  r1   = c(min(rs1),
           max(rs1),
           mean(rs1),
           median(rs1),
           sd(rs1),
           names(rs)[(which.min(rs1))],
           names(rs)[(which.max(rs1))])
  r2   = paste(sort(paste(names(rs)  , rs1, sep = '=')), collapse = '~', sep =
                 '~')
  return(c(r1, r2))
}

SummaryStatisticsOriginal = function(x,
                                     depVar,
                                     sex = 'sex',
                                     genotype = 'biological_sample_group',
                                     label = 'Raw data summary statistics',
                                     lower = FALSE,
                                     drop = TRUE,
                                     sep = '_',
                                     removeSpecialChars = FALSE,
                                     what = '[^0-9A-Za-z]',
                                     replace = '_') {
  r   = NULL
  # do not move me
  if (any(dim(x) == 0) ||
      nrow(x) < 1)
    return(list('Error in input data' = 'empty dataset'))

  if (sum(complete.cases0(x[, depVar])) < 1)
    return(list('Error in input data' = 'No value for the response - all NA'))

  x = x[complete.cases0(x[, depVar]), , drop = FALSE]
  #x      = droplevels0(x)
  v0 = varsInColsOrReturn(x, c(sex, genotype, depVar))
  if (is.numeric(x[, depVar])) {
    lvls   = interaction(x[, varsInColsOrReturn(x, c(sex, genotype))], sep = sep, drop = drop)
    vvls   = interaction(x[, varsInColsOrReturn(x, c(sex, genotype))], sep = sep, drop = drop)
  } else {
    if (length(v0) > 1) {
      lvls   = interaction(x[, varsInColsOrReturn(x, c(sex, genotype, depVar))], sep = sep, drop = drop)
      vvls   = interaction(x[, varsInColsOrReturn(x, c(sex, genotype))], sep = sep, drop = drop)

    } else{
      lvls   = interaction(x[, c(genotype, depVar)], sep = sep, drop = drop)
      vvls   = interaction(x[, c(genotype)], sep = sep, drop = drop)
    }
  }


  isNumeric = is.numeric(x[, depVar])
  summaryT    = as.list(tapply(x[, depVar], INDEX = lvls, function(xx) {
    c         = ifelse(length(na.omit(xx)) > 0, length(na.omit(xx)), 0)
    uc        = ifelse(length(na.omit(xx)) > 0, length(unique(na.omit(xx))), 0)
    ##############################
    if (isNumeric) {
      m        = ifelse(length(na.omit(xx)) > 0, mean(xx, na.rm = TRUE) , NA)
      sd       = ifelse(length(na.omit(xx)) > 0, sd0(xx, na.rm = TRUE)  , NA)
      NormTest = OpenStats:::normality.test0(xx)
      r = list(
        'Count' = c                                 ,
        'Unique N'     = uc                         ,
        'Mean'  = m                                 ,
        'SD'    = sd                                ,
        'Diversity in response' = uc > 1            ,
        'Normality test p-val'  = NormTest
      )
    } else{
      r = list('Count'                 = c  ,
               'Unique N'              = uc ,
               'Diversity in response' = uc > 1)
    }
    ##############################
    return(r)
  }, default = -999.991233210123))
  ##
  # fTmp = function(isNum) {
  #   if (isNum) {
  #     r = list('Count' = 0,
  #              'Mean' = NA,
  #              'SD' = NA)
  #   } else{
  #     r = list('Count' = 0)
  #   }
  #   return(r)
  # }
  # summaryT[summaryT %in% c(-999.991233210123)] = fTmp(isNum = isNumeric)
  if (lower)
    nnames = tolower(names(summaryT))
  else
    nnames =  names(summaryT)

  if (removeSpecialChars) {
    nnames = RemoveSpecialChars(nnames, replaceBy = replace, what = what)
  }
  vls = sort(unique(vvls))
  rl = summaryTDetailed = summaryT
  rl2 = NULL
  if (length(vls) > 0) {
    for (n in vls) {
      if (is.null(rl))
        next
      grp = grepl(pattern = n,
                  x = names(rl),
                  fixed = TRUE)
      if (sum(grp) < 1)
        next
      rl2[n] = list(rl[grp])
      rl = rl[!grp]
    }
    if (!is.null(rl) && length(rl)>0)
      rl2['others'] = rl
    summaryT = rl2
  }
  r = list(lbl = list('Collapsed' = summaryT, 'Detailed' = summaryTDetailed))
  names(r) = label

  return(r)
}


# Other columns from the data
OtherExtraColumns = function(obj,
                             ColNames,
                             names,
                             lower = TRUE,
                             what = '[^0-9A-Za-z]',
                             replaceBy = '_') {
  p4 = NULL
  ColNameInd = ColNames %in% colnames(obj)
  newNames = ColNames[ColNameInd]
  if (all(dim(obj) > 0) && length(newNames) > 0) {
    if (lower)
      newNames = tolower(newNames)
    p0 = RemoveSpecialChars(names, what = what, replaceBy = replaceBy)
    p4 = obj[, newNames, drop = FALSE]
    p4 = as.list(p4)
    names(p4) = p0[ColNameInd]
    p4 = RemoveNAFromList(p4)
  } else{
    p4$'Error in extra columns' = "Column names do not exist or the dataset is empty"
  }
  return(p4)
}

RemoveNAFromList = function(x = NULL) {
  if (is.null(x))
    return(x)

  if (length(x) < 1)
    return(x)

  if (length(names(x)) < 1)
    return(x)

  if (length(names(x)) == 1)
    return(x)

  ind = sapply(x, function(x) {
    return (sum(is.na(x))==length(x))
  })

  return(x[!ind])
}

# min including null/empty
min0 = function(x, ...) {
  if (!is.null(x) && length(na.omit(x)) > 0)
    r = min(x, ...)
  else
    r = 0

  return(r)
}

max0 = function(x, ...) {
  if (!is.null(x) && length(na.omit(x)) > 0)
    r = max(x, ...)
  else
    r = 0
  return(r)
}

# Get possible categories for the categorical variables
GetPossibleCategories = function(procedure = NULL, method = 'file') {
  requireNamespace('pingr')
  # Predefine list
  CatList = data.frame(parameter_stable_id = NA, categories = NA)
  message0('Loading the list of possible categories for categorical variables ...')
  if (method %in% c('file', 'update')) {
    message0('* We recommend you to update the file only *once* otherwise the whole process may take longer ...')
    if (method %in% 'update') {
      tmp = updateImpress(updateImpressFileInThePackage = TRUE)
      rm(tmp)
    }
    path = file.path(system.file("extdata", package = "DRrequiredAgeing"),
                     'AllCts.csv')
    message0('Reading the Category levels from the file : \n\t\t ====> ',
             path)
    CatList = read.csv(file = path)#[, -1]
  } else if (method %in% 'solr') {
    if (pingr::is_online()) {
      url = paste0(
        #'http://ves-ebi-d0:8986/solr/pipeline/select?q=procedure_stable_id:',
        'https://wwwdev.ebi.ac.uk/mi/impc/dev/solr/pipeline/select?q=procedure_stable_id:',
        ifelse(is.null(procedure), '*', paste0('"', procedure, '"')),
        '&fq=categories%3A*&fl=parameter_stable_id%2Ccategories&wt=csv&indent=true&rows=10000000'
      )
      message0('Reading the Category levels from : \n\t\t ====> ',
               url)
      CatList = read.csv(file = url)
      if (all(dim(na.omit(CatList)) > 0)) {
        CatList = CatList[!duplicated(CatList),]
        message0(nrow(CatList),
                 ' parameter(s) found ...')
      }
    } else{
      message0('Warning. Please check the internet connection ....')
    }
  } else {
    message0('please select either `file`, `solr`, `update` (latter needs admin permission) ....')
  }
  return(CatList)
}

# Find levels for the parameter
ReadFactorLevelsFromSolr = function(parameter,
                                    CatList,
                                    initReplacePattern = '~~~~',
                                    skip = '\\,',
                                    sep = ',',
                                    fixed = TRUE) {
  note = NULL
  ind = which(CatList$parameter_stable_id == parameter)
  if (is.null(CatList) || any(dim(CatList) < 1) || length(ind) < 1)
    return(list(levels = NULL, note = note))

  if (length(ind) > 1) {
    message0(
      'It seems there are two the same procedure with the same levels. Please check the dataset'
    )
  }

  catCorr = gsub(
    pattern =  skip,
    x = as.character(CatList$categories[ind]),
    replacement = initReplacePattern,
    fixed = fixed
  )

  Levels = unlist(strsplit(x = catCorr, split = sep, fixed = fixed))
  fLevels = gsub(
    pattern = initReplacePattern,
    replacement = skip,
    x = Levels,
    fixed = fixed
  )
  return(list(levels = unique(trimws(fLevels)), note = note))
}

# Fast replacement of nulls
replaceNull = function(x, replace = '-') {
  if (length(x) > 1) {
    x[is.null(x)] == replace
  } else{
    x = replace
  }
  return(x)
}

mapLevelsToFactor = function(levels, newlevels, name = 'response') {
  #############
  res  = NULL
  n1   = paste0('All levels in the ',
                name,
                collapse = '')
  res$n1 = newlevels
  ind    = which(newlevels %in% levels)
  if (length(ind) < 1 ||
      length(newlevels) < 0 ||
      length(levels) < 0 ||
      is.null(newlevels) ||
      is.null(levels) ||
      length(ind) == length(levels)) {
    if (!is.null(res$n1))
      names(res) = n1
    return(list(levels = levels, note = res))
  }
  ############
  n0 = paste0('Missing levels in ',
              name,
              collapse = '')
  res$n0 = newlevels[-ind]
  names(res) = c(n1, n0)

  flevels = c(levels, newlevels[-ind])
  return(list(levels = flevels, note = res))
}



RemoveNonExsistingFilesFromList = function(fileList) {
  logical = file.exists(fileList)
  if (sum(!logical) > 0) {
    message0('Non-exsisting file(s): \n\t\t ~> ',
             paste(fileList[!logical], collapse  = '\n\t\t ~> '))
  }
  return(fileList[logical])
}
# Compress a file
compressFiles = function(fileList,
                         dir,
                         filename,
                         overwrite = TRUE,
                         flags = '-D -j',
                         rmSource = TRUE,
                         quiet = TRUE) {
  if (length(fileList) < 1 ||
      length(RemoveNonExsistingFilesFromList(fileList)) < 1) {
    message0('No file to compress or files do not exists')
    return(NULL)
  }
  if (quiet)
    flags = paste0(flags, ' -q')
  fileList      = RemoveNonExsistingFilesFromList(fileList)
  fileListPaste = paste('\n\t\t ~~> ', fileList, sep = '')
  zipFile =  file.exists0(file.path(dir, paste0(filename,  '.zip')), overwrite = overwrite)
  MB      = 1024
  message0('Compressing file ...')
  result = zip(zipfile = zipFile,
               files = fileList,
               flags = flags)
  if (file.exists(zipFile) &&
      (file.size(zipFile) / MB) > 0 &&
      result == 0) {
    message0('Successfully compressed. file(s): ', fileListPaste)
    if (rmSource)
      unlink(fileList)
  } else{
    message0('Failed in compressing file. file(s): ', fileListPaste)
  }
  return(list(status = result, file = zipFile))
}

# Base64 encoder
base64 = function(x, active) {
  if (active) {
    message0('Based64 encryption as well as (gzip) compression in progress ...')
    r = base64enc::base64encode(memCompress(x, "g"))
  } else{
    r = x
  }
  return(r)
}

RR_thresholdCheck = function(data,
                             depVar,
                             parameter,
                             controllab = 'control',
                             Experimentallab = 'experimental',
                             genotypeColumn = 'biological_sample_group',
                             sexColumn = 'sex',
                             threshold = 60,
                             methodmap) {
  r = list('Criteria result' = TRUE           ,
           'Value'     = 'not applicable'     ,
           'Threshold' = 'not applicable')
  if (!is.null(data) && !is.null(depVar) && !is.null(parameter)) {
    method = getMethodi(
      var = parameter,
      type = ifelse(is.numeric(data[, depVar]),
                    'numeric',
                    'charachter'),
      methodMap = methodmap
    )

    if (method %in% 'RR') {
      data = droplevels0(data[complete.cases(data[, depVar]), , drop = FALSE])
      lv   = levels(droplevels0(subset(data, data[, genotypeColumn] %in% Experimentallab))[, sexColumn])
      tbl  = table(subset(data, (data[, genotypeColumn] %in% controllab) & (data[,sexColumn] %in% lv))[, sexColumn])
      if (!is.null(tbl)     &&
          sum(dim(tbl) > 0) &&
          all(tbl <= threshold)) {
        message0(
          'Not enough data for running RR, threshold = ',
          threshold,
          ': ',
          paste(names(tbl), tbl, sep = '=', collapse = ', ')
        )
        r = list(
          'Criteria result' = FALSE,
          'Value'           = paste(names(tbl), tbl, sep = '=', collapse = ', '),
          'Threshold'       = threshold
        )
      }
    }
  }
  return(r)
}

# rename if file exists (Danger! if overwrite==TRUE then it removes the files!)
file.exists0 = function(file, overwrite = FALSE, ...) {
  counter = 1
  ###
  fname = basename(file)
  dotLocation = which(strsplit(fname, "")[[1]] == ".")
  if (length(dotLocation) > 0) {
    name = substr(x = fname,
                  start = 0,
                  stop = dotLocation - 1)

    ext = substr(x = fname,
                 start =  dotLocation,
                 stop = nchar(fname))

  } else{
    ext  = ''
    name = fname
  }
  ###
  if (overwrite && file.exists(file, ...)) {
    unlink(file, recursive = FALSE, force = FALSE)
  }
  while (!overwrite && file.exists(file, ...)) {
    file = file.path(dirname(path = file),
                     paste0(name, '_', counter, ext))
    counter = counter + 1
  }
  if (counter > 1) {
    write(
      counter - 1,
      file = file.path0(
        dirname(file),
        'duplicates_detected.txt',
        IncludedFileName = TRUE,
        create = TRUE,
        check = FALSE
      )
    )
  }
  return(file)
}

# add bits to file name
addBitsToFileName = function(file,
                             bit = '',
                             sep = '_',
                             fullpath = TRUE) {
  if (fullpath)
    file = paste(dirname(file), bit, basename(file), sep = sep)
  else
    file = paste(bit, basename(file), sep = sep)
  return(file)
}
###
fillNAgaps <- function(x, firstBack=FALSE) {
  ## NA's in a vector or factor are replaced with last non-NA values
  ## If firstBack is TRUE, it will fill in leading NA's with the first
  ## non-NA value. If FALSE, it will not change leading NA's.

  # If it's a factor, store the level labels and convert to integer
  lvls <- NULL
  if (is.factor(x)) {
    lvls <- levels(x)
    x    <- as.integer(x)
  }

  goodIdx <- !is.na(x)

  # These are the non-NA values from x only
  # Add a leading NA or take the first good value, depending on firstBack
  if (firstBack)   goodVals <- c(x[goodIdx][1], x[goodIdx])
  else             goodVals <- c(NA,            x[goodIdx])

  # Fill the indices of the output vector with the indices pulled from
  # these offsets of goodVals. Add 1 to avoid indexing to zero.
  fillIdx <- cumsum(goodIdx)+1

  x <- goodVals[fillIdx]

  # If it was originally a factor, convert it back
  if (!is.null(lvls)) {
    x <- factor(x, levels=seq_along(lvls), labels=lvls)
  }

  x
}


# plot analised data with window
plot_win = function(phenlistObject, r, depVariable, check, ...) {
  requireNamespace('rlist')
  col = pch = as.factor(interaction(rlist::list.clean(
    list(
      phenlistObject@datasetPL$Genotype,
      phenlistObject@datasetPL$Sex,
      phenlistObject@datasetPL$LifeStage
    )
  )))
  if (!is.null(r) &&
      !is.null(r$finalModel$data$outlierWeight)) {
    r$finalModel$FullWeight[r$finalModel$data$outlierWeight < 1] = NA
    r$finalModel$FullWeight  =  fillNAgaps(r$finalModel$FullWeight, firstBack = TRUE)
    lwd = phenlistObject@datasetPL$outlierWeight
    lwd[is.na(lwd)] = 0
    lwd = (2 - lwd) ^ 2
  } else{
    lwd = 1
  }
  par(mar = c(5.1, 4.1, 4.1, 9.0))
  requireNamespace('SmoothWin')
  plot(r,
       col = as.integer(col),
       pch = as.integer(pch),
       cex.main = .7,
       lwd = lwd    ,
       ...)
  legend(
    'topright',
    legend = unique(col),
    col = as.integer(unique(col)),
    pch = as.integer(unique(pch)),
    inset = c(-.3, 0),
    xpd = TRUE,
    cex = .7
  )

  if (check == 1) {
    rt = table(phenlistObject@datasetPL$Batch)
    rt2 = rt[which(rt < 2)]
    rdata2 = phenlistObject@datasetPL[phenlistObject@datasetPL$Batch %in% names(rt2),]
    points(as.numeric(as.Date(rdata2$Batch)),
           rdata2[, depVariable],
           pch = 8,
           col = 'gray')
  }
}


detach_package <- function(pkg, character.only = FALSE)
{
  if (!character.only)
  {
    pkg <- deparse(substitute(pkg))
  }
  search_item <- paste("package", pkg, sep = ":")
  while (search_item %in% search())
  {
    detach(search_item,
           unload = TRUE,
           character.only = TRUE)
  }
}

# boxplot for window
boxplot_win = function(phenlistObject, we, threshold, ...) {
  boxplot(
    data_point ~ Genotype + Sex,
    data = phenlistObject@datasetPL[we > threshold,],
    col = rep(2:3, 2),
    cex.axis = .7,
    cex.main = .7,
    ...
  )
}

# make JSON vector (with [ or not)
SingleOrMore = function(name, value, comma = FALSE) {
  NewVal = ifelse(
    is.numeric(value),
    paste0(value, collapse = ', '),
    paste0('"', value, '"', collapse = ', ')
  )
  if (length(value) > 1) {
    outp = paste0('"', name, '":[', NewVal, ']', collapse = '')
  } else{
    outp = paste0('"', name, '":', NewVal, '', collapse = '')
  }
  if (comma)
    outp = paste(outp, ', ', sep = '', collapse = '')
  return(outp)
}


# cut the filename is it is too long
CutFileNameIfTooLong = function(fullPath,
                                max = 255,
                                checkRoot = TRUE) {
  file = basename(fullPath)
  dir  = dirname(fullPath)
  nch  = nchar(file)
  message0('The number of characters in file name: ', nch)
  k = 7
  if (nch >= max) {
    message0('File name too long. cut 255 characters from the RIGHT')
    file = paste0(
      RandomRegardSeed(
        n = 1,
        decimal = TRUE,
        round = 6
      ) * 10 ^ 6,
      '_',
      substr(x = file, start =  k + 2, nch),
      collapse = ''
    )
  }
  return(ifelse(checkRoot, file.path(dir, file), file.path(file)))
}



# plot windowing
PlotWindowingResult = function(args, overwrite = FALSE, ...) {
  if (args$plot      &&
      !is.null(args$we)) {
    message0('ploting the results ...')
    file = ifelse(
      is.null(args$PicDir),
      file.path(getwd(), args$filename),
      file.path(dir.create0(args$PicDir, recursive = TRUE), args$filename)
    )
    file = file.exists0(paste(file, '.png', sep = ''), overwrite = overwrite)
    file = CutFileNameIfTooLong(file)
    ##########
    if (args$storeplot)
      png(
        filename = file,
        width = 1000,
        height = 1000,
        res = 120
      )
    plot_win(
      phenlistObject = args$phenlistObject,
      main = args$main,
      r = args$r,
      check = args$check,
      depVariable = args$depVariable,
      xaxt = 'n'
    )
    tick =  round(seq(min(args$tt), max(args$tt), length.out = 5))
    axis(1, at = tick, labels = as.Date(tick, '1970-01-01'))
    if (args$storeplot)
      graphics.off()
    # if (args$storeplot && 1 == 0) {
    #   png(filename  = addBitsToFileName(file = file, bit = 'boxplot'))
    #   boxplot_win(
    #     phenlistObject = args$phenlistObject,
    #     we = args$we,
    #     threshold = args$threshold
    #   )
    #   if (args$storeplot)
    #     graphics.off()
    # }
    return(file)
  } else{
    return(NULL)
  }
}

## outlier detection
PhenListOutlierDetection = function(pl                       ,
                                    alpha   = .9             ,
                                    wgtFUN  = "sm2.adaptive" ,
                                    plot    = TRUE           ,
                                    active  = TRUE) {
  requireNamespace("robustbase")
  if (is.null(pl) || !is.numeric(pl@datasetPL$data_point))
    return(pl)
  if(!active){
    message0('\t  Outlier detection is disabled ...')
    pl@datasetPL$outlierWeight= 1
    return(pl)
  }
  #######
  message0('\tApplying the outlier detection algorithm ...')
  df   = pl@datasetPL
  cols = c('BatchInt', 'data_point','Weight')
  #######
  if (CheckIfNameExistInDataFrame(obj = df, 'Batch',checkLevels = TRUE))
    df$BatchInt = Date2Integer(df$Batch)

  # if (CheckIfNameExistInDataFrame(obj = df, 'LifeStage',checkLevels = TRUE))
  #   df$LifeStageInt  = as.integer(df$LifeStage)
  #
  # if (CheckIfNameExistInDataFrame(obj = df, 'Sex',checkLevels = TRUE))
  #   df$SexInt  = as.integer(df$Sex)

  plC      =  subset(df, df$Genotype %in% pl@refGenotype)
  message0(
    '\tOutlier detection based on the following variables:\n\t  ',
    paste(cols[cols %in% names(df)], sep = ', ', collapse = ', ')
  )
  plOutSet =	plC[, cols[cols %in% names(df)]]

  oat = tryCatch(
    expr = {
      ouObj = robustbase::covMcd(data.matrix(plOutSet), alpha = alpha, wgtFUN = wgtFUN)
      return(ouObj)
    },
    warning = function(war) {
      message0('Ops! outlier detection algorithm failed! no outlier will be marked ...')
      message0(war)
      return(NULL)
    },
    error   = function(err) {
      message0('Ops! outlier detection algorithm failed! no outlier will be marked ...')
      message0(err)
      return(NULL)
    }
  )
  #################
  if (!is.null(oat)   &&
      plot            &&
      is(oat, 'mcd')  &&
      active) {
    plot(
      plOutSet[, 1],
      plOutSet[, 2],
      col = oat$mcd.wt * 1 + 2,
      pch =  2 + (1 - oat$mcd.wt) * 10
    )
    plC$outlierWeight     = as.vector(oat$mcd.wt)
  } else{
    plC$outlierWeight     = 1
  }
  #################
  plC$UniqueId          = paste0(plC$external_sample_id, plC$discrete_point)
  pl@datasetPL$UniqueId = paste0(pl@datasetPL$external_sample_id,
                                 pl@datasetPL$discrete_point)
  #################
  pl@datasetPL = merge(
    x   = pl@datasetPL                          ,
    y   = plC[, c('outlierWeight', 'UniqueId')] ,
    by  = 'UniqueId'                            ,
    all = TRUE
  )
  pl@datasetPL$outlierWeight[is.na(pl@datasetPL$outlierWeight)] = 1
  pl@datasetPL$UniqueId = NULL
  return(pl)
}
## date to integer
Date2Integer = function(x) {
  dates = as.Date(x)
  dates = as.numeric(dates)
  return(dates)
}

# check if variable name exists in the data frame
CheckIfNameExistInDataFrame = function(obj, name, checkLevels = TRUE) {
  if (!NullOrError(obj) && !is.null(name) && length(name) > 0) {
    r = name %in% names(obj)
    if (r &&
        checkLevels &&
        !is.numeric(obj[, name]) &&
        length(unique(obj[, name])) == 1) {
      r = FALSE
    }
  } else{
    r = FALSE
  }
  return(r)
}


# file.path to check for space and special characters in the path
file.path0 = function(...,
                      check = TRUE,
                      what = '[^0-9A-Za-z/:]',
                      replace = '_',
                      create = TRUE,
                      IncludedFileName = FALSE) {
  path = file.path(...)
  if (check) {
    path = RemoveSpecialChars (path, what = what, replaceBy = replace)
  }
  if (create) {
    if (!IncludedFileName) {
      dir.create0(path, recursive = TRUE)
    } else{
      dir.create0(dirname(path), recursive = TRUE)
    }
  }
  return(path)
}




# exist in list
IsInList = function(item = NULL,
                    list = NULL ,
                    message = '') {
  r = FALSE
  if (length(list) < 1 || length(item) < 1)
    return(r)
  if (any(list %in% item)) {
    r = TRUE
    message0(message)
  }
  return(r)
}

sortList = function(x,...){
  x[order(names(x),...)]
}

subModelVectorOutput = function(object ,
                                subBreakColumns = NULL    ,
                                ExtraCols = NULL          ,
                                JSON      = FALSE         ,
                                ReportNullSchema = TRUE) {
  if (is.null(object)  ||
      is.null(ExtraCols))
    return(NULL)
  #
  message0('Submodel estimation in progress ...')
  o  = OpenStats:::OpenStatsComplementarySplit(object = object, variables = subBreakColumns)
  #
  if (is.null(o))
    return(NULL)
  #
  r = lapply(o, function(x) {
    res = OpenStatsReport(
      object            = x           ,
      othercolumns      = ExtraCols   ,
      JSON              = JSON        ,
      RemoveNullKeys    = TRUE        ,
      ReportNullSchema  = ReportNullSchema
    )
    if (!is.null(res)) {
      #res = vectorOutputDictionary2List(res)
      res$`Classification tag` = OpenStats:::classificationTag(x)
    }
    return(res)
  })
  if (!is.null(r))
    names(r) = names(o)
  return(r)
}

WhileListForDic2list=function(){
  r = c("Genotype contribution",
        "Genotype effect size",
        "Genotype estimate",
        "Genotype p-value",
        "Genotype percentage change",
        "Interactions p-value",
        "Intercept estimate",
        "Intercept p-value",
        "LifeStage EvKO effect size",
        "LifeStage EvKO estimate",
        "LifeStage EvKO p-value",
        "LifeStage EvKO standard error",
        "LifeStage LvKO effect size",
        "LifeStage LvKO estimate",
        "LifeStage LvKO p-value",
        "LifeStage LvKO standard error",
        "LifeStage effect size",
        "LifeStage estimate",
        "LifeStage p-value",
        "LifeStage standard error",
        "LifeStageSexGenotype FvEvKO effect size",
        "LifeStageSexGenotype FvEvKO estimate",
        "LifeStageSexGenotype FvEvKO p-value",
        "LifeStageSexGenotype FvEvKO standard error",
        "LifeStageSexGenotype FvLvKO effect size",
        "LifeStageSexGenotype FvLvKO estimate",
        "LifeStageSexGenotype FvLvKO p-value",
        "LifeStageSexGenotype FvLvKO standard error",
        "LifeStageSexGenotype MvEvKO effect size",
        "LifeStageSexGenotype MvEvKO estimate",
        "LifeStageSexGenotype MvEvKO p-value",
        "LifeStageSexGenotype MvEvKO standard error",
        "LifeStageSexGenotype MvLvKO effect size",
        "LifeStageSexGenotype MvLvKO estimate",
        "LifeStageSexGenotype MvLvKO p-value",
        "LifeStageSexGenotype MvLvKO standard error",
        "Sex FvKO effect size",
        "Sex FvKO estimate",
        "Sex FvKO p-value",
        "Sex FvKO standard error",
        "Sex MvKO effect size",
        "Sex MvKO estimate",
        "Sex MvKO p-value",
        "Sex MvKO standard error",
        "Sex effect size",
        "Sex estimate",
        "Sex p-value",
        "Sex standard error",
        "Weight effect size",
        "Weight estimate",
        "Weight p-value",
        "Weight standard error",
        "Genotype contribution",
        "Genotype effect size",
        "Genotype estimate",
        "Genotype p-value",
        "Genotype percentage change",
        "Interactions p-value",
        "Intercept estimate",
        "Intercept p-value",
        "LifeStage EvKO effect size",
        "LifeStage EvKO estimate",
        "LifeStage EvKO p-value",
        "LifeStage EvKO standard error",
        "LifeStage LvKO effect size",
        "LifeStage LvKO estimate",
        "LifeStage LvKO p-value",
        "LifeStage LvKO standard error",
        "LifeStage effect size",
        "LifeStage estimate",
        "LifeStage p-value",
        "LifeStage standard error",
        "LifeStageSexGenotype FvEvKO effect size",
        "LifeStageSexGenotype FvEvKO estimate",
        "LifeStageSexGenotype FvEvKO p-value",
        "LifeStageSexGenotype FvEvKO standard error",
        "LifeStageSexGenotype FvLvKO effect size",
        "LifeStageSexGenotype FvLvKO estimate",
        "LifeStageSexGenotype FvLvKO p-value",
        "LifeStageSexGenotype FvLvKO standard error",
        "LifeStageSexGenotype MvEvKO effect size",
        "LifeStageSexGenotype MvEvKO estimate",
        "LifeStageSexGenotype MvEvKO p-value",
        "LifeStageSexGenotype MvEvKO standard error",
        "LifeStageSexGenotype MvLvKO effect size",
        "LifeStageSexGenotype MvLvKO estimate",
        "LifeStageSexGenotype MvLvKO p-value",
        "LifeStageSexGenotype MvLvKO standard error",
        "Sex FvKO effect size",
        "Sex FvKO estimate",
        "Sex FvKO p-value",
        "Sex FvKO standard error",
        "Sex MvKO effect size",
        "Sex MvKO estimate",
        "Sex MvKO p-value",
        "Sex MvKO standard error",
        "Sex effect size",
        "Sex estimate",
        "Sex p-value",
        "Sex standard error",
        "Weight effect size",
        "Weight estimate",
        "Weight p-value",
        "Weight standard error")
  return(r)
}
# A new vector output for this package only
VectorOutput0 = function(c.ww0,
                         ExtraCols                 ,
                         activeWindowing           ,
                         na              = 'string',
                         null            = 'null',
                         subBreakColumns = NULL) {
  ####
  if (!NullOrError(c.ww0$NormalObj)) {
    p1 = OpenStats::OpenStatsReport(
      c.ww0$NormalObj,
      othercolumns = ExtraCols,
      JSON = FALSE,
      ReportNullSchema = TRUE
    )
    p1$`Classification tag` = OpenStats:::classificationTag(c.ww0$NormalObj)
    #####
    if (!is.null(subBreakColumns)) {
      p1$'Sub Models' = subModelVectorOutput(
        object = c.ww0$NormalObj,
        subBreakColumns = subBreakColumns,
        ExtraCols = ExtraCols,
        JSON = FALSE,
        ReportNullSchema = TRUE
      )
    }
    #####
  } else{
    p1 = NULL
  }
  # The second piece (Windowing results)
  if (activeWindowing && !NullOrError(c.ww0$WindowedObj)) {
    p2 =  OpenStats::OpenStatsReport(
      c.ww0$WindowedObj  ,
      othercolumns = ExtraCols,
      JSON = FALSE,
      ReportNullSchema = TRUE
    )
    p3 =  OpenStats::OpenStatsReport(
      c.ww0$FullObj  ,
      othercolumns = ExtraCols,
      JSON = FALSE,
      ReportNullSchema = TRUE
    )
    p4 =  OpenStats::OpenStatsReport(
      c.ww0$FullWindowedObj  ,
      othercolumns = ExtraCols,
      JSON = FALSE,
      ReportNullSchema = TRUE
    )
    #### Code for
    # OpenStats:::OpenStatsComplementarySplit()
    p2$`Classification tag` = OpenStats:::classificationTag(c.ww0$WindowedObj)
    p3$`Classification tag` = OpenStats:::classificationTag(c.ww0$FullObj)
    p4$`Classification tag` = OpenStats:::classificationTag(c.ww0$FullWindowedObj)
    #####
    if (!is.null(subBreakColumns)) {
      p2$'Sub Models' =
        subModelVectorOutput(
          object = c.ww0$WindowedObj        ,
          subBreakColumns = subBreakColumns ,
          ExtraCols = ExtraCols             ,
          JSON = FALSE                      ,
          ReportNullSchema = TRUE
        )
      ###
      p3$'Sub Models' =
        subModelVectorOutput(
          object = c.ww0$FullObj           ,
          subBreakColumns = subBreakColumns,
          ExtraCols = ExtraCols            ,
          JSON = FALSE                     ,
          ReportNullSchema = TRUE
        )
      ###
      p4$'Sub Models' =
        subModelVectorOutput(
          object = c.ww0$FullWindowedObj    ,
          subBreakColumns = subBreakColumns ,
          ExtraCols = ExtraCols             ,
          JSON = FALSE                      ,
          ReportNullSchema = TRUE
        )
    }
    #####
  } else{
    p2 = p3 = p4 = NULL
  }
  ###
  # p1 = vectorOutputDictionary2List(p1)
  # p2 = vectorOutputDictionary2List(p2)
  # p3 = vectorOutputDictionary2List(p3)
  # p4 = vectorOutputDictionary2List(p4)

  output = list(
    'Normal result'       = p1,
    'Windowed result'     = p2,
    'Full model result'   = p3,
    'Full model windowed result' = p4
  )

  return(list(json = NULL,
              list = output))
}

# vectorOutputDictionary2List = function(x, whitelist = WhileListForDic2list()) {
#   if (is.null(x))
#     return(x)
#   nx = names(x)
#   if (length(nx) < 1)
#     return(x)
#   for (name in nx) {
#     if (!name %in% whitelist)
#       next
#     if (is(x[[name]], 'list')) {
#       x$WILLBECHANGED = dictionary2listConvert(x[[name]])
#       names(x)[names(x) %in% 'WILLBECHANGED'] = paste0(name, '_2')
#     }
#
#   }
#   return(x)
# }

# Area under the curve
AUC = function(x, y) {
  lx = length(x)
  s = 0
  if (lx > 1) {
    for (i in 1:(lx - 1)) {
      s = s + abs(y[i] + y[i + 1]) * abs(x[i + 1] - x[i]) / 2
    }
  }
  return(s)
}

# Chunk a vector of in particular colonies into chunks (used mainly in Batch production)
chunkVector =  function(x,
                        n = 10,
                        min = 25,
                        activate = TRUE) {
  r = list()
  lx = length(x)
  if (lx > min && n > 1 && activate) {
    r = split(x, cut(seq_along(x), n, labels = FALSE))
  } else{
    r$'1' = x
  }
  return(r)
}

# Plan B for JSON not quouted keys!
checkQouteNAandNaN = function(pattern, replacement, x, ignoreCase = FALSE) {
  for (i in 1:length(pattern)) {
    pt  = pattern[i]
    rep = replacement[i]
    x = gsub(
      pattern = pt,
      replacement = rep,
      x = x,
      fixed = TRUE,
      ignore.case = ignoreCase
    )
  }
  return(x)
}

# Add windowing weights and store the data again
StoreRawDataAndWindowingWeights = function(storeRawData,
                                           activeWindowing        ,
                                           c.ww0                  ,
                                           ## 2
                                           orgData                ,
                                           RawoutputFile          ,
                                           files                  ,
                                           dir                    ,
                                           filename               ,
                                           ### 3
                                           ReadMeTxt              ,
                                           ReadMeFile             ,
                                           ### 4
                                           colnames = c('external_sample_id', 'AllModelWeights'),
                                           byid     = 'external_sample_id'                      ,
                                           compressRawData = TRUE,
                                           ### 5
                                           methodmap  ) {
  message0('Extra columns in the stored data: ',
           paste(colnames, sep = ', ', collapse = ', '))
  if (storeRawData    &&
      activeWindowing &&
      !NullOrError(c.ww0$WindowedObj) &&
      (c.ww0$method %in% 'MM')        &&
      all(colnames %in% names(c.ww0$InputObject@datasetPL))) {
    n3.5.w = merge(
      x = orgData,
      y = c.ww0$InputObject@datasetPL[, colnames],
      by = byid,
      all = TRUE
    )
    message0('writting the raw data + windowing weights to disk ... \n \t\t===> ',
             RawoutputFile)
    write.csv(x = n3.5.w,
              file = RawoutputFile,
              row.names = FALSE)
    write(x        = ReadMeTxt ,
          file     = ReadMeFile,
          ncolumns = 1)

    if (compressRawData) {
      comRes = compressFiles(
        fileList = files    ,
        dir = dir           ,
        filename = filename ,
        overwrite = TRUE
      )
      if (comRes$status == 0)
        message0('[data after adding the windowing weight column] Compression successful')
      else
        message0('[data after adding the windowing weight column] Compression failed')
    }
  }
}

# new message with timestamp in the begining
message0 = function(...,
                    ShowTimeTag = TRUE,
                    active = TRUE) {
  if (active) {
    if (ShowTimeTag)
      message(Sys.time(), '. ', ...)
    else
      message(...)
  }
}

## ReadMe
ReadMe = function(obj, URL = NULL, skip = NULL) {
  if('biological_sample_group' %in% names(obj)){
    obj = subset(obj, obj$biological_sample_group %in% 'experimental')
  }
  if (!is.null(obj) && nrow(obj) > 0) {
    ReadMeTxt = paste(
      c(
        'Gene_symbol',
        'Procedure',
        'Parameter',
        'Centre',
        'Strain',
        #'Metadata',
        'Zygosity',
        'Colony_id',
        'Metadata_hash',
        'URL'
      ),
      c(
        UniqueAndNNull(obj$gene_symbol, removeSpecials = FALSE)          ,
        UniqueAndNNull(obj$procedure_group, removeSpecials = FALSE)      ,
        UniqueAndNNull(obj$parameter_stable_id, removeSpecials = FALSE)  ,
        UniqueAndNNull(obj$phenotyping_center, removeSpecials = FALSE)   ,
        UniqueAndNNull(obj$strain_accession_id, removeSpecials = FALSE)  ,
        #UniqueAndNNull(obj$metadata,removeSpecials = FALSE)             ,
        UniqueAndNNull(obj$zygosity, removeSpecials = FALSE)             ,
        UniqueAndNNull(obj$colony_id[obj$biological_sample_group %in% 'experimental'],removeSpecials = FALSE),
        UniqueAndNNull(obj$metadata_group, removeSpecials = FALSE)       ,
        UniqueAndNNull(URL, removeSpecials = FALSE)
      ),
      sep = ' = '
    )
    message0('[',
             paste(if (is.null(skip)) {
               ReadMeTxt
             } else{
               ReadMeTxt[-skip]
             },
             collapse =  ']~>['),
             ']')
  } else{
    ReadMeTxt = 'Not found!'
  }

  return(ReadMeTxt)
}




# Replace in list
replaceInList <- function (x, FUN, ...)
{
  if (is.list(x)) {
    for (i in seq_along(x)) {
      x[i] <- list(replaceInList(x[[i]], FUN, ...))
    }
    x
  }
  else
    FUN(x, ...)
}


WriteToDB = function(df,
                     dbname    = 'db' ,
                     TableName = 'DR10',
                     maxtry    = 20000,
                     steptry   = 40,
                     maxdelay  = .30) {
  requireNamespace("DBI")
  dbname  = RemoveSpecialChars(dbname, what = '[^0-9A-Za-z/.]')
  message0('Writting to the SQLite ...')
  message0('\tDB name: ', dbname)
  dbtemp                = df
  names(dbtemp)[names(dbtemp) %in% '']  =  'Results'
  dbtemp                = as.data.frame(as.list(dbtemp))
  dbpath                = file.path0(
    'db'                    ,
    dbname                  ,
    check = FALSE           ,
    create = TRUE           ,
    IncludedFileName = TRUE
  )
  res     = NULL
  counter = 1
  for (i in 1:maxtry) {
    res = tryCatch(
      expr = {
        con =
          DBI::dbConnect(
            drv    = RSQLite::SQLite(),
            dbname = dbpath           ,
            synchronous       = NULL
          )
        DBI::dbWriteTable(con, TableName, dbtemp, append = TRUE)
        Sys.sleep(.1)
        DBI::dbDisconnect(conn     = con)
      },
      warning = function(war) {
        DBI::dbDisconnect(conn     = con)
        return(NULL)
      },
      error   = function(err) {
        DBI::dbDisconnect(conn     = con)
        return(NULL)
      }
    )

    if (is.null(res)) {
      sleepTime = RandomRegardSeed(1, max = maxdelay)
      message0('\t',i, '. Retrying ', sleepTime, 's ...')
      Sys.sleep(sleepTime)
      if (i %% steptry == 0) {
        newBase = sub(pattern = '.*(_Dup_)', '' , basename(dbpath), perl = TRUE)
        dbpath  = file.path(dirname(dbpath), paste0(counter, '_Dup_', newBase))
        counter = counter + 1
      }

    } else{
      message0('\tWritting to the database successful ...')
      break
    }
    if (i >= maxtry) {
      stop('\tsomething wrong with the data or the database')
    }
  }
}

FinalJson2ObjectCreator = function(FinalList,
                                   null = 'null',
                                   na = 'null' ,
                                   auto_unbox = TRUE,
                                   SpecialString = '==!!(:HAMED:)!!==',
                                   rep = 3,
                                   removeSpecialsFromNames = FALSE) {
  requireNamespace('jsonlite')
  message0('Forming the JSON object ...')
  FinalList = replaceInList(
    FinalList,
    FUN  = function(x) {
      if (is.null(x)) {
        x = SpecialString
      } else{
        x
      }
    }
  )

  FinalList = replaceInList(
    FinalList,
    FUN  = function(x) {
      if (!is.null(x) && any(is.infinite(x))) {
        x[is.infinite(x)] = sign(x[is.infinite(x)])*10^16
      } else{
        x
      }
    }
  )

  if (removeSpecialsFromNames)
    FinalList = LowerandRemoveSpecials(FinalList)
  for (i in 1:rep) {
    JsonObj  = jsonlite::toJSON(
      x = FinalList,
      null = null,
      na   = na,
      digits = NA,
      auto_unbox = ifelse(i == rep, auto_unbox, FALSE)
    )
    FinalList = jsonlite::fromJSON(txt =   JsonObj)
  }
  JsonObj = gsub(
    pattern = paste0('"', SpecialString, '"'),
    replacement = '{}',
    x = JsonObj,
    ignore.case = FALSE,
    fixed = TRUE
  )
  # Remove new lines! e.g. in Corneal Reflex
  JsonObj = gsub(pattern = '\n',
                 replacement = ' ',
                 x = JsonObj)
  return(JsonObj)
}



LowerandRemoveSpecials <- function(x)
{
  cnames <- names(x)
  if (is.null(cnames))
    return (x)
  x1 <- lapply(cnames, function(y)
    LowerandRemoveSpecials(x[[y]]))
  if ("data.frame" %in% class(x))
    x1 <- as.data.frame(x1)
  names(x1) <- gsub("[^[:alnum:]]", "_", tolower(cnames))
  return(x1)
}


UnzipAndfilePath = function(file, quiet = TRUE, order = TRUE) {
  if (!grepl(x = basename(file),
             pattern = '.zip',
             fixed = TRUE)) {
    message0('It is not a zip file:')
    message0('\t', file)
    return(file)
  }
  message0('Unziping file ...')
  # get the file url
  forecasturl = file
  # create a temporary directory
  td = file.path0(
    tempdir(),
    paste0(sample(LETTERS, 4, replace = TRUE), collapse = ''),
    check = FALSE,
    IncludedFileName = FALSE,
    create = TRUE
  )
  # create the placeholder file
  tf =  tempfile(
    pattern = paste0(sample(LETTERS, 4, replace = TRUE), collapse = ''),
    tmpdir = td,
    fileext = ".zip"
  )
  if (file.exists(forecasturl)) {
    file.copy(from = forecasturl,
              to =  tf,
              overwrite = TRUE)
  } else{
    # download into the placeholder file
    download.file(forecasturl, tf, quiet = quiet)
  }
  # get the name of the first file in the zip archive
  fname = unzip(tf, list = TRUE)$Name
  # unzip the file to the temporary directory
  unzip(tf,
        files = fname,
        exdir = td,
        overwrite = TRUE)
  # fpath is the full path to the extracted file
  fpath = file.path(td, fname)
  if (order) {
    fpath = fpath[order(nchar(fpath), fpath)]
  }
  message0(
    'Successfully unziped. List of files:\n',
    paste0('\t\t', 1:length(fpath), ' ==> ', fpath, collapse = '\n')
  )
  return(fpath)
}

# Create the relative path from the full path
relativePath = function(path, reference) {
  if (is.null(path))
    return(NULL)
  ######
  p = gsub(
    pattern = reference,
    replacement = '/',
    x = path,
    fixed = TRUE
  )
  p = gsub(
    pattern = '//',
    replacement = '/',
    x = p,
    fixed = TRUE
  )
  return(p)
}


# Closest points for the mutants
closest.time.value  = function(x,
                               time1,
                               y,
                               time2,
                               yind,
                               TimeOnly = TRUE,
                               maxIter  = 200) {
  requireNamespace('abind')
  message0('Resampling. Creating the index ...')
  x.df = data.frame(x.val = x, x.time = time1)
  y.df = data.frame(y.val = y,
                    y.time  = time2,
                    y.index = yind)
  message0('Resampling. Finding the distance of values from the mutants ...')
  output2  = ol = lapply(1:nrow(x.df), function(i) {
    tt <-
      cbind(x.df[i, ],
            lapply(x.df[i, ]$x.val, function(v) {
              diff <- abs(y.df$y.val - v)
              y.df$dist.V = diff
              out <- y.df
            }),
            ind = i,
            row.names = NULL)
    tt$dist.T <- abs(tt$x.time - tt$y.time)
    tt$totalD  = tt$dist.V + tt$dist.T
    if (TimeOnly) {
      tt = tt[sample(1:nrow(tt)), ]
    } else{
      tt = tt[order(tt$totalD)  , ]
      tt = tt[order(tt$dist.V)  , ]
    }
    tt   = tt[order(tt$dist.T)  , ]
  })
  dol = 1
  message0('Resampling. Refining ...')
  counter = 1
  while (sum(dol) > 0 && counter < maxIter) {
    ol  = lapply(
      X = output2,
      FUN = function(x) {
        if (!is.null(x)  && nrow(x) > 0) {
          x[1, ]
        } else{
          NULL
        }
      }
    )
    ol2  = abind(ol, along = 1)
    dol  = duplicated(ol2[, 'y.index'])
    if (sum(dol)) {
      output2[dol] = lapply(
        output2[dol],
        FUN = function(x) {
          x[-1, , drop = FALSE]
        }
      )
    }
    counter = counter + 1
  }
  ####################
  message0('Resampling. Creating output ...')
  return(as.data.frame(abind(ol[!unlist(lapply(
    ol,
    FUN = function(x) {
      is.null(x) || length(x) < 1
    }
  ))], along = 1)))
}

Factor2CharAndSubstitution = function(df,
                                      Ind,
                                      column,
                                      ArtifLabel,
                                      neutralise,
                                      mutInd,
                                      baselines) {
  df[, column]         = as.character(df[, column])
  df[Ind, column]      = ArtifLabel
  if (neutralise)      {
    df[mutInd, column] = baselines
  }
  df[, column]         = as.factor(df[, column])
  return(df)
}

mimicControls = function(df                             ,
                         depVariable = 'data_point'     ,
                         ArtifLabel = 'ArtificialMutant',
                         mutLabel = 'experimental'      ,
                         baselines = 'control'          ,
                         sex       = 'sex'              ,
                         neutralise = TRUE              ,
                         removeMutants = FALSE          ,
                         resample  = FALSE              ,
                         minSampRequired                ,
                         SexGenResLevels                ,
                         indicator       = NULL         ,
                         plot            = FALSE        ,
                         CutIfLongerThanMutants = TRUE) {
  if (!is.null(indicator))
    set.seed(indicator)
  requireNamespace('stringi')
  df.bckOrg    = df
  note         = list(
    'Shift' = NULL,
    'New mutant indices' = NULL,
    'Original mutant indices' = NULL
  )
  shift        = Indices = NULL
  #####
  # remove zero frequency categories
  df_rzeros = RemoveZeroFrequencyCategories(
    x = df,
    minSampRequired = minSampRequired,
    depVar = depVariable,
    totalLevels = SexGenResLevels
  )
  df                               = df_rzeros$x
  note$'Removed categories detail' = df_rzeros$note
  ###########
  if (is.null(df)                ||
      length(df) < 1             ||
      nrow  (df) < 1             ||
      nlevels(df[, sex]) < 1     ||
      !is.numeric(df[, depVariable])) {
    note$'Dataset status' = 'Empty dataset'
    return(list(
      df = df          ,
      ctv = NULL       ,
      input = df.bckOrg,
      note  = note
    ))
  }
  # Definind id is crutial otherwise it may end up final duplicated points
  # id must be numeric otherwise will change the entire output to factor! Stu.. R
  # Do not include 0 to have the same length always 021=21!
  df$id_d   = as.numeric(stri_rand_strings(nrow(df), length = 10, pattern = '[1-9]'))
  df$cTime  = Date2Integer(df$date_of_experiment)
  mutInd    = (df$biological_sample_group == mutLabel)
  conInd    = (df$biological_sample_group == baselines)
  if (!sum(mutInd)  ||
      !sum(conInd)) {
    note$'Dataset status'          = 'Missing controls or mutants'
    note$'Original mutant indices' =  mutInd
    return(list(
      df    = df.bckOrg,
      ctv   = NULL     ,
      input = df.bckOrg,
      note  = note
    ))
  }
  ##################
  for (l in levels(df[, sex])) {
    mutIndSex    = (df$biological_sample_group == mutLabel  &
                      df[, sex] == l)
    conIndSex    = (df$biological_sample_group == baselines &
                      df[, sex] == l)
    ###
    y     = df[conIndSex, depVariable]
    time2 = df[conIndSex, 'cTime']
    yid   = df[conIndSex, 'id_d']
    x     = df[mutIndSex, depVariable]
    lx    = length(x)
    if (is.null(shift)) {
      message0(l, ', Generating new SHIFT in progress ...')
      shift = ifelse(resample, round(runif(
        1                                             ,
        min(df$cTime)[1] - min(df$cTime[mutIndSex])[1],
        max(df$cTime)[1] - max(df$cTime[mutIndSex])[1]
      )), 0)
    }
    message0('\t ~> ', l, '\t ~~> Shift = ', round(shift))
    time1 = df$cTime[mutIndSex] + shift
    ####
    ctv = closest.time.value(
      x     = x      ,
      time1 = time1  ,
      y     = y      ,
      time2 = time2  ,
      yind  = yid
    )
    if (!is.null(ctv) && CutIfLongerThanMutants && nrow(ctv) > lx) {
      message0 ('The length of the simulated data is longer than the input. Cutting in progress ....')
      ctv = ctv[1:lx,]
    }
    Indices = c(Indices, ctv$y.index)
  }
  Ind                              =  df$id_d %in% Indices
  note$'New mutant indices'        =  which(Ind)
  note$'shift'                     = shift
  note$'Original mutant indices'   =  mutInd
  note$'Dataset status'            = 'No problem detected'
  # Replace the biological_sample_group and colony_id
  df$colony_id[Ind] = paste0(unique(na.omit(df$colony_id[df$biological_sample_group == mutLabel])), collapse = '')
  df = Factor2CharAndSubstitution(
    df = df,
    Ind = Ind,
    column = 'biological_sample_group',
    ArtifLabel = ArtifLabel,
    neutralise = neutralise,
    mutInd = mutInd,
    baselines = baselines
  )
  ####
  if (plot) {
    plot(
      time2,
      y,
      main = paste0('Shift = ', shift)     ,
      xlim = c(min(ctv$y.time, ctv$x.time) ,
               max(ctv$y.time, ctv$x.time)),
      type = 'p'
    )
    points(df$cTime     ,
           df[, depVariable],
           #col = as.integer(df$biological_sample_group) - 1,
           col = as.integer(as.factor(
             interaction(df$biological_sample_group, df$sex)
           )))
    points(
      df$cTime[mutInd],
      df[mutInd, depVariable],
      pch = paste0(1:length(mutInd)),
      cex = 1,
      lwd = 1,
      col = 2
    )
    points(
      time1,
      x,
      pch = paste0(1:length(x)),
      cex = 1.5,
      lwd = 2,
      col = 3
    )
    points(
      ctv$y.time,
      ctv$y.val,
      pch = paste0(1:length(ctv$y.val)),
      cex = 2,
      lwd = 3,
      col = 4
    )
    legend(
      'top',
      c('TRUE', 'Shifted Mut', 'Closest'),
      fill = c(2, 3, 4),
      horiz = TRUE
    )
  }
  # Must be in the end
  if (removeMutants) {
    df = df[!mutInd, ]
  }

  df$id_d = df$cTime = NULL
  return(list(
    df = droplevels(df),
    ctv = ctv,
    input = df.bckOrg,
    note = note
  ))
}
###
ignoromeGenes = function() {
  list = c(
    "MGI:3512453",
    "MGI:1914868",
    "MGI:3625331",
    "MGI:1347061",
    "MGI:1913332",
    "MGI:1913860",
    "MGI:1926116",
    "MGI:1914135",
    "MGI:2386323",
    "MGI:1919129",
    "MGI:2385289",
    "MGI:1890773",
    "MGI:1924748",
    "MGI:87911",
    "MGI:2661081",
    "MGI:1345162",
    "MGI:2179942",
    "MGI:2429637",
    "MGI:2442875",
    "MGI:1916320",
    "MGI:2151118",
    "MGI:99677",
    "MGI:108450",
    "MGI:1277167",
    "MGI:1925499",
    "MGI:2441950",
    "MGI:2675492",
    "MGI:106675",
    "MGI:2675256",
    "MGI:1197012",
    "MGI:2147658",
    "MGI:2448704",
    "MGI:1924809",
    "MGI:1932075",
    "MGI:87966",
    "MGI:1316648",
    "MGI:1338803",
    "MGI:3041226",
    "MGI:107796",
    "MGI:1919785",
    "MGI:1340024",
    "MGI:1914039",
    "MGI:1914731",
    "MGI:1924753",
    "MGI:1914917",
    "MGI:87997",
    "MGI:2151224",
    "MGI:2444854",
    "MGI:1096344",
    "MGI:104837",
    "MGI:1922680",
    "MGI:1098673",
    "MGI:1337008",
    "MGI:1915673",
    "MGI:1919865",
    "MGI:3052714",
    "MGI:1925726",
    "MGI:101919",
    "MGI:1919020",
    "MGI:2141861",
    "MGI:1929214",
    "MGI:1337060",
    "MGI:1336993",
    "MGI:1930124",
    "MGI:1915566",
    "MGI:88059",
    "MGI:1891066",
    "MGI:99595",
    "MGI:99433",
    "MGI:1917747",
    "MGI:2443687",
    "MGI:1922654",
    "MGI:2441869",
    "MGI:1924919",
    "MGI:1920591",
    "MGI:1860493",
    "MGI:2442308",
    "MGI:3028577",
    "MGI:1915496",
    "MGI:1923959",
    "MGI:1921442",
    "MGI:107511",
    "MGI:1334448",
    "MGI:1913845",
    "MGI:109384",
    "MGI:105121",
    "MGI:1924290",
    "MGI:1929492",
    "MGI:1354735",
    "MGI:104653",
    "MGI:1351597",
    "MGI:2153480",
    "MGI:1859660",
    "MGI:1330826",
    "MGI:108028",
    "MGI:1341628",
    "MGI:1099442",
    "MGI:1919772",
    "MGI:1338011",
    "MGI:894678",
    "MGI:3588200",
    "MGI:1270862",
    "MGI:2387643",
    "MGI:2652819",
    "MGI:1891372",
    "MGI:2143311",
    "MGI:1924210",
    "MGI:1922986",
    "MGI:1332238",
    "MGI:1891828",
    "MGI:2677212",
    "MGI:1338017",
    "MGI:101770",
    "MGI:88169",
    "MGI:2385271",
    "MGI:88180",
    "MGI:1101778",
    "MGI:1916418",
    "MGI:3045315",
    "MGI:1915082",
    "MGI:2146836",
    "MGI:1890651",
    "MGI:1338871",
    "MGI:2442001",
    "MGI:1920594",
    "MGI:1913751",
    "MGI:2141979",
    "MGI:1916433",
    "MGI:1925911",
    "MGI:1923029",
    "MGI:1914576",
    "MGI:88227",
    "MGI:1097680",
    "MGI:3607716",
    "MGI:1914181",
    "MGI:2668347",
    "MGI:88236",
    "MGI:1352750",
    "MGI:1920910",
    "MGI:2444177",
    "MGI:894644",
    "MGI:1316660",
    "MGI:1270839",
    "MGI:1914338",
    "MGI:2182269",
    "MGI:1913208",
    "MGI:2685431",
    "MGI:2179723",
    "MGI:1309469",
    "MGI:1196251",
    "MGI:1924106",
    "MGI:1914327",
    "MGI:107570",
    "MGI:1309992",
    "MGI:105369",
    "MGI:88289",
    "MGI:3512628",
    "MGI:1924122",
    "MGI:1923707",
    "MGI:1289263",
    "MGI:2685134",
    "MGI:1289168",
    "MGI:88340",
    "MGI:1334419",
    "MGI:1336885",
    "MGI:2442676",
    "MGI:106211",
    "MGI:1929745",
    "MGI:1915099",
    "MGI:3588198",
    "MGI:2685856",
    "MGI:894318",
    "MGI:88357",
    "MGI:1921765",
    "MGI:1278336",
    "MGI:1919641",
    "MGI:1332236",
    "MGI:3505689",
    "MGI:1098230",
    "MGI:1349448",
    "MGI:2684927",
    "MGI:1917704",
    "MGI:1923800",
    "MGI:1914244",
    "MGI:108084",
    "MGI:1915511",
    "MGI:1921451",
    "MGI:1919199",
    "MGI:2384581",
    "MGI:1927237",
    "MGI:1915817",
    "MGI:2135796",
    "MGI:99779",
    "MGI:1919386",
    "MGI:1931825",
    "MGI:1913761",
    "MGI:1914185",
    "MGI:2444926",
    "MGI:1346342",
    "MGI:2385186",
    "MGI:1930088",
    "MGI:1929288",
    "MGI:1917912",
    "MGI:88421",
    "MGI:2144529",
    "MGI:1914047",
    "MGI:3643623",
    "MGI:1923428",
    "MGI:1095396",
    "MGI:2155345",
    "MGI:104688",
    "MGI:1916706",
    "MGI:1349400",
    "MGI:1915164",
    "MGI:103226",
    "MGI:105959",
    "MGI:1923953",
    "MGI:2135874",
    "MGI:891996",
    "MGI:109176",
    "MGI:88513",
    "MGI:88490",
    "MGI:1347062",
    "MGI:1860086",
    "MGI:1298216",
    "MGI:1340053",
    "MGI:1351825",
    "MGI:88562",
    "MGI:1914535",
    "MGI:103556",
    "MGI:2387642",
    "MGI:2685586",
    "MGI:1316658",
    "MGI:88607",
    "MGI:1927669",
    "MGI:2183535",
    "MGI:2677061",
    "MGI:1915039",
    "MGI:1931838",
    "MGI:1915337",
    "MGI:1917890",
    "MGI:102563",
    "MGI:2444529",
    "MGI:1343154",
    "MGI:1919297",
    "MGI:1919240",
    "MGI:1346328",
    "MGI:3036254",
    "MGI:1196287",
    "MGI:2442474",
    "MGI:1916442",
    "MGI:1920081",
    "MGI:1918965",
    "MGI:1914737",
    "MGI:1354963",
    "MGI:2447771",
    "MGI:1918478",
    "MGI:1921379",
    "MGI:1914935",
    "MGI:1931881",
    "MGI:1921580",
    "MGI:1915848",
    "MGI:107384",
    "MGI:1261827",
    "MGI:3584043",
    "MGI:1890621",
    "MGI:1913882",
    "MGI:1330238",
    "MGI:1919357",
    "MGI:104627",
    "MGI:1922469",
    "MGI:2685183",
    "MGI:1922715",
    "MGI:1858208",
    "MGI:95281",
    "MGI:1343498",
    "MGI:95284",
    "MGI:1915293",
    "MGI:1924877",
    "MGI:107444",
    "MGI:99252",
    "MGI:1924933",
    "MGI:1916219",
    "MGI:1313286",
    "MGI:1890496",
    "MGI:1343095",
    "MGI:3576783",
    "MGI:1315195",
    "MGI:2142593",
    "MGI:2444896",
    "MGI:106645",
    "MGI:1202295",
    "MGI:1919340",
    "MGI:1201683",
    "MGI:97838",
    "MGI:103582",
    "MGI:1913321",
    "MGI:1890682",
    "MGI:3045306",
    "MGI:1351611",
    "MGI:1916889",
    "MGI:1915376",
    "MGI:2446163",
    "MGI:1921116",
    "MGI:1917613",
    "MGI:2657115",
    "MGI:1921192",
    "MGI:2445194",
    "MGI:3046463",
    "MGI:1923676",
    "MGI:1914000",
    "MGI:1925188",
    "MGI:2385126",
    "MGI:1922869",
    "MGI:1354738",
    "MGI:3039600",
    "MGI:1919429",
    "MGI:1926014",
    "MGI:1354708",
    "MGI:1920223",
    "MGI:95500",
    "MGI:2147790",
    "MGI:99501",
    "MGI:1914362",
    "MGI:95517",
    "MGI:1919764",
    "MGI:1913687",
    "MGI:1336879",
    "MGI:1932127",
    "MGI:1925642",
    "MGI:2443410",
    "MGI:1858193",
    "MGI:3028075",
    "MGI:102949",
    "MGI:1914004",
    "MGI:106315",
    "MGI:1194495",
    "MGI:2442579",
    "MGI:108076",
    "MGI:1890391",
    "MGI:1096879",
    "MGI:1888513",
    "MGI:108571",
    "MGI:108460",
    "MGI:2442040",
    "MGI:95626",
    "MGI:1926176",
    "MGI:2429943",
    "MGI:95678",
    "MGI:3055306",
    "MGI:1341724",
    "MGI:95709",
    "MGI:95710",
    "MGI:1923847",
    "MGI:95722",
    "MGI:95718",
    "MGI:2153041",
    "MGI:1891112",
    "MGI:2685452",
    "MGI:1921748",
    "MGI:107852",
    "MGI:1891703",
    "MGI:1919201",
    "MGI:1289257",
    "MGI:1934765",
    "MGI:2685211",
    "MGI:2139054",
    "MGI:1891463",
    "MGI:2685341",
    "MGI:2685858",
    "MGI:892973",
    "MGI:1346334",
    "MGI:2441763",
    "MGI:105102",
    "MGI:1914555",
    "MGI:102683",
    "MGI:95821",
    "MGI:1351343",
    "MGI:95862",
    "MGI:2385191",
    "MGI:2685307",
    "MGI:2446110",
    "MGI:1344360",
    "MGI:2158340",
    "MGI:1923858",
    "MGI:2444115",
    "MGI:2685817",
    "MGI:106209",
    "MGI:1196297",
    "MGI:1889802",
    "MGI:1352504",
    "MGI:1314872",
    "MGI:1314882",
    "MGI:95929",
    "MGI:95901",
    "MGI:96112",
    "MGI:2677838",
    "MGI:96120",
    "MGI:107159",
    "MGI:101947",
    "MGI:2678023",
    "MGI:96171",
    "MGI:107730",
    "MGI:96193",
    "MGI:96194",
    "MGI:1919862",
    "MGI:1859384",
    "MGI:2685814",
    "MGI:1333853",
    "MGI:96239",
    "MGI:1921627",
    "MGI:3036260",
    "MGI:1330288",
    "MGI:1924292",
    "MGI:1917625",
    "MGI:1858745",
    "MGI:96413",
    "MGI:1924183",
    "MGI:2683287",
    "MGI:2429859",
    "MGI:1915509",
    "MGI:1890359",
    "MGI:107973",
    "MGI:3655979",
    "MGI:99954",
    "MGI:1338071",
    "MGI:1342542",
    "MGI:96538",
    "MGI:109380",
    "MGI:96541",
    "MGI:103014",
    "MGI:1333800",
    "MGI:1890473",
    "MGI:1924375",
    "MGI:105304",
    "MGI:96560",
    "MGI:1917685",
    "MGI:1927753",
    "MGI:1917672",
    "MGI:2442377",
    "MGI:1922168",
    "MGI:1352757",
    "MGI:2677208",
    "MGI:107420",
    "MGI:2429603",
    "MGI:1926262",
    "MGI:1197515",
    "MGI:1197522",
    "MGI:2685110",
    "MGI:109442",
    "MGI:1096361",
    "MGI:96654",
    "MGI:3037820",
    "MGI:1336208",
    "MGI:2685627",
    "MGI:1918269",
    "MGI:2385276",
    "MGI:2145579",
    "MGI:99780",
    "MGI:2667167",
    "MGI:1919347",
    "MGI:2669829",
    "MGI:108426",
    "MGI:1921054",
    "MGI:107688",
    "MGI:109564",
    "MGI:1098269",
    "MGI:109187",
    "MGI:2444612",
    "MGI:2445185",
    "MGI:102849",
    "MGI:107540",
    "MGI:1923714",
    "MGI:103561",
    "MGI:2384899",
    "MGI:3629975",
    "MGI:96705",
    "MGI:3690448",
    "MGI:2143628",
    "MGI:109321",
    "MGI:1913758",
    "MGI:1890494",
    "MGI:104808",
    "MGI:104576",
    "MGI:96759",
    "MGI:2685031",
    "MGI:96778",
    "MGI:1891214",
    "MGI:1924819",
    "MGI:1915671",
    "MGI:1353635",
    "MGI:1921392",
    "MGI:108429",
    "MGI:108424",
    "MGI:96828",
    "MGI:1919666",
    "MGI:1916956",
    "MGI:1342770",
    "MGI:2389177",
    "MGI:104797",
    "MGI:96836",
    "MGI:109151",
    "MGI:99502",
    "MGI:1914113",
    "MGI:2443598",
    "MGI:1346867",
    "MGI:1346877",
    "MGI:2444554",
    "MGI:2444136",
    "MGI:2443731",
    "MGI:1333811",
    "MGI:1333813",
    "MGI:96924",
    "MGI:1920977",
    "MGI:1858420",
    "MGI:1916245",
    "MGI:98446",
    "MGI:1922863",
    "MGI:1914249",
    "MGI:1917967",
    "MGI:1924140",
    "MGI:106477",
    "MGI:1914819",
    "MGI:1913697",
    "MGI:1918127",
    "MGI:1918398",
    "MGI:1914277",
    "MGI:1924015",
    "MGI:1347361",
    "MGI:1919891",
    "MGI:1336894",
    "MGI:97052",
    "MGI:2146995",
    "MGI:1927340",
    "MGI:1924265",
    "MGI:1196612",
    "MGI:1913743",
    "MGI:1915822",
    "MGI:1915985",
    "MGI:1913542",
    "MGI:2684990",
    "MGI:1329021",
    "MGI:1342005",
    "MGI:1338850",
    "MGI:1923616",
    "MGI:1928394",
    "MGI:1915485",
    "MGI:1915364",
    "MGI:107624",
    "MGI:101785",
    "MGI:2183924",
    "MGI:1919192",
    "MGI:1341430",
    "MGI:1915241",
    "MGI:1915896",
    "MGI:2138939",
    "MGI:2180167",
    "MGI:1196326",
    "MGI:1289164",
    "MGI:109186",
    "MGI:1929915",
    "MGI:109166",
    "MGI:1341799",
    "MGI:1914523",
    "MGI:2443241",
    "MGI:1933754",
    "MGI:1920024",
    "MGI:101784",
    "MGI:3043305",
    "MGI:1340031",
    "MGI:1098547",
    "MGI:1352751",
    "MGI:1278343",
    "MGI:2449121",
    "MGI:1915074",
    "MGI:1341898",
    "MGI:108011",
    "MGI:97363",
    "MGI:99460",
    "MGI:1924833",
    "MGI:107605",
    "MGI:1860130",
    "MGI:1858233",
    "MGI:2444210",
    "MGI:97376",
    "MGI:2183436",
    "MGI:2385017",
    "MGI:104750",
    "MGI:2443642",
    "MGI:2142581",
    "MGI:1915243",
    "MGI:1351500",
    "MGI:2180849",
    "MGI:97402",
    "MGI:1919325",
    "MGI:2139018",
    "MGI:1328337",
    "MGI:1930265",
    "MGI:1916308",
    "MGI:2147616",
    "MGI:1918448",
    "MGI:1919451",
    "MGI:109147",
    "MGI:2137026",
    "MGI:1921441",
    "MGI:2385206",
    "MGI:1339975",
    "MGI:2679420",
    "MGI:1351623",
    "MGI:1913699",
    "MGI:1916207",
    "MGI:1347094",
    "MGI:1340806",
    "MGI:2143990",
    "MGI:2685692",
    "MGI:1891832",
    "MGI:1914933",
    "MGI:1923998",
    "MGI:1855700",
    "MGI:1861733",
    "MGI:1333863",
    "MGI:97524",
    "MGI:1916047",
    "MGI:104864",
    "MGI:2384308",
    "MGI:2140945",
    "MGI:97535",
    "MGI:2677270",
    "MGI:1344408",
    "MGI:1914495",
    "MGI:1924963",
    "MGI:2181202",
    "MGI:2441675",
    "MGI:2446138",
    "MGI:2443812",
    "MGI:2140712",
    "MGI:2156864",
    "MGI:2444341",
    "MGI:97577",
    "MGI:97578",
    "MGI:1916211",
    "MGI:891978",
    "MGI:1934659",
    "MGI:1261910",
    "MGI:3039607",
    "MGI:1203729",
    "MGI:1098772",
    "MGI:1916193",
    "MGI:2152214",
    "MGI:1298224",
    "MGI:1197006",
    "MGI:1916867",
    "MGI:1858231",
    "MGI:104747",
    "MGI:1351327",
    "MGI:2445289",
    "MGI:1333782",
    "MGI:2443207",
    "MGI:1270860",
    "MGI:2154240",
    "MGI:1859214",
    "MGI:1920328",
    "MGI:2180564",
    "MGI:97740",
    "MGI:1097163",
    "MGI:2675617",
    "MGI:97742",
    "MGI:106687",
    "MGI:1913411",
    "MGI:1926321",
    "MGI:105086",
    "MGI:104740",
    "MGI:1888712",
    "MGI:1917475",
    "MGI:104871",
    "MGI:1919362",
    "MGI:2685193",
    "MGI:2387581",
    "MGI:1321161",
    "MGI:1920949",
    "MGI:2442104",
    "MGI:107171",
    "MGI:102666",
    "MGI:97753",
    "MGI:2444067",
    "MGI:2685844",
    "MGI:103294",
    "MGI:97602",
    "MGI:2139971",
    "MGI:3043083",
    "MGI:1914479",
    "MGI:1914171",
    "MGI:3605764",
    "MGI:1923810",
    "MGI:2183441",
    "MGI:1914339",
    "MGI:1913284",
    "MGI:1914248",
    "MGI:1096367",
    "MGI:1095405",
    "MGI:102966",
    "MGI:97801",
    "MGI:97807",
    "MGI:102694",
    "MGI:1277956",
    "MGI:2448514",
    "MGI:2385107",
    "MGI:104606",
    "MGI:1915040",
    "MGI:1201692",
    "MGI:2442633",
    "MGI:3652039",
    "MGI:1922896",
    "MGI:1921584",
    "MGI:1921262",
    "MGI:1923596",
    "MGI:1918573",
    "MGI:1914060",
    "MGI:1099460",
    "MGI:1861774",
    "MGI:97879",
    "MGI:96522",
    "MGI:1931028",
    "MGI:109406",
    "MGI:97897",
    "MGI:1098222",
    "MGI:1914692",
    "MGI:1913464",
    "MGI:1918996",
    "MGI:2179276",
    "MGI:1921984",
    "MGI:1918325",
    "MGI:97930",
    "MGI:1346341",
    "MGI:1919206",
    "MGI:1261771",
    "MGI:2442859",
    "MGI:1915045",
    "MGI:104661",
    "MGI:2141142",
    "MGI:98038",
    "MGI:98105",
    "MGI:1921620",
    "MGI:1924467",
    "MGI:2136886",
    "MGI:1917682",
    "MGI:1858752",
    "MGI:1913771",
    "MGI:2153463",
    "MGI:1339467",
    "MGI:98223",
    "MGI:1933169",
    "MGI:2388100",
    "MGI:2679336",
    "MGI:2135937",
    "MGI:1919443",
    "MGI:103033",
    "MGI:98250",
    "MGI:1349165",
    "MGI:1914195",
    "MGI:2148802",
    "MGI:1349635",
    "MGI:1924621",
    "MGI:1916858",
    "MGI:1329016",
    "MGI:1916941",
    "MGI:1931466",
    "MGI:1100878",
    "MGI:98279",
    "MGI:1349457",
    "MGI:3041197",
    "MGI:1922997",
    "MGI:1933199",
    "MGI:2136890",
    "MGI:1935121",
    "MGI:2137677",
    "MGI:1350341",
    "MGI:1098703",
    "MGI:2446215",
    "MGI:1920973",
    "MGI:98297",
    "MGI:3605641",
    "MGI:98299",
    "MGI:1195268",
    "MGI:2445031",
    "MGI:108563",
    "MGI:1927664",
    "MGI:1915596",
    "MGI:1913390",
    "MGI:1099835",
    "MGI:1351663",
    "MGI:108402",
    "MGI:1201406",
    "MGI:3037150",
    "MGI:1929691",
    "MGI:2156052",
    "MGI:1353479",
    "MGI:104516",
    "MGI:1345283",
    "MGI:1353498",
    "MGI:1928369",
    "MGI:2442682",
    "MGI:1931249",
    "MGI:2385166",
    "MGI:2140361",
    "MGI:1919305",
    "MGI:2443383",
    "MGI:1915010",
    "MGI:2150150",
    "MGI:1890216",
    "MGI:1270850",
    "MGI:1336891",
    "MGI:94862",
    "MGI:1921337",
    "MGI:1351872",
    "MGI:2679449",
    "MGI:95453",
    "MGI:1859183",
    "MGI:1339795",
    "MGI:1919247",
    "MGI:1916186",
    "MGI:1096393",
    "MGI:109356",
    "MGI:1915076",
    "MGI:2139270",
    "MGI:1919232",
    "MGI:1916274",
    "MGI:1917729",
    "MGI:2685966",
    "MGI:1923992",
    "MGI:1921728",
    "MGI:1924574",
    "MGI:1927715",
    "MGI:2444120",
    "MGI:1927170",
    "MGI:1915196",
    "MGI:1923823",
    "MGI:1924834",
    "MGI:2446175",
    "MGI:1913433",
    "MGI:107931",
    "MGI:1344414",
    "MGI:1930252",
    "MGI:1891338",
    "MGI:105372",
    "MGI:2445190",
    "MGI:1341828",
    "MGI:2448556",
    "MGI:1915678",
    "MGI:1923396",
    "MGI:109355",
    "MGI:1890156",
    "MGI:1342296",
    "MGI:102928",
    "MGI:1923723",
    "MGI:1919445",
    "MGI:2441711",
    "MGI:1928849",
    "MGI:1921325",
    "MGI:99515",
    "MGI:1354961",
    "MGI:2153070",
    "MGI:2685071",
    "MGI:1196415",
    "MGI:2443028",
    "MGI:98483",
    "MGI:3576210",
    "MGI:1916222",
    "MGI:2144164",
    "MGI:1919488",
    "MGI:109567",
    "MGI:1925082",
    "MGI:98640",
    "MGI:1928486",
    "MGI:2144865",
    "MGI:1918968",
    "MGI:1914846",
    "MGI:103270",
    "MGI:3039623",
    "MGI:1913775",
    "MGI:2441683",
    "MGI:2669033",
    "MGI:1098686",
    "MGI:2444222",
    "MGI:1919150",
    "MGI:1932411",
    "MGI:2136977",
    "MGI:2147810",
    "MGI:2685030",
    "MGI:2142624",
    "MGI:2443597",
    "MGI:1919899",
    "MGI:106402",
    "MGI:2442082",
    "MGI:98775",
    "MGI:1921050",
    "MGI:894675",
    "MGI:2673064",
    "MGI:1930958",
    "MGI:2446193",
    "MGI:105070",
    "MGI:2450248",
    "MGI:1197527",
    "MGI:3582693",
    "MGI:2181659",
    "MGI:2138319",
    "MGI:2182472",
    "MGI:98797",
    "MGI:1920198",
    "MGI:1336209",
    "MGI:1913476",
    "MGI:3029307",
    "MGI:1930005",
    "MGI:1918576",
    "MGI:1328317",
    "MGI:106657",
    "MGI:1914199",
    "MGI:2442815",
    "MGI:2685973",
    "MGI:96270",
    "MGI:1330305",
    "MGI:1919383",
    "MGI:106244",
    "MGI:2384576",
    "MGI:2141418",
    "MGI:107848",
    "MGI:1329045",
    "MGI:1916092",
    "MGI:1933134",
    "MGI:1346078",
    "MGI:2145316",
    "MGI:98878",
    "MGI:1913355",
    "MGI:1914378",
    "MGI:1918957",
    "MGI:2443123",
    "MGI:1913405",
    "MGI:1915384",
    "MGI:1858178",
    "MGI:1888998",
    "MGI:1923429",
    "MGI:2444541",
    "MGI:1321389",
    "MGI:1916165",
    "MGI:102718",
    "MGI:2143698",
    "MGI:1913435",
    "MGI:2444304",
    "MGI:2443189",
    "MGI:1261847",
    "MGI:3642995",
    "MGI:104630",
    "MGI:107577",
    "MGI:1337100",
    "MGI:1927241",
    "MGI:1917819",
    "MGI:2685541",
    "MGI:2142282",
    "MGI:2442327",
    "MGI:1914258",
    "MGI:109484",
    "MGI:1922830",
    "MGI:2140248",
    "MGI:99187",
    "MGI:1918025",
    "MGI:1340045",
    "MGI:1917140",
    "MGI:2443465",
    "MGI:1929117",
    "MGI:1306812",
    "MGI:2442788",
    "MGI:1918381",
    "MGI:1914233",
    "MGI:3045312",
    "MGI:1919153",
    "MGI:2664358",
    "MGI:1351661",
    "MGI:1920597",
    "MGI:1922971",
    "MGI:101901",
    "MGI:3624119",
    "MGI:1919908",
    "MGI:1889619",
    "MGI:1920701",
    "MGI:1925378",
    "MGI:1915598",
    "MGI:2384969",
    "MGI:2388711",
    "MGI:2145458",
    "MGI:3043522",
    "MGI:1345189",
    "MGI:1860484",
    "MGI:2446235",
    "MGI:1309515",
    "MGI:2157522",
    "MGI:1917351",
    "MGI:3040701",
    "MGI:1270857",
    "MGI:2683295",
    "MGI:1915524",
    "MGI:2676828",
    "MGI:2676874",
    "MGI:3619440",
    "MGI:1919828",
    "MGI:3030226",
    "MGI:1913300",
    "MGI:2681306",
    "MGI:1921812",
    "MGI:2681253",
    "MGI:3619366",
    "MGI:3030932"
  )
  return(list)
}


requiredDataColumns = function(x){
  ColumnsList = c(
    'allele_accession_id',
    'metadata',
    'gene_accession_id',
    'project_name',
    'genetic_background',
    'strain_accession_id',
    'litter_id',
    'phenotyping_center',
    'time_point',
    'external_sample_id',
    'observation_id',
    'developmental_stage_name',
    'datasource_name',
    'procedure_group',
    'pipeline_stable_id',
    'parameter_stable_id',
    'age_in_days',
    'date_of_experiment',
    'weight',
    'pipeline_name',
    'procedure_stable_id',
    'observation_type',
    'developmental_stage_acc',
    'procedure_name',
    'date_of_birth',
    'gene_symbol',
    'metadata_group',
    'biological_sample_group',
    'discrete_point',
    'experiment_source_id',
    'data_point',
    'sex',
    'production_center',
    'colony_id',
    'parameter_name',
    'allele_symbol',
    'age_in_weeks',
    'zygosity',
    'strain_name',
    'data_type',
    'category'
  )
  return(ColumnsList)
}

updateEquationMap = function(updatePackage = TRUE) {
  # Extract Centers name from the skip list
  metapars = read.csv(file = file.path(local(), 'metadataParameters.csv'))
  centers = lapply(strsplit(
    metapars$parameter_stable_id,
    split = '_',
    fixed = TRUE
  ), function(x) {
    return (x[1])
  })
  centers = unique(unlist(centers))

  # Extract parameters from the method map
  methodmap = readConf('EquationMap.conf')
  orgparsColon = paste(names(methodmap), methodmap, sep = ':')
  orgpars = paste(names(methodmap), methodmap, sep = '_')
  parameters = lapply(strsplit(orgpars,
                               split = '_',
                               fixed = TRUE), function(x) {
                                 if (length(x) <= 5) {
                                   r = c()
                                   for (i in 1:6) {
                                     y = x
                                     y[length(y) - 1] = paste(paste0('00', i, collapse = ''), y[length(y)], sep=':',collapse = ':')
                                     r = c(r, paste(y[-c(1,length(y))], sep = '_', collapse = '_'))
                                   }
                                   return (r)
                                 } else{
                                   return(NULL)
                                 }
                               })
  parameters = unlist(parameters)

  #combine two lists
  #r = apply(expand.grid(centers, parameters), 1, paste, collapse = "_")
  r = as.vector(outer(centers, parameters, paste, sep="_"))
  r = unique(c(orgparsColon, sort(r)))
  if (updatePackage) {
    write.table(
      r,
      file = system.file("extdata", "EquationMap.conf", package = "DRrequiredAgeing"),
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  } else{
    write.table(
      r,
      file = file.path(getwd(), "EquationMap.conf"),
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  }
  gc()
  return(invisible(r))
}

updateMethodMap = function(updatePackage = TRUE) {
  # Extract Centers name from the skip list
  metapars = read.csv(file = file.path(local(), 'metadataParameters.csv'))
  centers = lapply(strsplit(
    metapars$parameter_stable_id,
    split = '_',
    fixed = TRUE
  ), function(x) {
    return (x[1])
  })
  centers = unique(unlist(centers))

  # Extract parameters from the method map
  methodmap = readConf('MethodMap.conf')
  orgparsColon = paste(names(methodmap), methodmap, sep = ':')
  orgpars = paste(names(methodmap), methodmap, sep = '_')
  parameters = lapply(strsplit(orgpars,
                               split = '_',
                               fixed = TRUE), function(x) {
                                 if (length(x) <= 5) {
                                   r = c()
                                   for (i in 1:6) {
                                     y = x
                                     y[length(y) - 1] = paste(paste0('00', i, collapse = ''), y[length(y)], sep=':',collapse = ':')
                                     r = c(r, paste(y[-c(1,length(y))], sep = '_', collapse = '_'))
                                   }
                                   return (r)
                                 } else{
                                   return(NULL)
                                 }
                               })
  parameters = unlist(parameters)

  #combine two lists
  #r = apply(expand.grid(centers, parameters), 1, paste, collapse = "_")
  r = as.vector(outer(centers, parameters, paste, sep="_"))
  r = unique(c(orgparsColon, sort(r)))
  if (updatePackage) {
    write.table(
      r,
      file = system.file("extdata", "MethodMap.conf", package = "DRrequiredAgeing"),
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  } else{
    write.table(
      r,
      file = file.path(getwd(), "MethodMap.conf"),
      row.names = FALSE,
      col.names = FALSE,
      quote = FALSE
    )
  }
  gc()
  return(invisible(r))
}

updateImpress = function(updateImpressFileInThePackage = FALSE,
                         updateOptionalParametersList  = FALSE,
                         updateTheSkipList             = FALSE,
                         saveRdata                     = NULL ,
                         updateMethodMapHuristic       = TRUE,
                         updateEquationMapHuristic     = TRUE
                         ) {
  outP = df = NULL
  requireNamespace('pingr')
  requireNamespace('jsonlite')
  if (!pingr::is_online()) {
    stop(
      'You must be connected to the internet to be able to update the categorical categories from the IMPReSS ...'
    )
  }
  ###############
  message0('Updating the IMPReSS categories in progress ...')
  message0('\t Step 1. Getting the list of parameters ...')
  counter   =  1
  startTime = Sys.time()
  df = data.frame('parameterKey' = character(),
                  optionCollection = character())
  pipelineList = jsonlite:::fromJSON(txt = 'http://api.mousephenotype.org/impress/pipeline/list')
  for (pipelineId in names(pipelineList)) {
    message0('Pipeline id: ', pipelineId)
    ProcedureList = jsonlite:::fromJSON(
      txt = paste0(
        'http://api.mousephenotype.org/impress/procedure/belongingtopipeline/keys/',
        pipelineId
      )
    )
    for (procedureId in names(ProcedureList)) {
      message0('\t Procedure id: ', procedureId)
      parameterOptions = jsonlite:::fromJSON(
        txt = paste0(
          'http://api.mousephenotype.org/impress/parameter/belongingtoprocedure/full/',
          procedureId
        )
      )
      if (is.null(dim(parameterOptions)))
        next

      if (counter < 2) {
        df =  parameterOptions
      } else{
        df = rbind(df , parameterOptions)
      }
      counter = counter + 1

    }
  }
  ###################################################
  if (!is.null(saveRdata) && saveRdata)
    save(df, file = paste0(Sys.Date(), '_', saveRdata, '_Impress.Rdata'))
  ###################################################

  if (updateTheSkipList)
    UpdateTheSkipListfromIMPReSSAPI(df)
  ###################################################
  message0('\t Step2. Fetching the category names from the category ids ...')
  dfSelected  = df[lapply(df$optionCollection, length) > 0,]
  #dfSelected  = dfSelected[dfSelected$isAnnotation, ]
  dfSelected  = dfSelected[dfSelected$type %in% 'simpleParameter', ]
  #dfSelected  = dfSelected[dfSelected$valueType %in% 'TEXT', ]
  dfSelected  = dfSelected[, c('parameterKey', 'optionCollection', 'parameterId')]
  dfSelected  = dfSelected[!duplicated(dfSelected$parameterKey),]

  message0('\t\t Total items to look up: ', nrow(dfSelected))
  dfSelected$categories = sapply(dfSelected$parameterId, function(x) {
    #message0('Pid = ', x)
    l = unlist(jsonlite:::fromJSON(
      paste0(
        'http://api.mousephenotype.org/impress/option/belongingtoparameter/names/',
        x
      )
    ))
    paste(trimws(l), collapse = ',', sep = ',')
  })

  ###################################################
  message0('Finished in ', round(difftime(Sys.time() , startTime, units = 'min'), 2), 'm')
  ###################################################
  if (updateImpressFileInThePackage) {
    fileName = system.file("extdata", "AllCts.csv", package = "DRrequiredAgeing")
  } else{
    fileName = file.path(getwd(), 'AllCts.csv')
  }
  ###################################################
  if(updateOptionalParametersList && updateImpressFileInThePackage){
    fileNameMeta = system.file("extdata", "metadataParameters.csv", package = "DRrequiredAgeing")
  } else{
    fileNameMeta = file.path(getwd(), 'metadataParameters.csv')
  }
  ###################################################
  outP = data.frame(
    parameter_stable_id = dfSelected$parameterKey,
    categories          = dfSelected$categories
  )
  ###################################################
  outM = data.frame(
    parameter_stable_id = unique(df$parameterKey[df$type %in% 'procedureMetadata'])
  )
  ###################################################
  message0('\tThe output file:\n\t  => ', fileName)
  write.csv(x         = outP    ,
            file      = fileName,
            row.names = FALSE)
  ###################################################
  message0('\tThe metadatafile file:\n\t  => ', fileNameMeta)
  write.csv(x         = outM    ,
            file      = fileNameMeta,
            row.names = FALSE)
  ###################################################
  if(updateMethodMapHuristic){
    updateMethodMap(updateImpressFileInThePackage)
  }
  ###################################################
  if(updateEquationMapHuristic){
    updateEquationMap(updateImpressFileInThePackage)
  }
  ###################################################
  return(invisible(list(
    categories = outP, dfObject = df
  )))

}

UpdateTheSkipListfromIMPReSSAPI = function(df) {
  if (is.null(df))
    message0('There is an error in the input data. Please check the IMPReSS website is on!')
  dfSkPar = subset(df, df$isAnnotation == FALSE)
  dfSkPar = dfSkPar[!duplicated(dfSkPar$parameterKey),]
  ###########
  fileName = system.file("extdata", "ExceptionMap.list", package = "DRrequiredAgeing")
  message0('Reading the current skiplist in the package ...\n\t Path = ',
           fileName)
  CurrentSkipList = readLines(con =  fileName)
  ###########
  NewSkipList  = c(CurrentSkipList, dfSkPar$parameterKey)
  NewSkipList  = NewSkipList[!duplicated(NewSkipList)]
  message0('Writting the output back to the path,\n\t', fileName)
  writeLines(text = NewSkipList, con = fileName)
}

CreateVirtualDrive = function(active = FALSE, currentwd = NULL) {
  wd = ifelse(is.null(currentwd), getwd(), currentwd)
  if (.Platform$OS.type != 'windows') {
    message0('Virtual drive only works for windows OS.')
    return(wd)
  }
  if (active) {
    message0('Creating a virtual drive ... ')
    system('subst U: /D', wait = TRUE)
    if (system(paste0('subst U: "', wd, '"'), wait = TRUE) < 1) {
      message0('Virtual directory successfully created.')
      wd = 'U:'
    } else{
      message0('Cannot create the virtual drive. It may already exist.')
    }
  }
  return(wd)
}



readInputDatafromFile = function(file = NULL,
                                 checkNamesForMissingColNames = TRUE,
                                 sep          = ',',
                                 na.strings   = 'NA') {
  message0('Reading the input file ...\n\t ~> ', file)
  if (!is.null(file))
    file = gsub(
      pattern = '//',
      replacement = '/',
      x = file,
      fixed = TRUE
    )
  if (!file.exists(file))
    message0('File is not local or does not exist!')
  message0('Reading the input data ...')
  if (!grepl(pattern = '.Rdata',
             x = head(file, 1),
             fixed = TRUE)) {
    rdata = read.csv(
      file = file                                    ,
      check.names      = checkNamesForMissingColNames,
      sep              = sep                         ,
      na.strings       = na.strings                  ,
      stringsAsFactors = TRUE
    )
  } else{
    message0('\tReading the input Rdata ...')
    loadfile = load(file = file)
    if (length(loadfile) < 1)
      stop('The loaded Rdata is blank ...')
    rdata = get(loadfile[1])
  }
  message0('Input file dimentions: ',
           paste0(dim(rdata), collapse  = ', '))
  rdata = removeLeadingSpaceFromDataFrameFactors(rdata)
  return(rdata)
}

removeLeadingSpaceFromDataFrameFactors = function(x) {
  if (is.null(x))
    return(x)
  message0('Remove leading/trailins space from a data frame ...')
  message0 ('\t Note that the factors get to the charachter and t/l space will be removed ...')
  message0 ('\t The order of the factors may not be preserved!')
  #x = as.data.frame(x)
  catV = unname(unlist(sapply(x, function(xx) {
    return(is.factor(xx))
  })))
  if (!sum(catV))
    return(x)

  for (i in which(catV == TRUE)) {
    x[, i] = as.factor(trimws(as.character(x[, i])))
  }
  return(x)
}

RecordSpentTime = function(timeSt               ,
                           dirName  = 'ParaTime',
                           fileName = 'file.txt',
                           rnd      = NULL      ,
                           active   = TRUE) {
  if(!active)
    return(NULL)
  message0('Recording time ...')
  if (!dir.exists(dirName))
    dir.create(dirName)
  r = round(difftime(Sys.time() , timeSt, units = 'sec'), 2)
  write(
    c(fileName, r)                                                ,
    file = file.path(dirName, paste(
      rnd,
      paste(fileName, sep = '_', collapse = '_'),
      '.txt',
      collapse = '_',
      sep = '_'
    ))                                                            ,
    ncolumns = 10 ^ 3                                             ,
    append = TRUE                                                 ,
    sep = '\t'
  )
}

TransformVariableByFunction = function(varType = NULL ,
                                       types    = c('unidimensional', 'time_series'),
                                       ###
                                       data    = NULL ,
                                       colName = NULL ,
                                       FUN     = as.numeric0,
                                       FUNData = is.factor) {
  if (is.null(data))
    return(data)
  if (varType %in% types       &&
      colName %in% names(data) &&
      FUNData(data[, colName])) {
    message0('Applying the function to the data [column = ', colName, ']')
    #data[, colName] = FUN(levels(data[, colName])[data[, colName]])
    data[, colName] = FUN(data[, colName])
  }
  return(data)
}

MissingPercent = function(var = NULL, data = NULL) {
  MissPercent = 1
  if (length(var) > 0  &&
      !is.null(data)   &&
      nrow(data) > 0   &&
      var %in% names(data)) {
    MissPercent = sum(is.na(data[, var])) / length(data[, var])
    message0('The percentage of missings in variable `',
             var,
             '`: ',
             round(MissPercent * 100),
             '%')
  }
  return(MissPercent)
}


EnoughWeightForTheSexGenInteraction = function(df,
                                               cols = c('Sex'          ,
                                                        'Genotype')    ,
                                               weightCol = 'Weight'    ,
                                               thresh    = 4) {
  if (is.null(df) || nrow(df) < 1)
    return(FALSE)
  if (!weightCol %in% names(df) || all(!cols %in% names(df)))
    return(FALSE)

  GSinter = tapply(X = df[, weightCol], INDEX = interaction(df[, cols[cols %in% names(df)]]), function(x) {
    length(na.omit(x))
  }, default = 0)

  if (min(GSinter) < thresh) {
    message0(
      'Less than ',
      thresh,
      ' sample in the ',
      paste(cols, collapse = '-'),
      ' table for ',
      weightCol,
      '\n\t ~~> min = ',
      min(GSinter)
    )
    return(FALSE)
  } else{
    return(TRUE)
  }
}

list.dirsDepth = function(path  = getwd(),
                          depth = 0      ,
                          cumulative = FALSE,
                          ...) {
  dirs = path
  if (depth > 0)
    for (i in 1:depth) {
      path   = list.dirs(path = path, recursive = FALSE, ...)
      dirs   = if(cumulative)  c(path, dirs) else  path
    }
  return(unique(dirs))
}


minijobsCreator = function(path  = getwd(),
                           depth = 2,
                           fname = 'minijobs.bch',
                           type = '*.tsv') {
  lf = list.dirsDepth(path = path, depth = depth)
  a = paste0(
    'sbatch --job-name=impc_stats_pipeline_job --mem=1G --time=2-00',
    ' -e ', basename(lf), '_error.err',
    ' -o ', basename(lf), '_output.log --wrap="find ',
    lf,
    ' -type f -name "',
    type,
    '" > ',
    path,
    '/FileIndex_',
    basename(lf),
    '_',
    round(runif(length(lf)), 6),
    '.Ind"'
  )
  write(
    x = a,
    file = fname,
    ncolumns = 10 ^ 5,
    append = FALSE,
    sep = '\n'
  )
}

dictionary2listConvert = function(x) {
  if (is.null(x) || !is(x, 'list'))
    return(x)
  x = as.list(x)
  r2 = lapply(names(x), function(name) {
    r = list(category = name, data = lapply(x[name], function(y) {
      return(y)
    }))
    return(r)
  })
  return(r2)
}


waitTillCommandFinish = function(checkcommand = 'squeue --format="%A %.30j"',
                                 WaitIfTheOutputContains = 'impc_stats_pipeline_job',
                                 WaitBeforeRetrySec = 60,
                                 ignoreline = 0,
                                 ...) {
  r = system(command = checkcommand,
             wait = TRUE,
             intern = TRUE,
             ...)
  totalSeconds = 0
  if (any(ignoreline > 0))
    r = r[-c(ignoreline)]
  while (length(r) > 0 && any(grepl(
    pattern = WaitIfTheOutputContains,
    x = r,
    perl = TRUE,
    ignore.case = FALSE
  ))) {
    totalSeconds = totalSeconds + WaitBeforeRetrySec
    message0('waiting for ',
             WaitBeforeRetrySec,
             's',
             '. Total waiting time: ',
             totalSeconds,
             's.')
    Sys.sleep(WaitBeforeRetrySec)
    r = system(command = checkcommand,
               wait = TRUE,
               intern = TRUE,
               ...)
    if (ignoreline > 0)
      r = r[-c(ignoreline)]

  }
  message0('continuing the pipeline ...')
}

submit_limit_jobs = function(bch_file,
                             job_id_logfile,
                             max_jobs=600) {
  message0("Start submit_limit_jobs")
  system(paste("echo > ", job_id_logfile))
  file <- file(bch_file, "r")
  while (length(command <- readLines(file, n = 1, warn = FALSE)) > 0) {
    while(TRUE) {
      num_running <- as.integer(system("squeue | wc -l", wait=TRUE, intern = TRUE))
      if(num_running <= max_jobs) {
        break
      }
      Sys.sleep(0.5)
    }
    job_id <- system(command, wait=TRUE, intern = TRUE)
    system(paste("echo '", job_id, "' >> ", job_id_logfile))
  }
  close(file)
  message0("End submit_limit_jobs")
}

filesContain = function(path = getwd(),
                        extension = NULL,
                        containWhat = 'Exit',
                        ...) {
  res = FALSE
  files = list.files(
    path = path,
    pattern = extension,
    all.files = TRUE,
    full.names = TRUE,
    include.dirs = FALSE,
    recursive = TRUE,
    ...
  )
  for (file in files) {
    message0('checking for term "', containWhat, '":', file)
    fcontent = readLines(file, warn = FALSE)
    for (l in fcontent) {
      if (is.null(l))
        next
      if (grepl(pattern = containWhat, x = l))
        return(TRUE)
    }
    message0('\t Passed ...')
  }

  return(res)
}

jobCreator = function(path = getwd(),
                      pattern = "\\.Rdata",
                      JobListFile = "DataGenerationJobList.bch",
                      mem = "45G",
                      time = "6-00") {
  if (!dir.exists(file.path("DataGeneratingLog")))
    dir.create(
      path = file.path("DataGeneratingLog"),
      recursive = TRUE,
      showWarnings = FALSE
    )
  if (file.exists(JobListFile))
    unlink(JobListFile)
  files = list.files(
    path = path,
    pattern = pattern,
    full.names = TRUE,
    recursive = FALSE,
    include.dirs = FALSE,
    all.files = TRUE
  )
  proc = tools::file_path_sans_ext(basename(files))
  write(
    paste0(
      "sbatch --job-name=impc_stats_pipeline_job --mem=", mem,
      " --time=", time,
      ' -e ',
      file.path('DataGeneratingLog', paste0(proc, '_errorlog.log')),
      ' -o ',
      file.path('DataGeneratingLog', paste0(proc, '_outputlog.log')),
      " --wrap='Rscript InputDataGenerator.R ",
      '"',
      files,
      '" "',
      proc,
      '"',
      "'"
    ),
    file = JobListFile
  )
}


packageBackup = function(package = NULL,
                         storepath = NULL,
                         flags = '-r9Xq') {
  if (is.null(package)) {
    message('please specify the package first ...')
    return(NULL)
  }
  for (path in .libPaths()) {
    ppath = file.path(path, package)

    if (dir.exists(ppath)) {
      if (!is.null(storepath) && !dir.exists(storepath))
        dir.create(storepath, recursive = TRUE)
      zipfile = file.path(ifelse(is.null(storepath), getwd(), storepath),
                          paste0(
                            package,
                            '_',
                            format(Sys.time(), '%Y-%m-%d %H_%M_%S'),
                            '.zip'
                          ))
      message('\t', package, ' -> backup created in: ', zipfile)
      zip(zipfile = zipfile  ,
          files = ppath,
          flags = flags)

    }
  }
}

ReplaceWordInFile = function(file,
                             pattern = '',
                             replaceBy = '',
                             fixed = TRUE,
                             ...) {
  if (!file.exists(file)) {
    message0('file does not exists. file \n\t =>', file)
    return (invisible(file))
  }
  tx  <- readLines(file)
  tx2  <-
    gsub(
      pattern = pattern,
      replace = replaceBy,
      x = tx,
      fixed = fixed,
      ...
    )
  writeLines(tx2, con = file)
  message0('Replace in file successful. Pattern: ',
           pattern,
           ', Replaced by: ',
           replaceBy)
  return (invisible(file))
}

StatsPipeline = function(path = getwd(),
                         SP.results = file.path(getwd(), 'SP'),
                         waitUntillSee = 'impc_stats_pipeline_job',
                         ignoreThisLineInWaitingCheck = 0,
                         windowingPipeline = TRUE,
                         DRversion = 'not_specified') {
  startTime = Sys.time()
  message0('Starting the IMPC statistical pipeline ... ')
  message0('\t Parquet files path:  ', path)
  message0('\t Output path: ', SP.results)
  ###############################################
  path0 = path
  path = SP.results
  if (!dir.exists(path))
    dir.create(path)
  setwd(path)
  currentuser = Sys.info()['user']

  ### Phase I: Preparing parquet files
  ###############################################
  message0('Phase I. Converting parquet files into Rdata ...')
  message0('Step 1. Reading the data from the parquets directory and creating jobs')
  source(file.path(local(),
                   'StatsPipeline/0-ETL/Step1MakePar2RdataJobs.R'))
  f(path0)
  rm0('f')
  ###############################################
  message0('Step 2. Reading the data from the parquets files and creating psudo Rdata')
  file.copy(
    from = file.path(local(),
                     'StatsPipeline/0-ETL/Step2Parquet2Rdata.R'),
    to = file.path(path, 'Step2Parquet2Rdata.R'),
    overwrite = TRUE
  )
  system('chmod 775 jobs_step2_Parquet2Rdata.bch', wait = TRUE)
  system('sbatch --job-name=impc_stats_pipeline_job --time=01:00:00 --mem=1G -o ../compressed_logs/step2_job_id.txt --wrap="bash ./jobs_step2_Parquet2Rdata.bch"', wait = TRUE)
  waitTillCommandFinish(
    WaitIfTheOutputContains = waitUntillSee,
    ignoreline = ignoreThisLineInWaitingCheck
  )
  file.remove(file.path(path, 'Step2Parquet2Rdata.R'))
  if (filesContain(path = path0,
                   extension = '\\.log',
                   containWhat = 'Exit'))
    stop('An error occured in step 2. Parquet2Rdata conversion')
  
  system(command = "find ../ -type f -name '*.log' -exec zip -m ../compressed_logs/step2_logs.zip {} +", wait = TRUE)
  system(command = "find ../ -type f -name '*.err' -exec zip -m ../compressed_logs/step2_logs.zip {} +", wait = TRUE)

  ###############################################
  message0('Step 3. Merging psudo Rdata files into single file for each procedure - jobs creator')
  source(file.path(local(),
                   'StatsPipeline/0-ETL/Step3MergeRdataFilesJobs.R'))
  f(file.path(path, 'ProcedureScatterRdata'))
  rm0('f')

  ###############################################
  message0('Step 4. Merging psudo Rdata files into single files per procedure')
  file.copy(
    from = file.path(local(),
                     'StatsPipeline/0-ETL/Step4MergingRdataFiles.R'),
    to = file.path(path, 'Step4MergingRdataFiles.R'),
    overwrite = TRUE
  )
  system('chmod 775 jobs_step4_MergeRdatas.bch', wait = TRUE)
  system('sbatch --job-name=impc_stats_pipeline_job --time=01:00:00 --mem=1G -o ../compressed_logs/step4_job_id.txt --wrap="bash ./jobs_step4_MergeRdatas.bch"', wait = TRUE)
  waitTillCommandFinish(
    WaitIfTheOutputContains = waitUntillSee,
    ignoreline = ignoreThisLineInWaitingCheck
  )
  file.remove(file.path(path, 'Step4MergingRdataFiles.R'))
  if (filesContain(path = path,
                   extension = '\\.log',
                   containWhat = 'Exit'))
    stop('An error occured in step 4. Merging Rdata files into one single Rdata file per procedure')
  
  system(command = "find . -type f -name '*.log' -exec zip -m ../compressed_logs/step4_logs.zip {} +", wait = TRUE)
  system(command = "find . -type f -name '*.err' -exec zip -m ../compressed_logs/step4_logs.zip {} +", wait = TRUE)

  ###############################################
  ## Compress logs
  message0('Phase I. Compressing the log files and house cleaning ...')
  system(command = 'zip -rm ../compressed_logs/phase1_jobs.zip *.bch', wait = TRUE)
  system(command = 'rm -rf ProcedureScatterRdata', wait = TRUE)
  ###########  END of Phase I ###################

  ##### Phase II. Reprocessing the data
  message0('Starting Phase II, packaging the big data into small packages ...')
  jobCreator(path = file.path(path, 'Rdata/'))
  file.copy(
    from = file.path(local(),
                     'StatsPipeline/jobs/InputDataGenerator.R'),
    to = file.path(path, 'InputDataGenerator.R'),
    overwrite = TRUE
  )
  system('chmod 775 DataGenerationJobList.bch', wait = TRUE)
  system('sbatch --job-name=impc_stats_pipeline_job --time=01:00:00 --mem=1G -o ../compressed_logs/phase2_job_id.txt --wrap="bash ./DataGenerationJobList.bch"', wait = TRUE)
  waitTillCommandFinish(
    WaitIfTheOutputContains = waitUntillSee,
    ignoreline = ignoreThisLineInWaitingCheck
  )
  file.remove(file.path(path, 'InputDataGenerator.R'))
  if (filesContain(
    path = file.path(path, 'DataGeneratingLog'),
    extension = '\\.log',
    containWhat = 'Exit'
  ))
    stop('An error occured in Phase II step 1. Packaging the big data into small packages')

  ## Compress logs
  message0('End of packaging data. ')
  message0('Phase II. Compressing the log files and house cleaning ... ')
  system(command = 'mv *.bch  DataGeneratingLog/', wait = TRUE)
  system(command = 'zip -rm phase2_logs.zip DataGeneratingLog/', wait = TRUE)
  system(command = 'mv phase2_logs.zip ../compressed_logs/', wait = TRUE)

  ## Add all single jobs into one single job
  message0('Appending all procedure based jobs into one single file ...')
  if (!dir.exists('jobs'))
    system(command = 'mkdir jobs', wait = TRUE)
  if (file.exists('jobs/AllJobs.bch'))
    system(command = 'rm jobs/AllJobs.bch', wait = TRUE)
  system(command = 'find ./*/*_RawData/*.bch -type f | xargs  cat >> jobs/AllJobs.bch', wait = TRUE)

  message0('Phase III. Initialising the statistical analysis ...')
  path = file.path(path, 'jobs')
  setwd(path)


  ## Update procedure/parameters from the IMPReSS
  message0('Updating the dynamic contents from the IMPReSS ...')
  updateImpress(
    updateImpressFileInThePackage = TRUE,
    updateOptionalParametersList = TRUE,
    updateTheSkipList = TRUE,
    saveRdata = FALSE
  )

  message0('Running the IMPC statistical pipeline by submitting jobs ...')
  # copy stats pipeline driver script
  if (windowingPipeline) {
    file.copy(
      from = file.path(local(),
                       'StatsPipeline/jobs/function_windowed.R'),
      to = file.path(path, 'function.R'),
      overwrite = TRUE
    )
  } else{
    file.copy(
      from = file.path(local(),
                       'StatsPipeline/jobs/function.R'),
      to = file.path(path, 'function.R'),
      overwrite = TRUE
    )
  }

  # Specify the DR version in the pipeline
  ReplaceWordInFile(file.path(path, 'function.R'),
                    'DRversionNotSpecified',
                    DRversion)

  system('chmod 775 AllJobs.bch', wait = TRUE)
  submit_limit_jobs(bch_file="AllJobs.bch", job_id_logfile="../../compressed_logs/phase3_job_id.txt")
  waitTillCommandFinish(
    WaitIfTheOutputContains = waitUntillSee,
    ignoreline = ignoreThisLineInWaitingCheck
  )

  message0('Postprocessing the IMPC statistical analysis results ...')
  setwd(SP.results)

  message0('Creating backups from the main R packages ...')
  if (dir.exists('RPackage_backup')) {
    system(command = 'rm -rf RPackage_backup', wait = TRUE)
  }
  packageBackup('DRrequiredAgeing', storepath = file.path(getwd(), 'RPackage_backup'))
  packageBackup('OpenStats', storepath = file.path(getwd(), 'RPackage_backup'))

  if (dir.exists('logs')) {
    system(command = 'rm -rf logs', wait = TRUE)
  }

  system(command = 'mkdir logs', wait = TRUE)


  system(command = "find . -type f -name '*.ClusterOut' -exec zip -m ../compressed_logs/phase3_logs.zip {} +", wait = TRUE)
  system(command = "find . -type f -name '*.ClusterErr' -exec zip -m ../compressed_logs/phase3_errs.zip {} +", wait = TRUE)

  message0(
    'This is the last step. If you see no file in the list below, the SP is successfully completed.'
  )
  setwd(file.path(SP.results, 'logs'))
  
  message0('SP finished in ', round(difftime(Sys.time(), startTime, units = "min"), 2))
}

editfile = function(file,
                    searchwhat = '',
                    replaceby = '') {
  tx  <- readLines(file)
  tx2  <-
    gsub(
      pattern = searchwhat,
      replace = replaceby,
      x = tx,
      fixed = TRUE
    )
  writeLines(tx2, con = file)
}


install.packages.auto <- function(x) {
  x <- as.character(substitute(x))
  if(isTRUE(x %in% .packages(all.available=TRUE))) {
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else {
    #update.packages(ask= FALSE) #update installed packages.
    eval(parse(text = sprintf("install.packages(\"%s\", dependencies = TRUE)", x)))
  }
  if(isTRUE(x %in% .packages(all.available=TRUE))) {
    eval(parse(text = sprintf("require(\"%s\")", x)))
  } else {
    source("http://bioconductor.org/biocLite.R")
    #biocLite(character(), ask=FALSE) #update installed packages.
    eval(parse(text = sprintf("biocLite(\"%s\")", x)))
    eval(parse(text = sprintf("require(\"%s\")", x)))
  }
}

##########################################
############## For the annotation pipeline
##########################################
GetMethodStPa <- function(x) {
  if (is.null(x)) {
    return("UNK")
  }
  if (grepl(pattern = "Fisher", x)) {
    r <- "FE"
  } else if (grepl(pattern = "(Mixed Model)|(Linear Model)", x)) {
    r <- "MM"
  } else if (grepl(pattern = "Reference Range", x)) {
    r <- "RR"
  } else {
    r <- "UNK"
  }
  return(r)
}

DirectionTagMM = function(x,
                          pvalue = NULL,
                          threshold = .0001) {
  tag = if (is.null(pvalue)) {
    'NoPValueAvailable'
  } else if (pvalue > threshold) {
    paste0('Test not significant at the level of ',
           threshold,
           ' (pvalue = ',
           pvalue,
           ')')
  }  else if (is.null(x) || length(x) < 1) {
    'NoEffectCalculated'
  } else if (x > 0 && pvalue < threshold) {
    'INCREASED'
  } else if (x < 0 && pvalue < threshold) {
    'DECREASED'
  } else if (x == 0 && pvalue < threshold) {
    'STEADY'
  } else{
    'Something is gone wrong. Please report this case to the IMPC team.'
  }
  return(tag)
}

DirectionTagFE = function(x,
                          threshold = .0001,
                          group = c('ABNORMAL', 'INCREASED', 'DECREASED')) {
  default = c('ABNORMAL', 'INCREASED', 'DECREASED')
  if (length(group) < 1) {
    message(
      'No group for the so called ABNORMAL category provided, swith to the backup plan:\n\t ',
      default
    )
    group = default
  }
  tag = if (is.null(x) || length(x) < 1) {
    'NoEffectCalculated'
  } else if (x < threshold) {
    group
  } else {
    paste0('NORMAL - Test not significant at the level of ',
           threshold,
           ' (pvalue = ',
           x,
           ')')
  }
  return(tag)
}

listM = function(x, What2Attach = list('OVERALL' = list('MPTERM' = NA))) {
  if (is.null(x))
    return(NULL)
  al = as.list(x)
  names(al) = al
  al[[names(al)]] = What2Attach
  return(al)
}


numbers_only <- function(x) {
  r = !grepl("\\D", x)
  return(r)
}

multiGrepl = function (x = NULL, pattern = NULL, ...)
{
  if (is.null(x)) {
    return(NULL)
  }
  if (is.null(pattern)) {
    return(TRUE)
  }
  r <- rep(TRUE, length(x))
  for (p in pattern) {
    r <- r & grepl(pattern = p, x = x, ...)
  }
  return(r)
}

returnWhatBasedOnThreshold = function(x = NULL,
                                      threshold = .0001,
                                      Return = 'ABNORMAL',
                                      ReturnNull = 'IgnoreThisCaseAtALL') {
  if (is.null(x)      ||
      is.null(Return) ||
      is.null(threshold) ||
      length(x) < 1      ||
      length(Return) < 1 ||
      length(threshold) < 1) {
    return(ReturnNull)
  } else if (is.numeric(x) && x < threshold) {
    return(Return)
  } else{
    return(ReturnNull)
  }
}

GenotypeTag = function(obj,
                       parameter_stable_id = NULL,
                       threshold = 10 ^ -4,
                       expDetailsForErrorOnly = NULL,
                       rrlevel = 10 ^ -4) {
  if (is.null(obj))
    return('NotAnalysed')
  method = GetMethodStPa(x = obj$`Applied method`)
  message('\t The analysis method = ', method)
  message('\t The decision threshold = ',
          threshold,
          ' and for the RR method it is ',
          rrlevel)
  if (method %in% 'MM') {
    tag = list(
      UNSPECIFIED = listM(
        DirectionTagMM(
          x      = obj$`Genotype estimate`$Value,
          pvalue = obj$`Genotype p-value`,
          threshold = threshold
        )
      ),
      UNSPECIFIED = listM(
        returnWhatBasedOnThreshold(
          x = obj$`Genotype p-value`,
          threshold = threshold,
          Return = 'ABNORMAL'
        )
      ),
      UNSPECIFIED = listM(
        returnWhatBasedOnThreshold(
          x = obj$`Genotype p-value`,
          threshold = threshold,
          Return = 'INFERRED'
        )
      ),
      FEMALE        = listM(
        DirectionTagMM(
          x      = obj$`Sex FvKO estimate`$Value,
          pvalue = obj$`Sex FvKO p-value`,
          threshold = threshold
        )
      ),
      FEMALE = listM(
        returnWhatBasedOnThreshold(
          x = obj$`Sex FvKO p-value`,
          threshold = threshold,
          Return = 'ABNORMAL'
        )
      ),
      FEMALE = listM(
        returnWhatBasedOnThreshold(
          x = obj$`Sex FvKO p-value`,
          threshold = threshold,
          Return = 'INFERRED'
        )
      ),
      MALE      = listM(
        DirectionTagMM(
          x      = obj$`Sex MvKO estimate`$Value,
          pvalue = obj$`Sex MvKO p-value`,
          threshold = threshold
        )
      ),
      MALE = listM(
        returnWhatBasedOnThreshold(
          x = obj$`Sex MvKO p-value`,
          threshold = threshold,
          Return = 'INFERRED'
        )
      ),
      MALE = listM(
        returnWhatBasedOnThreshold(
          x = obj$`Sex MvKO p-value`,
          threshold = threshold,
          Return = 'INFERRED'
        )
      )
    )
  } else if (method %in% 'FE') {
    ###########################################################################
    fmodels = obj$`Additional information`$Analysis$`Further models`$category
    if (is.null(fmodels))
      return(NULL)

    if (length(names(fmodels)) == 1 &&
        names(fmodels) == 'Complete table') {
      fmodels$Genotype$`Complete table` = fmodels$`Complete table`
      fmodels$`Complete table` = NULL
    }
    #fmodels$Genotype$`Complete table`
    AllCombinations = lapply(fmodels, function(x) {
      if (length(fmodels) > 0) {
        lapply(x, function(y) {
          DirectionTagFE(x = y$p.value, threshold = threshold)
        })
      } else{
        DirectionTagFE(x = x$p.value, threshold = threshold)
      }
    })
    if (is.null(AllCombinations))
      return(NULL)
    #### Make the list as sequence of names attached with dot (.)
    AllCombinations1 = unlist(AllCombinations)
    if (!is.null(AllCombinations1)){
      AllCombinations1 = AllCombinations1[grepl('Complete table',names(AllCombinations1))]
    }
    if (is.null(AllCombinations1))
      return(NULL)
    #### Keep only genotype analysis
    AllCombinations2 = AllCombinations1[grepl(pattern = 'Genotype', x =
                                                names(AllCombinations1))]
    #### The abnormal case is made from the data
    for (i in seq_along(AllCombinations2)) {
      ############# step 1
      names(AllCombinations2)[i] = gsub(
        pattern = 'Complete table',
        replacement = paste0(
          AllCombinations2[i],
          ifelse(
            AllCombinations2[i] == 'ABNORMAL',
            '.OVERALL.MPTERM',
            '.MPTERM'
          )
        ),
        x = names(AllCombinations2)[i]
      )
      ############# step 2
      gind = !multiGrepl(pattern = c('Genotype', 'OVERALL'),
                         x = names(AllCombinations2)[i])
      if (gind) {
        strs = unlist(strsplit(
          names(AllCombinations2)[i],
          split = '.',
          fixed = TRUE
        ))
        u = unique(c(strs[1],
                     AllCombinations2[i],
                     strs[-1]))
        names(AllCombinations2)[i] = paste(u,
                                           collapse = '.',
                                           sep = '.')
      }
      ############# step 3
      names(AllCombinations2)[i] = gsub(
        pattern = 'Genotype.',
        replacement = '',
        x = names(AllCombinations2)[i]
      )
      ############# step 4
      names(AllCombinations2)[i] = toupper(names(AllCombinations2)[i])
    }
    ############# step 5
    AllCombinations2 = AllCombinations2[grepl(
      pattern = '(MALE)|(FEMALE)|(OVERALL)',
      x = names(AllCombinations2),
      ignore.case = TRUE
    )]
    ############# step 6
    if (length(AllCombinations2) < 1)
      return(NULL)
    nam = names(AllCombinations2)
    for (i in 0:1000) {
      sb = substr(nam, start = nchar(nam) - i, nchar(nam))
      for (j in 1:length(sb)) {
        if (numbers_only(sb[j]) && i < nchar(sb[j]))
          nam[j] = substr(nam[j], start = 0, stop = nchar(nam[j]) - i -
                            1)
      }
    }
    ############# step 7
    gp7 = grepl(pattern = '.MPTERM',
                x = nam,
                fixed = TRUE)
    nam[!gp7] = paste0(nam[!gp7], '.MPTERM')

    ############# step 8
    if (!is.null(parameter_stable_id)) {
      controlCat = read.delim(file.path(local(),'annotation','CategoryRemapping.tsv'), sep = '\t')
      conCat = controlCat[controlCat$IMPRESS.ID %in% parameter_stable_id,]
      control = conCat[conCat$CATEGORY.1 == 0, 2]
      if (length(control) > 0) {
        for (c in control) {
          nam =   gsub(
            pattern = paste0('.', toupper(c), '.'),
            replacement = '.',
            x = nam,
            fixed = TRUE
          )
        }
        nam = gsub(
          pattern = '..',
          replacement = '.',
          x = nam,
          fixed = TRUE
        )
      }
    }
    ############# Finally!
    names(AllCombinations2) = nam
    tag = AllCombinations2

  } else if (method %in% 'RR') {
    ###########################################################################
    fmodels = obj$`Additional information`$Analysis$`Further models`
    if (is.null(fmodels))
      return(NULL)

    AllCombinations = lapply(fmodels, function(x) {
      lapply(x, function(y) {
        DirectionTagFE(x = y$Result$p.value, threshold = rrlevel)
      })
    })

    if (is.null(AllCombinations))
      return(NULL)
    #### Make the list as sequence of names attached with dot (.)
    AllCombinations1 = unlist(AllCombinations)
    ####
    AllCombinations1 = AllCombinations1[!(grepl(
      pattern = 'Low.',
      x = names(AllCombinations1),
      fixed = TRUE
    ) & AllCombinations1 %in% 'INCREASED')]
    AllCombinations1 = AllCombinations1[!(grepl(
      pattern = 'High.',
      x = names(AllCombinations1),
      fixed = TRUE
    ) & AllCombinations1 %in% 'DECREASED')]

    #### Keep only genotype analysis
    AllCombinations2 = AllCombinations1[grepl(pattern = '(Genotype.Genotype)|(Genotype.Male.Genotype)|(Genotype.Female.Genotype)', x =
                                                names(AllCombinations1))]
    #### The abnormal case is made from the data
    for (i in seq_along(AllCombinations2)) {
      names(AllCombinations2)[i]  = gsub(
        pattern = 'data_point.Genotype.Male.Genotype',
        paste0('MALE.', AllCombinations2[i], '.OVERALL.MPTERM'),
        x = names(AllCombinations2)[i]
      )
      names(AllCombinations2)[i]  = gsub(
        pattern = 'data_point.Genotype.Female.Genotype',
        paste0('FEMALE.', AllCombinations2[i], '.OVERALL.MPTERM'),
        x = names(AllCombinations2)[i]
      )
      names(AllCombinations2)[i]  = gsub(
        pattern = 'data_point.Genotype.Genotype',
        paste0('UNSPECIFIED.', AllCombinations2[i], '.OVERALL.MPTERM'),
        x = names(AllCombinations2)[i]
      )
      names(AllCombinations2)[i] = gsub(
        pattern = '..',
        replacement = '.',
        x = names(AllCombinations2)[i],
        fixed = TRUE
      )
      ############# step 4
      names(AllCombinations2)[i] = toupper(names(AllCombinations2)[i])
    }
    AllCombinations2 = AllCombinations2[!grepl('(_FEMALE)|(_MALE)',names(AllCombinations2),fixed = FALSE)]
    ############# step 6
    if (length(AllCombinations2) < 1)
      return(NULL)
    nam = names(AllCombinations2)
    for (i in 0:1000) {
      sb = substr(nam, start = nchar(nam) - i, nchar(nam))
      for (j in 1:length(sb)) {
        if (numbers_only(sb[j]) && i < nchar(sb[j]))
          nam[j] = substr(nam[j], start = 0, stop = nchar(nam[j]) - i -
                            1)
      }
    }

    ############# Finally!
    names(AllCombinations2) = nam #gsub(pattern = 'LOW.|HIGH.',replacement = '',x = nam)
    tag = AllCombinations2
  } else{
    tag = NULL
    write(
      x = paste(
        head(expDetailsForErrorOnly, 18),
        sep = '\t',
        collapse = '\t'
      ),
      file = 'ErrorneousCases.tsv.err',
      ncolumns = 5000
    )
  }
  return(tag)
}


unScrewProcedure = function(x) {
  r = unlist(strsplit(x = x, split = '~'))
  return(r)
}


merge.two = function(l1, l2) {
  if (is.null(l1) && is.null(l2))
    return(NULL)
  if (is.null(l1) && !is.null(l2))
    return(l1)
  if (!is.null(l1) && is.null(l2))
    return(l2)
  #l = c(l1, l2[!names(l2) %in% names(l1)])
  l = c(l1, l2)
  return(l)
}
###########################

removeAbnormalIfIncDecDetected = function(x = NULL,
                                          method = 'MM',
                                          active = 1) {
  if (!active)
    return(x)
  if (length(x) > 0 &&
      !method %in% 'RR' &&
      any(grepl(
        x = names(x),
        pattern = ('\\.INCREASED\\.|\\.DECREASED\\.')
      ))) {
    message('\t INCREASED/DECREASED detected. ABNORMAL will be ignored')
    x = x[!grepl(x = names(x),
                 pattern = ('.ABNORMAL.'),
                 fixed = TRUE)]
  }
  return(x)
}

MatchTheRestHalfWithTheFirstOne = function(x) {
  if (length(x) < 1)
    return(x)
  #x = x[!is.na(x)]
  mid = length(x) / 2
  x[(mid + 1):(2 * mid)] = x[1:mid]
  return(x)
}

DecIncDetector = function(x) {
  if (length(x) < 1)
    return(NA)
  r = c()
  if (any(grepl(pattern = 'DECREASED', x = names(x))))
    r = c(r, 'DECREASED')
  if (any(grepl(pattern = 'INCREASED', x = names(x))))
    r = c(r, 'INCREASED')
  if (any(grepl(pattern = 'ABNORMAL', x = names(x))))
    r = c(r, 'ABNORMAL')
  if (any(grepl(pattern = 'INFERRED', x = names(x))))
    r = c(r, 'INFERRED')

  if (length(r) < 1)
    r = NA
  if (length(r) > 1) {
    if (any(grepl(pattern = 'DECREASED', x = r)) && !any(grepl(pattern = 'INCREASED', x = r)))
      r = 'DECREASED'
    if (any(grepl(pattern = 'INCREASED', x = r)) && !any(grepl(pattern = 'DECREASED', x = r)))
      r = 'INCREASED'
    if (any(grepl(pattern = 'ABNORMAL', x = r)))
      r = 'ABNORMAL'
    if (any(grepl(pattern = 'INFERRED', x = r)))
      r = 'INFERRED'
  }
  return(r)
}

DecIncDetectorRR = function(x) {
  if (length(x) < 1)
    return(NA)
  r = c()
  ###################################### MALE
  ###################################### LOW
  r1 = multiGrepl(pattern = c('(DECREASED)|(INCREASED)', 'LOW', 'MALE'),
                  x = names(x))
  if (sum(r1) > 1) {
    rlowM = c('ABNORMAL')
  } else if (sum(r1) == 1) {
    if (grepl(pattern = 'INCREASED', x = names(x[r1])))
      rlowM = c('INCREASED')
    if (grepl(pattern = 'DECREASED', x = names(x[r1])))
      rlowM = c('DECREASED')
  } else{
    rlowM = NA
  }
  ###################################### High
  r2 = multiGrepl(pattern = c('(DECREASED)|(INCREASED)', 'HIGH', 'MALE'),
                  x = names(x))
  if (sum(r2) > 1) {
    rhighM = c('ABNORMAL')
  } else if (sum(r2) == 1) {
    if (grepl(pattern = 'INCREASED', x = names(x[r2])))
      rhighM = c('INCREASED')
    if (grepl(pattern = 'DECREASED', x = names(x[r2])))
      rhighM = c('DECREASED')
  } else{
    rhighM = NA
  }
  ######################################
  r3 = na.omit(c(rlowM, rhighM))
  if (length(r3) > 0) {
    if (length(r3) == 1) {
      r = c(r, na.omit(r3))
    } else{
      r = c(r, 'ABNORMAL')
    }
  }
  ###################################### FEMALE
  ###################################### LOW
  r4 = multiGrepl(pattern = c('(DECREASED)|(INCREASED)', 'LOW', 'FEMALE'),
                  x = names(x))
  if (sum(r4) > 1) {
    rlowF = c('ABNORMAL')
  } else if (sum(r4) == 1) {
    if (grepl(pattern = 'INCREASED', x = names(x[r4])))
      rlowF = c('INCREASED')
    if (grepl(pattern = 'DECREASED', x = names(x[r4])))
      rlowF = c('DECREASED')
  } else{
    rlowF = NA
  }
  ######################################
  r5 = multiGrepl(
    pattern = c('(DECREASED)|(INCREASED)', 'HIGH', 'FEMALE'),
    x = names(x)
  )
  if (sum(r5) > 1) {
    rhighF = c('ABNORMAL')
  } else if (sum(r5) == 1) {
    if (grepl(pattern = 'INCREASED', x = names(x[r5])))
      rhighF = c('INCREASED')
    if (grepl(pattern = 'DECREASED', x = names(x[r5])))
      rhighF = c('DECREASED')
  } else{
    rhighF = NA
  }
  ######################################
  r6 = na.omit(c(rlowF, rhighF))
  if (length(r6) > 0) {
    if (length(r6) == 1) {
      r = c(r, na.omit(r6))
    } else{
      r = c(r, 'ABNORMAL')
    }
  }
  ######################################
  ###################################### ABNORMAL
  r7 = multiGrepl(
    pattern = c('(DECREASED)|(INCREASED)', 'LOW', 'ABNORMAL'),
    x = names(x)
  )
  if (sum(r7) > 1) {
    rlowA = c('ABNORMAL')
  } else if (sum(r7) == 1) {
    if (grepl(pattern = 'INCREASED', x = names(x[r7])))
      rlowA = c('INCREASED')
    if (grepl(pattern = 'DECREASED', x = names(x[r7])))
      rlowA = c('DECREASED')
  } else{
    rlowA = NA
  }
  ######################################
  r8 = multiGrepl(
    pattern = c('(DECREASED)|(INCREASED)', 'HIGH', 'ABNORMAL'),
    x = names(x)
  )
  if (sum(r8) > 1) {
    rhighA = c('ABNORMAL')
  } else if (sum(r8) == 1) {
    if (grepl(pattern = 'INCREASED', x = names(x[r8])))
      rhighA = c('INCREASED')
    if (grepl(pattern = 'DECREASED', x = names(x[r8])))
      rhighA = c('DECREASED')
  } else{
    rhighA = NA
  }
  ######################################
  r9 = na.omit(c(rlowA, rhighA))
  if (length(r9) > 0) {
    if (length(r9) == 1) {
      r = c(r, na.omit(r9))
    } else{
      r = c(r, 'ABNORMAL')
    }
  }
  ######################################
  abpattern = inferpattern = FALSE

  if (length(x) > 0) {
    abpattern = grepl(pattern = 'ABNORMAL', x = names(x))
    inferpattern = grepl(pattern = 'INFERRED', x = names(x))
  }

  if (length(r) < 1 && !(any(abpattern) || any(inferpattern)))
    r = NA
  ######################################
  if (length(r) < 1 && (any(abpattern) || any(inferpattern)))
    if (any(grepl(pattern = 'ABNORMAL', x = names(x))))
      r = 'ABNORMAL'
  if (any(grepl(pattern = 'INFERRED', x = names(x))))
    r = 'INFERRED'
  ######################################
  if (length(r) > 1) {
    if (any(grepl(pattern = 'ABNORMAL', x = r)))
      r = 'ABNORMAL'
    if (any(grepl(pattern = 'INFERRED', x = r)))
      r = 'INFERRED'
  }
  ######################################
  return(r)
}

detectLevel = function(x, level) {
  if (is.null(x) || length(x) < 1)
    return(FALSE)
  r = grepl(pattern = level, x = x)
  return(any(r))
}

NullOrvalueReturn = function(x, list) {
  if (length(x) < 1)
    return(NULL)
  else
    return(list)
}

fA = function(x, pasteterms = TRUE) {
  if (length(x) > 0)
    x = x[!duplicated(x)]
  if (length(x) > 1) {
    x2 = x[!grepl(pattern = c('(MALE)|(FEMALE)'), names(x)) &
             grepl(pattern = c('ABNORMAL'), names(x))]
    x3 = x[grepl(pattern = c('ABNORMAL'), names(x))]
    if (length(x2) > 1 && length(x3) > 0 && pasteterms)
      x = paste(x3, sep = '~', collapse = '~')
    else if (length(x2) > 1 && length(x3) < 1 && pasteterms)
      x = paste(x2, sep = '~', collapse = '~')
    else if (pasteterms)
      x = paste(x, sep = '~', collapse = '~')
    else
      x = x
  }
  return(x)
}

fM = function(x, pasteterms = TRUE) {
  if (length(x) > 0)
    x = x[!duplicated(x)]
  if (length(x) > 1) {
    x2 = x[!grepl(pattern = c('ABNORMAL'), names(x)) &
             grepl(pattern = c('MALE'), names(x))]
    x3 = x[grepl(pattern = c('ABNORMAL'), names(x))]
    if (length(x2) > 1 && length(x3) > 0 && pasteterms)
      x = paste(x3, sep = '~', collapse = '~')
    else if (length(x2) > 1 && length(x3) < 1 && pasteterms)
      x = paste(x2, sep = '~', collapse = '~')
    else if (pasteterms)
      x = paste(x, sep = '~', collapse = '~')
    else
      x = x
  }
  return(x)
}

fF = function(x, pasteterms = TRUE) {
  if (length(x) > 0)
    x = x[!duplicated(x)]
  if (length(x) > 1) {
    x2 = x[!grepl(pattern = c('ABNORMAL'), names(x)) &
             grepl(pattern = c('FEMALE'), names(x))]
    x3 = x[grepl(pattern = c('ABNORMAL'), names(x))]
    if (length(x2) > 1 && length(x3) > 0 && pasteterms)
      x = paste(x3, sep = '~', collapse = '~')
    else if (length(x2) > 1 && length(x3) < 1 && pasteterms)
      x = paste(x2, sep = '~', collapse = '~')
    else if (pasteterms)
      x = paste(x, sep = '~', collapse = '~')
    else
      x = x
  }
  return(x)
}

bselect = function(x) {
  # order important, see MoreThan2Length
  if(length(x)<1)
    return(x)
  nx = names(x)
  xx = x
  if (length(nx) > 1 && any(grepl(pattern = 'INFERRED', x = nx))) {
    xx = x[grepl(pattern = 'INFERRED', x = nx)][1]
    return(xx)
  }
  if (length(nx) > 1 && any(grepl(pattern = 'ABNORMAL', x = nx))) {
    xx = x[grepl(pattern = 'ABNORMAL', x = nx)][1]
    return(xx)
  }

  if (length(nx) > 1 && any(grepl(pattern = 'OVERAL', x = nx))) {
    xx = x[grepl(pattern = 'OVERAL', x = nx)][1]
    return(xx)
  }

  if (length(xx) > 0)
    xx = na.omit(xx)
  return(xx)
}



MoreThan2Length = function(xx,
                           index,
                           MPTERMS = NULL,
                           general = TRUE) {
  # order important, see bselect
  if (length(xx) > 1 &&
      sum(index) > 1 &&
      !is.null(MPTERMS) &&
      length(MPTERMS) > 0) {
    if (general && any(grepl(pattern = 'INCREASED', xx)) &&
        any(grepl(pattern = 'INCREASED', names(MPTERMS)))) {
      return(grepl(pattern = 'INCREASED', xx))

    } else if (general && any(grepl(pattern = 'DECREASED', xx)) &&
               any(grepl(pattern = 'DECREASED', names(MPTERMS)))) {
      return(grepl(pattern = 'DECREASED', xx))

    } else if (any(grepl(pattern = 'INFERRED', xx)) &&
               any(grepl(pattern = 'INFERRED', names(MPTERMS)))) {
      return(grepl(pattern = 'INFERRED', xx))

    } else if (any(grepl(pattern = 'ABNORMAL', xx)) &&
               any(grepl(pattern = 'ABNORMAL', names(MPTERMS)))) {
      return(grepl(pattern = 'ABNORMAL', xx))

    } else if (any(grepl(pattern = 'OVERAL', xx)) &&
               any(grepl(pattern = 'OVERAL', names(MPTERMS)))) {
      return(grepl(pattern = 'OVERAL', xx))
    }

  }
  return(index)
}

MaleFemaleAbnormalCategories = function(x, method = 'AA', MPTERMS = NULL,json=NULL) {
  fgrep = grepl(pattern = 'FEMALE', names(x), fixed = TRUE)
  mgrep = grepl(pattern = 'MALE', names(x), fixed = TRUE) & !fgrep
  agrep = grepl(pattern = '(ABNORMAL)|(INFERRED)|(OVERAL)', names(x)) &
    !fgrep & !mgrep
  sexlevels = json$Result$`Vector output`$`Normal result`$`Classification tag`$`Active Sex levels`

  if (method %in% 'MM') {
    fgrep = MoreThan2Length(names(x), fgrep, MPTERMS)
    mgrep = MoreThan2Length(names(x), mgrep, MPTERMS)
    agrep = MoreThan2Length(names(x), agrep, MPTERMS)
  } else{
    # bug reported 30/10/2020 - (improvement)
    fgrep = MoreThan2Length(names(x), fgrep, MPTERMS, general = FALSE)
    mgrep = MoreThan2Length(names(x), mgrep, MPTERMS, general = FALSE)
    agrep = MoreThan2Length(names(x), agrep, MPTERMS, general = FALSE)
  }

  if (method %in% 'RR') {
    fevent = DecIncDetectorRR(x[fgrep])
    mevent = DecIncDetectorRR(x[mgrep])
    oevent = DecIncDetectorRR(x[agrep])
  } else{
    fevent = DecIncDetector(x[fgrep])
    mevent = DecIncDetector(x[mgrep])
    oevent = DecIncDetector(x[agrep])
  }
  # if males and females are in cross direction set abnormal otherwise select male/females
  fmevents = c(fevent, mevent)
  if (all(!is.na(fmevents)) && fmevents[1] == fmevents[2])
    oevent = mevent
  if (is.na(fmevents[1]) && !is.na(fmevents[2]))
    oevent = fmevents[2]
  if (!is.na(fmevents[1]) && is.na(fmevents[2]))
    oevent = fmevents[1]
  #############################
  # if male and/or female specific term found then ignore the combined sex mp terms
  if (length(x[fgrep]) > 0 || length(x[mgrep]) > 0)
    agrep = FALSE
  #############################
  if (method %in% 'RR') {
    MPTERM = list(
      NullOrvalueReturn(
        x[agrep],
        list(
          'term_id' = ifelse(length(fA(x[agrep], pasteterms = FALSE)) >
                               1, fA(x[agrep], pasteterms = FALSE)[1], fA(x[agrep])),
          event = oevent,
          sex = ifelse(
            length(unique(sexlevels)) > 1 ,
            "not_considered",
            unique(sexlevels)[1]
          ),
          otherPossibilities = ifelse(length(fA(x[agrep], pasteterms = FALSE)) >
                                        1, fA(x[agrep]), '')
        )
      ),
      NullOrvalueReturn(
        x[fgrep],
        list(
          'term_id' = ifelse(length(fF(x[fgrep], pasteterms = FALSE)) >
                               1, fF(x[fgrep], pasteterms = FALSE)[1], fF(x[fgrep])),
          event = fevent,
          sex = "female",
          otherPossibilities = ifelse(length(fF(x[fgrep], pasteterms = FALSE)) >
                                        1, fF(x[fgrep]), '')
        )
      ),
      NullOrvalueReturn(
        x[mgrep],
        list(
          'term_id' = ifelse(length(fM(x[mgrep], pasteterms = FALSE)) >
                               1, fM(x[mgrep], pasteterms = FALSE)[1], fM(x[mgrep])),
          event = mevent,
          sex = "male",
          otherPossibilities = ifelse(length(fM(x[mgrep], pasteterms = FALSE)) >
                                        1, fM(x[mgrep]), '')
        )
      )
    )
  } else{
    MPTERM = list(
      NullOrvalueReturn(x[agrep],
                        list(
                          'term_id' = bselect(x[agrep]),
                          event = oevent,
                          sex = ifelse(
                            length(unique(sexlevels)) > 1 ,
                            "not_considered",
                            unique(sexlevels)[1]
                          )
                        )),
      NullOrvalueReturn(x[fgrep], list(
        'term_id' = bselect(x[fgrep]),
        event = fevent,
        sex = "female"
      )),
      NullOrvalueReturn(x[mgrep], list(
        'term_id' = bselect(x[mgrep]),
        event = mevent,
        sex = "male"
      ))
    )
  }
  return(MPTERM)
}
##############################
annotationChooser = function(statpacket = NULL,
                             level = 10 ^ -4,
                             rrlevel = .005,
                             TermKey = 'MPTERM',
                             resultKey = 'Normal result',
                             mp_chooser_file = NULL) {
  requireNamespace('RPostgreSQL')
  requireNamespace('data.table')
  requireNamespace("data.table")
  requireNamespace("jsonlite")
  requireNamespace("rlist")
  requireNamespace("Tmisc")

  ulistTag3 = MPTERMS = NA
  message('Running the annotation pipeline')
  if (is.null(statpacket) ||
      length(statpacket) < 1 ||
      is.null(statpacket$V2) ||
      length(statpacket$V2)<1 ||
      !statpacket$V2 %in% 'Successful') {
    message('Not a successfull StatPackage!')
    return(invisible(list(
      MPTERM = ulistTag3, statpacket = statpacket
    )))
  }
  pipeline = statpacket$V15
  procedure = statpacket$V3
  parameter = statpacket$V6
  json      =  jsonlite::fromJSON(statpacket$V20)
  method   =   GetMethodStPa(json$Result$`Vector output`[[resultKey]]$`Applied method`)
  ##################################################################
  Gtag = GenotypeTag(
    obj = json$Result$`Vector output`[[resultKey]],
    parameter_stable_id = statpacket$V6,
    threshold = level,
    rrlevel = rrlevel
  )
  ##################################################################
  message('\t Reading the index file ...')
  load(mp_chooser_file)
  message(
    '\t~>',
    paste(
      pipeline,
      procedure,
      parameter,
      sep = '~>',
      collapse = '~>'
    ),
    '\n\t~>',
    json$Result$Details$`Gene page URL`
  )

  ################################
  b = a[[unScrewProcedure(pipeline)[1]]]
  c = b[[unScrewProcedure(procedure)[1]]]
  d = c[[unScrewProcedure(parameter)[1]]]
  if (is.null(d)) {
    message('No annotation available by IMPC. See https://www.mousephenotype.org/impress/')
    return(invisible(list(
      MPTERM = ulistTag3, statpacket = statpacket
    )))
  }
  ################################
  if (length(d) > 0 &&
      length(Gtag) > 0) {
    ################################
    ulistTag  = unlist(Gtag)
    ulistD    = unlist(d)
    names(ulistD)  = toupper(names(ulistD))
    ulistD         = removeAbnormalIfIncDecDetected(ulistD, method = method, active = FALSE)
    if (length(ulistD) < 1) {
      message('No annotation available by IMPC. See https://www.mousephenotype.org/impress/')
      return(invisible(list(
        MPTERM = ulistTag3, statpacket = statpacket
      )))
    }
    ################################
    ulistTag2 = ulistTag
    names(ulistTag2) = gsub(
      pattern = 'MALE|FEMALE',
      x = names(ulistTag2),
      replacement = 'UNSPECIFIED'
    )

    ulistTag3 = merge.two(ulistTag2, ulistTag)
    for (name in names(ulistTag3)) {
      # print(name)
      splN = unlist(strsplit(name, split = '.', fixed = TRUE))
      splN = splN[!splN %in% c('LOW', 'HIGH','DATA_POINT','GENOTYPE')]
      splnI = multiGrepl(pattern = splN,
                         x = names(ulistD),
                         fixed = TRUE)
      ulistTag3[names(ulistTag3) %in% name] = ifelse(is.null(ulistD[splnI]), 'CanNotFindMPTerm', head(ulistD[splnI], 1))
    }
    #print(ulistD)
    ################################
    if (length(ulistTag3) < 1)
      return(invisible(list(
        MPTERM = ulistTag3, statpacket = statpacket
      )))
    ##########################
    ulistTag3 =  MatchTheRestHalfWithTheFirstOne(ulistTag3)
    ulistTag3 =  ulistTag3[!duplicated(names(ulistTag3))]
    ##########################
    if (length(ulistTag3) < 1)
      return(invisible(list(
        MPTERM = ulistTag3, statpacket = statpacket
      )))
    ##########################
    ulistTag3 = ulistTag3[!is.na(ulistTag3)]
    ################################
    if (length(ulistTag3) > 0 &&
        method %in% 'RR') {
      lIncDec = multiGrepl(pattern = c('INCREASED|DECREASED'),
                           x = names(ulistTag3))
      lAbn    = multiGrepl(pattern = c('ABNORMAL'),
                           x = names(ulistTag3))
      if (sum(lIncDec) >= 2 && sum(lAbn) > 0) {
        ulistTag3 = ulistTag3[lAbn][1]
        names(ulistTag3) = gsub(pattern = '(LOW\\.)|(\\.HIGH)',
                                replacement = '',
                                names(ulistTag3))
      }
    }
    # Do not use unique!
    if (length(ulistTag3) > 0) {
      ulistTag3 =  ulistTag3[!duplicated(names(ulistTag3))]
    }
    ################################
    MPTERMS = MaleFemaleAbnormalCategories(x = ulistTag3,
                                           method = method,
                                           MPTERMS = ulistD,
                                           json = json)
    MPTERMS = rlist::list.clean(MPTERMS)
  }
  if (!is.null(MPTERMS)     &&
      length(ulistTag3) > 0 &&
      length(na.omit(ulistTag3)) > 0) {
    json$Result$Details[[TermKey]] = MPTERMS
    statpacket$V20 = FinalJson2ObjectCreator(FinalList = json)
  }else{
    message('No MP term found ...')
  }
  return(invisible(list(
    MPTERM = MPTERMS, statpacket = statpacket
  )))
}

#################################################################################
Write2Postg = function(df,
                       dbname = 'test',
                       host = "hh-yoda-05-01",
                       port = '5432',
                       user = 'impc',
                       password = 'impc',
                       tablename = paste0('db_',
                                         RemoveSpecialChars(
                                           format(Sys.time(), '%a%b%d %Y'),
                                           replaceBy = '_',
                                           what = ' '
                                         ))) {
  requireNamespace('RPostgreSQL')
  requireNamespace('data.table')
  requireNamespace('DBI')
  drv <- dbDriver("PostgreSQL")
  con <- dbConnect(
    drv,
    dbname = dbname,
    host = host,
    port = port,
    user = user,
    password = password
  )

  dbBegin(conn = con)
  r = dbWriteTable(
    conn = con,
    name = tablename,
    value = df,
    append = TRUE,
    row.names = FALSE
  )
  dbCommit(conn = con)
  dbDisconnect (conn = con)
  return(r)
}

randomIdGenerator = function(l = 10) {
  r=paste0(sample(
    x = 0:9,
    size = l,
    replace = TRUE
  ), collapse = '')
  return(r)
}

StratifiedMPTerms = function(object, name = 'MPTERM') {
  # overall, male, female
  r = c(NA, NA, NA)
  if (is.null(object))
    return(r)
  if (is.null(object) ||
      is.null(object[[name]]) ||
      length(object[[name]]) < 1 ||
      all(is.na(object[[name]])))
    return (r)

  for (i in seq_along(object[[name]])) {
    if (!is.null(object[[name]][[i]]$sex) && object[[name]][[i]]$sex == 'not_considered')
      r[1] = paste(object[[name]][[i]]$term_id,
                   collapse = '~',
                   sep = '~')
    if (!is.null(object[[name]][[i]]$sex) && object[[name]][[i]]$sex == 'male')
      r[2] = paste(object[[name]][[i]]$term_id,
                   collapse = '~',
                   sep = '~')
    if (!is.null(object[[name]][[i]]$sex) && object[[name]][[i]]$sex == 'female')
      r[3] = paste(object[[name]][[i]]$term_id,
                   collapse = '~',
                   sep = '~')
  }
  return(r)
}

annotationIndexCreator = function(path = getwd(),
                                  pattern = 'split_index',
                                  mem = "12G",
                                  time = "2-00",
                                  outputfile = 'annotation_jobs.bch') {
  lf = list.files(
    path = path,
    pattern = pattern,
    full.names = TRUE,
    recursive = FALSE,
    include.dirs = FALSE
  )

  if (!dir.exists(file.path(path, 'err')))
    dir.create(file.path(path, 'err'), recursive = TRUE)

  if (!dir.exists(file.path(path, 'out')))
    dir.create(file.path(path, 'out'), recursive = TRUE)

  if (!dir.exists(file.path(path, 'log')))
    dir.create(file.path(path, 'log'), recursive = TRUE)

  err = paste0(' -e ', dirname(lf), '/err/', basename(lf))
  out = paste0(' -o ', dirname(lf), '/out/', basename(lf))
  batch = paste0("sbatch --job-name=impc_stats_pipeline_job --mem=", mem,
                 " --time=", time,
                 ' ',
                 err,
                 out,
                 " --wrap='Rscript loader.R ",
                 '"',
                 basename(lf),
                 '"',
                 "'")
  write(batch, outputfile)
}
##########################################
############## End of annotation pipeline
##########################################
splitIndexFileIntoPiecesForPvalueExtraction = function(indexFilePath = NULL,
                                                       split = 1500,
                                                       mem = "9G",
                                                       time = "2-00",
                                                       outputfile = 'ExtractPValJobs.bch',
                                                       mp_chooser_file = NULL) {
  ldf = length(readLines(indexFilePath))
  ind = round(seq(1, ldf, length.out = split))

  if (file.exists(outputfile))
    unlink(outputfile)

  if (!dir.exists('error'))
    dir.create('error', recursive = TRUE)

  if (!dir.exists('output'))
    dir.create('output', recursive = TRUE)

  if (!dir.exists('resultF'))
    dir.create('resultF', recursive = TRUE)

  for (i in 2:(length(ind) + 1)) {
    write(
      paste(
        "sbatch --job-name=impc_stats_pipeline_job --mem=", mem,
        " --time=", time,
        ' -e "error/err',
        i,
        '" -o "output/out',
        i,
        '" ',
        "--wrap='Rscript ExtractPVals.R ",
        ind[i - 1],
        ' ',
        ind[i] - 1,
        ' "',
        indexFilePath,
        '" "',
        mp_chooser_file,
        '"',
        "'",
        sep = ''
      ),
      file = outputfile,
      ncolumns = 10 ^ 3,
      append = TRUE
    )
  }
}

changeRpackageDirectory = function(path = '~/DRs/R/packages') {
  #system('module load r-4.0.3-gcc-9.3.0-xiarbub',wait = TRUE)
  v = paste(
    R.version$major,
    gsub(
      pattern = '.',
      replacement = '_',
      R.version$minor,
      fixed = TRUE
    ),
    sep = '_',
    collapse = '_'
  )
  wdirc = file.path(path, v)
  if (!dir.exists(wdirc))
    dir.create(wdirc,
               showWarnings = FALSE,
               recursive = TRUE)
  .libPaths(new = wdirc)
  message(' => new package path set to: ', wdirc)
}


IMPC_HadoopLoad = function(SP.results = getwd(),
                           waitUntillSee = 'impc_stats_pipeline_job',
                           ignoreThisLineInWaitingCheck = 0,
                           ###
                           mp_chooser_file = NULL,
                           host =  "hh-hdp-master-02.ebi.ac.uk",
                           path = 'impc/statpackets',
                           prefix = 'DRXXX_',
                           port = '50070',
                           user =  Sys.info()['user'],
                           password = 'impc',
                           level = .0001,
                           rrlevel = .0001,
                           transfer = TRUE
) {

  configlist = list(
    mp_chooser_file = mp_chooser_file,
    host = host,
    path = path,
    prefix = prefix,
    port = port,
    user = user,
    password = password,
    level = level,
    rrlevel = rrlevel,
    transfer = transfer
  )

  DRrequiredAgeing:::message0('Step 1: Clean ups and creating the global index for results')
  DRrequiredAgeing:::message0('Indexing the results ...')
  setwd(file.path(SP.results))

  system('rm -f minijobs.bch', wait = TRUE)
  system('rm -f error.err', wait = TRUE)
  system('rm -f output.out', wait = TRUE)
  system('rm -f *.Ind', wait = TRUE)
  system('rm -rf AnnotationExtractorAndHadoopLoader/', wait = TRUE)

  DRrequiredAgeing:::minijobsCreator()
  system('chmod 775 minijobs.bch', wait = TRUE)
  DRrequiredAgeing:::submit_limit_jobs(bch_file="minijobs.bch", job_id_logfile="../../../compressed_logs/minijobs_job_id.txt")
  DRrequiredAgeing:::waitTillCommandFinish(
    WaitIfTheOutputContains = waitUntillSee,
    ignoreline = ignoreThisLineInWaitingCheck
  )
  
  system('mv minijobs.bch ../../../compressed_logs', wait = TRUE)
  system(command = "find . -type f -name '*_output.log' -exec zip -m ../../../compressed_logs/minijobs_logs.zip {} +", wait = TRUE)
  system(command = "find . -type f -name '*_error.err' -exec zip -m ../../../compressed_logs/minijobs_logs.zip {} +", wait = TRUE)

  DRrequiredAgeing:::message0('Moving single indeces into a separate directory called AnnotationExtractorAndHadoopLoader ...')
  system('rm -rf AnnotationExtractorAndHadoopLoader/', wait = TRUE)
  system('mkdir AnnotationExtractorAndHadoopLoader', wait = TRUE)
  system('chmod 775 AnnotationExtractorAndHadoopLoader/', wait = TRUE)
  system('mv *.Ind AnnotationExtractorAndHadoopLoader/', wait = TRUE)

  setwd(file.path(SP.results, 'AnnotationExtractorAndHadoopLoader'))
  DRrequiredAgeing:::message0('Concating single index files to create a global index for the results ...')
  system('cat *.Ind >> AllResultsIndeces.txt', wait = TRUE)

  DRrequiredAgeing:::message0('Zipping the single indeces ...')
  system('zip -rm allsingleindeces.zip *.Ind', wait = TRUE)

  system('rm -f split_index_*', wait = TRUE)
  system('rm -f splits.zip', wait = TRUE)
  system('split -50 AllResultsIndeces.txt split_index_', wait = TRUE)


  system('rm -f annotation_jobs.bch', wait = TRUE)
  system('rm -f config.Rdata', wait = TRUE)
  system('rm -f configHadoop.Rdata', wait = TRUE)
  system('rm -f logs.zip' , wait = TRUE)
  system('rm -f loader.R', wait = TRUE)

  system('rm -rf err', wait = TRUE)
  system('rm -rf out', wait = TRUE)
  system('rm -rf log', wait = TRUE)

  DRrequiredAgeing:::message0('Writing the config file ...')
  save(configlist,file = 'configHadoop.Rdata')

  if (is.null(mp_chooser_file) || !file.exists(mp_chooser_file)) {
    DRrequiredAgeing:::message0(
      'ERROR: mp_chooser not found at location\n\t',
      mp_chooser_file
    )
    quit()
  }

  DRrequiredAgeing:::annotationIndexCreator(
    path = getwd(),
    pattern = 'split_index',
    mem = "5G",
    time = "2-00",
    outputfile = 'annotation_jobs.bch'
  )
  system('chmod 775 annotation_jobs.bch', wait = TRUE)

  DRrequiredAgeing:::message0('Downloading the action script ...')
  system(
    'wget -O loader.R https://raw.githubusercontent.com/mpi2/impc_stats_pipeline/dev/IMPC%20annotation%20pipeline/loaderHadoop.R',
    wait = TRUE
  )

  DRrequiredAgeing:::submit_limit_jobs(bch_file="annotation_jobs.bch", job_id_logfile="../../../../compressed_logs/annotation_job_id.txt")
  DRrequiredAgeing:::waitTillCommandFinish(
    WaitIfTheOutputContains = waitUntillSee,
    ignoreline = ignoreThisLineInWaitingCheck
  )

  DRrequiredAgeing:::message0('Check for errors in the process ...')
  setwd(file.path(SP.results, 'AnnotationExtractorAndHadoopLoader','out'))
  if(system(command = 'grep "exit" * -lR', wait = TRUE)==0)
    stop ('Oh no! we found some errors in the process!')

  DRrequiredAgeing:::message0('Zipping logs ...')
  setwd(file.path(SP.results, 'AnnotationExtractorAndHadoopLoader'))
  system('mv annotation_jobs.bch ../../../../compressed_logs', wait = TRUE)
  system('zip -rm ../../../../compressed_logs/annotation_logs.zip log/* err/* out/*', wait = TRUE)
  system('zip -rm splits.zip split_index_*', wait = TRUE)

  DRrequiredAgeing:::message0('Job done.')

}

ClearReportsAfterCreation = function(dirpath = getwd(),
                                     pattern = 'csv',
                                     sep = ',') {
  files = list.files(
    path = dirpath,
    pattern = pattern,
    full.names = TRUE,
    recursive = FALSE
  )

  for (file in files) {
    message('Processing: ', file)
    df = data.table::fread(input = file, sep = sep)
    df = subset(df, df$status == 'Successful')
    maxmethod = table(df$`applied Method`)
    maxgroup = names(maxmethod)[which.max(maxmethod)]
    df = subset(df, df$`applied Method` == maxgroup)
    df = as.data.frame(df)
    write.csv(df,
              file = file.path(dirname(file), paste0('AllSuccessful_', basename(file))),
              row.names = FALSE)
  }

  return('Done!')
}


HadoopReTransferSCP = function(wd = getwd(),
                               pattern = '_.statpackets',
                               today = format(Sys.time(), '%d%m%Y'),
                               host =  "hadoop-login-02",
                               path = '/hadoop/user/mi_stats/impc/statpackets/',
                               prefix = 'DRXXX_',
                               user =  Sys.info()['user']) {
  library('data.table')
  ########################### Annotation pipeline #################################
  ##############################
  library(data.table)
  library(jsonlite)
  library(rlist)
  library(Tmisc)
  library(rwebhdfs)

  hadoopbase = file.path(path,
                         prefix,
                         today)

  files1 = list.files(pattern = '_.statpackets.gz', full.names = TRUE)
  message('Step1. Deleting .gz files ...')
  for (file1 in files1) {
    unlink(file1)
  }

  command = paste0('scp -r ', wd, ' ', host, ':', hadoopbase)
  message('Step2. Transfering files.  \n',
          command)

  system(paste0('ssh ', host, ' "mkdir -p ', hadoopbase, '"'), wait = TRUE)

  command = gsub(
    pattern = '//',
    replacement = '/',
    x = command,
    fixed = TRUE
  )

  system(command = command,
         wait = TRUE)

  return('Done!')

}

HadoopReTransfer = function(wd = getwd(),
                            today = format(Sys.time(), '%d%m%Y'),
                            host =  "hh-hdp-master-02.ebi.ac.uk",
                            path = 'impc/statpackets',
                            prefix = 'DRXXX_',
                            port = '50070',
                            user =  Sys.info()['user'],
                            password = 'impc') {
  library('data.table')
  ########################### Annotation pipeline #################################
  ##############################
  library(data.table)
  library(jsonlite)
  library(rlist)
  library(Tmisc)
  library(rwebhdfs)

  hadoopbase = file.path(path,
                         prefix,
                         today)

  hdfs <-
    webhdfs(
      namenode_host = host,
      namenode_port = port,
      hdfs_username =  user
    )
  rwebhdfs::mkdir(hdfs, hadoopbase)

  files = list.files(wd, full.names = TRUE,recursive = FALSE,include.dirs = FALSE)
  for (file in files) {
    message('processing:  ', file)
    hadoopPath = file.path(hadoopbase,
                           paste0(basename(file), '_.statpackets'))



    transfered = rwebhdfs::write_file(
      fs = hdfs,
      targetPath = hadoopPath,
      srcPath = file,
      sizeWarn = 10 ^ 30,
      append = FALSE,
      overwrite = TRUE
    )
    gc()

    if (transfered) {
      message('Transfer successful.')
      message('Compressing the temp statpacket file.')
      # system(command = paste0('gzip ', file),
      #        wait = TRUE)
      message('Done!')
    }  else{
      stop('Transfered not successful!')
    }
  }

  return('Done!')

}

#library(quantreg)
extractRiskyGenesFromDRs = function(newDRReportpath = '../DR19_Reports/DR19_AllSuccessful_WithQvaluesAndMPtermsdata_point_of_type_unidimensional.csv.gz',
                                    oldDRReportpath = '../DR18_Reports/DR18StatisticalResultsReportContinuous.csv.gz') {
  file = paste0('RiskyGenesToCheck_',
                DRrequiredAgeing:::RemoveSpecialChars(date()),
                '.txt')
  message(
    'please wait and the results will be displayed on screen (and a version will be stored in ',
    file,
    ')'
  )
  dfold = read.csv(file = oldDRReportpath)
  dfnew = read.csv(file = newDRReportpath)

  dfold = subset(dfold, dfold$status == 'Successful')
  dfnew = subset(dfnew, dfnew$status == 'Successful')

  dfold = subset(dfold, dfold$applied.Method == 'MM')
  dfnew = subset(dfnew, dfnew$applied.Method == 'MM')

  res10 = paste0('There are ',
                 nrow(dfnew) - nrow(dfold),
                 ' more statistical results in the new DR')


  res20 = paste0(
    'There are ',
    sum(!dfold$metadata_group %in% dfnew$metadata_group),
    ' metadata in the old DR that does not exist in the new DR'
  )

  # deep1 = dfold[!dfold$metadata_group %in% dfnew$metadata_group, ]
  # table(deep1$phenotyping_center)


  f1 = function(x) {
    if (is.numeric(x) && !all(is.na(x)))
      r = mean(x, na.rm = TRUE)
    else if ((is.factor(x) || is.character(x)) && !all(is.na(x)))
      r = paste(unique(as.character(x)),
                sep = ' | ',
                collapse = ' | ')
    else
      r = NA

    return (r)
  }

  f2 = function(x) {
    if (!all(is.na(x)))
      r = mean(as.numeric(x), na.rm = TRUE)
    else
      r = NA

    return (r)
  }

  f3 = function(x) {
    if (!all(is.na(x)))
      r = sum(as.numeric(x), na.rm = TRUE)
    else
      r = NA

    return (r)
  }

  f4 = function(df,
                col = 'diff',
                extractcol = 'Group.1',
                n = 10) {
    r = df[order(df[, col]),]
    r0 = c(head(r[, extractcol], n), tail(r[, extractcol], n))
    return (r0)
  }

  aoldmp = aggregate(dfold[, c('Total.KO.female',
                               'Total.KO.male',
                               'Total.WT.female',
                               'Total.WT.male')], by = list(dfold$gene_symbol), f3)

  anewmp = aggregate(dfnew[, c('Total.KO.female',
                               'Total.KO.male',
                               'Total.WT.female',
                               'Total.WT.male')], by = list(dfnew$gene_symbol), f3)


  aoldm = aggregate(dfold[, c('Total.KO.female',
                              'Total.KO.male',
                              'Total.WT.female',
                              'Total.WT.male')], by = list(dfold$gene_symbol), f2)

  anewm = aggregate(dfnew[, c('Total.KO.female',
                              'Total.KO.male',
                              'Total.WT.female',
                              'Total.WT.male')], by = list(dfnew$gene_symbol), f2)


  aoldm = merge(aoldm,
                data.frame('Group.1' = aoldmp$Group.1, 'total' = rowSums(aoldmp[,-1])),
                by = 'Group.1')
  anewm = merge(anewm,
                data.frame('Group.1' = anewmp$Group.1, 'total' = rowSums(anewmp[,-1])),
                by = 'Group.1')


  aold = aggregate(dfold$gene_symbol, by = list(dfold$gene_symbol), length)
  anew = aggregate(dfnew$gene_symbol, by = list(dfnew$gene_symbol), length)


  dfm0 = merge(anewm,
               aoldm,
               by = 'Group.1',
               suffixes = c('_dfnew', '_dfold'))

  dfm1 = merge(anew,
               aold,
               by = 'Group.1',
               suffixes = c('_dfnew', '_dfold'))
  dfm1$diff = dfm1$x_dfnew - dfm1$x_dfold
  dfm12 = dfm1[dfm1$diff != 0,]
  dfm12[order(dfm12$diff),]


  dfm  = merge(dfm0, dfm1, by = 'Group.1')
  dfm = dfm[, order(names(dfm))]

  dfm$KOFemaleDiff = dfm$Total.KO.female_dfnew < dfm$Total.KO.female_dfold
  dfm$KOMaleDiff   = dfm$Total.KO.male_dfnew  < dfm$Total.KO.male_dfold


  res30 = f4(dfm[dfm$KOFemaleDiff, c('Group.1', 'diff')])
  res40 = f4(dfm[dfm$KOMaleDiff, c('Group.1', 'diff')])

  res50 = f4(dfm[dfm$Total.WT.female_dfnew < dfm$Total.WT.female_dfold, c('Group.1', 'diff')])
  res60 = f4(dfm[dfm$Total.WT.male_dfnew < dfm$Total.WT.male_dfold, c('Group.1', 'diff')])


  # dif  = log(dfm1$x_dfnew) - log(dfm1$x_dfold)
  # hist(dif, breaks = 100)
  # summary(dif)
  # plot(dfm1$xdfnew,dfm1$xdfold,pch='.')
  # qline=rq(dfm1$xdfnew~dfm1$xdfold)
  # abline(qline)
  #
  # dfchange = dfm[dfm$x_dfnew < dfm$x_dfold &
  #                  dfm$total_dfnew < dfm$total_dfold,]
  # dfchange

  result = unique(c(
    res10,
    res20,
    "Check these risky genes: ",
    res30,
    res40,
    res50,
    res60
  ))

  print(result)
  print(paste0('Check: ./', file))
  write(result,
        file)
  return(invisible(result))
}
