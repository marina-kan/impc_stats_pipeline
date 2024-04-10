## You must check the 'check' parameter
mainAgeing = function(file = NULL                                    ,
                      sep = ','                                      ,
                      na.strings = 'NA'                              ,
                      normalisedPhenlist = FALSE                     ,
                      subdir = 'Results'                             ,
                      seed = 123456                                  ,
                      MethodOfReadingCategoricalCategories  = 'file' , # `file`, `solr` or `update`
                      OverwriteExistingFiles  = FALSE                ,
                      ignoreSkipList          = FALSE                ,
                      onlyFillNotExisitingResults = FALSE            ,
                      WhiteListMethods  = NULL                       ,
                      # OpenStats
                      MMOptimise              = c(1, 1, 1, 1, 1, 1) ,
                      FERROptimise            = c(TRUE,TRUE)        ,
                      FERRrep                 = 1500       ,
                      equation                = 'auto'     ,
                      # Only for simulations
                      simulation              = FALSE      ,
                      Simulation.iteration    = 1          ,
                      skiptimeseries          = TRUE       ,
                      # Only for windowing
                      activeWindowing = FALSE              ,
                      sensitivity = c(1, 1, 1, 0)          ,
                      pvalThreshold = c(0, 0, 0, 0)        ,
                      check = 1                            ,
                      direction = c(1, 1)                  ,
                      weightORthreshold = 'threshold'      ,
                      predFunction = function(m) {
                        predict(m)
                      }                                    ,
                      residFunction = function(m) {
                        resid(m)
                      }                                    ,
                      messages = TRUE                      ,
                      threshold = sqrt(.Machine$double.eps) * 10,
                      outdelim = '\t'                      ,
                      debug = TRUE                         ,
                      encode = FALSE                       ,
                      noSpaceAllowed = TRUE                ,
                      plotWindowing = TRUE                 ,
                      storeplot = TRUE                     ,
                      virtualDrive   = FALSE               ,
                      checkNamesForMissingColNames = TRUE  ,
                      # Just for the OpenStats manuscript
                      measureExecutionTime   = FALSE       ,
                      # Raw data
                      storeRawData           = TRUE        ,
                      compressRawData        = TRUE        ,
                      writeOutputToDB        = FALSE       ,
                      # Only for Batch generator
                      BatchProducer          =  FALSE      ,
                      cpu = 1                              ,
                      memory = "7G"                        ,
                      time = "08:00:00"                    ,
                      nMax                   = 10000       ,
                      ChunkSize              = 24          ,
                      MinColoniesInChunks    = 32          ,
                      controlSize            = 1500        ,
                      ### Just for outlier detection
                      outlierDetection       = FALSE       ,
                      ### Just for Ageing Batch GEnerator  ,
                      combineEAandLA        = FALSE        ,
                      solrBaseURL           = 'http://hx-noah-74-10:8090' ,
                      ### Just for debuging
                      superDebug             = FALSE                      ,
                      subBreakColumns        = NULL,
                      extraBatchParameters   = NULL,
                      DRversion              = 'not_specified',
                      ...) {
  message0('DRrequiredAgeing loaded')
  message0(
    Sys.time(),
    '  ############################################################\n',
    Sys.time(),
    '  % Make sure ExceptionList and CategoryList are updated!    %\n',
    Sys.time(),
    '  % Please check SENSITIVITY for the windowing algorithm     %\n',
    Sys.time(),
    '  % Sensitivy:',
    sensitivity,
    '                                                              \n',
    Sys.time(),
    '  ############################################################',
    ShowTimeTag = FALSE
  )
  message0('Process started ...')
  message0('Machine info:  ', paste(Sys.info(), collapse = ', '))
  message0('Loading dependent packages ...')
  requireNamespace('OpenStats'     )
  requireNamespace('foreach')
  requireNamespace('SmoothWin')
  requireNamespace('nlme')
  requireNamespace('base64enc')
  requireNamespace('RJSONIO'    )
  requireNamespace('jsonlite'   )
  requireNamespace('DBI'        )
  # Config files
  message0('Loading configuration ...')
  methodmap                      = readConf('MethodMap.conf')
  equationmap                    = readConf('EquationMap.conf')
  CategoryMap                    = readConf('CategoryMap.conf')
  MergeCategoryParameters        = read.csv(file = file.path(local(), 'MergeParameterList.txt'))
  initial                        = readConf('Initialize.conf')
  exceptionList                  = readFile(file = 'ExceptionMap.list')
  EA2LAMApping                   = read.csv(file = file.path(local(), 'EA2LA_parameter_mappings_2019-09-24.csv'))
  MetaDataList                   = read.csv(file = file.path(local(), 'metadataParameters.csv'))
  exceptionList                  = unique(c(exceptionList,MetaDataList$parameter_stable_id))
  
  # Main subdirectory/working directory
  message0('Preparing the working directory ...')
  cwd = getwd()
  wd  = file.path(cwd,
                  paste(subdir, sep = '_', collapse = '_'))
  dir.create0(wd, recursive = TRUE)
  wd = CreateVirtualDrive(active = virtualDrive,currentwd = wd)
  message0('Setting the working directory to: \n\t\t ===> ', wd)
  setwd(dir = wd)
  ################## The rnd must be above seed!
  initialRandomValue = round(runif(1,min = 1000,max = 9999))
  set.seed(seed)
  # Read file
  rdata = readInputDatafromFile(
    file = file,
    checkNamesForMissingColNames = checkNamesForMissingColNames,
    sep = sep,
    na.strings = na.strings
  )
  rdata                 = rdata[!is.na(rdata$phenotyping_center), ] # Just to remove NA centers
  if (!BatchProducer) {
    if ('DRversion' %in% names(rdata)) {
      message0('`DRversion` column already exists in the data and will be used!')
    } else  {
      rdata$DRversion       = DRversion
    }
  }
  new.data              = rdata
  new.data              = new.data[order(Date2Integer(new.data$date_of_experiment)), ]
  #########
  new.data$colony_id    = as.character(new.data$colony_id)
  new.data$external_sample_id = as.factor(new.data$external_sample_id)
  new.data$observation_id     = as.factor(new.data$observation_id)
  new.data$metadata = as.character(new.data$metadata)
  ################
  # Start analysis
  ################
  closeAllConnections()
  # Get possible categories for the categorical variables
  CatList = GetPossibleCategories (procedure = NULL, method = MethodOfReadingCategoricalCategories)
  message0('Filtering the dataset in progress ....')
  Strtime      = Sys.time()
  procedures   = as.character(unique(na.omit(new.data$procedure_stable_id)))
  for (procedure in procedures) {
    StrtimePro  = Sys.time()
    ###
    n2.9 = base::subset(new.data,  new.data$procedure_stable_id %in% procedure)
    parameters  = as.character(unique(na.omit(n2.9$parameter_stable_id)))
    for (parameter in parameters) {
      StrtimePar = Sys.time()
      FactorLevels = ReadFactorLevelsFromSolr(parameter = parameter, CatList = CatList)
      ### counter starts here ....
      counter   = 1
      outP      = list()
      n3.0 = base::subset(n2.9,  n2.9$parameter_stable_id %in% parameter)
      ############## Read The Ageing parameters from Solr
      if(BatchProducer && combineEAandLA)
        n3.0 = getEarlyAdultsFromParameterStableIds(
          LA_parameter_stable_id = parameter    ,
          map                    = EA2LAMApping ,
          LA_data                =  n3.0        ,
          solrBaseURL            = solrBaseURL
        )
      if(is.null(n3.0))
        next
      ##############
      centers   = as.character(unique(na.omit(n3.0$phenotyping_center)))
      for (center in centers) {
        n3.1     = base::subset(n3.0, n3.0$phenotyping_center %in% center)
        pipelines = as.character(unique(na.omit(n3.1$pipeline_stable_id)))
        for (pipeline in pipelines) {
          n3.11 = base::subset(n3.1, n3.1$pipeline_stable_id %in% pipeline)
          strains  = as.character(unique(na.omit(n3.11$strain_accession_id)))
          for (strain in strains) {
            n3.2  = base::subset(n3.11,  n3.11$strain_accession_id %in% strain)
            metas = as.character(unique(na.omit(n3.2$metadata_group)))
            for (meta in metas) {
              n3.3   = base::subset(n3.2,  n3.2$metadata_group %in% meta)
              n3.3.c = base::subset(n3.3,  n3.3$biological_sample_group %in% 'control')
              n3.3.m = base::subset(n3.3,!(n3.3$biological_sample_group %in% 'control'))
              zygositys = as.character(unique(na.omit(n3.3.m$zygosity)))
              for (zyg in zygositys) {
                n3.3.m_zyg = base::subset(n3.3.m, n3.3.m$zygosity %in% zyg)
                colonys    = as.character(unique(na.omit(n3.3.m_zyg$colony_id)))
                nColonies  = length(colonys)
                if (BatchProducer && nColonies > 0) {
                  #nMax                 = 10000
                  ChunkSizeFromNumbers = ((nrow(n3.3.c) < nMax) * max(1, round(nColonies /
                                                                                 ChunkSize)) +
                                            (nrow(n3.3.c) >= nMax) * nColonies)
                  minCol = ((nrow(n3.3.c) < nMax) * MinColoniesInChunks + (nrow(n3.3.c) >=
                                                                             nMax) * 1)
                  ColonyChunks         = chunkVector(
                    x = colonys,
                    n = ChunkSizeFromNumbers,
                    min = minCol,
                    activate = (nColonies >= MinColoniesInChunks) &&
                      (nrow(n3.3.c) >= controlSize)
                  )
                  outpDir  = file.path0(
                    wd,
                    paste0(Sys.Date(), '_', subdir, '_RawData/'),
                    check = FALSE,
                    create = TRUE,
                    IncludedFileName = FALSE
                  )
                  SubSubDirOrdFileName = RemoveSpecialChars(
                    paste(
                      Sys.Date()             ,
                      #RandomRegardSeed()    ,
                      #procedure             ,
                      parameter              ,
                      center                 ,
                      pipeline               ,
                      zyg                    ,
                      strain                 ,
                      meta                   ,
                      collapse = '_'
                    )
                  )
                  outpfile = file.path0(
                    outpDir,
                    SubSubDirOrdFileName,
                    check = FALSE,
                    create = TRUE,
                    IncludedFileName = TRUE
                  )
                  mess = paste0(
                    Sys.time(),
                    '. Processed file: ',
                    SubSubDirOrdFileName,
                    '. #Colonies = ',
                    nColonies,
                    ', #Controls = ',
                    nrow(n3.3.c),
                    ', Chunks = ',
                    length(ColonyChunks)
                  )
                  message0(mess, ShowTimeTag = FALSE)
                  write(
                    x = mess,
                    file = paste0(Sys.Date(),
                                  '_',
                                  subdir,
                                  '_DataGenerationLog.log'),
                    10 ^ 5,
                    append = TRUE
                  )
                  counter = 1
                  for (ChunkedColonies in ColonyChunks) {
                    BatchData     = rbind (
                      subset(
                        n3.3.m_zyg,
                        n3.3.m_zyg$colony_id %in% ChunkedColonies
                      ),
                      n3.3.c
                    )
                    BatchFileName = file.exists0(
                      paste0(
                        outpfile,
                        '_C',
                        length(ColonyChunks),
                        '_',
                        RandomRegardSeed(),
                        '_',
                        counter,
                        '.csv'
                      ),
                      overwrite = OverwriteExistingFiles
                    )
                    if (all(dim(BatchData) > 0)) {
                      write.csv(BatchData,
                                file = BatchFileName,
                                row.names = FALSE)

                      BatchFileNamezip  = compressFiles(
                        fileList   = BatchFileName,
                        dir        = dirname (BatchFileName) ,
                        filename   = basename(BatchFileName) ,
                        overwrite  = FALSE                   ,
                        rmSource   = TRUE
                      )
                      if (BatchFileNamezip$status == 0) {
                        BatchFileName = BatchFileNamezip$file
                      }
                      out = BatchGenerator(
                        file                 = BatchFileName    ,
                        dir                  = outpDir          ,
                        procedure            = procedure        ,
                        parameter            = parameter        ,
                        center               = center           ,
                        cpu                  = cpu              ,
                        memory               = memory           ,
                        time                 = time             ,
                        extraBatchParameters =  extraBatchParameters
                      )
                      write(
                        x = out,
                        file = paste0(outpDir, '/', subdir, '_Batch.bch'),
                        ncolumns = 10 ^ 5,
                        append = TRUE
                      )
                      counter = counter + 1
                    }
                    rm0(c('BatchData', 'BatchFileName'), silent = TRUE)
                    gc()
                  }
                } else{
                  message0(
                    ' [',
                    paste(
                      procedure,
                      parameter,
                      center   ,
                      pipeline ,
                      strain   ,
                      meta     ,
                      zyg      ,
                      length(colonys),
                      sep = ']~>['
                    ),
                    ']'
                  )
                  ### Single core
                  i = 1
                  MultiCoreRes = foreach::foreach (
                    i = 1:length(colonys),
                    .packages = c(
                      'SmoothWin'   ,
                      'base64enc'   ,
                      'nlme'        ,
                      'RJSONIO'     ,
                      'jsonlite'    ,
                      'OpenStats'  ,
                      'DRrequiredAgeing',
                      'DBI'
                    ),
                    .errorhandling = c('stop'),
                    .verbose = TRUE,
                    .inorder = FALSE
                  ) %do% {
                    # for (i in  1:length(colonys)){
                    message0('*~*~*~*~*~* ',
                             i,
                             '|',
                             length(colonys),
                             ' *~*~*~*~*~*')
                    for (sim.index in 1:ifelse(simulation, Simulation.iteration, 1)) {
                      # Removing the old objects if exist
                      ObjectsThatMustBeRemovedInEachIteration()
                      # Initialization before starting the analysis
                      note   = list()
                      note = c(note,
                               list('Random seed' = seed))
                      colony = colonys[i]
                      message0('Current colony: ', colony)

                      n3.4 = base::subset(n3.3.m_zyg,
                                          n3.3.m_zyg$colony_id %in% c(colony))
                      n3.5 = sortDataset(x = rbind (n3.4, n3.3.c),
                                         BatchCol = 'date_of_experiment')
                      note = c(
                        note,
                        list(
                          'Bodyweight included in the input data' = CheckIfNameExistInDataFrame(
                            obj   = n3.5,
                            name  = 'weight',
                            checkLevels = FALSE
                          )
                        )
                      )
                      # Imaginary URLs
                      note$'Gene page URL'        = GenePageURL       = GenePageURL(n3.5)
                      note$'Body weight page URL' = BodyWeightCurvURL = BodyWeightCurvURL(n3.5)

                      ReadMeTxt       = ReadMe (obj = n3.4, URL = GenePageURL)

                      # Define response column [do not move me!]
                      depVariable = getResponseColumn(n3.5$observation_type)
                      ################

                      n3.5 = TransformVariableByFunction(
                        varType = depVariable$lbl,
                        data = n3.5              ,
                        types = c('unidimensional', 'time_series'),
                        colName = depVariable$column              ,
                        FUN = as.numeric0
                      )
                      depVar      = depVariable$column
                      message0('Dependent variable: ', depVar)
                      note$'Response type'         = paste0(depVar,
                                                            '_of_type_',
                                                            paste(depVariable$lbl, sep = '.'))
                      note$'Observation type'    =
                        if (!is.null(unique(n3.5$observation_type))) {
                          paste(
                            unique(n3.5$observation_type),
                            sep = '~',
                            collapse = '~'
                          )
                        } else{
                          NULL
                        }
                      note$'Data type'           =
                        if (!is.null(unique(n3.5$data_type))) {
                          paste(unique(n3.5$data_type),
                                sep = '~',
                                collapse = '~')
                        } else{
                          NULL
                        }
                      # ABR min decreased to 1. 19/8/2020 meeting with Jeremy and Federico
                      minSampRequired = ifelse(
                        is.ABR(x = parameter),
                        as.numeric(initial$min_ABR_mut_each_sex),
                        ifelse(
                          nrow(n3.5) > 0            &&
                            depVar %in% names(n3.5) &&
                            is.numeric(n3.5[, depVar]),
                          as.numeric(initial$min_num_mut_each_sex),
                          0
                        )
                      )

                      # add missing levels to categorical variables
                      if (depVar %in% names(n3.5) &&
                          !is.numeric(n3.5[, depVar])) {
                        AllLevels = mapLevelsToFactor(levels = levels(n3.5[, depVar]),
                                                      newlevels = FactorLevels$levels)
                        levels(n3.5[, depVar]) = AllLevels$levels
                        #####
                        note = c(note,
                                 FactorLevels$note,
                                 AllLevels$note)
                        #####
                        SexGenResLevels = min(2 * 2 * length(AllLevels$levels), 4)
                      } else{
                        SexGenResLevels = 4
                      }
                      if (!depVariable$accepted) {
                        write(
                          paste(
                            ReadMeTxt,
                            sep = '\t',
                            collapse = '\t'
                          ),
                          file = 'NotProcessedFileImproperDataType.log',
                          append = TRUE
                        )
                        return('Not a proper dataset!')
                      }

                      if (simulation && is.numeric(n3.5[, depVar])) {
                        message0('Simulation in progress ... Round ',
                                 sim.index)
                        n3.5_tmp = mimicControls(
                          df = n3.5,
                          removeMutants = (sim.index == 1) ,
                          ArtifLabel  = 'experimental'     ,
                          mutLabel = 'experimental'        ,
                          baselines = 'control'            ,
                          neutralise = TRUE                ,
                          resample = (sim.index != 1)      ,
                          depVariable = depVar             ,
                          sex  = 'sex'                     ,
                          minSampRequired = minSampRequired,
                          SexGenResLevels = SexGenResLevels,
                          indicator       = sim.index      ,
                          plot = superDebug
                        )
                        n3.5 = n3.5_tmp$df
                        note = list(note , 'Simulation details' = n3.5_tmp$note)
                      }

                      # Summary statistics
                      n3.5_summary = SummaryStatisticsOriginal(x      = n3.5,
                                                               depVar = depVar,
                                                               label  = 'Raw data summary statistics')
                      note         = c(note, n3.5_summary)
                      if (CheckIfNameExistInDataFrame(n3.5, 'LifeStage')) {
                        LifeStageTable = table(n3.5$LifeStage)
                        message0(
                          'LifeStage Summary: ',
                          paste(
                            names(LifeStageTable),
                            LifeStageTable       ,
                            sep      = ':'       ,
                            collapse = ', '
                          )
                        )
                      }
                      # Remove zero frequency categories
                      n3.5.1_F_list = RemoveZeroFrequencyCategories(
                        x               = n3.5,
                        minSampRequired = minSampRequired,
                        depVar      = depVar,
                        totalLevels = SexGenResLevels
                      )
                      n3.5.1 = n3.5.1_F_list$x
                      note   = c(note, n3.5.1_F_list$note)

                      # Remove var categories
                      n3.5.1_v_list = RemoveZerovarCategories(
                        x = n3.5.1,
                        depVar = depVar,
                        minvar = 0,
                        method = getMethodi(
                          var = parameter,
                          type = ifelse(is.numeric(n3.5.1[, depVar]),
                                        'numeric',
                                        'charachter'),
                          methodMap = methodmap
                        )
                      )
                      n3.5.1 = n3.5.1_v_list$x
                      note   = c(note, n3.5.1_v_list$note)

                      OrgSpecIds      = OtherExtraColumns(
                        obj = n3.5,
                        ColNames = c(
                          'external_sample_id',
                          'observation_id',
                          'sex',
                          'biological_sample_group',
                          depVar,
                          'date_of_experiment',
                          'weight',
                          'metadata',
                          'discrete_point',
                          'time_point'
                        ),
                        names = c(
                          # all lower case
                          'Original external_sample_id',
                          'Original observation_id',
                          'Original sex',
                          'Original biological_sample_group',
                          'Original response',
                          'Original date_of_experiment',
                          'Original body weight',
                          'Original metadata',
                          'Original discrete_point',
                          'Original time_point'
                        )
                      )
                      note = c(note, OrgSpecIds)
                      message0('Creating output directory and file name ...')
                      SubSubDirOrdFileName = file.path0(
                        RemoveSpecialChars(center)    ,
                        RemoveSpecialChars(procedure) ,
                        RemoveSpecialChars(pipeline)  ,
                        RemoveSpecialChars(parameter) ,
                        RemoveSpecialChars(colony)    ,
                        RemoveSpecialChars(zyg)       ,
                        RemoveSpecialChars(meta)      ,
                        create = FALSE,
                        check = noSpaceAllowed
                      )
                      FileName = 'output'
                      outDir   = file.path0(
                        wd,
                        SubSubDirOrdFileName,
                        create           = TRUE,
                        check            = FALSE,
                        IncludedFileName = FALSE
                      )
                      outpfile = outpfile2 = paste0(outDir,
                                                    '/',
                                                    FileName,
                                                    collapse = '')
                      message0('Output directory: \n \t\t =>=>=> ',
                               outpfile)
                      if (onlyFillNotExisitingResults) {
                        if (any(file.exists(paste(
                          outpfile,
                          c('NotProcessed.tsv', 'Successful.tsv'),
                          sep = '_'
                        )))) {
                          message0('File already exists then skipped!')
                          return(NULL)
                        } else{
                          message0('Result does not exist! Adding in progress ...')
                          rmme = lapply(list.files(dirname(outpfile), full.names = TRUE), function(x) {
                            if (!is.null(x)    &&
                                file.exists(x) &&
                                (
                                  grepl(
                                    pattern = '.Rdata',
                                    x = x,
                                    fixed = TRUE
                                  ) ||
                                  grepl(
                                    pattern = 'Failed_critical_error',
                                    x = x,
                                    fixed = TRUE
                                  )
                                ))
                            file.remove(x)
                          })
                          write(outpfile,
                                file = 'DoesNotExists.log',
                                append = TRUE)
                        }
                      }

                      ####
                      if (storeRawData) {
                        # There is a second snippet for the rawdata + weights
                        RawoutputFile = RawoutputFile0 = file.exists0(paste(outpfile2,
                                                                            'rawData.csv',
                                                                            sep = '_'),
                                                                      overwrite  = OverwriteExistingFiles)
                        ReadMeFile  = file.exists0(file.path(dirname(RawoutputFile), 'ReadMe.txt'))
                        message0('writting the raw data file to disk ... \n \t\t ===> ',
                                 RawoutputFile)
                        write.csv(x = n3.5            ,
                                  row.names = FALSE   ,
                                  file = RawoutputFile)
                        write(x        = ReadMeTxt ,
                              file     = ReadMeFile,
                              ncolumns = 1)
                        if (compressRawData) {
                          comRes = compressFiles(
                            fileList  = c(ReadMeFile , RawoutputFile),
                            dir       = dirname (RawoutputFile)      ,
                            filename  = basename(RawoutputFile)      ,
                            overwrite = OverwriteExistingFiles
                          )
                          if (comRes$status == 0)
                            message0('Compression successful')
                          else
                            message0('Compression failed')
                          RawoutputFile0 = comRes$file
                        }
                      }
                      note$'Input file'             = relativePath(path = file, reference  = wd)
                      note$'Exported raw data file' = relativePath(path = if (storeRawData) {
                        RawoutputFile0
                      } else{
                        NULL
                      },
                      reference = wd)
                      note$'Readme file'         = relativePath(path = if (storeRawData &&
                                                                           !compressRawData) {
                        ReadMeFile
                      } else{
                        NULL
                      },
                      reference = wd)
                      ###'
                      isException    = IsInList(
                        item = c(parameter, procedure),
                        list = exceptionList,
                        message = 'Value found in the skip list'
                      ) && !ignoreSkipList
                      n3.5.2  = droplevels0(n3.5.1)
                      MergLev = MergeLevels(
                        x = n3.5.2[, depVar]                     ,
                        listOfLevelMaps = CategoryMap            ,
                        parameter_stable_id = parameter          ,
                        AllowedParametersList  = MergeCategoryParameters$parameter_stable_id,
                        report = TRUE
                      )
                      ###
                      n3.5.2[, depVar] = MergLev$x
                      n3.5.2           = droplevels0(n3.5.2[!is.na(n3.5.2[, depVar]), ])
                      n3.5.2OnlyKO     = subset(n3.5.2,
                                                n3.5.2$biological_sample_group %in% 'experimental')
                      note$'Relabeled levels for categorical variables'  = MergLev$note
                      if (!is.null(n3.5.2) &&
                          # data.frame is not zero
                          min0(dim(n3.5.2)) > 0 &&
                          # is it  really exist!
                          length(unique(n3.5.2$biological_sample_group)) > 1 &&
                          # include mut and cont
                          min0(table(n3.5.2$biological_sample_group)) >= minSampRequired &&
                          max0(table(
                            n3.5.2OnlyKO$biological_sample_group,
                            n3.5.2OnlyKO$sex
                          )) > 1 &&
                          # include at least 4/2 of each genotype
                          #length(unique(n3.5.2$colony_id)) > 1  &&
                          length(RepBlank(
                            unique(n3.5.2$colony_id),
                            match = c('', NA, 'NA')
                          )) > 1 &&
                          # include 2 colonies (cont & mut)
                          checkGenInCol(n3.5.2) &&
                          # each sex and genotype
                          depVariable$accepted  &&
                          length(na.omit(n3.5.2[, depVar])) > 0 &&
                          # response is not empty!
                          # there must be variation in data
                          NonZeroVariation(n3.5.2[, depVar]) &&
                          !isException &&
                          columnLevelsVariationRadio(dataset = n3.5.2, columnName = depVar) > 0.005 &&
                          RR_thresholdCheck(
                            data = n3.5.2,
                            depVar = depVar,
                            parameter = parameter,
                            methodmap = methodmap
                          )$'Criteria result' &&
                          ifelse (skiptimeseries &&
                                  depVariable$lbl %in% 'time_series',
                                  FALSE,
                                  TRUE)) {
                        message0('Analysing the dataset in progress ...')
                        message0('Creating OpenStats object ...')
                        a = OpenStats::OpenStatsList(
                          OpenStats:::RemoveSexWithZeroDataPointInGenSexTableOnlyStatsPipelinenotExposed(
                            n3.5.2,
                            cols = c('biological_sample_group', 'sex')
                          ),
                          testGenotype             = 'experimental',
                          refGenotype              = 'control',
                          dataset.colname.genotype = 'biological_sample_group',
                          dataset.colname.sex      = 'sex',
                          dataset.colname.weight   = 'weight',
                          dataset.colname.batch    = 'date_of_experiment'
                        )
                        a_summary_before_concurrent = SummaryStatisticsOriginal(
                          x = a@datasetPL,
                          depVar = depVar,
                          sex = 'Sex',
                          genotype = 'Genotype',
                          label = 'OpenStatsList object summary statistics'
                        )
                        note = c(note, a_summary_before_concurrent)
                        #
                        PhenListSpecIds = OtherExtraColumns (
                          obj = a@datasetPL,
                          ColNames = c('external_sample_id', 'observation_id'),
                          names    = c(
                            'OpenStatsList external_sample_id',
                            'OpenStatsList observation_id'
                          )
                        )
                        note = c(note, PhenListSpecIds)
                        ### Get method of analysis
                        method = getMethodi(
                          var = parameter,
                          type = ifelse(
                            is.numeric(a@datasetPL[, depVar]),
                            'numeric',
                            'charachter'
                          ),
                          methodMap = methodmap
                        )
                        # WhiteList methods
                        if (!is.null(WhiteListMethods) &&
                            !(method %in% WhiteListMethods)) {
                          message0('Black list applied. Method = ', method)
                          return(FALSE)
                        }
                        #### Check for concurrent control selection
                        concurrentControlSelect = FALSE
                        aTmp = concurrentContSelect(
                          activate = concurrentControlSelect &&
                            (method %in% 'MM') && !activeWindowing,
                          PhenListObj = a,
                          depVar = depVar,
                          minSampRequired = minSampRequired
                        )
                        a    = aTmp$obj
                        note = c(note, aTmp$note)
                        #### END OF concurrent

                        # Just before analysing data (because cuncurrent sampling)
                        if (concurrentControlSelect &&
                            (method %in% 'MM') && !activeWindowing) {
                          a_phenlist_concurrent_summary = SummaryStatisticsOriginal(
                            x = a@datasetPL,
                            depVar = depVar,
                            sex = 'Sex',
                            genotype = 'Genotype',
                            label = 'OpenStatList and cuncurrent data summary statistics'
                          )
                          note  =  c(note, a_phenlist_concurrent_summary)
                        }
                        # check the Weight column
                        message0('Checking whether Weight column exists in the raw data ...')
                        if (!CheckIfNameExistInDataFrame(a@datasetPL, 'Weight')) {
                          note$'Existence of the weight column in the OpenStatsList object' =
                            'Weight column does not exist in the raw data'
                        }
                        # Equation type
                        message0('equation is set to ', equation)
                        if (equation == 'auto') {
                          equationType = ifelse(
                            CheckIfNameExistInDataFrame(a@datasetPL, 'Weight')         &&
                              MissingPercent(
                                var = 'Weight',
                                data = a@datasetPL
                              ) <= .2 &&
                              EnoughWeightForTheSexGenInteraction(a@datasetPL)         ,
                            getEquation(var = parameter,
                                        equationMap = equationmap),
                            'withoutWeight'
                          )
                        } else{
                          equationType = equation
                        }
                        # This is the main engine!
                        note = c(
                          note,
                          list(
                            'Bodyweight initially included in the full model' = ifelse(method %in% 'MM', equationType, FALSE)
                          )
                        )
                        if (normalisedPhenlist) {
                          a = normalisePhenList(phenlist = a,
                                                colnames = c(depVar, 'Weight'))
                        }
                        message0('Fitting the model ...')
                        message0('Method: ',
                                 method,
                                 ', Equation:',
                                 equationType)
                        c.ww0 =	PhenStatWindow(
                          phenlistObject = a,
                          parameter = parameter,
                          minObs = minSampRequired,
                          method = method,
                          depVariable = depVar,
                          equation = equationType,
                          threshold = threshold,
                          pvalThreshold = pvalThreshold,
                          sensitivity = sensitivity,
                          messages = messages,
                          main = paste(
                            unique(n3.5.2$procedure_name)[1],
                            '\n',
                            unique(n3.5.2$parameter_name)[1],
                            '\n',
                            unique(center)[1],
                            unique(colony)[1],
                            unique(zyg)[1],
                            sep = '-',
                            collapse = ','
                          ),
                          seed = seed,
                          check = check,
                          storeplot = storeplot,
                          windowing = activeWindowing,
                          plot = plotWindowing,
                          PicDir =  dirname(outpfile),
                          filename = basename(outpfile),
                          OverwriteExistingFiles = OverwriteExistingFiles,
                          superDebug = superDebug,
                          predFunction = predFunction,
                          residFunction = residFunction,
                          weightORthreshold = weightORthreshold,
                          direction  = direction,
                          outlierDetection = outlierDetection,
                          FERROptimise = FERROptimise,
                          MMOptimise = MMOptimise,
                          FERRrep = FERRrep
                        )
                        note = c(
                          note                               ,
                          c.ww0$note                         ,
                          'Applied method' = c.ww0$method    ,
                          'Image URL'      = relativePath(
                            path = c.ww0$graphFileName ,
                            reference = wd
                          )
                        )
                        ExtraCols = c('external_sample_id', 'observation_id')
                        ####
                        message0('Preparing the output from VectorOutput function ...')
                        c.ww.vec       = VectorOutput0(
                          c.ww0        = c.ww0                             ,
                          ExtraCols    = ExtraCols                         ,
                          activeWindowing = activeWindowing                ,
                          subBreakColumns = subBreakColumns
                        )
                        ##
                        outP    = SuccessfulOutput(
                          args = list(
                            c.ww0 = c.ww0                                  ,
                            depVariable = depVariable                      ,
                            depVar      = depVar                           ,
                            c.ww.vec = c.ww.vec                            ,
                            procedure = procedure                          ,
                            parameter = parameter                          ,
                            center  = center                               ,
                            n3.5   = n3.5                                  ,
                            strain = strain                                ,
                            meta =   meta                                  ,
                            zyg = zyg                                      ,
                            colony = colony                                ,
                            note = note                                    ,
                            PhenListSpecIds = PhenListSpecIds              ,
                            OrgSpecIds = OrgSpecIds                        ,
                            BodyWeightCurvURL = BodyWeightCurvURL          ,
                            GenePageURL = GenePageURL                      ,
                            encode = encode                                ,
                            wd     = wd
                          ),
                          writeOutputToDB = writeOutputToDB
                        )# c(as.list(environment()), ls()))
                        if ((
                          NullOrError(c.ww0$NormalObj)          ||
                          !NullOrError(c.ww0$NormalObj$messages) ||
                          (
                            NullOrError(c.ww0$WindowedObj) &&
                            activeWindowing &&
                            (c.ww0$method %in% 'MM')
                          )
                        ) &&
                        debug)
                          save(
                            n3.5,
                            file = paste0(
                              outpfile,
                              RemoveSpecialChars(colony),
                              '_',
                              RandomRegardSeed(),
                              '.Rdata'
                            )
                          )
                        ##
                        StoreRawDataAndWindowingWeights(
                          storeRawData = storeRawData,
                          activeWindowing = activeWindowing,
                          c.ww0 = c.ww0,
                          RawoutputFile = RawoutputFile,
                          orgData = n3.5.2,
                          compressRawData = compressRawData,
                          files  = c(RawoutputFile, ReadMeFile),
                          dir        = dirname (RawoutputFile) ,
                          filename   = basename(RawoutputFile) ,
                          ReadMeTxt  = ReadMeTxt,
                          ReadMeFile = ReadMeFile,
                          methodmap  = methodmap

                        )
                        StatusSF = !NullOrError(c.ww0$NormalObj) &&
                          NullOrError(c.ww0$NormalObj$messages)
                        SucFaiFile = paste(
                          outpfile2,
                          ifelse(
                            StatusSF,
                            'Successful.tsv',
                            'Failed_critical_error.tsv'
                          ),
                          sep =
                            '_'
                        )
                        if (!StatusSF) {
                          write(
                            x    = SucFaiFile,
                            file = file.path(
                              wd,
                              paste0('Failed_analyses_', Sys.Date(), '.log')
                            ),
                            append = TRUE
                          )
                        }
                        write.table(
                          x = paste(outP,
                                    collapse =   outdelim),
                          file = file.exists0(SucFaiFile, overwrite = OverwriteExistingFiles),
                          append = TRUE,
                          row.names = FALSE,
                          col.names = FALSE,
                          quote = FALSE
                        )
                      } else {
                        message0('Dataset not processed ...')
                        optFail     = NotProcessedOutput(
                          args = list(
                            procedure         = procedure                           ,
                            parameter         = parameter                           ,
                            center            = center                              ,
                            n3.5.2            = n3.5.2                              ,
                            n3.5              = n3.5                                ,
                            strain            = strain                              ,
                            meta              = meta                                ,
                            zyg               = zyg                                 ,
                            colony            = colony                              ,
                            depVariable       = depVariable                         ,
                            #c.ww.vec         = c.ww.vec                            ,
                            note              = note                                ,
                            isException       = isException                         ,
                            minSampRequired   = minSampRequired                     ,
                            depVar            = depVar                              ,
                            GenePageURL       = GenePageURL                         ,
                            BodyWeightCurvURL = BodyWeightCurvURL                   ,
                            OrgSpecIds        = OrgSpecIds                          ,
                            encode 			      = encode                              ,
                            methodmap         = methodmap
                          ),
                          writeOutputToDB = writeOutputToDB
                        )
                        #optFail     = NotProcessedOutput(args = c(as.list(environment()), ls()))
                        NotProcFile = paste(outpfile2,
                                            'NotProcessed.tsv',
                                            sep =
                                              '_')
                        write.table(
                          paste(optFail,
                                collapse =   outdelim),
                          file =   file.exists0(NotProcFile, overwrite = OverwriteExistingFiles),
                          append = TRUE,
                          row.names = FALSE,
                          col.names = FALSE,
                          quote = FALSE
                        )
                        if (simulation)
                          break
                      }
                    }

                    message0(
                      'Finished in ',
                      round(difftime(
                        Sys.time() , Strtime, units = 'sec'
                      ), 2),
                      '(s).\n
                  -----------------------------------
                  \n\n '
                    )
                    counter  = counter  + 1
                    gc()
                  }
                }
              }
            }
          }
        }
      }
      RecordSpentTime(
        timeSt    = StrtimePar                ,
        dirName   = 'OpenStatsParameterTime'  ,
        fileName  = c(procedure, parameter)   ,
        rnd       = initialRandomValue        ,
        active    = measureExecutionTime
      )
    }
    RecordSpentTime(
      timeSt    = StrtimePro              ,
      dirName   = 'OpenStatsProcedureTime',
      fileName  = procedure               ,
      rnd       = initialRandomValue      ,
      active    = measureExecutionTime

    )
  }
  message0('Closing Connections ...')
  closeAllConnections()
  message0('Finished.')
  setwd(cwd)
  message0('Cleaning the memory ...')
  gc()
}


