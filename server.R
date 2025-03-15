library(ClusterGVis)
library(ggplot2)
library(org.Hs.eg.db)
library(org.Mm.eg.db)
library(SeuratData)
library(Seurat)
library(monocle)
library(monocle3)


load("data/pbmc.markers.rda")
load("data/diff_test_res.rda")
load("data/modulated_genes_ft.rda")

function(input, output, session) {
  # file upload limit
  options(shiny.maxRequestSize = 10240 * 1024^2)

  # ============================================================================
  # load user input data
  # ============================================================================
  mat <- reactive({
    if(input$example_data == 'NULL'){
      infile <- input$data$datapath

      if (is.null(infile)) return(NULL)

      read.csv(infile,header = T,sep = ',',check.names = F,row.names = 1)
    }else{
      read.csv("data/exps.csv",header = T,sep = ',',check.names = F,row.names = 1)
    }
  })

  # datatable show
  output$input <- DT::renderDT(options = list(scrollX = TRUE),{mat()})

  # ============================================================================
  # get clusters
  # ============================================================================
  output$get_cluster <- renderPlot({
    # check optimal cluster numbers
    if(!is.null(mat())){
      getClusters(obj = mat())
    }else{
      ggplot2::ggplot()
    }
  })

  # ============================================================================
  # load user single cell input data
  # ============================================================================
  scobj <- reactive({
    if(input$scexample_data == 'NULL'){
      infile <- input$sc_data$datapath

      if (is.null(infile)) return(NULL)

      readRDS(infile)
    }else{
      data("pbmc3k.final")
      UpdateSeuratObject(object = pbmc3k.final)
    }
  })

  # marker diff data
  diff <- reactive({
    if(input$scexample_data == 'NULL'){
      infile <- input$diff_data$datapath

      if (is.null(infile)) return(NULL)
      read.csv(infile,header = T,sep = ',',check.names = F)
    }else{
      pbmc.markers
    }
  })

  # prepare data
  prepare_data <- reactiveVal(NULL)

  observeEvent(input$parepare_sc_data, {
    if(!is.null(scobj())){
      # check diff data
      if(input$scexample_data == 'NULL'){
        diffData <- diff()
      }else{
        diffData <- pbmc.markers
      }

      st.data <- prepareDataFromscRNA(object = scobj(),
                                      diffData = diffData,
                                      showAverage = input$showAverage,
                                      cells = NULL,
                                      group.by = 'ident',
                                      assays = 'RNA',
                                      slot = 'data',
                                      scale.data = input$scale.data,
                                      cluster.order = NULL,
                                      keep.uniqGene = input$keep.uniqGene,
                                      sep = "_")

      prepare_data(st.data)
    }

  })

  # datatable show
  output$scprepareed_table <- DT::renderDT(options = list(scrollX = TRUE),{prepare_data()$wide.res})

  output$diff_table <- DT::renderDT(options = list(scrollX = TRUE),{diff()})

  # ============================================================================
  # load user monocle2 input data
  # ============================================================================
  m2obj <- reactive({
    if(input$m2example_data == 'NULL'){
      infile <- input$m2_data$datapath

      if(!is.null(infile)){
        readRDS(infile)
      }else{
        return(NULL)
      }

    }else{
      readRDS("data/HSMM.rds")
    }
  })

  # marker diff data
  m2diff <- reactiveVal(NULL)

  m2diff <- reactive({
    if(input$m2example_data == 'NULL'){
      infile <- input$m2diff_data$datapath

      if(is.null(infile)){
        return(NULL)
      }else{
        read.csv(infile,header = T,sep = ',',check.names = F,row.names = 1)
      }

    }else{
      diff_test_res
    }
  })


  # differentialGeneTest analysis
  m2difftest_data <- reactiveVal(NULL)

  observeEvent(input$run_m2_diff, {
    withProgress(message = 'Running differentialGeneTest analysis...', value = 0, {
      incProgress(0.1, detail = "Initializing...")

      if(!is.null(m2obj())){
        # check diff data
        if(is.null(m2diff())){

          # diff analysis
          expressed_genes <- row.names(subset(fData(m2obj()), num_cells_expressed >= input$num_cells_expressed))

          diff.data <- differentialGeneTest(cds = m2obj()[expressed_genes,],
                                            fullModelFormulaStr = as.character(input$fullModelFormulaStr),
                                            reducedModelFormulaStr = as.character(input$reducedModelFormulaStr),
                                            relative_expr = as.logical(input$relative_expr),
                                            cores = input$cores,
                                            verbose = FALSE)

          m2difftest_data(diff.data)

        }
      }

      incProgress(1, detail = "Done!")
    })

  })

  # check input data
  observe({
    if(!is.null(m2diff())){
      m2difftest_data(m2diff() %>% subset(pval < input$pval & status != "FAIL"))
    }else if(!is.null(m2difftest_data())){
      m2difftest_data(m2difftest_data %>% subset(pval < input$pval & status != "FAIL"))
    }
  })


  output$m2diff_table <- DT::renderDT(options = list(scrollX = TRUE),{m2difftest_data()})
  # ============================================================================
  # beam analysis
  beam_table <- reactiveVal(NULL)

  observeEvent(input$run_beam, {
    withProgress(message = 'Running BEAM analysis...', value = 0, {
      incProgress(0.1, detail = "Initializing...")

      if(!is.null(m2obj())){

        if(input$branch_labels == "NULL"){
          branch_labels <- NULL
        }else{
          branch_labels <- eval(parse(text = input$branch_labels))
        }

        bm <- BEAM(
          cds = m2obj(),
          fullModelFormulaStr = as.character(input$bfullModelFormulaStr),
          reducedModelFormulaStr = as.character(input$breducedModelFormulaStr),
          branch_states = NULL,
          branch_point = input$branch_point,
          relative_expr = as.logical(input$brelative_expr),
          branch_labels = branch_labels,
          verbose = FALSE,
          cores = input$bcores,
          progenitor_method = "duplicate"
        )

        beam_table(bm)
      }

      incProgress(1, detail = "Done!")
    })

  })

  output$beam_table <- DT::renderDT(options = list(scrollX = TRUE),{beam_table()})
  # ============================================================================
  # download table
  output$download_m2diff_table <- downloadHandler(
    filename = function() {
      paste('monocle2-diff-results','csv', sep = '.')
    },

    content = function(file) {
      write.csv(data.frame(m2difftest_data()),file,row.names = T)
    })

  output$download_beam_table <- downloadHandler(
    filename = function() {
      paste('monocle2-beam-results','csv', sep = '.')
    },

    content = function(file) {
      write.csv(data.frame(beam_table()),file,row.names = T)
    })
  # ============================================================================
  # monocle2 heatmap
  # ============================================================================
  m2htobj <- reactiveVal(NULL)

  m2ht <- eventReactive(input$m2plotst,{
    # check pseudotime heatmap type
    req(m2obj())
    req(m2difftest_data())

    withProgress(message = 'Generating plot...', value = 0, {
      incProgress(0.1, detail = "Initializing...")

      if(input$mpht_type == "pseudotime heatmap"){
        # parse parameters
        custom_args <- eval(parse(text = input$ph_params))

        cds_subset <- m2obj()[row.names(m2difftest_data()),]

        all_args <- c(list(cds_subset = cds_subset,
                           cores = 1,
                           show_rownames = as.logical(input$show_rownames),
                           return_heatmap = T),
                      custom_args)

        ht <- do.call(plot_pseudotime_heatmap2, all_args)

        # return list
        all_args2 <- c(list(cds_subset = cds_subset,
                            cores = 1,
                            show_rownames = as.logical(input$show_rownames)),
                       custom_args)

        htlst <- do.call(plot_pseudotime_heatmap2, all_args2)
        m2htobj(htlst)

        ht
      }else if(input$mpht_type == "genes branched heatmap"){
        # parse parameters
        custom_args <- eval(parse(text = input$pbh_params))

        cds_subset <- m2obj()[row.names(m2difftest_data()),]

        all_args <- c(list(cds_subset = cds_subset,
                           cores = 1,
                           show_rownames = as.logical(input$show_rownames),
                           return_heatmap = T),
                      custom_args)

        ht <- do.call(plot_genes_branched_heatmap2, all_args)

        # return list
        all_args2 <- c(list(cds_subset = cds_subset,
                            cores = 1,
                            show_rownames = as.logical(input$show_rownames)),
                       custom_args)

        htlst <- do.call(plot_genes_branched_heatmap2, all_args2)
        m2htobj(htlst)

        ht
      }else if(input$mpht_type == "multiple branches heatmap"){
        # parse parameters
        custom_args <- eval(parse(text = input$pmbh_params))

        cds_subset <- m2obj()[row.names(m2difftest_data()),]

        all_args <- c(list(cds_subset = cds_subset,
                           cores = 1,
                           show_rownames = as.logical(input$show_rownames),
                           return_heatmap = T),
                      custom_args)

        ht <- do.call(plot_multiple_branches_heatmap2, all_args)

        # return list
        all_args2 <- c(list(cds_subset = cds_subset,
                            cores = 1,
                            show_rownames = as.logical(input$show_rownames)),
                       custom_args)

        htlst <- do.call(plot_multiple_branches_heatmap2, all_args2)
        m2htobj(htlst)

        ht
      }
      incProgress(1, detail = "Done!")
    })

  })

  output$m2heatmap <- renderPlot({
    m2ht()
  })

  output$m2cl_table <- DT::renderDT(options = list(scrollX = TRUE),{m2htobj()})

  # download table
  output$download_m2cl_table <- downloadHandler(
    filename = function() {
      paste('monocle2-clustered-results','csv', sep = '.')
    },

    content = function(file) {
      write.csv(data.frame(m2htobj()$wide.res),file,row.names = T)
    })

  # download plot
  output$download_m2htplot <- downloadHandler(
    filename = function() {
      paste('monocle2_heatmap_plot', 'pdf',sep = '.')
    },

    content = function(file) {
      req(m2obj())
      req(m2difftest_data())

      pdf(file,width = input$htwidth, height = input$htheight)

      # plot
      if(input$mpht_type == "pseudotime heatmap"){
        # parse parameters
        custom_args <- eval(parse(text = input$ph_params))

        cds_subset <- m2obj()[row.names(m2difftest_data()),]

        all_args <- c(list(cds_subset = cds_subset,
                           cores = 1,
                           show_rownames = as.logical(input$show_rownames),
                           return_heatmap = T),
                      custom_args)

        ht <- do.call(plot_pseudotime_heatmap2, all_args)
      }else if(input$mpht_type == "genes branched heatmap"){
        # parse parameters
        custom_args <- eval(parse(text = input$pbh_params))

        cds_subset <- m2obj()[row.names(m2difftest_data()),]

        all_args <- c(list(cds_subset = cds_subset,
                           cores = 1,
                           show_rownames = as.logical(input$show_rownames),
                           return_heatmap = T),
                      custom_args)

        ht <- do.call(plot_genes_branched_heatmap2, all_args)
      }else if(input$mpht_type == "multiple branches heatmap"){
        # parse parameters
        custom_args <- eval(parse(text = input$pmbh_params))

        cds_subset <- m2obj()[row.names(m2difftest_data()),]

        all_args <- c(list(cds_subset = cds_subset,
                           cores = 1,
                           show_rownames = as.logical(input$show_rownames),
                           return_heatmap = T),
                      custom_args)

        ht <- do.call(plot_multiple_branches_heatmap2, all_args)
      }

      print(ht)
      dev.off()
    }
  )
  # ============================================================================
  # load user monocle3 input data
  # ============================================================================
  m3obj <- reactive({
    if(input$m3example_data == 'NULL'){
      infile <- input$m3_data$datapath

      if(!is.null(infile)){
        readRDS(infile)
      }else{
        return(NULL)
      }

    }else{
      readRDS("data/cds.rds")
    }
  })

  # marker diff data
  m3diff <- reactiveVal(NULL)

  m3diff <- reactive({
    if(input$m3example_data == 'NULL'){
      infile <- input$m3diff_data$datapath

      if(is.null(infile)){
        return(NULL)
      }else{
        read.csv(infile,header = T,sep = ',',check.names = F,row.names = 1)
      }

    }else{
      modulated_genes_ft
    }
  })


  # differentialGeneTest analysis
  m3difftest_data <- reactiveVal(NULL)

  observeEvent(input$run_graph_test, {
    withProgress(message = 'Running graph_test analysis...', value = 0, {
      incProgress(0.1, detail = "Initializing...")

      if(!is.null(m3obj())){
        # check diff data
        if(is.null(m3diff())){

          # diff analysis
          diff.data <- graph_test(
            cds = m3obj(),
            neighbor_graph = as.character(input$neighbor_graph),
            reduction_method = as.character(input$reduction_method),
            k = input$k,
            method = c("Moran_I"),
            alternative = as.character(input$alternative),
            expression_family = "quasipoisson",
            cores = input$m3cores,
            verbose = FALSE,
            nn_control = list()
          )

          m3difftest_data(diff.data)

        }
      }

      incProgress(1, detail = "Done!")
    })

  })

  # check input data
  observe({
    if(!is.null(m3diff())){
      m3difftest_data(m3diff() %>% subset(p_value < input$m3pval))
    }else if(is.null(m3diff()) & !is.null(m3difftest_data())){
      m3difftest_data(m3difftest_data() %>% subset(p_value < input$m3pval))
    }
  })


  output$m3diff_table <- DT::renderDT(options = list(scrollX = TRUE),{m3difftest_data()})

  # ============================================================================
  # download table
  output$download_m3diff_table <- downloadHandler(
    filename = function() {
      paste('monocle3-diff-results','csv', sep = '.')
    },

    content = function(file) {
      write.csv(data.frame(m3difftest_data()),file,row.names = T)
    })

  # ============================================================================
  # extract expression matrix

  m3_mat <- reactive({
    if(!is.null(m3obj())){
      mat <- pre_pseudotime_matrix(cds_obj = m3obj(),
                                   gene_list = rownames(m3difftest_data()))
      mat
    }else{
      return(NULL)
    }
  })

  output$m3mat_table <- DT::renderDT(options = list(scrollX = TRUE),{m3_mat()})

  # download table
  output$download_m3mat_table <- downloadHandler(
    filename = function() {
      paste('monocle3-pseudotime-matrix','csv', sep = '.')
    },

    content = function(file) {
      write.csv(data.frame(m3_mat()),file,row.names = T)
    })

  # ============================================================================
  # check input matrix data
  # ============================================================================
  input_mat <- reactiveVal(NULL)

  observe({
    # check input matrix
    if(!is.null(m3_mat()) & is.null(mat())){
      input_mat(m3_mat())
    }else if(is.null(m3_mat()) & !is.null(mat())){
      input_mat(mat())
    }
  })

  # ============================================================================
  # gene clustering
  # ============================================================================
  cluster_result <- reactiveVal(NULL)

  observeEvent(input$run_clustering, {
    req(input$cl_num, input$cl_method)

    cluster_num <- input$cl_num
    cluster_method <- input$cl_method

    result <- clusterData(obj = input_mat(),
                          cluster.method = cluster_method,
                          cluster.num = cluster_num)

    cluster_result(result)
  })


  output$clustered_table <- DT::renderDT(options = list(scrollX = TRUE),{cluster_result()$wide.res})

  # download table
  output$download_clustered_table <- downloadHandler(
    filename = function() {
      paste('clustered-results','csv', sep = '.')
    },

    content = function(file) {
      write.csv(cluster_result()$wide.res,file,row.names = T)
    })

  # ============================================================================
  # check input data for downstream analysis
  # ============================================================================

  observe({
    if(is.null(cluster_result()) & !is.null(prepare_data()) & is.null(m2htobj())){
      cluster_result(prepare_data())
    }else if(is.null(cluster_result()) & is.null(prepare_data()) & !is.null(m2htobj())){
      cluster_result(m2htobj())
    }
  })

  # ============================================================================
  # enrichment analysis for clusters
  # ============================================================================

  # load user defined enrichment data
  user_enrich_data <- reactive({
    infile <- input$user_enrich$datapath

    if (is.null(infile)){
      read.csv(infile,header = T,sep = ',',check.names = F)
    }else{
      return(NULL)
    }
  })

  # term2gene data
  term2gene_data <- reactive({
    infile <- input$TERM2GENE$datapath

    if (is.null(infile)){
      read.csv(infile,header = T,sep = ',',check.names = F)
    }else{
      return(NULL)
    }
  })

  if(is.null(term2gene_data)){
    TERM2GENE <- NULL
  }else{
    TERM2GENE <- term2gene_data
  }


  # ============================================================================
  # enrichment analysis
  enrichment_result <- reactiveVal(NULL)
  kegg_result <- reactiveVal(NULL)

  observeEvent(input$run_enrichment, {
    withProgress(message = 'Running enrichment analysis...', value = 0, {
      incProgress(0.1, detail = "Initializing...")

      if (input$OrgDb == "human") {
        OrgDb <- org.Hs.eg.db
      } else if (input$OrgDb == "mouse") {
        OrgDb <- org.Mm.eg.db
      }


      # enrichment analysis
      if(length(input$type) == 1){
        enrich <- enrichCluster(object = cluster_result(),
                                type = input$type,
                                TERM2GENE = TERM2GENE,
                                TERM2NAME = NULL,
                                OrgDb = OrgDb,
                                id.trans = input$id.trans,
                                fromType = input$fromType,
                                toType = input$toType,
                                readable = input$readable,
                                organism = input$organism,
                                pvalueCutoff  = input$pvalueCutoff,
                                topn = input$topn,
                                seed = 5201314,
                                add.gene = input$add.gene)

        enrichment_result(enrich)
      }else if(length(input$type) == 2){
        go <- enrichCluster(object = cluster_result(),
                            type = input$type[1],
                            TERM2GENE = TERM2GENE,
                            TERM2NAME = NULL,
                            OrgDb = OrgDb,
                            id.trans = input$id.trans,
                            fromType = input$fromType,
                            toType = input$toType,
                            readable = input$readable,
                            organism = input$organism,
                            pvalueCutoff  = input$pvalueCutoff,
                            topn = input$topn,
                            seed = 5201314,
                            add.gene = input$add.gene)

        enrichment_result(go)

        kegg <- enrichCluster(object = cluster_result(),
                              type = input$type[2],
                              TERM2GENE = TERM2GENE,
                              TERM2NAME = NULL,
                              OrgDb = OrgDb,
                              id.trans = input$id.trans,
                              fromType = input$fromType,
                              toType = input$toType,
                              readable = input$readable,
                              organism = input$organism,
                              pvalueCutoff  = input$pvalueCutoff,
                              topn = input$topn,
                              seed = 5201314,
                              add.gene = input$add.gene)

        kegg_result(kegg)
      }

      incProgress(1, detail = "Done!")
    })
  })

  # ============================================================================
  # show enriched table
  output$enriched_table <- DT::renderDT(options = list(scrollX = TRUE),{enrichment_result()})

  # download table
  output$download_enrich1 <- downloadHandler(
    filename = function() {
      paste('enrich1-results','csv', sep = '.')
    },

    content = function(file) {
      write.csv(data.frame(enrichment_result()),file,row.names = T)
    })

  output$enriched_table2 <- DT::renderDT(options = list(scrollX = TRUE),{kegg_result()})

  # download table
  output$download_enrich2 <- downloadHandler(
    filename = function() {
      paste('enrich2-results','csv', sep = '.')
    },

    content = function(file) {
      write.csv(data.frame(kegg_result()),file,row.names = T)
    })
  # ============================================================================
  # visualization
  # ============================================================================
  cmbht <- eventReactive(input$plotst,{
    if(!is.null(cluster_result())){
      # check enrichment data
      if(!is.null(enrichment_result)){
        annoTerm.data <- enrichment_result()
      }else{
        annoTerm.data <- NULL
      }

      if(!is.null(kegg_result)){
        annoKegg.data <- kegg_result()
      }else{
        annoKegg.data <- NULL
      }

      # parse parameters
      custom_args <- eval(parse(text = input$custom_params))

      all_args <- c(list(object = cluster_result(),
                         ht.col.list = list(col_range = c(-2, 0, 2),
                                            col_color = c(input$low, input$mid, input$high)),
                         plot.type = input$plot_type,
                         ms.col = c(input$line_low, input$line_mid, input$line_high),
                         mline.col = input$mline.col,
                         ncol = input$ncol,
                         annoTerm.data = annoTerm.data,
                         annoTerm.mside = "right",
                         add.bar = input$add.bar,
                         # KEGG term annotation
                         annoKegg.data = annoKegg.data,
                         annoKegg.mside = "right"),
                    custom_args)

      ht <- do.call(visCluster, all_args)

      ht
    }else{
      ggplot()
    }

  })

  output$cmb_plot <- renderPlot({
    cmbht()
  })

  # download plot
  output$download_plot <- downloadHandler(
    filename = function() {
      paste('ClusterGvis_plot', 'pdf',sep = '.')
    },

    content = function(file) {
      if(!is.null(cluster_result())){
        pdf(file,width = input$pwidth, height = input$pheight)

        # check enrichment data
        if(!is.null(enrichment_result)){
          annoTerm.data <- enrichment_result()
        }else{
          annoTerm.data <- NULL
        }

        if(!is.null(kegg_result)){
          annoKegg.data <- kegg_result()
        }else{
          annoKegg.data <- NULL
        }

        # plot
        # parse parameters
        custom_args <- eval(parse(text = input$custom_params))

        all_args <- c(list(object = cluster_result(),
                           ht.col.list = list(col_range = c(-2, 0, 2),
                                              col_color = c(input$low, input$mid, input$high)),
                           plot.type = input$plot_type,
                           ms.col = c(input$line_low, input$line_mid, input$line_high),
                           mline.col = input$mline.col,
                           ncol = input$ncol,
                           annoTerm.data = annoTerm.data,
                           annoTerm.mside = "right",
                           add.bar = input$add.bar,
                           # KEGG term annotation
                           annoKegg.data = annoKegg.data,
                           annoKegg.mside = "right"),
                      custom_args)

        ht <- do.call(visCluster, all_args)

        print(ht)
        dev.off()
      }

    }
  )

  # ============================================================================
  # viscluster doc
  # ============================================================================
  # save html var
  help_content <- reactiveVal("")

  # click
  observeEvent(input$load_doc, {
    help_path <- utils::help("visCluster")
    if (!is.null(help_path)) {
      tryCatch({
        rd_content <- capture.output(
          tools::Rd2HTML(utils:::.getHelpFile(as.character(help_path)))
        )
        clean_content <- paste(rd_content, collapse = "\n")
        help_content(clean_content)
      }, error = function(e) {
        help_content("<p>fail to extract documentation.</p>")
      })
    } else {
      help_content("<p>visCluster not found！</p>")
    }
  })


  output$prams_doc <- renderUI({
    req(help_content())
    HTML(help_content())
  })

}


