library(shiny)
library(DT)
library(shinyBS)
library(ClusterGVis)
library(shinydashboard)
library(shinydashboardPlus)
library(colourpicker)
library(shinycssloaders)

# Define UI for application that draws a histogram
dashboardPage(
  dashboardHeader(title = "ClusterGvis"),

  # ============================================================================
  # dashboardSidebar
  # ============================================================================
  dashboardSidebar(
    sidebarMenu(
      menuItem("Home", tabName = "Home", icon = icon("house")),
      menuItem("Workflow", tabName = "Workflow", icon = icon("desktop"))
    )
  ),

  # ============================================================================
  # dashboardBody
  # ============================================================================
  dashboardBody(
    tabItems(
      # First tab content
      tabItem(tabName = "Home",
              tags$h1("Welcome to the ClusterGvis Home Page!", style = "color: #003366; text-align: center; font-weight: bold;"),
              tags$h5(
                HTML("<p><strong>ClusterGvis</strong> is an integrated R package designed to streamline
                the analysis and visualization of gene expression data. It combines clustering, functional
                enrichment, and publication-quality visualization into a single, user-friendly workflow.
                Key features include:</p>
              <ul><li><strong>Clustering Analysis</strong>: Supports fuzzy c-means and k-means clustering for time-series or bulk transcriptomics data.</li>
                  <li><strong>Functional Enrichment</strong>: Performs pathway enrichment (GO, KEGG, custom databases) for identified gene clusters.</li>
                  <li><strong>Visualization</strong>: Generates high-quality plots, including heatmaps, trend lines, and enrichment barplots, using <strong>ComplexHeatmap</strong> as the core engine.</li>
                  <li><strong>Single-Cell Integration</strong>: Seamlessly integrates with <strong>Seurat</strong> and <strong>monocle</strong> for single-cell RNA-Seq data analysis.</li>
              </ul><p><strong>ClusterGvis</strong> simplifies complex workflows, reduces coding requirements,
                     and enables researchers to generate comprehensive, publication-ready figures
                     with minimal effort. It is particularly valuable for users with limited programming
                     experience, offering a powerful yet accessible tool for transcriptomics data exploration.</p>"),
                style = "font-family: Arial; line-height: 1.5;"
              ),

              tags$div(style = "height: 5px;"),
              HTML("<b>Please send an email to 3219030654@stu.cpu.edu.cn or leave an issue on
                   https://github.com/junjunlab/ClusterGVis/issues if you have any questions or advice.</b>"),

              HTML(
                paste0(
                  "<b>More details about ClusterGVis can be reached at ",
                  "<a href='https://github.com/junjunlab/ClusterGVis' target='_blank'>GitHub Repository</a>",
                  " and ",
                  "<a href='https://junjunlab.github.io/ClusterGvis-manual/' target='_blank'>User Manual</a>",
                  ".</b>"
                )
              ),


              tags$div(style = "height: 10px;"),

              # picture
              tags$img(
                src = "clusterGvis_graphic_abstract.jpg",
                style = "width: 100%; height: auto;"
              )
      ),

      tabItem(tabName = "Workflow",
              fluidPage(
                tags$style(HTML(".nav-tabs { border-bottom: 2px solid #333; }.tab-content { margin-top: 20px; }")),
                tabsetPanel(
                  tabPanel(tags$b("1.input table"),
                           tabsetPanel(
                             tabPanel(HTML("<b><i>matrix input</i></b>"),
                                      fluidRow(
                                        column(4,
                                               box(title = "params",status = "orange",width = 12,solidHeader = T,

                                                   selectInput('example_data','Example data',
                                                               choices = c('exps','NULL'),selected = 'NULL'),
                                                   fileInput("data","Input data (gene expression matrix with csv format)")
                                               )
                                        ),
                                        column(8,
                                               box(title = "table",status = "orange",width = 12,solidHeader = T,
                                                   DT::DTOutput("input")
                                               ))
                                      )
                             ),
                             tabPanel(HTML("<b><i>seurat input</i></b>"),
                                      fluidRow(
                                        column(4,
                                               box(title = "single cell params",status = "success",width = 12,solidHeader = T,
                                                   selectInput('scexample_data','Example data',
                                                               choices = c('pbmc3k.final','NULL'),selected = 'NULL'),
                                                   div(style = "margin-bottom: -25px;",
                                                       fileInput("sc_data","Input data (seurat object with rds format)")
                                                   ),
                                                   div(style = "margin-bottom: -25px;",
                                                       fileInput("diff_data","marker genes diff data (csv format)")
                                                   ),
                                                   fluidRow(
                                                     column(6,
                                                            selectInput("showAverage","showAverage",
                                                                        choices = c("TRUE","FALSE"),
                                                                        selected = "TRUE")
                                                     ),
                                                     column(6,
                                                            selectInput("scale.data","scale.data",
                                                                        choices = c("TRUE","FALSE"),
                                                                        selected = "TRUE")
                                                     )
                                                   ),
                                                   selectInput("keep.uniqGene","keep.uniqGene",
                                                               choices = c("TRUE","FALSE"),
                                                               selected = "TRUE"),

                                                   actionButton("parepare_sc_data", HTML('<span style="color: red; font-weight: bold;">Preparing data</span>'))
                                               )
                                        ),
                                        column(8,
                                               box(title = "table",status = "success",width = 12,solidHeader = T,
                                                   tabsetPanel(
                                                     tabPanel("single cell data",DT::DTOutput("scprepareed_table")),
                                                     tabPanel("diff marker genes table",DT::DTOutput("diff_table"))
                                                   )
                                               )
                                        )
                                      )
                             ),
                             tabPanel(HTML("<b><i>monocle2 input</i></b>"),
                                      fluidRow(
                                        column(4,
                                               box(title = "monocle2 params",status = "primary",width = 12,solidHeader = T,

                                                   selectInput('m2example_data','Example data',
                                                               choices = c('HSMM','NULL'),selected = 'NULL'),
                                                   div(style = "margin-bottom: -10px;",
                                                       fileInput("m2_data","Input data (monocle2 CellDataSet object with rds format)")
                                                   ),

                                                   bsCollapsePanel(title = 'monocle differential gene test',value = '',style = 'primary',
                                                                   div(style = "margin-bottom: -20px;",
                                                                       fileInput("m2diff_data","monocle2 differentialGeneTest output with csv format")
                                                                   ),

                                                                   div(style = "border-top: 2px solid #333; margin: 20px 0;"),

                                                                   numericInput("num_cells_expressed","num_cells_expressed",value = 10,min = 1),
                                                                   textInput("fullModelFormulaStr","fullModelFormulaStr",value = "~sm.ns(Pseudotime, df=3)"),
                                                                   textInput("reducedModelFormulaStr","reducedModelFormulaStr",value = "~1"),
                                                                   fluidRow(
                                                                     column(6,
                                                                            selectInput("relative_expr","relative_expr",
                                                                                        choices = c("TRUE","FALSE"),
                                                                                        selected = "TRUE")
                                                                     ),
                                                                     column(6,
                                                                            numericInput("cores","cores",value = 1,min = 1)
                                                                     )
                                                                   ),
                                                                   actionButton("run_m2_diff", HTML('<span style="color: red; font-weight: bold;">Runing differentialGeneTest</span>'))
                                                   ),
                                                   selectizeInput("mpht_type","heatmap type",
                                                                  choices = c("pseudotime heatmap","genes branched heatmap","multiple branches heatmap"),
                                                                  selected = "pseudotime heatmap"),
                                                   numericInput("pval","pvalue threshold",value = 0.0001,min = 0,max = 1),
                                                   selectInput("show_rownames","show_rownames",
                                                               choices = c("TRUE","FALSE"),
                                                               selected = "TRUE"),
                                                   bsCollapsePanel(title = 'plot pseudotime heatmap',value = '',style = 'primary',
                                                                   textAreaInput("ph_params", "custom parameters:",
                                                                                 value = 'list(num_clusters = 4)',
                                                                                 resize = "vertical")
                                                   ),
                                                   bsCollapsePanel(title = 'plot genes branched heatmap',value = '',style = 'primary',
                                                                   textAreaInput("pbh_params", "custom parameters:",
                                                                                 value = 'list(num_clusters = 4,branch_point = 1)',
                                                                                 resize = "vertical")
                                                   ),
                                                   bsCollapsePanel(title = 'BEAM analysis',value = '',style = 'primary',
                                                                   textInput("bfullModelFormulaStr","fullModelFormulaStr",value = "~sm.ns(Pseudotime, df = 3)*Branch"),
                                                                   textInput("breducedModelFormulaStr","reducedModelFormulaStr",value = "~sm.ns(Pseudotime, df = 3)"),
                                                                   numericInput("branch_point","branch_point",value = 1,min = 1),
                                                                   textInput("branch_labels","branch_labels",value = "NULL"),
                                                                   fluidRow(
                                                                     column(6,
                                                                            selectInput("brelative_expr","relative_expr",
                                                                                        choices = c("TRUE","FALSE"),
                                                                                        selected = "TRUE")
                                                                     ),
                                                                     column(6,
                                                                            numericInput("bcores","cores",value = 1,min = 1)
                                                                     )
                                                                   ),
                                                                   textAreaInput("pmbh_params", "custom parameters for plot_multiple_branches_heatmap2:",
                                                                                 value = 'list(num_clusters = 4,branches = c(1,3,4,5))',
                                                                                 resize = "vertical"),
                                                                   actionButton("run_beam", HTML('<span style="color: red; font-weight: bold;">Runing BEAM</span>'))
                                                   ),
                                                   div(style = "margin-bottom: 15px;",
                                                       actionButton("m2plotst",HTML('<span style="color: red; font-weight: bold;">Submit Plot</span>'),width = "100%"),

                                                   ),
                                                   # ==============================download plot
                                                   bsCollapsePanel(title = 'plot download',value = '',style = 'primary',
                                                                   fluidRow(
                                                                     column(6,numericInput("htwidth",label = "plot width",min = 0,max = 50,value = 8)),
                                                                     column(6,numericInput("htheight",label = "plot height",min = 0,max = 50,value = 10))
                                                                   ),
                                                                   downloadButton('download_m2htplot','download plot')
                                                   )

                                               )
                                        ),
                                        column(8,
                                               box(title = "table",status = "primary",width = 12,solidHeader = T,
                                                   tabsetPanel(
                                                     tabPanel("heatmap",
                                                              plotOutput("m2heatmap")
                                                     ),
                                                     tabPanel("diff genes table",
                                                              DT::DTOutput("m2diff_table") %>% withSpinner(type = 8),
                                                              downloadButton('download_m2diff_table','download table')
                                                     ),
                                                     tabPanel("BEAM genes table",
                                                              DT::DTOutput("beam_table") %>% withSpinner(type = 8),
                                                              downloadButton('download_beam_table','download table')
                                                     ),
                                                     tabPanel("clustered table",
                                                              DT::DTOutput("m2cl_table") %>% withSpinner(type = 8),
                                                              downloadButton('download_m2cl_table','download table')
                                                     )
                                                   )

                                               ))
                                      )
                             ),
                             tabPanel(HTML("<b><i>monocle3 input</i></b>"),
                                      fluidRow(
                                        column(4,
                                               box(title = "params",status = "navy",width = 12,solidHeader = T,
                                                   selectInput('m3example_data','Example data',
                                                               choices = c('worm embryo','NULL'),selected = 'NULL'),
                                                   div(style = "margin-bottom: -10px;",
                                                       fileInput("m3_data","Input data (monocle3 CellDataSet object with rds format)")
                                                   ),
                                                   numericInput("m3pval","pvalue threshold",value = 0.0001,min = 0,max = 1),

                                                   bsCollapsePanel(title = 'graph test analysis',value = '',style = 'primary',
                                                                   div(style = "margin-bottom: -20px;",
                                                                       fileInput("m3diff_data","monocle3 graph_test output with csv format")
                                                                   ),

                                                                   div(style = "border-top: 2px solid #333; margin: 20px 0;"),

                                                                   selectInput("neighbor_graph","neighbor_graph",choices = c("knn", "principal_graph"),selected = "principal_graph"),
                                                                   textInput("reduction_method","reduction_method",value = "UMAP"),

                                                                   fluidRow(
                                                                     column(6,
                                                                            numericInput("k","k",value = 25,min = 1,max = 100)
                                                                     ),
                                                                     column(6,
                                                                            numericInput("m3cores","cores",value = 1,min = 1,max = 100)
                                                                     )
                                                                   ),
                                                                   selectInput("alternative","alternative",choices = c("greater","two.sided"),selected = "greater"),

                                                                   actionButton("run_graph_test", HTML('<span style="color: red; font-weight: bold;">Runing graph_test analysis</span>'))

                                                   )
                                               )
                                        ),
                                        column(8,
                                               box(title = "monocle3 diff table",status = "navy",width = 12,solidHeader = T,

                                                   tabsetPanel(
                                                     tabPanel("diff genes table",
                                                              DT::DTOutput("m3diff_table") %>% withSpinner(type = 8),
                                                              downloadButton('download_m3diff_table','download table')
                                                     ),
                                                     tabPanel("matrix table",
                                                              DT::DTOutput("m3mat_table") %>% withSpinner(type = 8),
                                                              downloadButton('download_m3mat_table','download table')
                                                     )
                                                   )
                                               ))
                                      )
                             )
                           )

                  ),
                  tabPanel(tags$b("2.cluster numbers"),
                           column(4,
                                  box(title = "params",status = "info",width = 12,solidHeader = T)
                           ),
                           column(8,
                                  box(title = "plot",status = "purple",width = 12,solidHeader = T,
                                      plotOutput("get_cluster") %>% withSpinner(type = 8))
                           )
                  ),
                  tabPanel(tags$b("3.gene clustering"),
                           column(4,
                                  box(title = "params",status = "info",width = 12,solidHeader = T,
                                      numericInput("cl_num","cluster numbers",value = 1,min = 1,max = 20),
                                      selectInput("cl_method","cluster methods",
                                                  choices = c("mfuzz","TCseq","kmeans"),
                                                  selected = "kmeans"),
                                      actionButton("run_clustering", HTML('<span style="color: red; font-weight: bold;">Clustering Started</span>'))
                                  )
                           ),
                           column(8,
                                  box(title = "clustered table",status = "purple",width = 12,solidHeader = T,
                                      DT::DTOutput("clustered_table"),
                                      downloadButton('download_clustered_table','download table'))
                           )
                  ),
                  tabPanel(tags$b("4.pathway enrichment"),
                           column(4,
                                  box(title = "params",status = "primary",width = 12,solidHeader = T,
                                      div(style = "margin-bottom: -25px;",
                                          fileInput("user_enrich", "Use own enrichment data")
                                      ),
                                      selectInput("type","type (one or two options surpported)",
                                                  choices = c("BP","MF","CC","KEGG","ownSet"),
                                                  selected = "BP",multiple = T),
                                      div(style = "margin-bottom: -25px;",
                                          fileInput("TERM2GENE","gene-term mapping data")
                                      ),
                                      fluidRow(
                                        column(6,
                                               selectInput("OrgDb","OrgDb",
                                                           choices = c("human","mouse"),
                                                           selected = "human")
                                        ),
                                        column(6,
                                               selectInput("id.trans","id.trans",
                                                           choices = c("TRUE","FALSE"),
                                                           selected = "TRUE")
                                        )
                                      ),
                                      fluidRow(
                                        column(6,
                                               textInput("fromType","fromType",value = "SYMBOL")
                                        ),
                                        column(6,
                                               textInput("toType","toType",value = "ENTREZID")
                                        )
                                      ),
                                      fluidRow(
                                        column(6,
                                               selectInput("readable","readable",
                                                           choices = c("TRUE","FALSE"),
                                                           selected = "TRUE")
                                        ),
                                        column(6,
                                               textInput("organism","organism",value = "hsa")
                                        )
                                      ),
                                      sliderInput("pvalueCutoff","pvalueCutoff",min = 0,max = 1,step = 0.01,value = 0.05),
                                      fluidRow(
                                        column(6,
                                               numericInput("topn","topn",min = 0,max = 50,value = 5)
                                        ),
                                        column(6,
                                               selectInput("add.gene","add.gene",
                                                           choices = c("TRUE","FALSE"),
                                                           selected = "FALSE")
                                        )
                                      ),


                                      actionButton("run_enrichment", HTML('<span style="color: red; font-weight: bold;">Run Enrichment for Clusters</span>'))
                                  )

                           ),
                           column(8,
                                  fluidPage(
                                    box(title = "enrichment table",status = "maroon",width = 12,solidHeader = T,
                                        DT::DTOutput("enriched_table"),
                                        downloadButton('download_enrich1','download table')
                                    )
                                  ),
                                  fluidPage(
                                    box(title = "enrichment table2",status = "maroon",width = 12,solidHeader = T,
                                        DT::DTOutput("enriched_table2"),
                                        downloadButton('download_enrich2','download table')
                                    )
                                  )
                           )
                  ),
                  tabPanel(tags$b("5.visualization"),
                           column(4,
                                  box(title = "params",status = "primary",width = 12,solidHeader = T,
                                      fluidRow(
                                        column(12,tags$p("Select colors for heatmap:",
                                                         style = "font-weight: bold; text-decoration: underline;")),
                                        column(4,colourInput("low","low","#08519C")),
                                        column(4,colourInput("mid","mid","white")),
                                        column(4,colourInput("high","high","#A50F15"))
                                      ),
                                      fluidRow(
                                        column(12,tags$p("Select colors for membership lines:",
                                                         style = "font-weight: bold; text-decoration: underline;")),
                                        column(4,colourInput("line_low","low","#0099CC")),
                                        column(4,colourInput("line_mid","mid","grey90")),
                                        column(4,colourInput("line_high","high","#CC3333"))
                                      ),
                                      fluidRow(
                                        column(6,
                                               selectInput("plot_type","plot type",choices = c("line","heatmap","both"),selected = "line")
                                        ),
                                        column(6,selectInput("add.bar","add.bar",
                                                             choices = c("TRUE","FALSE"),
                                                             selected = "FALSE"))
                                      ),
                                      textAreaInput("marker_gene","gene symbols to be marked ",
                                                    value = "supply gene symbols to be marked on heatmap with ‘,' seprated.",
                                                    resize = "vertical"),
                                      fluidRow(
                                        column(6,
                                               selectInput("line.side","line.side",choices = c("left","right"),selected = "right")
                                        ),
                                        column(6,selectInput("markGenes.side","markGenes.side",choices = c("left","right"),selected = "right"))
                                      ),
                                      fluidRow(
                                        column(6,
                                               numericInput("ncol","ncol",value = 4,min = 1,max = 30)
                                        ),
                                        column(6,colourInput("mline.col","mline.col","#CC3333"))
                                      ),



                                      # ==============================
                                      bsCollapsePanel(title = 'custom parameter setting',value = '',style = 'primary',
                                                      textAreaInput("custom_params", "custom parameters:",
                                                                    value = 'list(column_names_rot = 45)',
                                                                    resize = "vertical")

                                      ),
                                      div(style = "margin-bottom: 15px;",
                                          actionButton("plotst",HTML('<span style="color: red; font-weight: bold;">Submit Plot</span>'),width = "100%"),

                                      ),


                                      # ==============================download plot
                                      bsCollapsePanel(title = 'plot download',value = '',style = 'primary',
                                                      fluidRow(
                                                        column(6,numericInput("pwidth",label = "plot width",min = 0,max = 50,value = 8)),
                                                        column(6,numericInput("pheight",label = "plot height",min = 0,max = 50,value = 10))
                                                      ),
                                                      downloadButton('download_plot','download plot')
                                      )

                                  )
                           ),
                           column(8,
                                  fluidPage(
                                    box(title = "plot",status = "danger",width = 12,solidHeader = T,
                                        plotOutput("cmb_plot") %>% withSpinner(type = 8)

                                    )
                                  ),
                                  fluidPage(
                                    box(title = "documentation",status = "danger",width = 12,solidHeader = T,
                                        actionButton("load_doc", "load viscluster documentation"),
                                        htmlOutput("prams_doc", class = "doc-box")
                                    ),
                                    tags$style(HTML(".doc-box {
                                                      max-height: 500px;
                                                      overflow-y: auto;
                                                      padding: 10px;
                                                      border: 1px solid #ddd;
                                                      background-color: #fff;
                                                      font-family: Arial, sans-serif;
                                                      font-size: 14px;
                                                      line-height: 1.5;
                                                      color: #333;
                                                    }
                                                  "))
                                  ),

                           )
                  )
                )
              )
      )
    )
  )
)
