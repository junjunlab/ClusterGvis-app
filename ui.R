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
                  tabPanel(tags$b("1.Data input"),
                           tabsetPanel(
                             tabPanel(HTML("<b><i>Matrix input</i></b>"),
                                      fluidRow(
                                        column(4,
                                               box(title = "Params",status = "orange",width = 12,solidHeader = T,

                                                   selectInput('example_data','Example data',
                                                               choices = c('exps','NULL'),selected = 'NULL'),
                                                   fileInput("data","Input data (gene expression matrix with csv format)")
                                               )
                                        ),
                                        column(8,
                                               box(title = "Table",status = "orange",width = 12,solidHeader = T,
                                                   DT::DTOutput("input")
                                               ))
                                      )
                             ),
                             tabPanel(HTML("<b><i>Seurat input</i></b>"),
                                      fluidRow(
                                        column(4,
                                               box(title = "Single cell params",status = "success",width = 12,solidHeader = T,
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
                                                            selectInput("showAverage","ShowAverage",
                                                                        choices = c("TRUE","FALSE"),
                                                                        selected = "TRUE")
                                                     ),
                                                     column(6,
                                                            selectInput("scale.data","Scale.data",
                                                                        choices = c("TRUE","FALSE"),
                                                                        selected = "TRUE")
                                                     )
                                                   ),
                                                   selectInput("keep.uniqGene","Keep.uniqGene",
                                                               choices = c("TRUE","FALSE"),
                                                               selected = "TRUE"),

                                                   actionButton("parepare_sc_data", HTML('<span style="color: red; font-weight: bold;">Preparing data</span>'))
                                               )
                                        ),
                                        column(8,
                                               box(title = "Table",status = "success",width = 12,solidHeader = T,
                                                   tabsetPanel(
                                                     tabPanel("Single cell data",DT::DTOutput("scprepareed_table")),
                                                     tabPanel("Diff marker genes table",DT::DTOutput("diff_table"))
                                                   )
                                               )
                                        )
                                      )
                             ),
                             tabPanel(HTML("<b><i>Monocle2 input</i></b>"),
                                      fluidRow(
                                        column(4,
                                               box(title = "Monocle2 params",status = "primary",width = 12,solidHeader = T,

                                                   selectInput('m2example_data','Example data',
                                                               choices = c('HSMM','NULL'),selected = 'NULL'),
                                                   div(style = "margin-bottom: -10px;",
                                                       fileInput("m2_data","Input data (monocle2 CellDataSet object with rds format)")
                                                   ),

                                                   bsCollapsePanel(title = 'Monocle differential gene test',value = '',style = 'primary',
                                                                   div(style = "margin-bottom: -20px;",
                                                                       fileInput("m2diff_data","Monocle2 differentialGeneTest output with csv format")
                                                                   ),

                                                                   div(style = "border-top: 2px solid #333; margin: 20px 0;"),

                                                                   numericInput("num_cells_expressed","Num_cells_expressed",value = 10,min = 1),
                                                                   textInput("fullModelFormulaStr","FullModelFormulaStr",value = "~sm.ns(Pseudotime, df=3)"),
                                                                   textInput("reducedModelFormulaStr","ReducedModelFormulaStr",value = "~1"),
                                                                   fluidRow(
                                                                     column(6,
                                                                            selectInput("relative_expr","Relative_expr",
                                                                                        choices = c("TRUE","FALSE"),
                                                                                        selected = "TRUE")
                                                                     ),
                                                                     column(6,
                                                                            numericInput("cores","Cores",value = 1,min = 1)
                                                                     )
                                                                   ),
                                                                   actionButton("run_m2_diff", HTML('<span style="color: red; font-weight: bold;">Runing differentialGeneTest</span>'))
                                                   ),
                                                   selectizeInput("mpht_type","Heatmap type",
                                                                  choices = c("pseudotime heatmap","genes branched heatmap","multiple branches heatmap"),
                                                                  selected = "pseudotime heatmap"),
                                                   numericInput("pval","P value threshold",value = 0.0001,min = 0,max = 1),
                                                   selectInput("show_rownames","Show rownames",
                                                               choices = c("TRUE","FALSE"),
                                                               selected = "TRUE"),
                                                   bsCollapsePanel(title = 'Plot pseudotime heatmap',value = '',style = 'primary',
                                                                   textAreaInput("ph_params", "Custom parameters:",
                                                                                 value = 'list(num_clusters = 4)',
                                                                                 resize = "vertical")
                                                   ),
                                                   bsCollapsePanel(title = 'Plot genes branched heatmap',value = '',style = 'primary',
                                                                   textAreaInput("pbh_params", "Custom parameters:",
                                                                                 value = 'list(num_clusters = 4,branch_point = 1)',
                                                                                 resize = "vertical")
                                                   ),
                                                   bsCollapsePanel(title = 'BEAM analysis',value = '',style = 'primary',
                                                                   textInput("bfullModelFormulaStr","FullModelFormulaStr",value = "~sm.ns(Pseudotime, df = 3)*Branch"),
                                                                   textInput("breducedModelFormulaStr","ReducedModelFormulaStr",value = "~sm.ns(Pseudotime, df = 3)"),
                                                                   numericInput("branch_point","Branch_point",value = 1,min = 1),
                                                                   textInput("branch_labels","Branch_labels",value = "NULL"),
                                                                   fluidRow(
                                                                     column(6,
                                                                            selectInput("brelative_expr","Relative_expr",
                                                                                        choices = c("TRUE","FALSE"),
                                                                                        selected = "TRUE")
                                                                     ),
                                                                     column(6,
                                                                            numericInput("bcores","Cores",value = 1,min = 1)
                                                                     )
                                                                   ),
                                                                   textAreaInput("pmbh_params", "Custom parameters for plot_multiple_branches_heatmap2:",
                                                                                 value = 'list(num_clusters = 4,branches = c(1,3,4,5))',
                                                                                 resize = "vertical"),
                                                                   actionButton("run_beam", HTML('<span style="color: red; font-weight: bold;">Runing BEAM</span>'))
                                                   ),
                                                   div(style = "margin-bottom: 15px;",
                                                       actionButton("m2plotst",HTML('<span style="color: red; font-weight: bold;">Submit Plot</span>'),width = "100%"),

                                                   ),
                                                   # ==============================download plot
                                                   bsCollapsePanel(title = 'Plot download',value = '',style = 'primary',
                                                                   fluidRow(
                                                                     column(6,numericInput("htwidth",label = "Plot width",min = 0,max = 50,value = 8)),
                                                                     column(6,numericInput("htheight",label = "Plot height",min = 0,max = 50,value = 10))
                                                                   ),
                                                                   downloadButton('download_m2htplot','Download plot')
                                                   )

                                               )
                                        ),
                                        column(8,
                                               box(title = "Table",status = "primary",width = 12,solidHeader = T,
                                                   tabsetPanel(
                                                     tabPanel("Heatmap",
                                                              plotOutput("m2heatmap")
                                                     ),
                                                     tabPanel("Diff genes table",
                                                              DT::DTOutput("m2diff_table") %>% withSpinner(type = 8),
                                                              downloadButton('download_m2diff_table','Download table')
                                                     ),
                                                     tabPanel("BEAM genes table",
                                                              DT::DTOutput("beam_table") %>% withSpinner(type = 8),
                                                              downloadButton('download_beam_table','Download table')
                                                     ),
                                                     tabPanel("Clustered table",
                                                              DT::DTOutput("m2cl_table") %>% withSpinner(type = 8),
                                                              downloadButton('download_m2cl_table','Download table')
                                                     )
                                                   )

                                               ))
                                      )
                             ),
                             tabPanel(HTML("<b><i>Monocle3 input</i></b>"),
                                      fluidRow(
                                        column(4,
                                               box(title = "Params",status = "navy",width = 12,solidHeader = T,
                                                   selectInput('m3example_data','Example data',
                                                               choices = c('worm embryo','NULL'),selected = 'NULL'),
                                                   div(style = "margin-bottom: -10px;",
                                                       fileInput("m3_data","Input data (monocle3 CellDataSet object with rds format)")
                                                   ),
                                                   numericInput("m3pval","P value threshold",value = 0.0001,min = 0,max = 1),

                                                   bsCollapsePanel(title = 'Graph test analysis',value = '',style = 'primary',
                                                                   div(style = "margin-bottom: -20px;",
                                                                       fileInput("m3diff_data","Monocle3 graph_test output with csv format")
                                                                   ),

                                                                   div(style = "border-top: 2px solid #333; margin: 20px 0;"),

                                                                   selectInput("neighbor_graph","Neighbor graph",choices = c("knn", "principal_graph"),selected = "principal_graph"),
                                                                   textInput("reduction_method","Reduction method",value = "UMAP"),

                                                                   fluidRow(
                                                                     column(6,
                                                                            numericInput("k","k",value = 25,min = 1,max = 100)
                                                                     ),
                                                                     column(6,
                                                                            numericInput("m3cores","Cores",value = 1,min = 1,max = 100)
                                                                     )
                                                                   ),
                                                                   selectInput("alternative","Alternative",choices = c("greater","two.sided"),selected = "greater"),

                                                                   actionButton("run_graph_test", HTML('<span style="color: red; font-weight: bold;">Runing graph_test analysis</span>'))

                                                   )
                                               )
                                        ),
                                        column(8,
                                               box(title = "Monocle3 diff table",status = "navy",width = 12,solidHeader = T,

                                                   tabsetPanel(
                                                     tabPanel("Diff genes table",
                                                              DT::DTOutput("m3diff_table") %>% withSpinner(type = 8),
                                                              downloadButton('download_m3diff_table','download table')
                                                     ),
                                                     tabPanel("Matrix table",
                                                              DT::DTOutput("m3mat_table") %>% withSpinner(type = 8),
                                                              downloadButton('download_m3mat_table','Download table')
                                                     )
                                                   )
                                               ))
                                      )
                             )
                           )

                  ),
                  tabPanel(tags$b("2.Cluster number estimation"),
                           column(4,
                                  box(title = "Params",status = "info",width = 12,solidHeader = T)
                           ),
                           column(8,
                                  box(title = "Plot",status = "purple",width = 12,solidHeader = T,
                                      plotOutput("get_cluster") %>% withSpinner(type = 8))
                           )
                  ),
                  tabPanel(tags$b("3.Data clustering"),
                           column(4,
                                  box(title = "Params",status = "info",width = 12,solidHeader = T,
                                      numericInput("cl_num","Cluster numbers",value = 1,min = 1,max = 20),
                                      selectInput("cl_method","Cluster methods",
                                                  choices = c("mfuzz","TCseq","kmeans"),
                                                  selected = "kmeans"),
                                      actionButton("run_clustering", HTML('<span style="color: red; font-weight: bold;">Clustering Started</span>'))
                                  )
                           ),
                           column(8,
                                  box(title = "Clustered table",status = "purple",width = 12,solidHeader = T,
                                      DT::DTOutput("clustered_table"),
                                      downloadButton('download_clustered_table','Download table'))
                           )
                  ),
                  tabPanel(tags$b("4.Enrichment analysis"),
                           column(4,
                                  box(title = "Params",status = "primary",width = 12,solidHeader = T,
                                      div(style = "margin-bottom: -25px;",
                                          fileInput("user_enrich", "Use own enrichment data")
                                      ),
                                      selectInput("type","Type (one or two options surpported)",
                                                  choices = c("BP","MF","CC","KEGG","ownSet"),
                                                  selected = "BP",multiple = T),
                                      div(style = "margin-bottom: -25px;",
                                          fileInput("TERM2GENE","Gene-term mapping data")
                                      ),
                                      fluidRow(
                                        column(6,
                                               selectInput("OrgDb","OrgDb",
                                                           choices = c("human","mouse"),
                                                           selected = "human")
                                        ),
                                        column(6,
                                               selectInput("id.trans","Id trans",
                                                           choices = c("TRUE","FALSE"),
                                                           selected = "TRUE")
                                        )
                                      ),
                                      fluidRow(
                                        column(6,
                                               textInput("fromType","From type",value = "SYMBOL")
                                        ),
                                        column(6,
                                               textInput("toType","To type",value = "ENTREZID")
                                        )
                                      ),
                                      fluidRow(
                                        column(6,
                                               selectInput("readable","Readable",
                                                           choices = c("TRUE","FALSE"),
                                                           selected = "TRUE")
                                        ),
                                        column(6,
                                               textInput("organism","Organism",value = "hsa")
                                        )
                                      ),
                                      sliderInput("pvalueCutoff","P value cutoff",min = 0,max = 1,step = 0.01,value = 0.05),
                                      fluidRow(
                                        column(6,
                                               numericInput("topn","Top number",min = 0,max = 50,value = 5)
                                        ),
                                        column(6,
                                               selectInput("add.gene","Add gene",
                                                           choices = c("TRUE","FALSE"),
                                                           selected = "FALSE")
                                        )
                                      ),


                                      actionButton("run_enrichment", HTML('<span style="color: red; font-weight: bold;">Run Enrichment for Clusters</span>'))
                                  )

                           ),
                           column(8,
                                  fluidPage(
                                    box(title = "Enrichment table",status = "maroon",width = 12,solidHeader = T,
                                        DT::DTOutput("enriched_table"),
                                        downloadButton('download_enrich1','Download table')
                                    )
                                  ),
                                  fluidPage(
                                    box(title = "Enrichment table2",status = "maroon",width = 12,solidHeader = T,
                                        DT::DTOutput("enriched_table2"),
                                        downloadButton('download_enrich2','Download table')
                                    )
                                  )
                           )
                  ),
                  tabPanel(tags$b("5.Integrated visualization"),
                           column(4,
                                  box(title = "Params",status = "primary",width = 12,solidHeader = T,
                                      fluidRow(
                                        column(12,tags$p("Select colors for heatmap:",
                                                         style = "font-weight: bold; text-decoration: underline;")),
                                        column(4,colourInput("low","Low","#08519C")),
                                        column(4,colourInput("mid","Mid","white")),
                                        column(4,colourInput("high","High","#A50F15"))
                                      ),
                                      fluidRow(
                                        column(12,tags$p("Select colors for membership lines:",
                                                         style = "font-weight: bold; text-decoration: underline;")),
                                        column(4,colourInput("line_low","Low","#0099CC")),
                                        column(4,colourInput("line_mid","Mid","grey90")),
                                        column(4,colourInput("line_high","High","#CC3333"))
                                      ),
                                      fluidRow(
                                        column(6,
                                               selectInput("plot_type","Plot type",choices = c("line","heatmap","both"),selected = "line")
                                        ),
                                        column(6,selectInput("add.bar","Add bar",
                                                             choices = c("TRUE","FALSE"),
                                                             selected = "FALSE"))
                                      ),
                                      textAreaInput("marker_gene","Gene symbols to be marked ",
                                                    value = "Supply gene symbols to be marked on heatmap with ‘,' seprated.",
                                                    resize = "vertical"),
                                      fluidRow(
                                        column(6,
                                               selectInput("line.side","Line side",choices = c("left","right"),selected = "right")
                                        ),
                                        column(6,selectInput("markGenes.side","MarkGenes side",choices = c("left","right"),selected = "right"))
                                      ),
                                      fluidRow(
                                        column(6,
                                               numericInput("ncol","Ncol",value = 4,min = 1,max = 30)
                                        ),
                                        column(6,colourInput("mline.col","Middle line color","#CC3333"))
                                      ),



                                      # ==============================
                                      bsCollapsePanel(title = 'Custom parameter setting',value = '',style = 'primary',
                                                      textAreaInput("custom_params", "Custom parameters:",
                                                                    value = 'list(column_names_rot = 45)',
                                                                    resize = "vertical")

                                      ),
                                      div(style = "margin-bottom: 15px;",
                                          actionButton("plotst",HTML('<span style="color: red; font-weight: bold;">Submit Plot</span>'),width = "100%"),

                                      ),


                                      # ==============================download plot
                                      bsCollapsePanel(title = 'Plot download',value = '',style = 'primary',
                                                      fluidRow(
                                                        column(6,numericInput("pwidth",label = "Plot width",min = 0,max = 50,value = 8)),
                                                        column(6,numericInput("pheight",label = "Plot height",min = 0,max = 50,value = 10))
                                                      ),
                                                      downloadButton('download_plot','Download plot')
                                      )

                                  )
                           ),
                           column(8,
                                  fluidPage(
                                    box(title = "Plot",status = "danger",width = 12,solidHeader = T,
                                        plotOutput("cmb_plot") %>% withSpinner(type = 8)

                                    )
                                  ),
                                  fluidPage(
                                    box(title = "Documentation",status = "danger",width = 12,solidHeader = T,
                                        actionButton("load_doc", "Load viscluster documentation"),
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
                                  )

                           )
                  )
                )
              )
      )
    )
  )
)
