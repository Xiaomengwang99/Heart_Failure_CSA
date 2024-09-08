library(shiny)
library(shinythemes)

# Define UI
ui <- fluidPage(theme = shinytheme("cosmo"),
                tags$head(
                  tags$style(HTML("
            p, ul, li, input, .form-group label {
                font-size: 22px;
            }
            select, .checkbox label{
                font-size: 23px;
            }
            h4, h3{
                font-size: 22px;
                font-weight: bold;
            }
            .title-panel h1 {
                font-size: 40px;
            }
            .action-button {
                font-size: 22px;
            }
            .shiny-text-output {
                font-size: 24px;
                font-weight: bold;
            }
        "))
                ),
                
                
    titlePanel("ICU Survival Lower Bound Estimator"),
    fluidRow(
      column(12,
             tags$p(HTML("This Shiny app provides an interactive platform to predict survival lower bound using the model in <em>Predicting survival time for critically ill patients with heart failure using conformalized survival analysis</em> by Wang, et al.(2024). 
                         The app is designed to predict the lower bound of survival time after intensive care unit (ICU) discharge for patients with heart failure diagnosis that have been admitted to the ICU.
                         Given a set of patient's information and medical history, this app calculates the lower bound of survival time for a specific alpha level.
                         We also provide a plot of predicted survival curve for alpha from 0.02 to 0.5 to help users visualize the prediction.")),
             tags$p("Instructions:"),
             tags$ul(
               tags$li("Step 1: Input patient's information and medical history above the line."),
               tags$li("Step 2: Input a desired alpha level between 0 and 1."),
               tags$li("Step 3: Click 'Calculate' to view the results.")
             )
      )
    ),
    fluidRow(
      column(8, fluidRow(
      column(6,
             selectInput("insurance", "Insurance:", choices = c("Medicaid", "Medicare", "Other")),
             selectInput("language", "Language Spoken:", choices = c("ENGLISH", "OTHER")),
             selectInput("marital_status", "Marital Status:", choices = c("OTHER","DIVORCED", "MARRIED", "SINGLE", "WIDOWED")),
             numericInput("los", "Length of ICU Stay (days):", 10),
             selectInput("gender", "Gender:", choices = c("M", "F")),
             numericInput("anchor_age", "Age:", 50, min = 0, max = 120),
             selectInput("race", "Race:", choices = c("Black", "White",
                          "Other", "Unknown", "Asian", "Hispanic"))
            
      ),
      column(4,
             h4("Does the patient have (or ever had) any of the following diagnosis?"),
             checkboxInput("has_type1_diabetes", "Type 1 diabetes", FALSE),
             checkboxInput("has_type2_diabetes", "Type 2 diabetes", FALSE),
             checkboxInput("has_hypertension", "Hypertension", FALSE),
             checkboxInput("has_atrial_fibrillation", "Atrial Fibrillation", FALSE),
             checkboxInput("has_copd", "Chronic Obstructive Pulmonary Disease (COPD)", FALSE),
             checkboxInput("has_asthma", "Asthma", FALSE),
             checkboxInput("has_liver_disease", "Liver Disease", FALSE),
             checkboxInput("has_chronic_kidney_disease", "Chronic Kidney Disease", FALSE),
             checkboxInput("has_cancer", "Cancer", FALSE),
             checkboxInput("has_depression", "Depression", FALSE),
             checkboxInput("has_anemia", "Anemia", FALSE),
      )
      ),
      fluidRow(
        tags$hr(style = "border-top: 2px solid black;"),
        column(6, 
          h3("Calculate predicted lower bound for specific alpha level"),
          numericInput("alpha", "Alpha Level:", 0.1, min = 0, max = 1),
          actionButton("calcBtn", "Calculate")
        ),
        column(4,
        h4("The lower bound for survival time (days) is:"),
        textOutput("riskOutput")
        )
      )
      ),
      column(4,
             h4("Does the patient take (or ever taken) any of the following medication?"),
             checkboxInput("ACE.Inhibitor", "ACE Inhibitor", FALSE),
             checkboxInput("Beta.Blocker", "Beta Blocker", FALSE),
             checkboxInput("Diuretic", "Diuretic", FALSE),
             checkboxInput("Anticoagulant", "Anticoagulant", FALSE),
             h4("Has the patient had previous heart failure diagnosis?"),
             checkboxInput("PreviousHF", "Yes", FALSE),
             tags$hr(style = "border-top: 2px solid black;"),
             h4("Plot of predicted survival curve for alpha from 0.02 to 0.5"),
             plotOutput("alpha_plot")
      )
    )
)

# Define server logic
server <- function(input, output, session){
    observeEvent(input$calcBtn, {
        # Here you would implement the risk calculation logic
        # For demonstration, I'm just returning a static message
        suppressPackageStartupMessages(library(conTree))
        suppressPackageStartupMessages(library(ggplot2))
        suppressPackageStartupMessages(library(survival))
        suppressPackageStartupMessages(library(gbm))
        #source("cfsurv_paper-main/utils/source.R")
        source("HF_shiny_allfunc.r")
        res <- readRDS("trained_model2.rds")
        X = data.frame(insurance = input$insurance,
                       language = input$language,
                       marital_status = input$marital_status,
                       los = input$los, gender = input$gender,
                       anchor_age = input$anchor_age, race = input$race,
                       has_type1_diabetes=input$has_type1_diabetes,
                       has_type2_diabetes=input$has_type2_diabetes,
                       has_hypertension=input$has_hypertension,
                       has_atrial_fibrillation=input$has_atrial_fibrillation,
                       has_copd=input$has_copd,
                       has_asthma=input$has_asthma,
                       has_liver_disease=input$has_liver_disease,
                       has_chronic_kidney_disease=input$has_chronic_kidney_disease,
                       has_cancer=input$has_cancer,
                       has_depression=input$has_depression,
                       has_anemia=input$has_anemia,
                       PreviousHF=input$PreviousHF,
                       ACE.Inhibitor=input$ACE.Inhibitor,
                       Beta.Blocker=input$Beta.Blocker,
                       Diuretic=input$Diuretic,
                       Anticoagulant=input$Anticoagulant
                       )
        X[, "insurance"] = as.factor(X[, "insurance"])
        X[, "language"] = as.factor(X[, "language"])
        X[, "marital_status"] = as.factor(X[, "marital_status"])
        X[, "gender"] = as.factor(X[, "gender"])
        X[, "race"] = as.factor(X[, "race"])
        
        
        df_train = readRDS("df_train.rds")
        
        p = dim(df_train)[2]
        
        x_new = X
        colnames(x_new) = colnames(df_train)[1:p]
        x_new[, 8:p] <- lapply(x_new[, 8:p], function(x) ifelse(x == FALSE, "False", ifelse(x == TRUE, "True", x)))
        x_new[, 3] <- lapply(x_new[, 3], function(x) ifelse(x == "OTHER", "", x))
        x_new[, 2] <- lapply(x_new[, 2], function(x) ifelse(x == "OTHER", "?", x))
        x_new = rbind(x_new, df_train[, 1:p])
        x_new[, c(2:3, 8:p)] <- lapply(x_new[, c(2:3,8:p)], as.factor)

        seed = 123
        
        ret <- cfsurv_predict(x = x_new[1, ], c = res$c,
                                 alpha= input$alpha,
                                 type="quantile",
                                 seed = seed,
                                 model = "distBoost_predict",
                                 dist= "weibull",
                                 I_fit = NULL,
                                 ftol=.1,tol=.1,
                                 n.tree=100, 
                                 weight_calib = res$weight_calib, 
                                 rf_model = res$rf_model,
                                 DB_gbm_mdl = res$DB_gbm_mdl, 
                                 DB_mdlrb = res$DB_mdlrb, 
                                 score = res$score, 
                                 res_T = res$res_T
        )
        
        alpha_lst = seq(0.02, 0.5, 0.02)
        l = length(alpha_lst)
        ret_lst = rep(NA, l)
        for (i in 1:l){
          ret_lst[i] <- cfsurv_predict(x = x_new[1, ], c = res$c,
                                       alpha= alpha_lst[i],
                                       type="quantile",
                                       seed = seed,
                                       model = "distBoost_predict",
                                       dist= "weibull",
                                       I_fit = NULL,
                                       ftol=.1,tol=.1,
                                       n.tree=100, 
                                       weight_calib = res$weight_calib, 
                                       rf_model = res$rf_model,
                                       DB_gbm_mdl = res$DB_gbm_mdl, 
                                       DB_mdlrb = res$DB_mdlrb, 
                                       score = res$score, 
                                       res_T = res$res_T
          )
        }
        
        output$riskOutput <- renderText({
            ret
        })
        output$alpha_plot <- renderPlot({
          temp = ret_lst < 364
          plot(alpha_lst[temp], ret_lst[temp], 
               ylab = "Predicted lower bound for survival", 
               xlab = "Alpha Level", 
               type = "l", 
               xlim = c(0, 0.5))
        })
    })
}

# Run the application
shinyApp(ui = ui, server = server)
