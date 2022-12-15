#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(markdown)

# Define UI for application that draws a histogram
shinyUI(fluidPage(
    tags$head(includeHTML(("google-analytics.html"))),
    includeCSS("style.css"),
    
    verticalLayout(
        titlePanel("ampir: Antimicrobial Peptide Prediction in R"),
        
        wellPanel(
            div(class="header",
                includeMarkdown("about.md")
            )
        ),
        
        conditionalPanel(
            condition = ("input.go == 0"),
            wellPanel(
                h3("Usage Instructions"),
                h4("Select the classification model"),
                includeHTML("instructions.html"),
                selectInput("ampir_model","Classification model",choices = c("Full length precursor proteins (precursor)","Mature peptides (mature)")),
                h4("Enter protein sequences in fasta format"),
                p("You may either paste sequences into the text box below (suitable for small searches) or upload a file containing your sequences"),
                br(),
                textAreaInput("text_sequences", "Paste protein sequences", 
                              value = "", width = "100%",height=200, placeholder = ""),
                h5("..or.."),
                fileInput("uploaded_sequences","Upload sequences in a file (Max 100Mb)"),
                
                actionButton("go", "Submit Sequences"),
                )
        ),
        uiOutput("final_results")
    )


))
