#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(ampir)
library(tidyverse)

options(shiny.maxRequestSize = 100*1024^2)

parse_fasta <- function(lines) {

    ### get sequence names
    seq_name_index <- grep(">", lines)
    seq_name <- gsub(">", "", lines[seq_name_index])
    
    ### get sequence
    seq_aa_start_index <- seq_name_index + 1
    seq_aa_end_index <- c(seq_name_index, length(lines)+1)[-1]-1
    
    seq_aa <- rep(NA, length(seq_name_index))

    withProgress(message = 'Reading sequences', value=0, {
    
        ### replace NA content with actual sequence content, and concatenate the lines
        for(i in seq_along(seq_name_index)){
            seq_aa_start <- seq_aa_start_index[i]
            seq_aa_end   <- seq_aa_end_index[i]
            seq_aa[i] <- gsub("[[:space:]]", "",
                              paste(lines[seq_aa_start:seq_aa_end],
                                    collapse = ""))
            
            if ( i %% 100 == 0){
                incProgress(i/length(lines))
            }
        }
        
        
        
    })
    
    # TODO: Add validation code for the sequences and report problems to the user
    
    data.frame(seq_name, seq_aa, stringsAsFactors = FALSE)
}

shinyServer(function(input, output) {

    valid_inputs <- eventReactive(input$go,{
        
        
        
        if (is.null(input$uploaded_sequences)){
            lines <- stringi::stri_split_lines(input$text_sequences)[[1]]
            return(parse_fasta(lines))
        } else {
            con <- file(input$uploaded_sequences$datapath)
            
            upload_lines <- readLines(con)
            close(con)
            return(parse_fasta(upload_lines))
        }
    })
    
    predictions <- reactive({
        ampir_model <- ifelse( input$ampir_model == "Full length precursor proteins (precursor)","precursor","mature")
        
        seqs <- valid_inputs()
        if ( nrow(seqs)==0){
            return(list(dt=NULL,info=NULL))
        }
        
        n_chunks <- ceiling(nrow(seqs)/3000)
        
        chunking_group <- rep(1:n_chunks,3000)[1:nrow(seqs)]
        
        chunked_seqs <- cbind(seqs,factor(chunking_group))
        
        progress <- shiny::Progress$new()
        on.exit(progress$close())
        
        progress$set(message = paste("Running prediction on ",nrow(seqs),"sequences"), value = 0)
        
        predict_with_progress <- function(df){
            
            p <- predict_amps(df,n_cores = 4, model = ampir_model)
            progress$inc(1/n_chunks)
            p
        }
        
        preds <- split(chunked_seqs,chunking_group) %>% map_dfr(predict_with_progress)
        
        dt <- preds %>% select(Name=seq_name,`AMP Probability`=prob_AMP,`Sequence`=seq_aa)
        info <- c(model=ampir_model)
        list(dt=dt,info=info)
    })
    
    output$final_results <- renderUI({
        ampir_version <- packageVersion("ampir")
        preds <- predictions()$dt
        ampir_model <- predictions()$info['model']
        
        
        has_outputs=(!is.null(preds) && (nrow(preds) > 0) )

        if ( has_outputs){
            valid_preds <- preds %>% filter(!is.na(`AMP Probability`))
            
            if ( nrow(valid_preds)==0){
                retval <- tagList(p(paste("You provided ",nrow(preds)," proteins as input but ampir could not run predictions for any of these because they were too short or contained invalid amino acids")),
                                 p("Reload this page to try a different set of sequences"))
            } else {
            
                retval <- tagList(
                    p(paste("Calculated AMP probabilities for ",nrow(valid_preds)," out of ",nrow(preds)," input protein sequences")),
                    p(paste("This search was performed with ampir version ",ampir_version," using the ",ampir_model," model. The R version was ",R.version.string)),
                    p("Click download button to retrieve results in csv format"),
                    p("Reload this page to perform a new search"),
                    downloadButton("downloadData", "Download"),
                
                    plotOutput('phistogram'),
        
            
                    dataTableOutput('prediction_table')
                )
            }
        } else {
            retval <- tagList(p("Your input was empty. You must provide at least one valid protein sequence."),
                              p("Reload this page to try a different set of sequences"))  
        }
        retval
    })
    
    output$prediction_table <- renderDataTable(
        predictions()$dt
    )
    
    output$phistogram <- renderPlot({
        preds <- predictions()

        preds$dt %>% ggplot(aes(x=`AMP Probability`)) + geom_histogram() + 
            xlab("Probability of antimicrobial activity") +
            ylab("Number of sequences") + ggtitle("Histogram of predicted AMP probabilities")
        
    })
    
    output$downloadData <- downloadHandler(
        filename = function() {
            paste("ampir_predictions", ".csv", sep = "")
        },
        content = function(file) {
            write_csv(predictions()$dt, file)
        }
    )
    
})
