library(shiny)
library(shinyWidgets)
library(DT)

shinyUI(
  tagList(
    
    tags$style(".container{
                    width: 100%;
                    margin: 0 auto;
                    padding: 0;
                }
               @media screen and (min-width: 1300px){
                .container{
                    width: 1300px;
                }
               }"),
    
    tags$div(class="container",
             fluidPage(
               
               tags$head(
                 tags$style(HTML("hr {border-top: 1px solid #000000;}"))
               ),
               
               HTML('<br>'),
               
               titlePanel("COVID-19 Epidemic Modeling (Version 2 - Testing)"),
               
               actionLink('howtouse', 'Learn more about this tool.'),
               
               HTML('<br><br>'),
               
               sidebarLayout(
                 sidebarPanel(
                   HTML('<h4><b>Location Information</b></h4>'),
                   
                   dateInput(inputId = 'curr_date', 
                             label = 'Set Day 0 Date',
                   ),
                   
                   numericInput(inputId = 'num_people', 
                                label = 'Number of People in Area', 
                                value = 883305),
                   
                   # TODO: add in prediction fields for ICU patients, deaths and ventilators
                   uiOutput(outputId = 'prediction_fld'),
                   
                   uiOutput(outputId = 'prior_val'),
                   
                   hr(),
                   
                   HTML('<h4><b>Add Interventions</b></h4>'),
                   
                   checkboxInput(inputId = 'showint', 
                                 label = 'Add Intervention'),
                   
                   uiOutput(outputId = 'intervention_ui'),
                   
                   dataTableOutput(outputId = 'int_table'),
                   
                   hr(),
                   
                   HTML('<h4><b>Parameter Selection</b></h4>'),
                   
                   sliderInput('per.hosp', '% of Infections that Result in Hospitalizations', min = 0, max = 1, value = 0.06),
                   sliderInput('illness.length', 'Length of Illness if Not Hospitalized', min = 1, max = 20, value = 14), 
                   sliderInput('hosp.delay.time', 'Duration from Infection to Hospitalization', min = 1, max = 20, value = 10),
                   
                   HTML('<br><h4>If 100 people are in the non-ICU hospital on a given day, how many - the next day -
                        would:</h4>'),
                   
                   numericInput('p.g_icu', 'Go to ICU?', min = 0, max = 100, value = 10),
                   numericInput('p.g_d', 'Be discharged?', min = 0, max = 100, value = 5),

                   HTML('<br><h4>If 100 people are in the ICU unit on a given day, how many - the next day -
                        would:</h4>'),
                   
                   numericInput('p.icu_g', 'Move to the non-ICU hospital?', min = 0, max = 100, value = 15),
                   numericInput('p.icu_v', 'Go on a ventilator?', min = 0, max = 100, value = 50),


                   HTML('<br><h4>If 100 people are on a ventilator on a given day, how many - the next day -
                        would:</h4>'),
                   
                   numericInput('p.v_icu', 'Move to the non-ventilated ICU state?', min = 0, max = 100, value = 10),
                   numericInput('p.v_m', 'Die?', min = 0, max = 100, value = 7),
                   
                   hr(),
                   
                   HTML('<h4><b>Settings</b></h4>'),
                   
                   sliderInput(inputId = 'proj_num_days', 
                               label = 'Number of Days to Project', 
                               min = 10, 
                               max = 730, 
                               step = 5, 
                               value = 365),
                   
                   materialSwitch(inputId = "usedouble", 
                                  label = 'Use doubling time instead of Re', 
                                  status = 'primary'), 
                   
                   
                   HTML('<br><br><b>Notes</b>: This app is a modified version of the <a href="http://penn-chime.phl.io/">Penn Chime app</a>.
                 This is a beta version - the projections may or may not be accurate.
                        
                        <br><br> The code for this tool is on <a href="https://github.com/jpspeng/CovidShinyModel">Github</a>.'),
                   
                   tags$script("$(document).on('click', '#int_table button', function () {
                  Shiny.onInputChange('lastClickId',this.id);
                                             Shiny.onInputChange('lastClick', Math.random())
                                             });")
                 ),
                 
                 mainPanel(
                   wellPanel(
                     HTML('<h3><b>Day 0 Estimates</b></h3>'),
                     htmlOutput(outputId = 'infected_ct')
                   ),
                   
                   wellPanel(
                     HTML('<h3><b>Mean Hitting Times</b></h3>'),
                     htmlOutput(outputId = 'derived_stats')
                   ),
                   
                   wellPanel(
                     style = "background: white",
                     HTML('<h3><b>Projections</b></h3>'),
                     
                     radioGroupButtons(inputId = 'selected_graph', 
                                       label = '', 
                                       choices = c('Cases', 'Hospitalization', 'Hospital Resources'),
                                       justified = TRUE, 
                                       status = "primary"),
                     
                     uiOutput(outputId = 'plot_output'),
                     
                     HTML('<br>'),
                     
                     fluidRow(
                       align = 'center',
                       column(width = 1),
                       column(width = 1, 
                              HTML('<br>'),
                              actionButton("goleft", "", icon = icon("arrow-left"), width = '100%')),
                       column(width = 8,  
                              uiOutput(outputId = 'description')),
                       column(width = 1, 
                              HTML('<br>'),
                              actionButton("goright", "", icon = icon("arrow-right"), width = '100%')),
                       column(width = 1)
                     ),
                     
                     HTML('<br><br>'),
                     
                     div(
                       dataTableOutput(outputId = 'rendered.table'),
                       style = "font-size:110%"),
                     downloadButton(outputId = 'downloadData', 
                                    label = "Download Raw Outputs as CSV"),
                     HTML('<br><br>')
                   )
                 )
               )
             )
    )
  )
)