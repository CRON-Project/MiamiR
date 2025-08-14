library(shiny)

ui <- fluidPage(
  actionButton("go", "Run Task"),
  tags$pre(id = "plotProgressBar", "[                              ] 0% (Starting...)"),
  tags$script(HTML("
    Shiny.addCustomMessageHandler('plot_progress', function(message) {
      const el = document.getElementById('plotProgressBar');
      if (el) {
        const totalBars = 30;
        const filled = Math.round(message.pct / 100 * totalBars);
        const bar = '[' + '='.repeat(filled) + ' '.repeat(totalBars - filled) + '] ' + message.pct + '% (' + message.msg + ')';
        el.textContent = bar;
      }
    });
  "))
)

server <- function(input, output, session) {
  observeEvent(input$go, {
    for (i in 1:10) {
      session$sendCustomMessage("plot_progress", list(
        pct = round(i * 10),
        msg = paste("Line", i, "of 10")
      ))
      Sys.sleep(0.2)
    }
  })
}

shinyApp(ui, server)
