

CWT.mainPanel <- buildMainPanel(  id="CWT",render="plotOutput",
                                  title="Continuous Wavelet Transform (CWT)",
                                  title.AT="Acceleration (AT)",
                                  title.VT="Velocity (VT)",
                                  title.DT="Displacement (DT)",
                                  height="700px")

# CWT.tabPanel <- buildTabPanel(id="CWT")
CWT.tabPanel <- tabPanel(
  title="CWT",
  sidebarLayout(
    sidebarPanel(
      CWT.sidebarPanel(id="CWT")
    ),
    mainPanel(
      CWT.mainPanel
    )
  )
)



