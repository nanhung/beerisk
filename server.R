library(shiny); library(deSolve); library(ggplot2)

### Bee Risk model 
bee.risk <- function(t, x, parms) {
  with(as.list(c(parms, x)), 
  {
    # compute time-dependent parameter
    c1 <- c0*(0.448-0.090*sin(2*pi*(t-1)/365)-0.386*cos(2*pi*(t-1)/365))
    l1 <- l0*(0.588+0.149*sin(2*pi*(t-1)/365)-0.422*cos(2*pi*(t-1)/365))
    m1 <- m0*(0.584-0.139*sin(2*pi*(t-1)/365)-0.248*cos(2*pi*(t-1)/365))
    
    # compute toxicity
    md <- exp(0.0007767*d-4.487)
    
    #Function for survival of brood
    S <- (f^2/(f^2+b^2))*(H/(H+v))
    #recruitment function
    R <- amin-s*F/(H+F)+amax*(b^2/(b^2+f^2))
    
    # compute derivatives
    df <- c1*F-ga*(H+F)-gb*B
    dB <- l1*S-fi*B
    dH <- fi*B-R*H
    dF <- R*H-(m1+md)*F
    
    # combine results
    result <- c(df, dB, dH, dF)
    list(result)
  })
}


## Makes the data in a more ggplot-friendly format
stacker <- function(df){
  df.stack <- stack(df[, -1])
  df.stack$time <- rep(seq_len(nrow(df)), length(table(df.stack$ind)))
  names(df.stack)[1:2] <- c("population", "compartments")
  df.stack$compartments <- factor(df.stack$compartments, 
                                  levels=c("f", "B","H","F"), ordered=TRUE)
  return(df.stack)
  
}

## Start the Shiny stuff
shinyServer(function(input, output) { 
  output$guessPlot <- renderPlot(function() {
    
    timemax <- input$tmaxday  
    dt <- seq(0, timemax, by=1)
    
    y_max <- input$ymax
    
    # First, making sure we use the right model.
    ## SIR without seasonal effects
    c0 <- input$c0
    m0 <- input$m0
    l0 <- input$l0
    
    amin <- input$amin
    amax <- input$amax
    s <- input$s
    b <- input$b
    v <- input$v
    ga <- input$ga
    gb <- input$gb
    fi <- input$fi
    
    d<-input$d
    
    f <- input$f  
    B <- input$B
    H <- input$H
    F <- input$F
   
      inits <- c(f=f, B=B, H=H, F=F)
      guess_pars <- c(c0=c0, l0=l0, amin=amin, amax=amax, s=s,m0=m0, 
                      ga=ga, gb=gb, fi=fi, b=b, v=v, d=d)
      guess <- stacker(as.data.frame(ode(inits, dt, bee.risk, 
                                         parms=guess_pars)))
      
    # Base plot
      base.plot <- ggplot(guess, aes(x=time, y=population, 
                                     group=compartments, 
                                     color=compartments)) + 
        xlab("Time (day)") + ylab("Food Storagrs and Bee Population")
      
    # Label plot  
      label.plot <- base.plot + geom_line(alpha=.8) + 
      theme(legend.position="top") + theme(legend.title=element_blank()) +
        theme(axis.text=element_text(size=12), 
              axis.title=element_text(size=14,face="bold"))
      
    # Final plot  
      final.plot <- label.plot + 
      geom_line(aes(colour = compartments), size=1.5) + 
      scale_y_continuous("Food Storagrs and Bee Population \n", lim=c(0, y_max))+
      ggtitle("State Variables")
      
      print(final.plot)
    
  })  # closes reactivePlot
})  # Closes shinyServer