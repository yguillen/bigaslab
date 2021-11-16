library(shiny)
library(Biobase)
library(gridExtra)


ui <- fluidPage(
  
  # Application title
  titlePanel("Interrogating the T-ALL expression dataset"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("thegene","Gene to Analyse",
                  choices=unique(melt_vastout$NAME),
                  selected  = "NFKBIA")
      ),
    
    # Show a plot of the generated distribution
    mainPanel(
      tableOutput("table"),
      plotOutput("boxplot"),
      #plotOutput("boxplot2"),
      plotOutput("boxplot4"),
      plotOutput("boxplot3"),
      plotOutput("boxplotclust"),
      plotOutput("boxplot5")
    )
  )
  
)

server <- function(input, output) {
  
  # If your data are stored in a csv or txt file, you could add the read.csv, read.delim commands here instead

output$boxplot <- renderPlot({
    
    gene <- input$thegene
    ggplot(data=melt_vastout[melt_vastout$NAME==gene,],aes(x=state,y=value))+
      geom_point(color="black",size=6)+
      geom_point(data=melt_vastout[melt_vastout$GSEA!="line" & melt_vastout$NAME==gene,],aes(color=Type),size=5)+
      geom_boxplot(data=melt_vastout[melt_vastout$NAME==gene,],alpha=0)+
      #geom_point(data=melt_vastout[melt_vastout$GSEA=="GSE69239" & melt_vastout$NAME==gene,],size=5)+
      geom_point(data=melt_vastout[melt_vastout$GSEA=="line" & melt_vastout$NAME==gene,],aes(color=Pheno),size=5)+
      #scale_shape_manual(values=c(0,8,2,5,6,7,10,11,12,9))+
      scale_color_brewer(palette="Spectral")+
      facet_wrap(~GSEA,scales="free_x",nrow=1)+
      theme_bw()+
      theme(axis.text.x = element_text(angle=45, hjust = 1),
            strip.text = element_text(size=5),
            legend.position = "bottom")
  })
  
output$boxplot2 <- renderPlot({
      gene <- input$thegene
      ggplot(data=melt_vastout_noncor[melt_vastout_noncor$NAME==gene,],aes(x=state,y=value))+
      geom_point(color="black",size=6)+
      geom_point(data=melt_vastout_noncor[melt_vastout_noncor$GSEA!="line" & melt_vastout_noncor$NAME==gene,],aes(color=Type),size=5)+
      geom_boxplot(data=melt_vastout_noncor[melt_vastout_noncor$NAME==gene,],alpha=0)+
      #geom_point(data=melt_vastout[melt_vastout$GSEA=="GSE69239" & melt_vastout$NAME==gene,],size=5)+
      geom_point(data=melt_vastout_noncor[melt_vastout_noncor$GSEA=="line" & melt_vastout_noncor$NAME==gene,],aes(color=Pheno),size=5)+
      #scale_shape_manual(values=c(0,8,2,5,6,7,10,11,12,9))+
      scale_color_brewer(palette="Spectral")+
      facet_wrap(~GSEA,scales="free",nrow=1)+
      theme_bw()+
      theme(axis.text.x = element_text(angle=45, hjust = 1),
            strip.text = element_text(size=5),
            legend.position = "bottom")
    
  })


output$boxplot3 <- renderPlot({
  gene <- input$thegene
  ggplot(data=melt_vastout[melt_vastout$NAME==gene,],aes(x=First_event,y=value))+
    geom_point(color="black",size=6)+
    geom_point(data=melt_vastout[melt_vastout$GSEA!="line" & melt_vastout$NAME==gene,],aes(color=First_event),size=5)+
    geom_boxplot(data=melt_vastout[melt_vastout$NAME==gene,],alpha=0)+
    #geom_point(data=melt_vastout[melt_vastout$GSEA=="GSE69239" & melt_vastout$NAME=="RBM39",],size=5)+
    geom_point(data=melt_vastout[melt_vastout$GSEA=="line" & melt_vastout$NAME==gene,],aes(color=Pheno),size=5)+
    #scale_shape_manual(values=c(0,8,2,5,6,7,10,11,12,9))+
    scale_color_brewer(palette="Spectral")+
    #facet_wrap(~GSEA,scales="free",nrow=1)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          legend.position = "bottom")
  
  
  
})

#melt_vastout_clust<-merge(melt_vastout,cuns,by.x="sample",by.y="SampID")
#melt_vastout_clust$uns_clust<-as.factor(melt_vastout_clust$uns_clust)

output$boxplotclust <- renderPlot({
  gene <- input$thegene
  ggplot(data=melt_vastout_clust[melt_vastout_clust$NAME==gene,],aes(x=as.factor(uns_clust),y=value))+
    geom_point(color="black",size=6)+
    geom_point(data=melt_vastout_clust[melt_vastout_clust$NAME==gene,],aes(color=uns_clust),size=5)+
    geom_boxplot(data=melt_vastout_clust[melt_vastout_clust$NAME==gene,],aes(color=uns_clust),alpha=0)+
    #geom_point(data=melt_vastout[melt_vastout$GSEA=="GSE69239" & melt_vastout$NAME=="RBM39",],size=5)+
    #scale_shape_manual(values=c(0,8,2,5,6,7,10,11,12,9))+
    scale_color_brewer(palette="Spectral")+
    #facet_wrap(~GSEA,scales="free",nrow=1)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          legend.position = "bottom")
  
  
  
})

output$boxplot4 <- renderPlot({
  gene <- input$thegene
  ggplot(data=melt_vastout[melt_vastout$NAME==gene & melt_vastout$GSEA=="GSE69239",],aes(x=Source,y=value))+
    geom_point(color="black",size=6)+
    geom_point(data=melt_vastout[melt_vastout$GSEA!="line" & melt_vastout$NAME==gene & melt_vastout$GSEA=="GSE69239",],aes(color=Type),size=5)+
    #geom_boxplot(data=melt_vastout[melt_vastout$NAME=="NOTCH1" & melt_vastout$GSEA=="GSE69239",],alpha=0)+
    geom_point(data=melt_vastout[melt_vastout$GSEA=="line" & melt_vastout$NAME=="NOTCH1" & melt_vastout$GSEA=="GSE69239",],aes(color=Pheno),size=5)+
    geom_ribbon(aes(group=GSEA,ymin=0,ymax=value),color="black",alpha=0.1)+
    scale_color_brewer(palette="Spectral")+
    facet_wrap(~GSEA,scales="free_x",nrow=1)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=45, hjust = 1),
          strip.text = element_text(size=9),
          legend.position = "bottom")
  
  
  
})

output$boxplot5 <- renderPlot({
  gene <- input$thegene
ggplot(difASm[difASm$GENE==gene,],aes(x=Type,y=value))+
  geom_point(color="black",size=2)+
  geom_point(aes(color=Type),size=1)+
  scale_color_brewer(palette="Spectral")+
  geom_boxplot(alpha=0)+
  new_scale_color()+
  geom_point(data=difASm[difASm$GENE==gene & difASm$GSEA=="line",],aes(color=Pheno),size=5)+
  scale_color_manual(values=c("yellow","pink"))+
  facet_wrap(~EVENT,scales="free_x",ncol = 4)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=45, hjust = 1))

})



output$table <- renderDataTable({
 subset(melt_vastout,melt_vastout$NAME==input$thegene)
})

}


# Run the application 
shinyApp(ui = ui, server = server)
