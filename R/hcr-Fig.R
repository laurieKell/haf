#' HCR Plot with Enhanced Regions
#' 
#' @title HCR Plot with Enhanced Regions
#' 
#' @description Creates a comprehensive harvest control rule plot showing different management zones with enhanced visual styling.
#' 
#' @param btarget \code{numeric} Target biomass reference point (default: 1.0)
#' @param bmax \code{numeric} Maximum biomass reference point (default: 1.5)
#' @param blim \code{numeric} Limit biomass reference point (default: 0.1)
#' @param bthresh \code{numeric} Threshold biomass reference point (default: 0.8)
#' @param btrigger \code{numeric} Trigger biomass reference point (default: 0.9)
#' @param ftarget \code{numeric} Target fishing mortality (default: 0.8)
#' @param fmax \code{numeric} Maximum fishing mortality (default: 1.5)
#' @param fmsy \code{numeric} MSY fishing mortality (default: 1.0)
#' @param textSize \code{numeric} Text size for annotations (default: 4)
#' 
#' @return \code{ggplot} object showing the HCR with colored management zones
#' 
#' @details
#' Creates a plot showing the hockey-stick HCR with five management zones:
#' \itemize{
#'   \item \strong{Collapsed} (red): B < Blim
#'   \item \strong{Rebuilding} (yellow): Blim ≤ B < Bthreshold
#'   \item \strong{Sustainable} (green): Bthreshold ≤ B < Btarget, F ≤ FMSY
#'   \item \strong{Overfishing} (orange): Bthreshold ≤ B < Btarget, F > FMSY
#'   \item \strong{Overfished} (red): Btarget ≤ B < Bmax, F > FMSY
#' }
#' 
#' @examples
#' \dontrun{
#' # Create standard HCR plot
#' hcrPlot2()
#' 
#' # Customize reference points
#' hcrPlot2(btarget=1.2, blim=0.2, ftarget=0.6)
#' }
#' 
#' @export
hcrPlot2<-function(btarget=1.0,bmax=1.5,blim=0.1,bthresh=0.8,btrigger=0.9,
                   ftarget=0.8,fmax=1.5,fmsy=1.0,
                   textSize=4){
  
  bref=c(target=btarget,max=bmax,lim=blim,thresh=bthresh,trigger=btrigger)
  fref=c(target=ftarget,max=fmax,msy=fmsy)
  
  regions=data.frame(
    xmin=c( 0.00, bref["lim"], bref["thresh"], bref["thresh"],  bref["thresh"]),
    xmax=c( bref["lim"], bref["thresh"], bref["max"], bref["max"], bref["lim"]),
    ymin=c(-0.03,-0.03,-0.03, fref["msy"], fref["msy"]),
    ymax=c( fref["max"], fref["max"], fref["msy"], fref["max"], fref["max"]),
    fill=c("Collapsed","Rebuilding", "Sustainable", "Overfishing","Overfished"))
  
  region_colors=c(
    "Collapsed"  ="#b03a3a",
    "Overfished" ="#b03a3a",
    "Overfishing"="orange",
    "Rebuilding" ="yellow",
    "Sustainable"="#2e6f40")
  
  p=ggplot()+
    geom_rect(data=regions,
              aes(xmin=xmin, xmax=xmax,
                  ymin=ymin, ymax=ymax,
                  fill=fill), alpha=1)+
    scale_fill_manual(values=region_colors)+
    scale_x_continuous(expand=c(0,0))+
    scale_y_continuous(expand=c(0,0))+
    geom_vline(xintercept=c(bref["lim"]), linetype="dashed")+
    geom_hline(yintercept=c(1),   linetype="solid")+
    geom_line(aes(x=rep(bref["target"],2),y=c(-0.03,fref["target"])))+
    geom_line(aes(x=c(0,bref["lim"],bref["thresh"],bref["max"]), y=c(0,0,rep(fref["target"],2))),col="steelblue",linewidth=1)+
    geom_line(aes(x=c(0,0.4,bref["thresh"]), y=c(0,0.0,fref["target"])),col="steelblue",linewidth=1)+
    geom_line(aes(x=c(0,0.6,bref["thresh"]), y=c(0,0.0,fref["target"])),col="steelblue",linewidth=1)+
    labs(x=expression(B/B[target]),
         y=expression(F/F[target]))+
    theme(legend.position="none")+
    annotate("text", label="Collapsed",        x=bref["lim"]*0.5,   y=fref["msy"]*0.5,vjust=0.5,hjust=1,  angle=90, size=textSize)+
    annotate("text", label="Overfished",       x=bref["thresh"]*0.5,y=(fref["msy"]+fref["max"])*0.5,vjust=0,  hjust=0,  angle=0,  size=textSize,col="grey95")+
    annotate("text", label="Overfishing",      x=btarget,           y=(fref["msy"]+fref["max"])*0.5,vjust=0,  hjust=0,            size=textSize)+
    annotate("text", label="Rebuilding",       x=bref["thresh"]*0.5,y=fref["msy"]*0.5,vjust=0,  hjust=0,  angle=0,  size=textSize)+
    annotate("text", label=expression(B[lim]), x=bref["lim"],       y=fref["target"], vjust=0,  hjust=1,            size=textSize)+
    annotate("text", label=expression(B[thr]), x=bref["thresh"],    y=fref["target"], vjust=0,  hjust=1,            size=textSize)+
    annotate("text", label=expression(B[tgt]), x=bref["target"],    y=fref["target"], vjust=0,  hjust=1,            size=textSize)+
    annotate("text", label="Sustainable",      x=bref["target"],    y=fref["msy"]*0.5,vjust=0,  hjust=0,    angle=0,size=textSize,col="grey95")+ 
    annotate("text", label=expression(F[tgt]), x=bref["max"],       y=fref["target"], vjust=0,  hjust=1,            size=textSize)
  
    p}

#' Standard HCR Plot
#' 
#' @title Standard HCR Plot
#' 
#' @description Creates a standard harvest control rule plot showing the hockey-stick relationship between biomass and fishing mortality.
#' 
#' @param btarget \code{numeric} Target biomass reference point (default: 1.0)
#' @param bmax \code{numeric} Maximum biomass reference point (default: 1.5)
#' @param blim \code{numeric} Limit biomass reference point (default: 0.1)
#' @param bthresh \code{numeric} Threshold biomass reference point (default: 0.8)
#' @param btrigger \code{numeric} Trigger biomass reference point (default: 0.9)
#' @param ftarget \code{numeric} Target fishing mortality (default: 1.0)
#' @param fmax \code{numeric} Maximum fishing mortality (default: 1.5)
#' @param fadv \code{numeric} Advice fishing mortality (default: ftarget)
#' @param textSize \code{numeric} Text size for annotations (default: 4)
#' 
#' @return \code{ggplot} object showing the standard HCR
#' 
#' @details
#' Creates a plot showing the standard hockey-stick HCR with three main zones:
#' \itemize{
#'   \item \strong{Critical} (brown): B < Blim
#'   \item \strong{Rebuilding} (yellow): Blim ≤ B < Btrigger
#'   \item \strong{Sustainable} (green): Btrigger ≤ B < Bmax
#' }
#' 
#' @examples
#' \dontrun{
#' # Create standard HCR plot
#' hcrPlot()
#' 
#' # Customize for specific stock
#' hcrPlot(btarget=1.2, blim=0.15, ftarget=0.8)
#' }
#' 
#' @export
hcrPlot<-function(btarget=1.0,bmax=1.5,blim=0.1,bthresh=0.8,btrigger=0.9,
                  ftarget=1.0,fmax=1.5,fadv=ftarget,
                  textSize=4){
  
  bref=c(target=btarget,max=bmax,lim=blim,thresh=bthresh,trigger=btrigger)
  fref=c(target=ftarget,max=fmax,adv=fadv)
  
  btarget=bref["target"]
  if (is.na(btarget)) btarget=1
  
  regions=data.frame(
    xmin=c(0.00, bref["lim"], bref["thresh"]),
    xmax=c(      bref["lim"], bref["thresh"], bref["max"]),
    ymin=c(0.00, 0.00, 0.00),
    ymax=rep(fref["max"], 3),
    fill=c("Critical", "Overfished", "Overfishing"))
  
  region1=data.frame(
    x=c(0, bref["lim"],bref["lim"],0),
    y=c(0, fref["target"]*bref["lim"]/bref["trigger"], 0,0))
  region2=data.frame(
    x=c(0, bref["thresh"],bref["thresh"],0),
    y=c(0, fref["target"]/bref["trigger"]*bref["thresh"], 0,0))
  region3=data.frame(
    x=c(bref["thresh"],bref["thresh"],bref["trigger"],bref["max"],bref["max"],bref["thresh"]),
    y=c(0,             fref["target"]*bref["thresh"]/bref["trigger"],fref["target"], fref["target"], 0,0))
    
  region_colors=c(
    "Critical"   ="brown",
    "Overfished" ="#b03a3a",
    "Overfishing"="orange",
    "Rebuilding" ="yellow",
    "Sustainable"="#2e6f40")
  
  # Create the plot
  p=ggplot()+
    geom_rect(data=regions,
              aes(xmin=xmin, xmax=xmax,
                  ymin=ymin, ymax=ymax,
                  fill=fill), alpha=1)+
    geom_polygon(data=region2,
                  aes(x,y),fill="yellow")+
    geom_polygon(data=region3,
                  aes(x,y),fill="#2e6f40")+
    geom_polygon(data=region1,
                 aes(x,y),fill="brown3")+
    scale_fill_manual(values=region_colors)+
    scale_x_continuous(expand=c(0, 0), limits=c(0, bref["max"]))+
    scale_y_continuous(expand=c(0, 0), limits=c(0, fref["max"]))+
    geom_line(aes(x=c(0,bref["trigger"],bref["trigger"]),
                  y=c(  fref["target"], fref["target"],0)), color="black",linetype="dashed")+
    geom_line(aes(x=c(0,bref["trigger"],bref["max"]),
                   y=c(0,fref["target"], fref["target"])), color="black")+
    geom_line(aes(x=rep(bref["lim"],2),
                  y=c(0,bref["lim"]*fref["target"]/bref["trigger"])), color="black")+
    geom_line(aes(x=rep(bref["thresh"],2),
                  y=c(0,bref["thresh"]*fref["target"]/bref["trigger"])), color="black")
  
  p=p+
    annotate("text", label=expression(B[lim]),      x=bref["lim"],   y=0.1, hjust=1, vjust=0,size=textSize*0.5,col="grey95")+
    annotate("text", label=expression(B[threshold]),x=bref["thresh"],y=fref["target"]/bref["trigger"]*bref["thresh"], 
                                                                                       hjust=1, vjust=0,size=textSize*0.5,col="grey95")+
    annotate("text", label=expression(B[trigger]), x=bref["trigger"],y=fref["target"], hjust=1, vjust=0,size=textSize*0.5)+
    annotate("text", label=expression(F[target]),  x=bref["max"],    y=fref["target"], hjust=1, vjust=0,size=textSize*0.5)+
    annotate("text", label="Critical",             x=0.05, y=bref["thresh"], hjust=0.5, vjust=0.5,angle=90,size=textSize*0.5,col="grey95")+
    annotate("text", label="Rebuilding",           x=0.50, y=0.25, hjust=0.5, vjust=1,size=textSize*0.5)+
    annotate("text", label="Sustainable",          x=(btarget+bref["max"])/2,  y=0.25, hjust=0.5, vjust=1,size=textSize*0.5,col="grey95")+
    annotate("text", label="Overfished",           x=0.50, y=1.35, hjust=0.5, vjust=1,size=textSize*0.5,col="grey95")+
    annotate("text", label="Overfishing",          x=(btarget+bref["max"])/2,  y=1.35, hjust=0.5, vjust=1,size=textSize*0.5)+
    labs(x=expression(B/B[target]),
         y=expression(F/F[target]))+
    theme_minimal()+theme(
      legend.position="none",
      axis.text.x =element_text(size=textSize),
      axis.text.y =element_text(size=textSize),
      axis.title.x=element_text(size=textSize*1.2),
      axis.title.y=element_text(size=textSize*1.2))
  
    if (!is.na(bref["target"])){ 
       p=p+geom_line(aes(x=rep(btarget,2),
                  y=c(0,fref["target"])), color="blue")+
       annotate("text", label=expression(B[target]),  x=btarget, y=fref["target"], hjust=1, vjust=1,size=textSize*0.5)}
    
     suppressWarnings(p)}
