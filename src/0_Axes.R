Add_Left_Axes <- function(
    Plot, 
    Meta_Data, 
    X_Col = NULL, 
    Y_Col = NULL,
    X_Label = NULL, 
    Y_Label = NULL,
    Arrow_Color = "black",
    Arrow_Width = 0.8,
    Text_Size = 4,
    Arrow_Len = 0.25,
    X_Label_Offset = 2.5,
    Y_Label_Offset = 3) {
  if (is.null(X_Label)) X_Label <- X_Col
  if (is.null(Y_Label)) Y_Label <- Y_Col
  
  X_Min <- min(Meta_Data[[X_Col]])
  Y_Min <- min(Meta_Data[[Y_Col]])
  Lower_Left <- floor(min(X_Min, Y_Min)) - 2
  Line_Len <- abs(0.3 * Lower_Left) + Lower_Left
  Mid_Point <- abs(0.3 * Lower_Left) / 2 + Lower_Left
  
  # Y Axes
  Plot <- Plot + 
    geom_segment(data = data.frame(X = Lower_Left, Y = Lower_Left,
                                   X_End = Lower_Left, Y_End = Line_Len),
                 mapping = aes(x = X, y = Y, xend = X_End, yend = Y_End),
                 arrow = arrow(length = unit(Arrow_Len, "cm"), ends = "last"),
                 color = Arrow_Color, linewidth = Arrow_Width, inherit.aes = FALSE)
  # X Axes
  Plot <- Plot + 
    geom_segment(data = data.frame(X = Lower_Left, Y = Lower_Left,
                                   X_End = Line_Len, Y_End = Lower_Left),
                 mapping = aes(x = X, y = Y, xend = X_End, yend = Y_End),
                 arrow = arrow(length = unit(Arrow_Len, "cm"), ends = "last"),
                 color = Arrow_Color, linewidth = Arrow_Width, inherit.aes = FALSE)
  # Axes Label
  Plot <- Plot + 
    geom_text(data = data.frame(Label = c(Y_Label, X_Label),
                                Angle = c(90, 0),
                                X = c(Lower_Left - Y_Label_Offset, Mid_Point),
                                Y = c(Mid_Point, Lower_Left - X_Label_Offset)),
              mapping = aes(x = X, y = Y, label = Label, angle = Angle),
              size = Text_Size, fontface = "bold", hjust = 0.5, vjust = 0.5, color = Arrow_Color, inherit.aes = FALSE)
  return(Plot)}
