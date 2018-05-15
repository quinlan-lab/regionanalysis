library("UpSetR")
expressionInput <- c(CCR=3707, MPC=671, pLI=506, `CCR&MPC`=609, `CCR&pLI`=1668, `MPC&pLI`=131, `CCR&MPC&pLI`=925) # change to have numbers from venndiagram.py
upset(fromExpression(expressionInput), order.by = "freq", mainbar.y.label = "Number of Genes", sets.x.label = "Total Number of Genes Per Metric", text.scale = c(2.1, 1.9, 2.1, 1.9, 2.1, 1.8), point.size = 3.5, line.size = 2)
dev.copy(pdf,"upset.pdf", width=18, height=12)
dev.off()
