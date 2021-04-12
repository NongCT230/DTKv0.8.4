

# Basic statistical data  --------------------------------------------


#v0.7.5

#class contains MFA data and MFA analysis calculations.
# - buildMFAtable MFA and MFAplot

#the design goal with this is that this is part of a composition where the plotting classes need their data as an attribute.
#this will be faster because it erases the need to copy paste. This does not mean that a class can't have their specific data.
#this also enables us to add stuff onto MFA easier but makes for a really big document. I am changing this to a class since it seems to be the design choice for the other analysis but it doesn't have to be a class.
#I am just following the design choices here.

Class.Analysis.MFA <- R6Class("Class.Analysis.MFA",
  
	private = list(
		chem_data = NULL,
		chem_data_mfa = NULL,
		
		#updated FactoMineR::plot.MFA function which records the plots it makes so they can be printed individually later
		plot.MFA = function (x, axes = c(1, 2), choix = c("ind", "var", "group", 
                                       "axes", "freq"), ellipse = NULL, ellipse.par = NULL, lab.grpe = TRUE, 
          lab.var = TRUE, lab.ind = TRUE, lab.par = FALSE, lab.col = TRUE, 
          habillage = "group", col.hab = NULL, invisible = c("none", 
                                                             "ind", "ind.sup", "quanti", "quanti.sup", "quali", "quali.sup", 
                                                             "row", "row.sup", "col", "col.sup"), partial = NULL, 
          lim.cos2.var = 0, chrono = FALSE, xlim = NULL, ylim = NULL, 
          title = NULL, palette = NULL, autoLab = c("auto", "yes", 
                                                    "no"), new.plot = FALSE, select = NULL, unselect = 0.7, 
          shadowtext = FALSE, legend = list(bty = "y", x = "topleft"), 
          ...) 
				{
				  
				  #ADDITION - plot objects
				  plotlist <- list();
				  #
				  
				  res.mfa <- x
				  if (!inherits(res.mfa, "MFA")) 
					stop("non convenient data")
				  if (is.numeric(unselect)) 
					if ((unselect > 1) | (unselect < 0)) 
					  stop("unselect should be betwwen 0 and 1")
				  autoLab <- match.arg(autoLab, c("auto", "yes", "no"))
				  if (autoLab == "yes") 
					autoLab = TRUE
				  if (autoLab == "no") 
					autoLab = FALSE
				  choix <- match.arg(choix, c("ind", "var", "group", "axes", 
											  "freq"))
				  invisible <- match.arg(invisible, c("none", "ind", "ind.sup", 
													  "quanti", "quanti.sup", "quali", "row", "row.sup", "col", 
													  "col.sup"), several.ok = TRUE)
				  if ("none" %in% invisible) 
					invisible = NULL
				  lab.x <- paste("Dim ", axes[1], " (", format(res.mfa$eig[axes[1], 
																		   2], nsmall = 2, digits = 2), "%)", sep = "")
				  lab.y <- paste("Dim ", axes[2], " (", format(res.mfa$eig[axes[2], 
																		   2], nsmall = 2, digits = 2), "%)", sep = "")
				  group <- res.mfa$call$group
				  nbre.grpe <- length(group)
				  type <- res.mfa$call$type
				  num.group.sup = NULL
				  if (!is.null(res.mfa$call$num.group.sup)) {
					num.group.sup <- res.mfa$call$num.group.sup
					nbre.grpe.sup <- length(num.group.sup)
					type.sup <- type[num.group.sup]
					type.act <- type[-num.group.sup]
					nbre.grpe <- nbre.grpe - length(num.group.sup)
				  }
				  if (choix == "axes") {
					if ((new.plot) & !nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) 
					  dev.new()
					if (is.null(palette)) 
					  palette(c("black", "red", "green3", "blue", "magenta", 
								"darkgoldenrod", "darkgreen", "darkgray", "cyan", 
								"violet", "turquoise", "orange", "lightpink", 
								"lavender", "yellow", "lightgreen", "lightgrey", 
								"lightblue", "darkkhaki", "darkmagenta", "darkolivegreen", 
								"lightcyan", "darkorange", "darkorchid", "darkred", 
								"darksalmon", "darkseagreen", "darkslateblue", 
								"darkslategray", "darkslategrey", "darkturquoise", 
								"darkviolet", "lightgray", "lightsalmon", "lightyellow", 
								"maroon"))
					if (is.null(title)) 
					  title <- "Partial axes"
					plot(0, 0, xlab = lab.x, ylab = lab.y, xlim = c(-1.1, 
																	1.1), ylim = c(-1.1, 1.1), col = "white", asp = 1, 
						 main = title, ...)
					x.cercle <- seq(-1, 1, by = 0.01)
					y.cercle <- sqrt(1 - x.cercle^2)
					lines(x.cercle, y = y.cercle, ...)
					lines(x.cercle, y = -y.cercle, ...)
					abline(v = 0, lty = 2, ...)
					abline(h = 0, lty = 2, ...)
					coord.axes <- res.mfa$partial.axes$coord[, axes, drop = FALSE]
					if (!is.null(select)) {
					  if (mode(select) == "numeric") 
						selection <- select
					  else {
						if (sum(rownames(res.mfa$partial.axes$coord) %in% 
								select) != 0) 
						  selection <- which(rownames(res.mfa$partial.axes$coord) %in% 
											   select)
						else {
						  if (grepl("contrib", select)) 
							selection <- (rev(order(res.mfa$partial.axes$contrib[, 
																				 axes[1], drop = FALSE] * res.mfa$eig[axes[1], 
																													  1] + res.mfa$partial.axes$contrib[, axes[2], 
																																						drop = FALSE] * res.mfa$eig[axes[2], 1])))[1:min(nrow(res.mfa$partial.axes$coord), 
																																																		 sum(as.integer(unlist(strsplit(select, 
																																																										"contrib"))), na.rm = T))]
						  if (grepl("coord", select)) 
							selection <- (rev(order(apply(res.mfa$partial.axes$coord[, 
																					 axes]^2, 1, sum))))[1:min(nrow(res.mfa$partial.axes$coord), 
																											   sum(as.integer(unlist(strsplit(select, 
																																			  "coord"))), na.rm = T))]
						  if (is.integer(select)) 
							selection <- select
						}
					  }
					}
					if (habillage == "group") {
					  if (is.null(col.hab) | length(col.hab) < length(group)) {
						if (is.null(res.mfa$call$num.group.sup)) 
						  col.hab <- 2:(length(group) + 1)
						else {
						  col.hab[which(!(1:length(group)) %in% (res.mfa$call$num.group.sup))] <- 2:(1 + 
																									   length(group) - length(res.mfa$call$num.group.sup))
						  col.hab[res.mfa$call$num.group.sup] <- length(group) - 
							length(res.mfa$call$num.group.sup) + 1 + 
							(1:length(res.mfa$call$num.group.sup))
						}
					  }
					  i = 1
					  couleur.axes <- col.hab[i]
					  auxil = strsplit(rownames(res.mfa$partial.axes$coord)[1], 
									   ".", fixed = TRUE)[[1]]
					  auxil2 = auxil[length(auxil)]
					  for (j in 2:nrow(res.mfa$partial.axes$coord)) {
						auxil = strsplit(rownames(res.mfa$partial.axes$coord)[j], 
										 ".", fixed = TRUE)[[1]]
						if (auxil2 != auxil[length(auxil)]) {
						  i = i + 1
						  auxil2 = auxil[length(auxil)]
						}
						couleur.axes <- c(couleur.axes, col.hab[i])
					  }
					}
					else {
					  couleur.axes <- NULL
					  for (i in 1:length(group)) couleur.axes <- c(couleur.axes, 
																   rep("black", ncol(res.mfa$partial.axes$coord)))
					}
					posi <- coll <- NULL
					col.legend = unique(couleur.axes)
					if (!is.null(select)) {
					  coord.axes <- coord.axes[selection, , drop = FALSE]
					  couleur.axes <- couleur.axes[selection]
					}
					for (v in 1:nrow(coord.axes)) {
					  arrows(0, 0, coord.axes[v, 1], coord.axes[v, 2], 
							 length = 0.1, angle = 15, code = 2, col = couleur.axes[v], 
							 ...)
					  if (abs(coord.axes[v, 1]) > abs(coord.axes[v, 2])) {
						if (coord.axes[v, 1] >= 0) 
						  posi <- c(posi, 4)
						else posi <- c(posi, 2)
					  }
					  else {
						if (coord.axes[v, 2] >= 0) 
						  posi <- c(posi, 3)
						else posi <- c(posi, 1)
					  }
					  labe <- rownames(coord.axes)
					}
					if (autoLab == "auto") 
					  autoLab = (length(labe) < 50)
					if (autoLab == FALSE) 
					  text(coord.axes[, 1], y = coord.axes[, 2], labels = labe, 
						   pos = posi, col = couleur.axes, ...)
					if (autoLab == TRUE) 
					  autoLab(coord.axes[, 1], y = coord.axes[, 2], labels = labe, 
							  col = couleur.axes, shadotext = shadowtext, 
							  ...)
					if (habillage == "group") {
					  L <- list(x = "topleft", legend = rownames(res.mfa$group$Lg)[-length(rownames(res.mfa$group$Lg))], 
								text.col = col.legend)
					  L <- modifyList(L, legend)
					  do.call(graphics::legend, L)
					}
					
					#ADDITION
					plotlist[["axes"]] <- recordPlot();
					#
					
				  }
				  if (choix == "group") {
					coord.actif <- res.mfa$group$coord[, axes, drop = FALSE]
					if (!is.null(res.mfa$group$coord.sup)) 
					  coord.illu <- res.mfa$group$coord.sup[, axes, drop = FALSE]
					selection <- selectionS <- NULL
					if (!is.null(select)) {
					  if (mode(select) == "numeric") 
						selection <- select
					  else {
						if (sum(rownames(res.mfa$group$coord) %in% select) + 
							sum(rownames(res.mfa$group$coord.sup) %in% 
								select) != 0) 
						  selection <- which(rownames(res.mfa$group$coord) %in% 
											   select)
						else {
						  if (grepl("contrib", select)) 
							selection <- (rev(order(res.mfa$group$contrib[, 
																		  axes[1]] * res.mfa$eig[axes[1], 1] + res.mfa$group$contrib[, 
																																	 axes[2]] * res.mfa$eig[axes[2], 1])))[1:min(nrow(res.mfa$group$coord), 
																																												 sum(as.integer(unlist(strsplit(select, 
																																																				"contrib"))), na.rm = T))]
						  if (grepl("coord", select)) 
							selection <- (rev(order(apply(res.mfa$group$coord[, 
																			  axes]^2, 1, sum))))[1:min(nrow(res.mfa$group$coord), 
																										sum(as.integer(unlist(strsplit(select, 
																																	   "coord"))), na.rm = T))]
						  if (grepl("cos2", select)) {
							if (sum(as.numeric(unlist(strsplit(select, 
															   "cos2"))), na.rm = T) >= 1) 
							  selection <- (rev(order(apply(res.mfa$group$cos2[, 
																			   axes], 1, sum))))[1:min(nrow(res.mfa$group$coord), 
																									   sum(as.numeric(unlist(strsplit(select, 
																																	  "cos2"))), na.rm = T))]
							else selection <- which(apply(res.mfa$group$cos2[, 
																			 axes], 1, sum) > sum(as.numeric(unlist(strsplit(select, 
																															 "cos2"))), na.rm = T))
						  }
						  if (is.integer(select)) 
							selection <- select
						}
					  }
					}
					if ((!is.null(select)) & (!is.null(res.mfa$group$coord.sup))) {
					  if (mode(select) == "numeric") 
						selectionS <- select
					  else {
						if (sum(rownames(res.mfa$group$coord) %in% select) + 
							sum(rownames(res.mfa$group$coord.sup) %in% 
								select) != 0) 
						  selectionS <- which(rownames(res.mfa$group$coord.sup) %in% 
												select)
						else {
						  if (grepl("contrib", select)) 
							selectionS <- NULL
						  if (grepl("coord", select)) 
							selectionS <- (rev(order(apply(res.mfa$group$coord.sup[, 
																				   axes]^2, 1, sum))))[1:min(nrow(res.mfa$group$coord.sup), 
																											 sum(as.integer(unlist(strsplit(select, 
																																			"coord"))), na.rm = T))]
						  if (grepl("cos2", select)) {
							if (sum(as.numeric(unlist(strsplit(select, 
															   "cos2"))), na.rm = T) >= 1) 
							  selectionS <- (rev(order(apply(res.mfa$group$cos2.sup[, 
																					axes], 1, sum))))[1:min(nrow(res.mfa$group$coord.sup), 
																											sum(as.numeric(unlist(strsplit(select, 
																																		   "cos2"))), na.rm = T))]
							else selectionS <- which(apply(res.mfa$group$cos2.sup[, 
																				  axes], 1, sum) > sum(as.numeric(unlist(strsplit(select, 
																																  "cos2"))), na.rm = T))
						  }
						  if (is.integer(select)) 
							selectionS <- select
						}
					  }
					}
					if (length(col.hab) == 1) 
					  col.hab = rep(col.hab, length(group))
					if (is.null(col.hab)) {
					  col.hab = rep("darkred", nrow(coord.actif))
					  if (!is.null(res.mfa$group$coord.sup)) 
						col.hab = c(col.hab, rep("darkolivegreen", nrow(coord.illu)))
					}
					if (habillage == "group") 
					  col.hab <- (2:(length(group) + 1))
					if ((new.plot) & !nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) 
					  dev.new()
					if (is.null(palette)) 
					  palette(c("black", "red", "green3", "blue", "magenta", 
								"darkgoldenrod", "darkgreen", "darkgray", "cyan", 
								"violet", "turquoise", "orange", "lightpink", 
								"lavender", "yellow", "lightgreen", "lightgrey", 
								"lightblue", "darkkhaki", "darkmagenta", "darkolivegreen", 
								"lightcyan", "darkorange", "darkorchid", "darkred", 
								"darksalmon", "darkseagreen", "darkslateblue", 
								"darkslategray", "darkslategrey", "darkturquoise", 
								"darkviolet", "lightgray", "lightsalmon", "lightyellow", 
								"maroon"))
					coo <- labe <- coll <- ipch <- fonte <- NULL
					if (is.null(xlim)) 
					  xlim <- c(0, 1)
					if (is.null(ylim)) 
					  ylim <- c(0, 1)
					if (is.null(title)) 
					  title <- "Groups representation"
					plot(0, 0, main = title, xlab = lab.x, ylab = lab.y, 
						 xlim = xlim, ylim = ylim, asp = 1, col = "white", 
						 ...)
					abline(v = 0, lty = 2, ...)
					abline(h = 0, lty = 2, ...)
					coo <- rbind(coo, coord.actif)
					if (lab.grpe) {
					  labe <- c(labe, rownames(coord.actif))
					}
					else labe <- c(labe, rep("", nrow(coord.actif)))
					coll <- c(coll, col.hab[1:nrow(coord.actif)])
					ipch <- c(ipch, rep(17, nrow(coord.actif)))
					fonte <- c(fonte, rep(1, nrow(coord.actif)))
					if (!is.null(selection)) {
					  if (is.numeric(unselect)) 
						coll[!((1:length(coll)) %in% selection)] = rgb(t(col2rgb(coll[!((1:length(coll)) %in% 
																						  selection)])), alpha = 255 * (1 - unselect), 
																	   maxColorValue = 255)
					  else coll[!((1:length(coll)) %in% selection)] = unselect
					  labe[!((1:length(coll)) %in% selection)] <- ""
					}
					if (!is.null(res.mfa$group$coord.sup)) {
					  coo <- rbind(coo, coord.illu)
					  if (lab.grpe) {
						labe2 <- rownames(coord.illu)
					  }
					  else labe2 <- rep("", nrow(coord.illu))
					  coll2 <- col.hab[(nrow(coord.actif) + 1):(nrow(coord.actif) + 
																  nrow(coord.illu))]
					  ipch2 <- rep(2, nrow(coord.illu))
					  fonte2 <- rep(3, nrow(coord.illu))
					  if (length(select) == 1) {
						if (grepl("contrib", select)) {
						  if (is.numeric(unselect)) 
							coll2[1:length(coll2)] = rgb(t(col2rgb(coll2[1:length(coll2)])), 
														 alpha = 255 * (1 - unselect), maxColorValue = 255)
						  else coll2[1:length(coll2)] = unselect
						  labe2[1:length(coll2)] <- ""
						}
					  }
					  if (!is.null(selectionS)) {
						if (is.numeric(unselect)) 
						  coll2[!((1:length(coll2)) %in% selectionS)] = rgb(t(col2rgb(coll2[!((1:length(coll2)) %in% 
																								selectionS)])), alpha = 255 * (1 - unselect), 
																			maxColorValue = 255)
						else coll2[!((1:length(coll2)) %in% selectionS)] = unselect
						labe2[!((1:length(coll2)) %in% selectionS)] <- ""
					  }
					  coll = c(coll, coll2)
					  labe = c(labe, labe2)
					  fonte = c(fonte, fonte2)
					  ipch = c(ipch, ipch2)
					}
					if (shadowtext) 
					  points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, 
							 ...)
					if (autoLab == "auto") 
					  autoLab = (length(labe) < 50)
					if (autoLab == TRUE) 
					  autoLab(coo[, 1], y = coo[, 2], labels = labe, col = coll, 
							  font = fonte, shadotext = shadowtext, ...)
					if (autoLab == FALSE) 
					  text(coo[, 1], y = coo[, 2], labels = labe, col = coll, 
						   font = fonte, pos = 3, ...)
					if (!shadowtext) 
					  points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, 
							 ...)
				  
					#ADDITION
					plotlist[["group"]] <- recordPlot();
					#
				  
				  }
				  if (choix == "var") {
					if ((new.plot) & !nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) 
					  dev.new()
					if (is.null(palette)) 
					  palette(c("black", "red", "green3", "blue", "magenta", 
								"darkgoldenrod", "darkgreen", "darkgray", "cyan", 
								"violet", "turquoise", "orange", "lightpink", 
								"lavender", "yellow", "lightgreen", "lightgrey", 
								"lightblue", "darkkhaki", "darkmagenta", "darkolivegreen", 
								"lightcyan", "darkorange", "darkorchid", "darkred", 
								"darksalmon", "darkseagreen", "darkslateblue", 
								"darkslategray", "darkslategrey", "darkturquoise", 
								"darkviolet", "lightgray", "lightsalmon", "lightyellow", 
								"maroon"))
					test.invisible <- vector(length = 2)
					if (!is.null(invisible)) {
					  test.invisible[1] <- match("quanti", invisible)
					  test.invisible[2] <- match("quanti.sup", invisible)
					}
					else test.invisible <- rep(NA, 2)
					col <- NULL
					if (habillage == "group") {
					  if (is.null(col.hab) | length(col.hab) < length(group[type == 
																			"c"])) {
						if (!is.null(res.mfa$call$num.group.sup)) {
						  col.hab[which(!(1:length(group)) %in% (res.mfa$call$num.group.sup))] <- 2:(1 + 
																									   length(group) - length(res.mfa$call$num.group.sup))
						  col.hab[res.mfa$call$num.group.sup] <- length(group) - 
							length(res.mfa$call$num.group.sup) + 1 + 
							(1:length(res.mfa$call$num.group.sup))
						  col <- c(1 + rep(which(res.mfa$call$nature.group[-res.mfa$call$num.group.sup] == 
												   "quanti"), times = group[which(res.mfa$call$nature.group == 
																					"quanti")]), length(group) - length(res.mfa$call$num.group.sup) + 
									 1 + rep(which((res.mfa$call$nature.group[res.mfa$call$num.group.sup]) == 
													 "quanti.sup"), times = group[which(res.mfa$call$nature.group == 
																						  "quanti.sup")]))
						}
						else {
						  col.hab <- 2:(length(group) + 1)
						  col <- 1 + rep(which(type == "c"), times = group[type == 
																			 "c"])
						}
					  }
					}
					else {
					  if (is.null(col.hab) | length(col.hab) < sum(group[type == 
																		 "c"])) 
						col <- rep(1, sum(group[type == "c"]))
					  else col <- col.hab
					}
					if (is.null(title)) 
					  title <- "Correlation circle"
					plot(0, 0, main = title, xlab = lab.x, ylab = lab.y, 
						 xlim = c(-1.1, 1.1), ylim = c(-1.1, 1.1), col = "white", 
						 asp = 1, ...)
					x.cercle <- seq(-1, 1, by = 0.01)
					y.cercle <- sqrt(1 - x.cercle^2)
					lines(x.cercle, y = y.cercle, ...)
					lines(x.cercle, y = -y.cercle, ...)
					abline(v = 0, lty = 2, ...)
					abline(h = 0, lty = 2, ...)
					if ((!is.null(select)) & (!is.null(res.mfa["quanti.var"]$quanti.var))) {
					  if (mode(select) == "numeric") 
						selection <- select
					  else {
						if (sum(rownames(res.mfa$quanti.var$coord) %in% 
								select) + sum(rownames(res.mfa$quanti.var.sup$coord) %in% 
											  select) != 0) 
						  selection <- which(rownames(res.mfa$quanti.var$coord) %in% 
											   select)
						else {
						  if (grepl("contrib", select)) 
							selection <- (rev(order(res.mfa$quanti.var$contrib[, 
																			   axes[1], drop = FALSE] * res.mfa$eig[axes[1], 
																													1] + res.mfa$quanti.var$contrib[, axes[2], 
																																					drop = FALSE] * res.mfa$eig[axes[2], 1])))[1:min(nrow(res.mfa$quanti.var$coord), 
																																																	 sum(as.integer(unlist(strsplit(select, 
																																																									"contrib"))), na.rm = T))]
						  if (grepl("coord", select)) 
							selection <- (rev(order(apply(res.mfa$quanti.var$coord[, 
																				   axes]^2, 1, sum))))[1:min(nrow(res.mfa$quanti.var$coord), 
																											 sum(as.integer(unlist(strsplit(select, 
																																			"coord"))), na.rm = T))]
						  if (grepl("cos2", select)) {
							if (sum(as.numeric(unlist(strsplit(select, 
															   "cos2"))), na.rm = T) >= 1) 
							  selection <- (rev(order(apply(res.mfa$quanti.var$cos2[, 
																					axes], 1, sum))))[1:min(nrow(res.mfa$quanti.var$coord), 
																											sum(as.numeric(unlist(strsplit(select, 
																																		   "cos2"))), na.rm = T))]
							else selection <- which(apply(res.mfa$quanti.var$cos2[, 
																				  axes], 1, sum) > sum(as.numeric(unlist(strsplit(select, 
																																  "cos2"))), na.rm = T))
						  }
						  if (is.integer(select)) 
							selection <- select
						}
					  }
					}
					if ((!is.null(select)) & (!is.null(res.mfa$quanti.var.sup))) {
					  if (mode(select) == "numeric") 
						selectionS <- select
					  else {
						if (sum(rownames(res.mfa$quanti.var$coord) %in% 
								select) + sum(rownames(res.mfa$quanti.var.sup$coord) %in% 
											  select) != 0) 
						  selectionS <- which(rownames(res.mfa$quanti.var.sup$coord) %in% 
												select)
						else {
						  if (grepl("contrib", select)) 
							selectionS <- NULL
						  if (grepl("coord", select)) 
							selectionS <- (rev(order(apply(res.mfa$quanti.var.sup$coord[, 
																						axes]^2, 1, sum))))[1:min(nrow(res.mfa$quanti.var.sup$coord), 
																												  sum(as.integer(unlist(strsplit(select, 
																																				 "coord"))), na.rm = T))]
						  if (grepl("cos2", select)) {
							if (sum(as.numeric(unlist(strsplit(select, 
															   "cos2"))), na.rm = T) >= 1) 
							  selectionS <- (rev(order(apply(res.mfa$quanti.var.sup$cos2[, 
																						 axes], 1, sum))))[1:min(nrow(res.mfa$quanti.var.sup$coord), 
																												 sum(as.numeric(unlist(strsplit(select, 
																																				"cos2"))), na.rm = T))]
							else selectionS <- which(apply(res.mfa$quanti.var.sup$cos2[, 
																					   axes], 1, sum) > sum(as.numeric(unlist(strsplit(select, 
																																	   "cos2"))), na.rm = T))
						  }
						  if (is.integer(select)) 
							selectionS <- select
						}
					  }
					}
					labe <- labe2 <- coll <- coll2 <- NULL
					if (!is.null(res.mfa["quanti.var"]$quanti.var)) {
					  coll <- col[1:nrow(res.mfa["quanti.var"]$quanti.var$coord)]
					  if (lab.var) 
						labe <- rownames(res.mfa["quanti.var"]$quanti.var$coord)
					  else labe <- rep("", nrow(res.mfa["quanti.var"]$quanti.var$coord))
					}
					if (!is.null(res.mfa$quanti.var.sup)) {
					  if (lab.var) 
						labe2 <- rownames(res.mfa$quanti.var.sup$coord)
					  else labe2 <- rep("", nrow(res.mfa$quanti.var.sup$coord))
					  coll2 <- col[(length(coll) + 1):length(col)]
					}
					if (!is.null(select)) {
					  if (!is.null(res.mfa["quanti.var"]$quanti.var) & 
						  is.na(test.invisible[1])) {
						if (is.numeric(unselect)) 
						  coll[!((1:length(coll)) %in% selection)] = rgb(t(col2rgb(coll[!((1:length(coll)) %in% 
																							selection)])), alpha = 255 * (1 - unselect), 
																		 maxColorValue = 255)
						else coll[!((1:length(coll)) %in% selection)] = unselect
						labe[!((1:length(coll)) %in% selection)] <- ""
					  }
					  if (!is.null(res.mfa$quanti.var.sup) & is.na(test.invisible[2])) {
						if (is.numeric(unselect)) 
						  coll2[!((1:length(coll2)) %in% selectionS)] = rgb(t(col2rgb(coll2[!((1:length(coll2)) %in% 
																								selectionS)])), alpha = 255 * (1 - unselect), 
																			maxColorValue = 255)
						else coll2[!((1:length(coll2)) %in% selectionS)] = unselect
						labe2[!((1:length(coll2)) %in% selectionS)] <- ""
					  }
					}
					col <- c(coll, coll2)
					labe <- c(labe, labe2)
					if (habillage == "group" & is.na(test.invisible[1]) & 
						is.na(test.invisible[2])) {
					  L <- list(x = "topleft", legend = rownames(res.mfa$group$Lg[-nrow(res.mfa$group$Lg), 
																				  , drop = FALSE])[type == "c"], text.col = col.hab[type == 
																																	  "c"], cex = 0.8 * par("cex"))
					  L <- modifyList(L, legend)
					  do.call(graphics::legend, L)
					}
					if (habillage == "group" & is.na(test.invisible[1]) & 
						!is.na(test.invisible[2])) {
					  if ("quanti.sup" %in% res.mfa$call$nature.var) {
						L <- list(x = "topleft", legend = rownames(res.mfa$group$Lg[-c(num.group.sup, 
																					   nrow(res.mfa$group$Lg)), , drop = FALSE])[type.act == 
																																   "c"], text.col = col.hab[which(!((1:length(group)) %in% 
																																									  res.mfa$call$num.group.sup))[type.act == "c"]], 
								  cex = 0.8 * par("cex"))
					  }
					  else {
						L <- list(x = "topleft", legend = rownames(res.mfa$group$Lg[-nrow(res.mfa$group$Lg), 
																					])[type == "c"], text.col = col.hab[type == 
																														  "c"], cex = 0.8 * par("cex"))
					  }
					  L <- modifyList(L, legend)
					  do.call(graphics::legend, L)
					}
					if (habillage == "group" & !is.na(test.invisible[1]) & 
						is.na(test.invisible[2])) {
					  if ("quanti" %in% res.mfa$call$nature.var) 
						L <- list(x = "topleft", legend = rownames(res.mfa$group$Lg[num.group.sup, 
																					, drop = FALSE])[type.sup == "c"], text.col = col.hab[res.mfa$call$num.group.sup[type.sup == 
																																									   "c"]], cex = 0.8 * par("cex"))
					  else L <- list(x = "topleft", legend = rownames(res.mfa$group$Lg[num.group.sup, 
																					   , drop = FALSE])[type.sup == "c"], text.col = col.hab[res.mfa$call$num.group.sup[type.sup == 
																																										  "c"]], cex = 0.8 * par("cex"))
					  L <- modifyList(L, legend)
					  do.call(graphics::legend, L)
					}
					nrow.coord.var <- 0
					coo <- posi <- NULL
					if ((!is.null(res.mfa["quanti.var"]$quanti.var)) & (!is.na(test.invisible[1]))) {
					  col[1:nrow(res.mfa$quanti.var$cor[, axes, drop = FALSE])] <- "transparent"
					}
					if (!is.null(res.mfa["quanti.var"]$quanti.var)) {
					  coord.var <- res.mfa$quanti.var$cor[, axes, drop = FALSE]
					  coo <- coord.var
					  if (length(which(apply(res.mfa$quanti.var$cos2[, 
																	 axes, drop = FALSE], 1, sum, na.rm = TRUE) < 
									   lim.cos2.var)) > 0) 
						col[which(apply(res.mfa$quanti.var$cos2[, axes, 
																drop = FALSE], 1, sum, na.rm = TRUE) < lim.cos2.var)] <- "transparent"
					  for (v in 1:nrow(coord.var)) {
						arrows(0, 0, coord.var[v, 1], coord.var[v, 2], 
							   length = 0.1, angle = 15, code = 2, col = col[v])
						if (lab.var) {
						  if (abs(coord.var[v, 1]) > abs(coord.var[v, 
																   2])) {
							if (coord.var[v, 1] >= 0) 
							  posi <- c(posi, 4)
							else posi <- c(posi, 2)
						  }
						  else {
							if (coord.var[v, 2] >= 0) 
							  posi <- c(posi, 3)
							else posi <- c(posi, 1)
						  }
						}
					  }
					}
					if ((!is.null(res.mfa$quanti.var.sup$coord)) & (!is.na(test.invisible[2]))) {
					  col[nrow(coo) + (1:nrow(res.mfa$quanti.var.sup$cor[, 
																		 axes, drop = FALSE]))] <- "transparent"
					}
					if (!is.null(res.mfa$quanti.var.sup$coord)) {
					  coord.quanti <- res.mfa$quanti.var.sup$cor[, axes, 
																 drop = FALSE]
					  coo <- rbind(coo, coord.quanti)
					  if (length(which(apply(res.mfa$quanti.var.sup$cos2[, 
																		 axes, drop = FALSE], 1, sum, na.rm = TRUE) < 
									   lim.cos2.var)) > 0) 
						col[nrow(coo) - nrow(coord.quanti) + which(apply(res.mfa$quanti.var.sup$cos2[, 
																									 axes, drop = FALSE], 1, sum, na.rm = TRUE) < 
																	 lim.cos2.var)] <- "transparent"
					  for (q in 1:nrow(coord.quanti)) {
						arrows(0, 0, coord.quanti[q, 1], coord.quanti[q, 
																	  2], length = 0.1, angle = 15, code = 2, lty = 2, 
							   col = col[nrow(coo) - nrow(coord.quanti) + 
										   q], ...)
						if (lab.var) {
						  if (abs(coord.quanti[q, 1]) > abs(coord.quanti[q, 
																		 2])) {
							if (coord.quanti[q, 1] >= 0) 
							  posi <- c(posi, 4)
							else posi <- c(posi, 2)
						  }
						  else {
							if (coord.quanti[q, 2] >= 0) 
							  posi <- c(posi, 3)
							else posi <- c(posi, 1)
						  }
						}
					  }
					}
					if (autoLab == "auto") 
					  autoLab = (length(labe) - sum(col == "transparent") < 
								   50)
					if (autoLab == FALSE) 
					  text(coo[, 1], y = coo[, 2], labels = labe, pos = posi, 
						   col = col, ...)
					if (autoLab == TRUE) 
					  autoLab(coo[which(col != "transparent"), 1], y = coo[which(col != 
																				   "transparent"), 2], labels = labe[which(col != 
																															 "transparent")], col = col[which(col != "transparent")], 
							  shadotext = shadowtext, ...)
					par(mar = c(5, 4, 4, 2) + 0.1)
					
					#ADDITION
					plotlist[["var"]] <- recordPlot();
					#
					
				  }
				  if (choix == "freq") {
					if ((new.plot) & !nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) 
					  dev.new()
					if (is.null(palette)) 
					  palette(c("black", "red", "green3", "blue", "magenta", 
								"darkgoldenrod", "darkgreen", "darkgray", "cyan", 
								"violet", "turquoise", "orange", "lightpink", 
								"lavender", "yellow", "lightgreen", "lightgrey", 
								"lightblue", "darkkhaki", "darkmagenta", "darkolivegreen", 
								"lightcyan", "darkorange", "darkorchid", "darkred", 
								"darksalmon", "darkseagreen", "darkslateblue", 
								"darkslategray", "darkslategrey", "darkturquoise", 
								"darkviolet", "lightgray", "lightsalmon", "lightyellow", 
								"maroon"))
					if (is.null(col.hab)) 
					  col.hab = c("black", "grey60", "darkblue", "blue")
					col.row = col.hab[1]
					col.row.sup = col.hab[2]
					col.col = col.hab[3]
					col.col.sup = col.hab[4]
					coord.col <- res.mfa$freq$coord[, axes, drop = FALSE]
					coord.row <- res.mfa$ind$coord[, axes]
					coord.row.sup <- coord.col.sup <- NULL
					if (!is.null(res.mfa$ind.sup)) 
					  coord.row.sup <- res.mfa$ind.sup$coord[, axes, drop = FALSE]
					if (!is.null(res.mfa$freq.sup)) 
					  coord.col.sup <- res.mfa$freq.sup$coord[, axes, 
															  drop = FALSE]
					test.invisible <- vector(length = 4)
					if (!is.null(invisible)) {
					  test.invisible[1] <- match("row", invisible)
					  test.invisible[2] <- match("col", invisible)
					  test.invisible[3] <- match("row.sup", invisible)
					  test.invisible[4] <- match("col.sup", invisible)
					}
					else test.invisible <- rep(NA, 4)
					if (is.null(xlim)) {
					  xmin <- xmax <- 0
					  if (is.na(test.invisible[1])) 
						xmin <- min(xmin, coord.row[, 1])
					  if (is.na(test.invisible[1])) 
						xmax <- max(xmax, coord.row[, 1])
					  if (is.na(test.invisible[2])) 
						xmin <- min(xmin, coord.col[, 1])
					  if (is.na(test.invisible[2])) 
						xmax <- max(xmax, coord.col[, 1])
					  if (is.na(test.invisible[3])) 
						xmin <- min(xmin, coord.row.sup[, 1])
					  if (is.na(test.invisible[3])) 
						xmax <- max(xmax, coord.row.sup[, 1])
					  if (is.na(test.invisible[4])) 
						xmin <- min(xmin, coord.col.sup[, 1])
					  if (is.na(test.invisible[4])) 
						xmax <- max(xmax, coord.col.sup[, 1])
					  xlim <- c(xmin, xmax) * 1.2
					}
					else {
					  xmin = xlim[1]
					  xmax = xlim[2]
					}
					if (is.null(ylim)) {
					  ymin <- ymax <- 0
					  if (is.na(test.invisible[1])) 
						ymin <- min(ymin, coord.row[, 2])
					  if (is.na(test.invisible[1])) 
						ymax <- max(ymax, coord.row[, 2])
					  if (is.na(test.invisible[2])) 
						ymin <- min(ymin, coord.col[, 2])
					  if (is.na(test.invisible[2])) 
						ymax <- max(ymax, coord.col[, 2])
					  if (is.na(test.invisible[3])) 
						ymin <- min(ymin, coord.row.sup[, 2])
					  if (is.na(test.invisible[3])) 
						ymax <- max(ymax, coord.row.sup[, 2])
					  if (is.na(test.invisible[4])) 
						ymin <- min(ymin, coord.col.sup[, 2])
					  if (is.na(test.invisible[4])) 
						ymax <- max(ymax, coord.col.sup[, 2])
					  ylim <- c(ymin, ymax) * 1.2
					}
					else {
					  ymin = ylim[1]
					  ymax = ylim[2]
					}
					col <- NULL
					if (habillage == "group") {
					  if (is.null(col.hab) | length(col.hab) < length(group[type == 
																			"f"])) 
						col.hab <- 2:(length(group[type == "f"]) + 1)
					  for (i in 1:length(group[type == "f"])) col <- c(col, 
																	   rep(col.hab[i], group[type == "f"][i]))
					}
					else {
					  if (is.null(col.hab) | length(col.hab) < sum(group[type == 
																		 "f"])) 
						col <- rep(1, sum(group[type == "f"]))
					  else col <- col.hab
					}
					if (is.null(title)) 
					  titre <- "Factor map for the contingency table(s)"
					else titre <- title
					plot(0, 0, main = titre, xlab = lab.x, ylab = lab.y, 
						 xlim = xlim, ylim = ylim, col = "white", asp = 1, 
						 ...)
					abline(h = 0, lty = 2, ...)
					abline(v = 0, lty = 2, ...)
					selection <- selectionC <- selectionS <- selectionCS <- NULL
					if (!is.null(select)) {
					  if (mode(select) == "numeric") 
						selection <- select
					  else {
						if (sum(rownames(res.mfa$freq.sup$coord) %in% 
								select) + sum(rownames(res.mfa$freq$coord) %in% 
											  select) + sum(rownames(res.mfa$ind$coord) %in% 
															select) + sum(rownames(res.mfa$ind.sup$coord) %in% 
																		  select) != 0) 
						  selection <- which(rownames(res.mfa$ind$coord) %in% 
											   select)
						else {
						  if (grepl("contrib", select)) 
							selection <- (rev(order(res.mfa$ind$contrib[, 
																		axes[1], drop = FALSE] * res.mfa$eig[axes[1], 
																											 1] + res.mfa$ind$contrib[, axes[2], drop = FALSE] * 
													  res.mfa$eig[axes[2], 1])))[1:min(nrow(res.mfa$ind$coord), 
																					   sum(as.integer(unlist(strsplit(select, 
																													  "contrib"))), na.rm = T))]
						  if (grepl("inertia", select)) 
							selection <- (rev(order(apply(res.mfa$ind$within.inertia[, 
																					 axes], 1, sum))))[1:min(nrow(res.mfa$ind$coord), 
																											 sum(as.integer(unlist(strsplit(select, 
																																			"inertia"))), na.rm = T))]
						  if (grepl("coord", select)) 
							selection <- (rev(order(apply(res.mfa$ind$coord[, 
																			axes]^2, 1, sum))))[1:min(nrow(res.mfa$ind$coord), 
																									  sum(as.integer(unlist(strsplit(select, 
																																	 "coord"))), na.rm = T))]
						  if (grepl("cos2", select)) {
							if (sum(as.numeric(unlist(strsplit(select, 
															   "cos2"))), na.rm = T) >= 1) 
							  selection <- (rev(order(apply(res.mfa$ind$cos2[, 
																			 axes], 1, sum))))[1:min(nrow(res.mfa$ind$coord), 
																									 sum(as.numeric(unlist(strsplit(select, 
																																	"cos2"))), na.rm = T))]
							else selection <- which(apply(res.mfa$ind$cos2[, 
																		   axes], 1, sum) > sum(as.numeric(unlist(strsplit(select, 
																														   "cos2"))), na.rm = T))
						  }
						  if (is.integer(select)) 
							selection <- select
						}
					  }
					}
					if ((!is.null(select)) & (!is.null(res.mfa$ind.sup$coord))) {
					  if (mode(select) == "numeric") 
						selectionS <- select
					  else {
						if (sum(rownames(res.mfa$freq.sup$coord) %in% 
								select) + sum(rownames(res.mfa$freq$coord) %in% 
											  select) + sum(rownames(res.mfa$ind$coord) %in% 
															select) + sum(rownames(res.mfa$ind.sup$coord) %in% 
																		  select) != 0) 
						  selectionS <- which(rownames(res.mfa$ind.sup$coord) %in% 
												select)
						else {
						  if (grepl("contrib", select)) 
							selectionS <- NULL
						  if (grepl("inertia", select)) 
							selectionS <- (rev(order(apply(res.mfa$ind.sup$within.inertia[, 
																						  axes]^2, 1, sum))))[1:min(nrow(res.mfa$ind.sup$coord), 
																													sum(as.integer(unlist(strsplit(select, 
																																				   "inertia"))), na.rm = T))]
						  if (grepl("coord", select)) 
							selectionS <- (rev(order(apply(res.mfa$ind.sup$coord[, 
																				 axes]^2, 1, sum))))[1:min(nrow(res.mfa$ind.sup$coord), 
																										   sum(as.integer(unlist(strsplit(select, 
																																		  "coord"))), na.rm = T))]
						  if (grepl("cos2", select)) {
							if (sum(as.numeric(unlist(strsplit(select, 
															   "cos2"))), na.rm = T) >= 1) 
							  selectionS <- (rev(order(apply(res.mfa$ind.sup$cos2[, 
																				  axes], 1, sum))))[1:min(nrow(res.mfa$ind.sup$coord), 
																										  sum(as.numeric(unlist(strsplit(select, 
																																		 "cos2"))), na.rm = T))]
							else selectionS <- which(apply(res.mfa$ind.sup$cos2[, 
																				axes], 1, sum) > sum(as.numeric(unlist(strsplit(select, 
																																"cos2"))), na.rm = T))
						  }
						  if (is.integer(select)) 
							selectionS <- select
						}
					  }
					}
					if ((!is.null(select)) & (!is.null(res.mfa$freq$coord))) {
					  if (mode(select) == "numeric") 
						selectionC <- select
					  else {
						if (sum(rownames(res.mfa$freq.sup$coord) %in% 
								select) + sum(rownames(res.mfa$freq$coord) %in% 
											  select) + sum(rownames(res.mfa$ind$coord) %in% 
															select) + sum(rownames(res.mfa$ind.sup$coord) %in% 
																		  select) != 0) 
						  selectionC <- which(rownames(res.mfa$freq$coord) %in% 
												select)
						else {
						  if (grepl("contrib", select)) 
							selectionC <- (rev(order(res.mfa$freq$contrib[, 
																		  axes[1], drop = FALSE] * res.mfa$eig[axes[1], 
																											   1] + res.mfa$freq$contrib[, axes[2], drop = FALSE] * 
													   res.mfa$eig[axes[2], 1])))[1:min(nrow(res.mfa$freq$coord), 
																						sum(as.integer(unlist(strsplit(select, 
																													   "contrib"))), na.rm = T))]
						  if (grepl("coord", select)) 
							selectionC <- (rev(order(apply(res.mfa$freq$coord[, 
																			  axes]^2, 1, sum))))[1:min(nrow(res.mfa$freq$coord), 
																										sum(as.integer(unlist(strsplit(select, 
																																	   "coord"))), na.rm = T))]
						  if (grepl("cos2", select)) {
							if (sum(as.numeric(unlist(strsplit(select, 
															   "cos2"))), na.rm = T) >= 1) 
							  selectionC <- (rev(order(apply(res.mfa$freq$cos2[, 
																			   axes], 1, sum))))[1:min(nrow(res.mfa$freq$coord), 
																									   sum(as.numeric(unlist(strsplit(select, 
																																	  "cos2"))), na.rm = T))]
							else selectionC <- which(apply(res.mfa$freq$cos2[, 
																			 axes], 1, sum) > sum(as.numeric(unlist(strsplit(select, 
																															 "cos2"))), na.rm = T))
						  }
						  if (is.integer(select)) 
							selectionC <- select
						}
					  }
					}
					if ((!is.null(select)) & (!is.null(res.mfa$freq.sup$coord))) {
					  if (mode(select) == "numeric") 
						selectionCS <- select
					  else {
						if (sum(rownames(res.mfa$freq.sup$coord) %in% 
								select) + sum(rownames(res.mfa$freq$coord) %in% 
											  select) + sum(rownames(res.mfa$ind$coord) %in% 
															select) + sum(rownames(res.mfa$ind.sup$coord) %in% 
																		  select) != 0) 
						  selectionCS <- which(rownames(res.mfa$freq.sup$coord) %in% 
												 select)
						else {
						  if (grepl("contrib", select)) 
							selectionCS <- NULL
						  if (grepl("coord", select)) 
							selectionCS <- (rev(order(apply(res.mfa$freq.sup$coord[, 
																				   axes]^2, 1, sum))))[1:min(nrow(res.mfa$freq.sup$coord), 
																											 sum(as.integer(unlist(strsplit(select, 
																																			"coord"))), na.rm = T))]
						  if (grepl("cos2", select)) {
							if (sum(as.numeric(unlist(strsplit(select, 
															   "cos2"))), na.rm = T) >= 1) 
							  selectionCS <- (rev(order(apply(res.mfa$freq.sup$cos2[, 
																					axes], 1, sum))))[1:min(nrow(res.mfa$freq.sup$coord), 
																											sum(as.numeric(unlist(strsplit(select, 
																																		   "cos2"))), na.rm = T))]
							else selectionCS <- which(apply(res.mfa$freq.sup$cos2[, 
																				  axes], 1, sum) > sum(as.numeric(unlist(strsplit(select, 
																																  "cos2"))), na.rm = T))
						  }
						  if (is.integer(select)) 
							selectionCS <- select
						}
					  }
					}
					coo <- labe <- coll <- ipch <- fonte <- NULL
					if (is.na(test.invisible[1])) {
					  coo <- rbind(coo, coord.row)
					  if (lab.ind) {
						labe <- rownames(coord.row)
					  }
					  else labe <- rep("", nrow(coord.row))
					  coll <- rep(col.row, nrow(coord.row))
					  ipch <- c(ipch, rep(20, nrow(coord.row)))
					  fonte <- c(fonte, rep(1, nrow(coord.row)))
					  if (!is.null(selection)) {
						if (is.numeric(unselect)) 
						  coll[!((1:length(coll)) %in% selection)] = rgb(t(col2rgb(coll[!((1:length(coll)) %in% 
																							selection)])), alpha = 255 * (1 - unselect), 
																		 maxColorValue = 255)
						else coll[!((1:length(coll)) %in% selection)] = unselect
						labe[!((1:length(coll)) %in% selection)] <- ""
					  }
					}
					if (is.na(test.invisible[2])) {
					  coo <- rbind(coo, coord.col)
					  if (lab.ind) {
						labe2 <- rownames(coord.col)
					  }
					  else labe2 <- rep("", nrow(coord.col))
					  coll2 <- rep(col.col, nrow(coord.col))
					  ipch <- c(ipch, rep(17, nrow(coord.col)))
					  fonte <- c(fonte, rep(1, nrow(coord.col)))
					  if (!is.null(selectionC)) {
						if (is.numeric(unselect)) 
						  coll2[!((1:length(coll2)) %in% selectionC)] = rgb(t(col2rgb(coll2[!((1:length(coll2)) %in% 
																								selectionC)])), alpha = 255 * (1 - unselect), 
																			maxColorValue = 255)
						else coll2[!((1:length(coll2)) %in% selectionC)] = unselect
						labe2[!((1:length(coll2)) %in% selectionC)] <- ""
					  }
					  coll <- c(coll, coll2)
					  labe <- c(labe, labe2)
					}
					if (!is.null(res.mfa$freq.sup) & is.na(test.invisible[4])) {
					  coo <- rbind(coo, coord.col.sup)
					  if (lab.ind) {
						labe2 <- rownames(coord.col.sup)
					  }
					  else labe2 <- rep("", nrow(coord.col.sup))
					  coll2 <- rep(col.col.sup, nrow(coord.col.sup))
					  ipch <- c(ipch, rep(17, nrow(coord.col.sup)))
					  fonte <- c(fonte, rep(1, nrow(coord.col.sup)))
					  if (!is.null(selectionCS)) {
						if (is.numeric(unselect)) 
						  coll2[!((1:length(coll2)) %in% selectionCS)] = rgb(t(col2rgb(coll2[!((1:length(coll2)) %in% 
																								 selectionCS)])), alpha = 255 * (1 - unselect), 
																			 maxColorValue = 255)
						else coll2[!((1:length(coll2)) %in% selectionCS)] = unselect
						labe2[!((1:length(coll2)) %in% selectionCS)] <- ""
					  }
					  coll <- c(coll, coll2)
					  labe <- c(labe, labe2)
					}
					if (!is.null(res.mfa$ind.sup) & is.na(test.invisible[3])) {
					  coo <- rbind(coo, coord.row.sup)
					  if (lab.ind) {
						labe2 <- rownames(coord.row.sup)
					  }
					  else labe2 <- rep("", nrow(coord.row.sup))
					  coll2 <- rep(col.row.sup, nrow(coord.row.sup))
					  ipch <- c(ipch, rep(17, nrow(coord.row.sup)))
					  fonte <- c(fonte, rep(1, nrow(coord.row.sup)))
					  if (!is.null(selectionS)) {
						if (is.numeric(unselect)) 
						  coll2[!((1:length(coll2)) %in% selectionS)] = rgb(t(col2rgb(coll2[!((1:length(coll2)) %in% 
																								selectionS)])), alpha = 255 * (1 - unselect), 
																			maxColorValue = 255)
						else coll2[!((1:length(coll2)) %in% selectionS)] = unselect
						labe2[!((1:length(coll2)) %in% selectionS)] <- ""
					  }
					  if (length(select) == 1) {
						if (grepl("contrib", select)) {
						  if (is.numeric(unselect)) 
							coll2[1:length(coll2)] = rgb(t(col2rgb(coll2[1:length(coll2)])), 
														 alpha = 255 * (1 - unselect), maxColorValue = 255)
						  else coll2[1:length(coll2)] = unselect
						  labe2[1:length(coll2)] <- ""
						}
					  }
					  coll <- c(coll, coll2)
					  labe <- c(labe, labe2)
					}
					if (shadowtext) 
					  points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, 
							 ...)
					if (any(labe != "")) {
					  if (autoLab == "auto") 
						autoLab = (length(which(labe != "")) < 50)
					  if (autoLab == TRUE) 
						autoLab(coo[labe != "", 1], y = coo[labe != 
															  "", 2], labels = labe[labe != ""], col = coll[labe != 
																											  ""], font = fonte[labe != ""], shadotext = shadowtext, 
								...)
					  if (autoLab == FALSE) 
						text(coo[labe != "", 1], y = coo[labe != "", 
														 2], labels = labe[labe != ""], col = coll[labe != 
																									 ""], font = fonte[labe != ""], pos = 3, ...)
					}
					if (!shadowtext) 
					  points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, 
							 ...)
					if (habillage == "group") {
					  L <- list(x = "topleft", legend = rownames(res.mfa$group$Lg[-nrow(res.mfa$group$Lg), 
																				  ])[type == "f"], text.col = col.hab, cex = 0.8 * 
								  par("cex"))
					  L <- modifyList(L, legend)
					  do.call(graphics::legend, L)
					}
					
					#ADDITION
					plotlist[["freq"]] <- recordPlot();
					#
					
				  }
				  if (choix == "ind") {
					test.invisible <- vector(length = 3)
					if (!is.null(invisible)) {
					  test.invisible[1] <- match("ind", invisible)
					  test.invisible[2] <- match("ind.sup", invisible)
					  test.invisible[3] <- match("quali", invisible)
					}
					else test.invisible <- rep(NA, 3)
					nb.ind.actif <- nrow(res.mfa$ind$coord)
					nb.ind.illu <- 0
					if (!is.null(res.mfa$ind.sup)) 
					  nb.ind.illu <- nrow(res.mfa$ind.sup$coord)
					nb.ind <- nb.ind.actif + nb.ind.illu
					coord.ind <- res.mfa$ind$coord[, axes, drop = FALSE]
					coord.ind.partiel <- res.mfa$ind$coord.partiel[, axes, 
																   drop = FALSE]
					coord.ind.sup <- NULL
					if (!is.null(res.mfa$ind.sup)) {
					  coord.ind.sup <- res.mfa$ind.sup$coord[, axes, drop = FALSE]
					  coord.ind.partiel.sup <- res.mfa$ind.sup$coord.partiel[, 
																			 axes, drop = FALSE]
					}
					coord.quali <- coord.quali.sup <- coord.quali.partiel <- coord.quali.sup.partiel <- NULL
					nrow.coord.quali <- 0
					if (!is.null(res.mfa["quali.var"]$quali.var)) {
					  coord.quali <- res.mfa$quali.var$coord[, axes, drop = FALSE]
					  coord.quali.partiel <- res.mfa$quali.var$coord.partiel[, 
																			 axes, drop = FALSE]
					  nrow.coord.quali <- nrow(coord.quali)
					}
					if (!is.null(res.mfa["quali.var.sup"])) {
					  coord.quali.sup <- res.mfa$quali.var.sup$coord[, 
																	 axes, drop = FALSE]
					  coord.quali.partiel.sup <- res.mfa$quali.var.sup$coord.partiel[, 
																					 axes, drop = FALSE]
					}
					group.ind.actif <- group.ind.sup <- group.quali <- group.quali.sup <- NULL
					if (!is.null(partial)) {
					  if (length(partial) == 1) {
						if (partial == "all") {
						  group.ind.actif <- 1:nrow(coord.ind)
						  if (!is.null(res.mfa$ind.sup)) 
							group.ind.sup <- 1:nrow(coord.ind.sup)
						  if (!is.null(res.mfa["quali.var"]$quali.var)) 
							group.quali <- 1:nrow(coord.quali)
						  if (!is.null(res.mfa["quali.var.sup"]$quali.var.sup)) 
							group.quali.sup <- 1:nrow(coord.quali.sup)
						}
						else {
						  for (i in 1:length(partial)) {
							if (partial[i] %in% rownames(coord.ind)) 
							  group.ind.actif <- c(group.ind.actif, 
												   match(partial[i], rownames(coord.ind)))
							if (partial[i] %in% rownames(coord.ind.sup)) 
							  group.ind.sup <- c(group.ind.sup, match(partial[i], 
																	  rownames(coord.ind.sup)))
							if (partial[i] %in% rownames(coord.quali)) 
							  group.quali <- c(group.quali, match(partial[i], 
																  rownames(coord.quali)))
							if (partial[i] %in% rownames(coord.quali.sup)) 
							  group.quali.sup <- c(group.quali.sup, 
												   match(partial[i], rownames(coord.quali.sup)))
						  }
						}
					  }
					  else {
						for (i in 1:length(partial)) {
						  if (partial[i] %in% rownames(coord.ind)) 
							group.ind.actif <- c(group.ind.actif, match(partial[i], 
																		rownames(coord.ind)))
						  if (partial[i] %in% rownames(coord.ind.sup)) 
							group.ind.sup <- c(group.ind.sup, match(partial[i], 
																	rownames(coord.ind.sup)))
						  if (partial[i] %in% rownames(coord.quali)) 
							group.quali <- c(group.quali, match(partial[i], 
																rownames(coord.quali)))
						  if (partial[i] %in% rownames(coord.quali.sup)) 
							group.quali.sup <- c(group.quali.sup, match(partial[i], 
																		rownames(coord.quali.sup)))
						}
					  }
					}
					if (!is.null(ellipse)) {
					  coord.ellipse <- ellipse$res
					  npoint.ellipse <- ellipse$call
					}
					else coord.ellipse <- NULL
					if (!is.null(ellipse.par)) {
					  coord.ellipse.par <- ellipse.par$res
					  npoint.ellipse.par <- ellipse.par$call
					}
					else coord.ellipse.par <- NULL
					if (is.null(xlim)) {
					  xmin <- xmax <- 0
					  if (is.na(test.invisible[1])) 
						xmin <- min(xmin, coord.ind[, 1])
					  if (is.na(test.invisible[1])) 
						xmax <- max(xmax, coord.ind[, 1])
					  if (!is.null(res.mfa$ind.sup) & is.na(test.invisible[2])) 
						xmin <- min(xmin, coord.ind.sup[, 1])
					  if (!is.null(res.mfa$ind.sup) & is.na(test.invisible[2])) 
						xmax <- max(xmax, coord.ind.sup[, 1])
					  if (is.na(test.invisible[1])) 
						xmin <- min(xmin, coord.ind.partiel[unlist(lapply(group.ind.actif, 
																		  function(k) seq(nbre.grpe * (k - 1) + 1, length = nbre.grpe))), 
															1])
					  if (is.na(test.invisible[1])) 
						xmax <- max(xmax, coord.ind.partiel[unlist(lapply(group.ind.actif, 
																		  function(k) seq(nbre.grpe * (k - 1) + 1, length = nbre.grpe))), 
															1])
					  if (!is.null(res.mfa$ind.sup) & is.na(test.invisible[2])) 
						xmin <- min(xmin, coord.ind.partiel.sup[unlist(lapply(group.ind.sup, 
																			  function(k) seq(nbre.grpe * (k - 1) + 1, length = nbre.grpe))), 
																1])
					  if (!is.null(res.mfa$ind.sup) & is.na(test.invisible[2])) 
						xmax <- max(xmax, coord.ind.partiel.sup[unlist(lapply(group.ind.sup, 
																			  function(k) seq(nbre.grpe * (k - 1) + 1, length = nbre.grpe))), 
																1])
					  if (!is.null(res.mfa["quali.var"]$quali.var) & is.na(test.invisible[3])) 
						xmin <- min(xmin, coord.quali[, 1])
					  if (!is.null(res.mfa["quali.var"]$quali.var) & is.na(test.invisible[3])) 
						xmax <- max(xmax, coord.quali[, 1])
					  if (!is.null(res.mfa["quali.var"]$quali.var) & is.na(test.invisible[3])) 
						xmin <- min(xmin, coord.quali.partiel[unlist(lapply(group.quali, 
																			function(k) seq(nbre.grpe * (k - 1) + 1, length = nbre.grpe))), 
															  1])
					  if (!is.null(res.mfa["quali.var"]$quali.var) & is.na(test.invisible[3])) 
						xmax <- max(xmax, coord.quali.partiel[unlist(lapply(group.quali, 
																			function(k) seq(nbre.grpe * (k - 1) + 1, length = nbre.grpe))), 
															  1])
					  if (!is.null(res.mfa$quali.var.sup) & is.na(test.invisible[3])) 
						xmin <- min(xmin, coord.quali[, 1], coord.quali.sup[, 
																			1])
					  if (!is.null(res.mfa$quali.var.sup) & is.na(test.invisible[3])) 
						xmax <- max(xmax, coord.quali[, 1], coord.quali.sup[, 
																			1])
					  if (!is.null(res.mfa$quali.var.sup) & is.na(test.invisible[3])) 
						xmin <- min(xmin, coord.quali.partiel.sup[unlist(lapply(group.quali.sup, 
																				function(k) seq(nbre.grpe * (k - 1) + 1, length = nbre.grpe))), 
																  1])
					  if (!is.null(res.mfa$quali.var.sup) & is.na(test.invisible[3])) 
						xmax <- max(xmax, coord.quali.partiel.sup[unlist(lapply(group.quali.sup, 
																				function(k) seq(nbre.grpe * (k - 1) + 1, length = nbre.grpe))), 
																  1])
					  xlim <- c(xmin, xmax) * 1.1
					}
					else {
					  xmin = xlim[1]
					  xmax = xlim[2]
					}
					if (is.null(ylim)) {
					  ymin <- ymax <- 0
					  if (is.na(test.invisible[1])) 
						ymin <- min(ymin, coord.ind[, 2])
					  if (is.na(test.invisible[1])) 
						ymax <- max(ymax, coord.ind[, 2])
					  if (!is.null(res.mfa$ind.sup) & is.na(test.invisible[2])) 
						ymin <- min(ymin, coord.ind.sup[, 2])
					  if (!is.null(res.mfa$ind.sup) & is.na(test.invisible[2])) 
						ymax <- max(ymax, coord.ind.sup[, 2])
					  if (is.na(test.invisible[1])) 
						ymin <- min(ymin, coord.ind.partiel[unlist(lapply(group.ind.actif, 
																		  function(k) seq(nbre.grpe * (k - 1) + 1, length = nbre.grpe))), 
															2])
					  if (is.na(test.invisible[1])) 
						ymax <- max(ymax, coord.ind.partiel[unlist(lapply(group.ind.actif, 
																		  function(k) seq(nbre.grpe * (k - 1) + 1, length = nbre.grpe))), 
															2])
					  if (!is.null(res.mfa$ind.sup) & is.na(test.invisible[2])) 
						ymin <- min(ymin, coord.ind.partiel.sup[unlist(lapply(group.ind.sup, 
																			  function(k) seq(nbre.grpe * (k - 1) + 1, length = nbre.grpe))), 
																2])
					  if (!is.null(res.mfa$ind.sup) & is.na(test.invisible[2])) 
						ymax <- max(ymax, coord.ind.partiel.sup[unlist(lapply(group.ind.sup, 
																			  function(k) seq(nbre.grpe * (k - 1) + 1, length = nbre.grpe))), 
																2])
					  if (!is.null(res.mfa["quali.var"]$quali.var) & is.na(test.invisible[3])) 
						ymin <- min(ymin, coord.quali[, 2])
					  if (!is.null(res.mfa["quali.var"]$quali.var) & is.na(test.invisible[3])) 
						ymax <- max(ymax, coord.quali[, 2])
					  if (!is.null(res.mfa["quali.var"]$quali.var) & is.na(test.invisible[3])) 
						ymin <- min(ymin, coord.quali.partiel[unlist(lapply(group.quali, 
																			function(k) seq(nbre.grpe * (k - 1) + 1, length = nbre.grpe))), 
															  2])
					  if (!is.null(res.mfa["quali.var"]$quali.var) & is.na(test.invisible[3])) 
						ymax <- max(ymax, coord.quali.partiel[unlist(lapply(group.quali, 
																			function(k) seq(nbre.grpe * (k - 1) + 1, length = nbre.grpe))), 
															  2])
					  if (!is.null(res.mfa$quali.var.sup) & is.na(test.invisible[3])) 
						ymin <- min(ymin, coord.quali[, 1], coord.quali.sup[, 
																			2])
					  if (!is.null(res.mfa$quali.var.sup) & is.na(test.invisible[3])) 
						ymax <- max(ymax, coord.quali[, 1], coord.quali.sup[, 
																			2])
					  if (!is.null(res.mfa$quali.var.sup) & is.na(test.invisible[3])) 
						ymin <- min(ymin, coord.quali.partiel.sup[unlist(lapply(group.quali.sup, 
																				function(k) seq(nbre.grpe * (k - 1) + 1, length = nbre.grpe))), 
																  2])
					  if (!is.null(res.mfa$quali.var.sup) & is.na(test.invisible[3])) 
						ymax <- max(ymax, coord.quali.partiel.sup[unlist(lapply(group.quali.sup, 
																				function(k) seq(nbre.grpe * (k - 1) + 1, length = nbre.grpe))), 
																  2])
					  ylim <- c(ymin, ymax) * 1.1
					}
					else {
					  ymin = ylim[1]
					  ymax = ylim[2]
					}
					selection <- NULL
					if (!is.null(select)) {
					  if (mode(select) == "numeric") 
						selection <- select
					  else {
						if (sum(rownames(res.mfa$ind$coord) %in% select) != 
							0) 
						  selection <- which(rownames(res.mfa$ind$coord) %in% 
											   select)
						else {
						  if (grepl("contrib", select)) 
							selection <- (rev(order(res.mfa$ind$contrib[, 
																		axes[1], drop = FALSE] * res.mfa$eig[axes[1], 
																											 1] + res.mfa$ind$contrib[, axes[2], drop = FALSE] * 
													  res.mfa$eig[axes[2], 1])))[1:min(nrow(res.mfa$ind$coord), 
																					   sum(as.integer(unlist(strsplit(select, 
																													  "contrib"))), na.rm = T))]
						  if (grepl("dist", select)) 
							selection <- (rev(order(res.mfa$ind$dist)))[1:min(nrow(res.mfa$ind$coord), 
																			  sum(as.integer(unlist(strsplit(select, 
																											 "dist"))), na.rm = T))]
						  if (grepl("coord", select)) 
							selection <- (rev(order(apply(res.mfa$ind$coord[, 
																			axes]^2, 1, sum))))[1:min(nrow(res.mfa$ind$coord), 
																									  sum(as.integer(unlist(strsplit(select, 
																																	 "coord"))), na.rm = T))]
						  if (grepl("cos2", select)) {
							if (sum(as.numeric(unlist(strsplit(select, 
															   "cos2"))), na.rm = T) >= 1) 
							  selection <- (rev(order(apply(res.mfa$ind$cos2[, 
																			 axes], 1, sum))))[1:min(nrow(res.mfa$ind$coord), 
																									 sum(as.numeric(unlist(strsplit(select, 
																																	"cos2"))), na.rm = T))]
							else selection <- which(apply(res.mfa$ind$cos2[, 
																		   axes], 1, sum) > sum(as.numeric(unlist(strsplit(select, 
																														   "cos2"))), na.rm = T))
						  }
						  if (is.integer(select)) 
							selection <- select
						}
					  }
					}
					if (habillage == "group") {
					  if (is.null(col.hab) | length(col.hab) != (nbre.grpe)) 
						col.hab <- 2:(nbre.grpe + 1)
					  col.ind <- c(rep(1, nb.ind.actif), rep(col.hab, 
															 nb.ind.actif))
					  if (!is.null(res.mfa$ind.sup)) 
						col.ind.sup <- c(rep(1, nb.ind - nb.ind.actif), 
										 rep(col.hab, nb.ind - nb.ind.actif))
					  if (length(group[type == "n"]) != 0) 
						col.quali <- c(rep(1, sum(res.mfa$call$group.mod[type == 
																		   "n"])), rep(col.hab, sum(res.mfa$call$group.mod[type == 
																															 "n"])))
					  if (!is.null(res.mfa$quali.var.sup)) 
						col.quali.sup <- c(rep(1, sum(res.mfa$call$group.mod[num.group.sup][type.sup == 
																							  "n"])), rep(col.hab, sum(res.mfa$call$group.mod[num.group.sup][type.sup == 
																																							   "n"])))
					  if (!is.null(ellipse)) 
						col.ellipse <- rep(1, nb.ind.actif)
					  if (!is.null(ellipse.par)) 
						col.ellipse.par <- rep(col.hab, nb.ind.actif)
					}
					if (habillage == "ind") {
					  if (is.null(col.hab) | length(col.hab) != nb.ind) 
						col.hab <- 1:nb.ind
					  col.ind <- c(col.hab[1:nb.ind.actif], rep(col.hab[1:nb.ind.actif], 
																each = nbre.grpe))
					  if (!is.null(res.mfa$ind.sup)) 
						col.ind.sup <- c(col.hab[(nb.ind.actif + 1):nb.ind], 
										 rep(col.hab[(nb.ind.actif + 1):nb.ind], each = nbre.grpe))
					  if (length(group[type == "n"]) != 0) 
						col.quali <- col.quali.sup <- rep("black", (1 + 
																	  nbre.grpe) * sum(res.mfa$call$group.mod[type == 
																												"n"]))
					  if (!is.null(ellipse)) 
						col.ellipse <- col.hab[1:nb.ind.actif]
					  if (!is.null(ellipse.par)) 
						col.ellipse.par <- rep(col.hab[1:nb.ind.actif], 
											   each = nbre.grpe)
					}
					if ((habillage != "none") & (habillage != "ind") & (habillage != 
																		"group")) {
					  group.act <- (1:length(group))
					  if (!is.null(num.group.sup)) 
						group.act <- group.act[-num.group.sup]
					  nbre.modalite <- nbre.modalite.sup <- NULL
					  liste.quali <- liste.quali.sup <- NULL
					  for (i in group.act) {
						if (type[i] == "n") {
						  for (k in 1:ncol(res.mfa$separate.analyses[[i]]$call$X)) nbre.modalite <- c(nbre.modalite, 
																									  nlevels(res.mfa$separate.analyses[[i]]$call$X[, 
																																					k]))
						  if (i == 1) 
							liste.quali <- c(liste.quali, colnames(res.mfa$call$X[1:group[1]]))
						  else liste.quali <- c(liste.quali, colnames(res.mfa$call$X[(sum(group[1:(i - 
																									 1)]) + 1):sum(group[1:i])]))
						}
					  }
					  if (!is.null(num.group.sup)) {
						for (i in num.group.sup) {
						  if (type[i] == "n") {
							if (i == 1) 
							  liste.quali.sup <- c(liste.quali.sup, 
												   colnames(res.mfa$call$X[1:group[1]]))
							else liste.quali.sup <- c(liste.quali.sup, 
													  colnames(res.mfa$call$X[(sum(group[1:(i - 
																							  1)]) + 1):sum(group[1:i])]))
							for (k in 1:ncol(res.mfa$separate.analyses[[i]]$call$X)) nbre.modalite.sup <- c(nbre.modalite.sup, 
																											nlevels(res.mfa$separate.analyses[[i]]$call$X[, 
																																						  k]))
						  }
						}
					  }
					  if (is.double(habillage)) 
						nom.quali <- colnames(res.mfa$call$X)[habillage]
					  else nom.quali = habillage
					  if (!(nom.quali %in% c(liste.quali, liste.quali.sup))) 
						stop("The variable ", habillage, " is not qualitative")
					  modalite <- levels(as.factor(res.mfa$call$X[, nom.quali]))
					  col.ind <- as.numeric(as.factor(res.mfa$call$X[, 
																	 nom.quali]))
					  if (is.null(col.hab) | length(col.hab) != length(modalite)) 
						col.hab <- 2:(1 + length(modalite))
					  col.ind <- col.hab[col.ind]
					  if (!is.null(res.mfa$call$ind.sup)) {
						col.ind.sup <- col.ind[res.mfa$call$ind.sup]
						col.ind <- col.ind[-res.mfa$call$ind.sup]
						col.ind.sup <- c(col.ind.sup, rep(col.ind.sup, 
														  each = nbre.grpe))
					  }
					  col.ind <- c(col.ind, rep(col.ind, each = nbre.grpe))
					  col.ellipse <- col.ind[1:nb.ind.actif]
					  col.ellipse.par <- col.ind[-c(1:nb.ind.actif)]
					  col.quali <- rep("black", sum(res.mfa$call$group.mod[type == 
																			 "n"]))
					  if (nom.quali %in% liste.quali) {
						indice.inf <- sum(nbre.modalite[0:(match(nom.quali, 
																 liste.quali) - 1)]) + 1
						indice.sup <- indice.inf + length(modalite) - 
						  1
						if (length(group[type == "n"]) != 0) {
						  for (i in 1:length(liste.quali)) {
							if (liste.quali[i] == nom.quali) 
							  col.quali[indice.inf:indice.sup] <- col.hab
						  }
						}
					  }
					  col.quali <- c(col.quali, rep(col.quali, each = nbre.grpe))
					  col.quali.sup <- rep("black", sum(res.mfa$call$group.mod[(type == 
																				  "n") %in% num.group.sup]))
					  if (nom.quali %in% liste.quali.sup) {
						indice.inf.sup <- sum(nbre.modalite.sup[0:(match(nom.quali, 
																		 liste.quali.sup) - 1)]) + 1
						indice.sup.sup <- indice.inf.sup + length(modalite) - 
						  1
						if (length(group[type == "n"]) != 0) {
						  for (i in 1:length(liste.quali.sup)) {
							if (liste.quali.sup[i] == nom.quali) 
							  col.quali.sup[indice.inf.sup:indice.sup.sup] <- col.hab
						  }
						}
					  }
					  col.quali.sup <- c(col.quali.sup, rep(col.quali.sup, 
															each = nbre.grpe))
					}
					if (habillage == "none") 
					  col.ind <- col.ind.sup <- col.quali.sup <- col.quali <- col.ellipse <- col.ellipse.par <- rep("black", 
																													nb.ind * (nbre.grpe + 1))
					if ((new.plot) & !nzchar(Sys.getenv("RSTUDIO_USER_IDENTITY"))) 
					  dev.new(width = min(14, max(8, 8 * (xmax - xmin)/(ymax - 
																		  ymin))), height = 8)
					if (is.null(palette)) 
					  palette(c("black", "red", "green3", "blue", "magenta", 
								"darkgoldenrod", "darkgreen", "darkgray", "cyan", 
								"violet", "turquoise", "orange", "lightpink", 
								"lavender", "yellow", "lightgreen", "lightgrey", 
								"lightblue", "darkkhaki", "darkmagenta", "darkolivegreen", 
								"lightcyan", "darkorange", "darkorchid", "darkred", 
								"darksalmon", "darkseagreen", "darkslateblue", 
								"darkslategray", "darkslategrey", "darkturquoise", 
								"darkviolet", "lightgray", "lightsalmon", "lightyellow", 
								"maroon"))
					if (is.null(title)) 
					  title <- "Individual factor map"
					plot(0, 0, main = title, xlab = lab.x, ylab = lab.y, 
						 xlim = xlim, ylim = ylim, col = "white", asp = 1, 
						 ...)
					abline(v = 0, lty = 2, ...)
					abline(h = 0, lty = 2, ...)
					coo <- labe <- coll <- ipch <- fonte <- NULL
					if (is.na(test.invisible[1])) {
					  coo <- rbind(coo, coord.ind)
					  if (lab.ind) {
						labe <- c(labe, rownames(coord.ind))
					  }
					  else labe <- c(labe, rep("", nrow(coord.ind)))
					  coll <- c(coll, col.ind[1:nb.ind.actif])
					  ipch <- c(ipch, rep(20, nrow(coord.ind)))
					  fonte <- c(fonte, rep(1, nrow(coord.ind)))
					  if (!is.null(selection)) {
						if (is.numeric(unselect)) 
						  coll[!((1:length(coll)) %in% selection)] = rgb(t(col2rgb(coll[!((1:length(coll)) %in% 
																							selection)])), alpha = 255 * (1 - unselect), 
																		 maxColorValue = 255)
						else coll[!((1:length(coll)) %in% selection)] = unselect
						labe[!((1:length(coll)) %in% selection)] <- ""
					  }
					  for (i in group.ind.actif) {
						if (col2rgb(coll[i], alpha = TRUE)[4] == 255) {
						  for (j in 1:nbre.grpe) {
							points(coord.ind.partiel[(i - 1) * nbre.grpe + 
													   j, 1], coord.ind.partiel[(i - 1) * nbre.grpe + 
																				  j, 2], cex = 0.8 * par("cex"), col = col.ind[nb.ind.actif + 
																																 (i - 1) * nbre.grpe + j], pch = 20)
							if (lab.par) 
							  text(coord.ind.partiel[(i - 1) * nbre.grpe + 
													   j, 1], y = coord.ind.partiel[(i - 1) * 
																					  nbre.grpe + j, 2], labels = rownames(coord.ind.partiel)[(i - 
																																				 1) * nbre.grpe + j], pos = 3, col = col.ind[nb.ind.actif + 
																																															   (i - 1) * nbre.grpe + j], ...)
							if (chrono) {
							  if (j > 1) 
								lines(c(coord.ind.partiel[(i - 1) * 
															nbre.grpe + (j - 1), 1], coord.ind.partiel[(i - 
																										  1) * nbre.grpe + j, 1]), c(coord.ind.partiel[(i - 
																																						  1) * nbre.grpe + (j - 1), 2], coord.ind.partiel[(i - 
																																																			 1) * nbre.grpe + j, 2]), col = col.ind[i], 
									  ...)
							}
							else lines(c(coord.ind[i, 1], coord.ind.partiel[(i - 
																			   1) * nbre.grpe + j, 1]), c(coord.ind[i, 
																													2], coord.ind.partiel[(i - 1) * nbre.grpe + 
																																			j, 2]), col = col.ind[nb.ind.actif + (i - 
																																													1) * nbre.grpe + j], lty = j, ...)
						  }
						}
					  }
					}
					if (!is.null(res.mfa$ind.sup) & is.na(test.invisible[2])) {
					  coo <- rbind(coo, coord.ind.sup)
					  if (lab.ind) {
						labe <- c(labe, rownames(coord.ind.sup))
					  }
					  else labe <- c(labe, rep("", nrow(coord.ind.sup)))
					  coll <- c(coll, col.ind.sup[1:(nb.ind - nb.ind.actif)])
					  ipch <- c(ipch, rep(21, nrow(coord.ind.sup)))
					  fonte <- c(fonte, rep(3, nrow(coord.ind.sup)))
					  for (i in group.ind.sup) {
						for (j in 1:nbre.grpe) {
						  points(coord.ind.partiel.sup[(i - 1) * nbre.grpe + 
														 j, 1], coord.ind.partiel.sup[(i - 1) * nbre.grpe + 
																						j, 2], cex = 0.8 * par("cex"), col = col.ind.sup[nb.ind - 
																																		   nb.ind.actif + (i - 1) * nbre.grpe + j], 
								 pch = 21)
						  if (lab.par) 
							text(coord.ind.partiel.sup[(i - 1) * nbre.grpe + 
														 j, 1], y = coord.ind.partiel.sup[nb.ind + 
																							(i - 1) * nbre.grpe + j, 2], labels = rownames(coord.ind.partiel.sup)[(i - 
																																									 1) * nbre.grpe + j], pos = 3, col = col.ind.sup[nb.ind - 
																																																					   nb.ind.actif + (i - 1) * nbre.grpe + j], 
								 cex = par("cex") * 0.8)
						  if (chrono) {
							if (j > 1) 
							  lines(c(coord.ind.partiel.sup[(i - 1) * 
															  nbre.grpe + (j - 1), 1], coord.ind.partiel.sup[(i - 
																												1) * nbre.grpe + j, 1]), c(coord.ind.partiel.sup[(i - 
																																									1) * nbre.grpe + (j - 1), 2], coord.ind.partiel.sup[(i - 
																																																						   1) * nbre.grpe + j, 2]), col = col.ind[nb.ind.actif + 
																																																																	i])
						  }
						  else lines(c(coord.ind.sup[i, 1], coord.ind.partiel.sup[(i - 
																					 1) * nbre.grpe + j, 1]), c(coord.ind.sup[i, 
																															  2], coord.ind.partiel.sup[(i - 1) * nbre.grpe + 
																																						  j, 2]), col = col.ind.sup[nb.ind - nb.ind.actif + 
																																													  (i - 1) * nbre.grpe + j], lty = j)
						}
					  }
					}
					if (!is.null(coord.quali) & is.na(test.invisible[3])) {
					  coo <- rbind(coo, coord.quali)
					  if (lab.var) {
						labe <- c(labe, rownames(coord.quali))
					  }
					  else labe <- c(labe, rep("", nrow(coord.quali)))
					  coll <- c(coll, col.quali[1:nrow.coord.quali])
					  ipch <- c(ipch, rep(15, nrow(coord.quali)))
					  fonte <- c(fonte, rep(2, nrow(coord.quali)))
					  for (i in group.quali) {
						for (j in 1:nbre.grpe) {
						  points(coord.quali.partiel[(i - 1) * nbre.grpe + 
													   j, 1], coord.quali.partiel[(i - 1) * nbre.grpe + 
																					j, 2], pch = 15, col = col.quali[nrow.coord.quali + 
																													   (i - 1) * nbre.grpe + j], cex = par("cex") * 
								   0.8)
						  if (lab.var & lab.par) 
							text(coord.quali.partiel[(i - 1) * nbre.grpe + 
													   j, 1], y = coord.quali.partiel[(i - 1) * 
																						nbre.grpe + j, 2], labels = rownames(coord.quali.partiel)[(i - 
																																					 1) * nbre.grpe + j], pos = 3, col = col.quali[nrow.coord.quali + 
																																																	 (i - 1) * nbre.grpe + j], ...)
						  if (chrono) {
							if (j > 1) 
							  lines(c(coord.quali.partiel[(i - 1) * 
															nbre.grpe + (j - 1), 1], coord.quali.partiel[(i - 
																											1) * nbre.grpe + j, 1]), c(coord.quali.partiel[(i - 
																																							  1) * nbre.grpe + (j - 1), 2], coord.quali.partiel[(i - 
																																																				   1) * nbre.grpe + j, 2]), col = col.quali[i])
						  }
						  else lines(c(coord.quali[i, 1], coord.quali.partiel[(i - 
																				 1) * nbre.grpe + j, 1]), c(coord.quali[i, 
																														2], coord.quali.partiel[(i - 1) * nbre.grpe + 
																																				  j, 2]), col = col.quali[nrow.coord.quali + 
																																											(i - 1) * nbre.grpe + j], lty = j)
						}
					  }
					}
					if (!is.null(coord.quali.sup) & is.na(test.invisible[3])) {
					  coo <- rbind(coo, coord.quali.sup)
					  if (lab.var) {
						labe <- c(labe, rownames(coord.quali.sup))
					  }
					  else labe <- c(labe, rep("", nrow(coord.quali.sup)))
					  coll <- c(coll, col.quali.sup[1:nrow(coord.quali.sup)])
					  ipch <- c(ipch, rep(22, nrow(coord.quali.sup)))
					  fonte <- c(fonte, rep(4, nrow(coord.quali.sup)))
					  for (i in group.quali.sup) {
						for (j in 1:nbre.grpe) {
						  points(coord.quali.partiel.sup[(i - 1) * nbre.grpe + 
														   j, 1], coord.quali.partiel.sup[(i - 1) * 
																							nbre.grpe + j, 2], pch = 22, col = col.quali.sup[nrow(coord.quali.sup) + 
																																			   (i - 1) * nbre.grpe + j], cex = par("cex") * 
								   0.8)
						  if (lab.var & lab.par) 
							text(coord.quali.partiel.sup[(i - 1) * nbre.grpe + 
														   j, 1], y = coord.quali.partiel.sup[(i - 
																								 1) * nbre.grpe + j, 2], labels = rownames(coord.quali.partiel.sup)[(i - 
																																									   1) * nbre.grpe + j], pos = 3, col = col.quali.sup[nrow(coord.quali.sup) + 
																																																						   (i - 1) * nbre.grpe + j], ...)
						  if (chrono) {
							if (j > 1) 
							  lines(c(coord.quali.partiel.sup[(i - 1) * 
																nbre.grpe + (j - 1), 1], coord.quali.partiel.sup[(i - 
																													1) * nbre.grpe + j, 1]), c(coord.quali.partiel.sup[(i - 
																																										  1) * nbre.grpe + (j - 1), 2], coord.quali.partiel.sup[(i - 
																																																								   1) * nbre.grpe + j, 2]), col = col.quali[nrow.coord.quali + 
																																																																			  i])
						  }
						  else lines(c(coord.quali.sup[i, 1], coord.quali.partiel.sup[(i - 
																						 1) * nbre.grpe + j, 1]), c(coord.quali.sup[i, 
																																	2], coord.quali.partiel.sup[(i - 1) * nbre.grpe + 
																																								  j, 2]), col = col.quali.sup[nrow(coord.quali.sup) + 
																																																(i - 1) * nbre.grpe + j], lty = j)
						}
					  }
					}
					if (shadowtext) 
					  points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, 
							 ...)
					if (any(labe != "")) {
					  if (autoLab == "auto") 
						autoLab = (length(which(labe != "")) < 50)
					  if (autoLab == TRUE) 
						autoLab(coo[labe != "", 1], y = coo[labe != 
															  "", 2], labels = labe[labe != ""], col = coll[labe != 
																											  ""], font = fonte[labe != ""], shadotext = shadowtext, 
								...)
					  if (autoLab == FALSE) 
						text(coo[labe != "", 1], y = coo[labe != "", 
														 2], labels = labe[labe != ""], col = coll[labe != 
																									 ""], font = fonte[labe != ""], pos = 3, ...)
					}
					if (!shadowtext) 
					  points(coo[, 1], y = coo[, 2], pch = ipch, col = coll, 
							 ...)
					if ((!is.null(partial)) & (habillage == "group")) {
					  L <- list(x = "topleft", legend = rownames(res.mfa$group$Lg)[-c(num.group.sup, 
																					  length(rownames(res.mfa$group$Lg)))], lty = 1:length(rownames(res.mfa$group$Lg)[-c(num.group.sup, 
																																										 length(rownames(res.mfa$group$Lg)))]), text.col = col.hab, 
								col = col.hab, cex = par("cex") * 0.8)
					  L <- modifyList(L, legend)
					  do.call(graphics::legend, L)
					}
					if ((!is.null(partial)) & (habillage != "group")) {
					  L <- list(x = "topleft", legend = rownames(res.mfa$group$Lg)[-c(num.group.sup, 
																					  length(rownames(res.mfa$group$Lg)))], lty = 1:length(rownames(res.mfa$group$Lg)[-c(num.group.sup, 
																																										 length(rownames(res.mfa$group$Lg)))]), cex = par("cex") * 
								  0.8)
					  L <- modifyList(L, legend)
					  do.call(graphics::legend, L)
					}
					if ((habillage != "none") & (habillage != "ind") & (habillage != 
																		"group")) {
					  L <- list(x = "topleft", legend = levels(res.mfa$call$X[, 
																			  habillage]), text.col = col.hab, cex = par("cex") * 
								  0.8)
					  L <- modifyList(L, legend)
					  do.call(graphics::legend, L)
					}
					if (!is.null(coord.ellipse) & is.na(test.invisible[2])) {
					  for (e in 1:nb.ind.actif) {
						debut <- ((nb.ind.actif - 1) * npoint.ellipse) + 
						  1
						fin <- debut + npoint.ellipse - 1
						data.elli <- coord.ellipse[debut:fin, -1]
						lines(data.elli[, 1], y = data.elli[, 2], col = col.ellipse[e])
					  }
					}
					if (!is.null(coord.ellipse)) {
					  for (e in 1:nlevels(coord.ellipse[, 1])) {
						data.elli <- coord.ellipse[(npoint.ellipse * 
													  (e - 1) + 1):(npoint.ellipse * e), -1]
						lines(data.elli[, 1], y = data.elli[, 2], col = col.ellipse[e])
					  }
					}
					if (!is.null(coord.ellipse.par)) {
					  for (i in group.ind.actif) {
						for (j in 1:nbre.grpe) {
						  ind.e <- (i - 1) * nbre.grpe + j
						  data.elli <- coord.ellipse.par[(npoint.ellipse.par * 
															(ind.e - 1) + 1):(npoint.ellipse.par * ind.e), 
														 -1]
						  lines(data.elli[, 1], y = data.elli[, 2], 
								col = col.ellipse.par[ind.e], lty = 2)
						}
					  }
					}
				  }
				  
				  #ADDITION
				  plotlist[["ind"]] <- recordPlot();
				  #
				  
				  return (plotlist);
				  
				}


		
	),
	
	public =  list(
		
		initialize = function(){
			self$basicData <- Class.Analysis.Data$new(); #Done here because it will be shared across all instances of the object otherwise.
		},
    
		#finalizer
		finalize = function() {},
		
		getChemData = function(){return (private$chem_data);},
		getChemDataMFA = function(){return (private$chem_data_mfa);},
		
		basicData = NULL, #main problem with putting this here is that R and R6 by extention is not a strongly typed language.
		
		#builds the MFA table using the basic statistics tables generated by the basicstats class

		buildMFATable = function( chemicals, label_by = "casn" ){
		  
			loginfo("Building MFA table");
		  
			private$chem_data = tribble(~casn, ~name, ~aeid, ~ac50, ~scalar_top, ~intended_target_family);
			#Gabriel could this all be done in one go and not need to go thru every chemical? sure it can.
			#add all required stats for each chemical to the table
		  
			
			self$basicData$buildBasicStatsTable( chemicals );
			self$basicData$computeScalarTop();

			tmp = self$basicData$getScalarTopTable();

			
			private$chem_data <- inner_join(self$basicData$getBasicStatsTable(), tmp ) %>% select(casn, name, aeid, intended_target_family, ac50,aeid, scalar_top);
			#format data correctly for the MFA script
			if(nrow(private$chem_data)>0){

			if(label_by == "casn"){
			  rnames <- paste0(private$chem_data$casn, "-", private$chem_data$aeid);
			} else if (label_by == "name"){
			  rnames <- paste0(private$chem_data$name, "-", private$chem_data$aeid);
			}
			private$chem_data = private$chem_data %>% select(ac50, scalar_top, intended_target_family);
			row.names(private$chem_data) <- rnames;
			private$chem_data$intended_target_family = as.factor(private$chem_data$intended_target_family);
			} else {
			
			private$chem_data = private$chem_data %>% select(ac50, scalar_top, intended_target_family);
			
			}
		  
		  #handle scalar_top data that may be an issue due to infinite values (x/0 scenario during calculation)
		  #also MFA seems to throw infinite value errors for superlarge numbers
		  # so convert all Inf and superhigh values into 10^100 which is still a ridiculously large number but works
		  private$chem_data$scalar_top <- sapply(private$chem_data$scalar_top, function(x){
			if (is.infinite(x)) {
			  return (10^100);
			} else if (x > 10^100 ){
			  return (10^100);
			} else {
			  return (x);
			}
		  });
		  #filter out NA values since MFA complains about those
		  private$chem_data = filter(private$chem_data, !is.na(ac50) & !is.na(scalar_top) & !is.na(intended_target_family));

		  return (private$chem_data);
		},
		#updated FactoMineR::MFA function which records the plots it makes so they can be printed individually later
		#updated FactoMineR::MFA function which records the plots it makes so they can be printed individually later
		#putting base = private$chem_data as a test here.
MFA = function (base = private$chem_data, group, type = rep("s", length(group)), excl = NULL, 
          ind.sup = NULL, ncp = 5, name.group = NULL, num.group.sup = NULL, 
          graph = TRUE, weight.col.mfa = NULL, row.w = NULL, axes = c(1, 
                                                                      2), tab.comp = NULL) 
{
  
  #ADDITION - plot objects
  plotlist <- list();
  #
  
  moy.p <- function(V, poids) {
    res <- sum(V * poids, na.rm = TRUE)/sum(poids[!is.na(V)])
  }
  ec <- function(V, poids) {
    res <- sqrt(sum(V^2 * poids, na.rm = TRUE)/sum(poids[!is.na(V)]))
  }
  funcLg <- function(x, y, ponderation.x, ponderation.y, wt = rep(1/nrow(x), 
                                                                  nrow(x)), cor = FALSE) {
    if (is.data.frame(x)) 
      x <- as.matrix(x)
    else if (!is.matrix(x)) 
      stop("'x' must be a matrix or a data frame")
    if (is.data.frame(y)) 
      y <- as.matrix(y)
    else if (!is.matrix(y)) 
      stop("'y' must be a matrix or a data frame")
    if (!all(is.finite(x))) 
      stop("'x' must contain finite values only")
    if (!all(is.finite(y))) 
      stop("'y' must contain finite values only")
    s <- sum(wt)
    wt <- wt/s
    center <- colSums(wt * x)
    x <- sqrt(wt) * sweep(t(t(x) * sqrt(ponderation.x)), 
                          2, center, check.margin = FALSE)
    center <- colSums(wt * y)
    y <- sqrt(wt) * sweep(t(t(y) * sqrt(ponderation.y)), 
                          2, center, check.margin = FALSE)
    Lg <- 0
    for (i in 1:ncol(x)) Lg <- Lg + sum(crossprod(x[, i], 
                                                  y)^2)
    Lg
  }
  if (!is.null(tab.comp)) {
    if (!is.null(weight.col.mfa)) 
      stop("Weightings on the variables are not allowed with the tab.comp argument")
    if (!is.null(ind.sup)) 
      stop("Supplementary individuals are not allowed with tab.comp")
    if (!is.null(num.group.sup)) 
      stop("Supplementary groups are not allowed with tab.comp")
  }
  nature.group <- NULL
  for (i in 1:length(group)) {
    if ((type[i] == "n") && (!(i %in% num.group.sup))) 
      nature.group <- c(nature.group, "quali")
    if ((type[i] == "n") && (i %in% num.group.sup)) 
      nature.group <- c(nature.group, "quali.sup")
    if (((type[i] == "s") || (type[i] == "c")) && (!(i %in% 
                                                     num.group.sup))) 
      nature.group <- c(nature.group, "quanti")
    if (((type[i] == "s") || (type[i] == "c")) && (i %in% 
                                                   num.group.sup)) 
      nature.group <- c(nature.group, "quanti.sup")
    if (((type[i] == "f") || (type[i] == "f2")) && (!(i %in% 
                                                      num.group.sup))) 
      nature.group <- c(nature.group, "contingency")
    if (((type[i] == "f") || (type[i] == "f2")) && (i %in% 
                                                    num.group.sup)) 
      nature.group <- c(nature.group, "contingency.sup")
  }
  nature.var <- rep(nature.group, times = group)
  type.var <- NULL
  for (i in 1:length(group)) {
    if ((type[i] == "n") && (!(i %in% num.group.sup))) 
      type.var <- c(type.var, rep("quali", group[i]))
    if ((type[i] == "n") && (i %in% num.group.sup)) 
      type.var <- c(type.var, rep("quali.sup", group[i]))
    if (((type[i] == "s") || (type[i] == "c")) && (!(i %in% 
                                                     num.group.sup))) 
      type.var <- c(type.var, rep("quanti", group[i]))
    if (((type[i] == "s") || (type[i] == "c")) && (i %in% 
                                                   num.group.sup)) 
      type.var <- c(type.var, rep("quanti.sup", group[i]))
    if (((type[i] == "f") || (type[i] == "f2") || (type[i] == 
                                                   "f3")) && (!(i %in% num.group.sup))) 
      type.var <- c(type.var, rep(type[i], group[i]))
    if (((type[i] == "f") || (type[i] == "f2") || (type[i] == 
                                                   "f3")) && (i %in% num.group.sup)) 
      type.var <- c(type.var, rep(paste(type[i], "sup", 
                                        sep = "_"), group[i]))
  }
  if (is.null(rownames(base))) 
    rownames(base) = 1:nrow(base)
  if (is.null(colnames(base))) 
    colnames(base) = paste("V", 1:ncol(base), sep = "")
  base <- as.data.frame(base)
  base <- droplevels(base)
  if (!is.null(ind.sup)) {
    base <- rbind.data.frame(base[-ind.sup, ], base[ind.sup, 
                                                    , drop = FALSE])
    ind.sup <- (nrow(base) - length(ind.sup) + 1):nrow(base)
  }
  nbre.var <- ncol(base)
  nbre.group <- length(group)
  group.actif <- NULL
  if ("n" %in% type) {
    niveau = NULL
    for (j in 1:ncol(base)) {
      if (!is.numeric(base[, j])) 
        niveau = c(niveau, levels(base[, j]))
    }
    for (j in 1:ncol(base)) {
      if (!is.numeric(base[, j])) {
        if (sum(niveau %in% levels(base[, j])) != nlevels(base[, 
                                                               j])) 
          levels(base[, j]) = paste(colnames(base)[j], 
                                    levels(base[, j]), sep = "_")
      }
    }
  }
  for (i in 1:nbre.group) if (!(i %in% num.group.sup)) 
    group.actif <- c(group.actif, i)
  group.mod <- group
  nbre.ind <- nrow(base)
  nb.actif <- nbre.ind - length(ind.sup)
  if (nbre.var != sum(group)) 
    stop("not convenient group definition")
  if (nbre.group != length(type)) 
    stop("not convenient type definition")
  base <- as.data.frame(base)
  if (!inherits(base, "data.frame")) 
    stop("base should be a data.frame")
  res.separe <- vector(mode = "list", length = nbre.group)
  if (is.null(name.group)) 
    name.group <- paste("group", 1:nbre.group, sep = ".")
  names(res.separe) <- name.group
  ind.grpe <- 0
  if (any(is.na(base))) {
    if (!("n" %in% type)) 
      for (j in 1:ncol(base)) base[, j] <- replace(base[, 
                                                        j], is.na(base[, j]), mean(base[, j], na.rm = TRUE))
      else {
        if (type[1] != "n") 
          for (j in 1:group[1]) base[, j] <- replace(base[, 
                                                          j], is.na(base[, j]), mean(base[, j], na.rm = TRUE))
          for (g in 2:nbre.group) {
            if (type[g] != "n") 
              for (j in (sum(group[1:(g - 1)]) + 1):sum(group[1:g])) base[, 
                                                                          j] <- replace(base[, j], is.na(base[, j]), 
                                                                                        mean(base[, j], na.rm = TRUE))
          }
          if (is.null(tab.comp)) {
            if (type[1] == "n") 
              for (j in 1:group[1]) base[, j] <- as.factor(replace(as.character(base[, 
                                                                                     j]), is.na(base[, j]), paste(colnames(base)[j], 
                                                                                                                  ".NA", sep = "")))
              for (g in 2:nbre.group) {
                if (type[g] == "n") 
                  for (j in (sum(group[1:(g - 1)]) + 1):sum(group[1:g])) base[, 
                                                                              j] <- as.factor(replace(as.character(base[, 
                                                                                                                        j]), is.na(base[, j]), paste(colnames(base)[j], 
                                                                                                                                                     ".NA", sep = "")))
              }
          }
      }
  }
  if (is.null(row.w)) 
    row.w <- rep(1, nb.actif)
  if (any("f" %in% type) + any("f2" %in% type) + any("f3" %in% 
                                                     type) > 1) 
    stop("For the contingency tables, the type must the the same")
  if (("f" %in% type) || ("f2" %in% type) || ("f3" %in% type)) {
    grfrec <- c(which(type == "f"), which(type == "f2"), 
                which(type == "f3"))
    for (i in grfrec) {
      if ((type[i] == "f2") || (type[i] == "f3") || (i %in% 
                                                     num.group.sup)) {
        if (i == 1) 
          base[, 1:group[1]] <- base[, 1:group[1]]/sum(base[1:nb.actif, 
                                                            1:group[1]])
        else base[, (sum(group[1:(i - 1)]) + 1):sum(group[1:i])] <- base[, 
                                                                         (sum(group[1:(i - 1)]) + 1):sum(group[1:i])]/sum(base[1:nb.actif, 
                                                                                                                               (sum(group[1:(i - 1)]) + 1):sum(group[1:i])])
      }
    }
    type.var == "f"
    if (!any(type.var == "f")) 
      sumT <- 1
    else sumT <- sum(base[1:nb.actif, as.logical((type.var == 
                                                    "f") + (type.var == "f2") + (type.var == "f3"))])
    if (sumT == 0) 
      sumT <- 1
    base[, as.logical((type.var == "f") + (type.var == "f_sup") + 
                        (type.var == "f2") + (type.var == "f2_sup") + (type.var == 
                                                                         "f3") + (type.var == "f3_sup"))] <- base[, as.logical((type.var == 
                                                                                                                                  "f") + (type.var == "f_sup") + (type.var == "f2") + 
                                                                                                                                 (type.var == "f2_sup") + (type.var == "f3") + (type.var == 
                                                                                                                                                                                  "f3_sup"))]/sumT
    F.jt <- list()
    Fi.t <- list()
    for (j in grfrec) {
      if (j == 1) {
        F.jt[[j]] <- apply(base[1:nb.actif, 1:group[1]], 
                           2, sum)
        Fi.t[[j]] <- apply(base[, 1:group[1]], 1, sum)
      }
      else {
        F.jt[[j]] <- apply(base[1:nb.actif, (sum(group[1:(j - 
                                                            1)]) + 1):(sum(group[1:(j - 1)]) + group[j])], 
                           2, sum)
        Fi.t[[j]] <- apply(base[, (sum(group[1:(j - 
                                                  1)]) + 1):(sum(group[1:(j - 1)]) + group[j])], 
                           1, sum)
      }
    }
    if ("f" %in% type.var) {
      row.w[1:nrow(base)] <- 0
      for (j in grfrec) {
        if (!(j %in% num.group.sup)) 
          row.w <- row.w + Fi.t[[j]]
      }
    }
    F..t <- numeric()
    for (j in grfrec) F..t[j] <- sum(Fi.t[[j]][1:nb.actif])
    for (t in grfrec) {
      if (t == 1) {
        base[, 1:group[t]] <- sweep(base[, 1:group[t]], 
                                    2, F.jt[[t]], FUN = "/")
        base[, 1:group[t]] = sweep(base[, 1:group[t]], 
                                   1, Fi.t[[t]]/F..t[t], FUN = "-")
        base[, 1:group[t]] = sweep(base[, 1:group[t]], 
                                   1, row.w, FUN = "/")
      }
      else {
        base[, (sum(group[1:(t - 1)]) + 1):sum(group[1:t])] <- sweep(base[, 
                                                                          (sum(group[1:(t - 1)]) + 1):sum(group[1:t])], 
                                                                     2, F.jt[[t]], FUN = "/")
        base[, (sum(group[1:(t - 1)]) + 1):sum(group[1:t])] <- sweep(base[, 
                                                                          (sum(group[1:(t - 1)]) + 1):sum(group[1:t])], 
                                                                     1, Fi.t[[t]]/F..t[t], FUN = "-")
        base[, (sum(group[1:(t - 1)]) + 1):sum(group[1:t])] <- sweep(base[, 
                                                                          (sum(group[1:(t - 1)]) + 1):sum(group[1:t])], 
                                                                     1, row.w, FUN = "/")
      }
    }
    row.w <- row.w[1:nb.actif]
  }
  if (!is.null(ind.sup)) 
    row.w.moy.ec <- c(row.w, rep(0, length(ind.sup)))
  else row.w.moy.ec <- row.w
  if (is.null(weight.col.mfa)) 
    weight.col.mfa <- rep(1, sum(group.mod))
  if (!is.null(tab.comp)) {
    group.mod <- tab.comp$call$group.mod
    ind.var.group <- tab.comp$call$ind.var
    tab.comp <- tab.comp$tab.disj
  }
  for (g in 1:nbre.group) {
    aux.base <- as.data.frame(base[, (ind.grpe + 1):(ind.grpe + 
                                                       group[g])])
    dimnames(aux.base) <- list(rownames(base), colnames(base)[(ind.grpe + 
                                                                 1):(ind.grpe + group[g])])
    if (type[g] == "s") 
      res.separe[[g]] <- PCA(aux.base, ind.sup = ind.sup, 
                             scale.unit = TRUE, ncp = ncp, row.w = row.w, 
                             graph = FALSE, col.w = weight.col.mfa[(ind.grpe + 
                                                                      1):(ind.grpe + group[g])])
    if (type[g] == "c") 
      res.separe[[g]] <- PCA(aux.base, ind.sup = ind.sup, 
                             scale.unit = FALSE, ncp = ncp, row.w = row.w, 
                             graph = FALSE, col.w = weight.col.mfa[(ind.grpe + 
                                                                      1):(ind.grpe + group[g])])
    if (type[g] == "f" || type[g] == "f2" || type[g] == 
        "f3") 
      res.separe[[g]] <- PCA(aux.base, ind.sup = ind.sup, 
                             scale.unit = FALSE, ncp = ncp, row.w = row.w, 
                             graph = FALSE, col.w = F.jt[[g]] * weight.col.mfa[(ind.grpe + 
                                                                                  1):(ind.grpe + group[g])])
    if (type[g] == "n") {
      for (v in (ind.grpe + 1):(ind.grpe + group[g])) {
        if (!is.factor(base[, v])) 
          stop("factors are not defined in the qualitative groups")
      }
      res.separe[[g]] <- MCA(aux.base, excl = excl[[g]], 
                             ind.sup = ind.sup, ncp = ncp, graph = FALSE, 
                             row.w = row.w)
    }
    if (!is.null(tab.comp)) {
      if (type[g] == "s") 
        res.separe[[g]] <- PCA(tab.comp[, ind.var.group[[g]]], 
                               scale.unit = TRUE, row.w = row.w, ind.sup = ind.sup, 
                               col.w = weight.col.mfa[(ind.grpe + 1):(ind.grpe + 
                                                                        group[g])], graph = FALSE)
      if (type[g] == "c") 
        res.separe[[g]] <- PCA(tab.comp[, ind.var.group[[g]]], 
                               scale.unit = FALSE, row.w = row.w, ind.sup = ind.sup, 
                               col.w = weight.col.mfa[(ind.grpe + 1):(ind.grpe + 
                                                                        group[g])], graph = FALSE)
      if (type[g] == "n") 
        res.separe[[g]] <- MCA(aux.base, ind.sup = ind.sup, 
                               ncp = ncp, graph = FALSE, row.w = row.w, tab.disj = tab.comp[, 
                                                                                            ind.var.group[[g]]])
    }
    ind.grpe <- ind.grpe + group[g]
  }
  data <- matrix(0, nbre.ind, 0)
  ind.grpe <- ind.grpe.mod <- 0
  ponderation <- vector(length = sum(group.mod))
  poids.bary <- NULL
  for (g in 1:nbre.group) {
    aux.base <- base[, (ind.grpe + 1):(ind.grpe + group[g]), 
                     drop = FALSE]
    if (!is.null(tab.comp)) {
      if (g == 1) 
        aux.base <- tab.comp[, 1:group.mod[1]]
      else aux.base <- tab.comp[, (cumsum(group.mod)[g - 
                                                       1] + 1):cumsum(group.mod)[g], drop = FALSE]
    }
    aux.base <- as.data.frame(aux.base)
    colnames(aux.base) <- colnames(base)[(ind.grpe + 1):(ind.grpe + 
                                                           group[g])]
    if (type[g] == "s") {
      centre.aux.base <- apply(as.data.frame(aux.base), 
                               2, moy.p, row.w.moy.ec)
      aux.base <- t(t(as.matrix(aux.base)) - centre.aux.base)
      ecart.type.aux.base <- apply(as.data.frame(aux.base), 
                                   2, ec, row.w.moy.ec)
      ecart.type.aux.base[ecart.type.aux.base <= 1e-08] <- 1
      aux.base <- t(t(aux.base)/ecart.type.aux.base)
      type[g] = "c"
    }
    if (type[g] == "c") {
      data <- cbind.data.frame(data, aux.base)
      ponderation[(ind.grpe.mod + 1):(ind.grpe.mod + group.mod[g])] <- 1/res.separe[[g]]$eig[1, 
                                                                                             1]
    }
    if (type[g] == "f" || type[g] == "f2") {
      data <- cbind.data.frame(data, aux.base)
      ponderation[(ind.grpe.mod + 1):(ind.grpe.mod + group[g])] <- F.jt[[g]]/res.separe[[g]]$eig[1, 
                                                                                                 1]
    }
    if (type[g] == "n") {
      if (!is.null(tab.comp)) {
        if (g == 1) 
          tmp <- tab.comp[, 1:group.mod[1]]
        else tmp <- tab.comp[, (cumsum(group.mod)[g - 
                                                    1] + 1):cumsum(group.mod)[g], drop = FALSE]
      }
      else {
        tmp <- tab.disjonctif(aux.base)
        group.mod[g] <- ncol(tmp)
      }
      centre.tmp <- apply(tmp, 2, moy.p, row.w.moy.ec)
      centre.tmp <- centre.tmp/sum(row.w.moy.ec)
      tmp2 <- tmp * (row.w.moy.ec/sum(row.w.moy.ec))
      poids.bary <- c(poids.bary, colSums(tmp2))
      poids.tmp <- 1 - apply(tmp2, 2, sum)
      if (!is.null(excl[[g]])) 
        poids.tmp[excl[[g]]] <- 0
      ponderation[(ind.grpe.mod + 1):(ind.grpe.mod + group.mod[g])] <- poids.tmp/(res.separe[[g]]$eig[1, 
                                                                                                      1] * group[g])
      tmp <- tmp/sum(row.w.moy.ec)
      tmp <- t(t(as.matrix(tmp)) - centre.tmp)
      ecart.type.tmp <- apply(tmp, 2, ec, row.w.moy.ec)
      if (!is.null(tab.comp)) 
        ecart.type.tmp <- sqrt(centre.tmp * sum(row.w.moy.ec) * 
                                 (1 - centre.tmp * sum(row.w.moy.ec)))/sum(row.w.moy.ec)
      ecart.type.tmp[ecart.type.tmp <= 1e-08] <- 1
      tmp <- t(t(as.matrix(tmp))/ecart.type.tmp)
      data <- cbind.data.frame(data, as.data.frame(tmp))
    }
    if (!is.null(excl)) 
      ponderation[ponderation == 0] <- 1e-15
    ind.grpe <- ind.grpe + group[g]
    ind.grpe.mod <- ind.grpe.mod + group.mod[g]
  }
  data.group.sup.indice <- data.group.sup <- NULL
  data.pca <- data
  rownames(data.pca) <- rownames(base)
  if (!is.null(num.group.sup)) {
    ponderation.tot <- ponderation
    ponderation.group.sup <- NULL
    nb.of.var <- 0
    supp.quanti <- supp.quali <- NULL
    colnames.data.group.sup <- NULL
    for (i in 1:nbre.group) {
      if (i %in% num.group.sup) {
        if ((type[i] == "c") || (type[i] == "f")) 
          supp.quanti <- c(supp.quanti, (1 + nb.of.var):(nb.of.var + 
                                                           group.mod[i]))
        if (type[i] == "n") 
          supp.quali <- c(supp.quali, (1 + nb.of.var):(nb.of.var + 
                                                         group.mod[i]))
        if (is.null(data.group.sup)) 
          data.group.sup <- as.data.frame(data[, (1 + 
                                                    nb.of.var):(nb.of.var + group.mod[i])])
        else data.group.sup <- cbind.data.frame(data.group.sup, 
                                                data[, (1 + nb.of.var):(nb.of.var + group.mod[i])])
        if (ncol(data.group.sup) > 1) 
          colnames.data.group.sup <- c(colnames.data.group.sup, 
                                       colnames(data)[(1 + nb.of.var):(nb.of.var + 
                                                                         group.mod[i])])
        else colnames.data.group.sup <- colnames(data)[1 + 
                                                         nb.of.var]
        ponderation.group.sup <- c(ponderation.group.sup, 
                                   ponderation[(1 + nb.of.var):(nb.of.var + group.mod[i])])
      }
      nb.of.var <- nb.of.var + group.mod[i]
    }
    colnames(data.group.sup) <- colnames.data.group.sup
    ponderation <- ponderation.tot[-c(supp.quanti, supp.quali)]
    data <- data[, -c(supp.quanti, supp.quali)]
    data.group.sup.indice <- (ncol(data) + 1):(ncol(data) + 
                                                 ncol(data.group.sup))
    data.pca <- cbind.data.frame(data, data.group.sup)
  }
  ncp.tmp <- min(nb.actif - 1, ncol(data))
  ind.var <- 0
  ind.quali <- NULL
  for (g in 1:nbre.group) {
    if (type[g] == "n") 
      ind.quali <- c(ind.quali, c((ind.var + 1):(ind.var + 
                                                   group[g])))
    ind.var <- ind.var + group[g]
  }
  aux.quali.sup.indice <- aux.quali.sup <- data.sup <- NULL
  if (!is.null(ind.quali)) {
    aux.quali.sup <- as.data.frame(base[, ind.quali, drop = FALSE])
    if (is.null(data.group.sup)) 
      aux.quali.sup.indice <- (ncol(data) + 1):(ncol(data) + 
                                                  ncol(aux.quali.sup))
    else aux.quali.sup.indice <- (ncol(data) + ncol(data.group.sup) + 
                                    1):(ncol(data) + ncol(data.group.sup) + ncol(aux.quali.sup))
    data.pca <- cbind.data.frame(data.pca, aux.quali.sup)
  }
  row.w = row.w[1:nb.actif]
  if ((!is.null(tab.comp)) & (any("n" %in% type))) {
    data.pca <- data.pca[, -aux.quali.sup.indice]
    aux.quali.sup.indice <- NULL
  }
  res.globale <- PCA(data.pca, scale.unit = FALSE, col.w = ponderation, 
                     row.w = row.w, ncp = ncp, ind.sup = ind.sup, quali.sup = aux.quali.sup.indice, 
                     quanti.sup = data.group.sup.indice, graph = FALSE)
  if ((!is.null(tab.comp)) & (any("n" %in% type))) {
    res.globale$quali.var$coord <- res.globale$var$coord[unlist(ind.var.group[type %in% 
                                                                                "n"]), ]
    res.globale$quali.var$contrib <- res.globale$var$contrib[unlist(ind.var.group[type %in% 
                                                                                    "n"]), ]
    res.globale$quali.var$cos2 <- res.globale$var$cos2[unlist(ind.var.group[type %in% 
                                                                              "n"]), ]
    res.globale$call$quali.sup$barycentre <- sweep(crossprod(tab.comp[, 
                                                                      unlist(ind.var.group[type %in% "n"])], as.matrix(data.pca)), 
                                                   1, apply(tab.comp[, unlist(ind.var.group[type %in% 
                                                                                              "n"])], 2, sum), FUN = "/")
    res.globale$quali.sup$coord <- sweep(crossprod(tab.comp[, 
                                                            unlist(ind.var.group[type %in% "n"])], res.globale$ind$coord), 
                                         1, apply(tab.comp[, unlist(ind.var.group[type %in% 
                                                                                    "n"])], 2, sum), FUN = "/")
  }
  ncp <- min(ncp, nrow(res.globale$eig))
  call <- res.globale$call
  call$group <- group
  call$type <- type
  call$ncp <- ncp
  call$group.mod <- group.mod
  call$num.group.sup <- num.group.sup
  call$name.group <- name.group
  call$X <- base
  call$XTDC <- data
  call$nature.group <- nature.group
  call$nature.var <- nature.var
  contrib.group <- matrix(NA, length(group.actif), ncp)
  dimnames(contrib.group) <- list(name.group[group.actif], 
                                  paste("Dim", c(1:ncp), sep = "."))
  dist2.group <- vector(length = length(group.actif))
  ind.var <- ind.var.sup <- 0
  for (g in 1:length(group.actif)) {
    if (group.mod[group.actif[g]] != 1) 
      contrib.group[g, ] <- apply(res.globale$var$contrib[(ind.var + 
                                                             1):(ind.var + group.mod[group.actif[g]]), 1:ncp]/100, 
                                  2, sum)
    else contrib.group[g, ] <- res.globale$var$contrib[ind.var + 
                                                         1, 1:ncp]/100
    ind.var <- ind.var + group.mod[group.actif[g]]
    dist2.group[g] <- sum((res.separe[[group.actif[g]]]$eig[, 
                                                            1]/res.separe[[group.actif[g]]]$eig[1, 1])^2)
  }
  coord.group <- t(t(contrib.group) * res.globale$eig[1:ncol(contrib.group), 
                                                      1])
  cos2.group <- coord.group^2/dist2.group
  if (!is.null(num.group.sup)) {
    coord.group.sup <- matrix(NA, length(num.group.sup), 
                              ncp)
    dimnames(coord.group.sup) <- list(name.group[num.group.sup], 
                                      paste("Dim", c(1:ncp), sep = "."))
    ind.gc <- 0
    for (gc in 1:length(num.group.sup)) {
      for (k in 1:ncp) {
        if (is.null(ind.sup)) 
          coord.group.sup[gc, k] <- funcLg(res.globale$ind$coord[, 
                                                                 k, drop = FALSE], data.group.sup[, (ind.gc + 
                                                                                                       1):(ind.gc + group.mod[num.group.sup[gc]]), 
                                                                                                  drop = FALSE], ponderation.x = 1/res.globale$eig[k, 
                                                                                                                                                   1], ponderation.y = ponderation.group.sup[(ind.gc + 
                                                                                                                                                                                                1):(ind.gc + group.mod[num.group.sup[gc]])], 
                                           wt = row.w/sum(row.w))
        else coord.group.sup[gc, k] <- funcLg(res.globale$ind$coord[-ind.sup, 
                                                                    k, drop = FALSE], data.group.sup[-ind.sup, 
                                                                                                     (ind.gc + 1):(ind.gc + group.mod[num.group.sup[gc]]), 
                                                                                                     drop = FALSE], ponderation.x = 1/res.globale$eig[k, 
                                                                                                                                                      1], ponderation.y = ponderation.group.sup[(ind.gc + 
                                                                                                                                                                                                   1):(ind.gc + group.mod[num.group.sup[gc]])], 
                                              wt = row.w/sum(row.w))
      }
      ind.gc <- ind.gc + group.mod[num.group.sup[gc]]
    }
  }
  Lg <- matrix(0, nbre.group + 1, nbre.group + 1)
  ind.gl <- 0
  for (gl in c(group.actif, num.group.sup)) {
    ind.gc <- 0
    for (gc in c(group.actif, num.group.sup)) {
      if (gc >= gl) {
        if (is.null(num.group.sup)) {
          if (is.null(ind.sup)) 
            Lg[gl, gc] <- Lg[gc, gl] <- funcLg(x = data[, 
                                                        ind.gl + (1:group.mod[gl]), drop = FALSE], 
                                               y = data[, ind.gc + (1:group.mod[gc]), 
                                                        drop = FALSE], ponderation.x = ponderation[ind.gl + 
                                                                                                     (1:group.mod[gl])], ponderation.y = ponderation[ind.gc + 
                                                                                                                                                       (1:group.mod[gc])], wt = row.w/sum(row.w))
          else Lg[gl, gc] <- Lg[gc, gl] <- funcLg(x = data[-ind.sup, 
                                                           ind.gl + (1:group.mod[gl]), drop = FALSE], 
                                                  y = data[-ind.sup, ind.gc + (1:group.mod[gc]), 
                                                           drop = FALSE], ponderation.x = ponderation[ind.gl + 
                                                                                                        (1:group.mod[gl])], ponderation.y = ponderation[ind.gc + 
                                                                                                                                                          (1:group.mod[gc])], wt = row.w/sum(row.w))
        }
        else {
          if (is.null(ind.sup)) 
            Lg[gl, gc] <- Lg[gc, gl] <- funcLg(x = cbind.data.frame(data, 
                                                                    data.group.sup)[, ind.gl + (1:group.mod[gl]), 
                                                                                    drop = FALSE], y = cbind.data.frame(data, 
                                                                                                                        data.group.sup)[, ind.gc + (1:group.mod[gc]), 
                                                                                                                                        drop = FALSE], ponderation.x = c(ponderation, 
                                                                                                                                                                         ponderation.group.sup)[ind.gl + (1:group.mod[gl])], 
                                               ponderation.y = c(ponderation, ponderation.group.sup)[ind.gc + 
                                                                                                       (1:group.mod[gc])], wt = row.w/sum(row.w))
          else Lg[gl, gc] <- Lg[gc, gl] <- funcLg(x = cbind.data.frame(data, 
                                                                       data.group.sup)[-ind.sup, ind.gl + (1:group.mod[gl]), 
                                                                                       drop = FALSE], y = cbind.data.frame(data, 
                                                                                                                           data.group.sup)[-ind.sup, ind.gc + (1:group.mod[gc]), 
                                                                                                                                           drop = FALSE], ponderation.x = c(ponderation, 
                                                                                                                                                                            ponderation.group.sup)[ind.gl + (1:group.mod[gl])], 
                                                  ponderation.y = c(ponderation, ponderation.group.sup)[ind.gc + 
                                                                                                          (1:group.mod[gc])], wt = row.w/sum(row.w))
        }
      }
      ind.gc <- ind.gc + group.mod[gc]
    }
    ind.gl <- ind.gl + group.mod[gl]
  }
  Lg[nbre.group + 1, ] <- Lg[, nbre.group + 1] <- apply(Lg[group.actif, 
                                                           ], 2, sum)/res.globale$eig[1, 1]
  Lg[nbre.group + 1, nbre.group + 1] <- sum(Lg[group.actif, 
                                               nbre.group + 1])/res.globale$eig[1, 1]
  dist2.group <- diag(Lg)
  if (!is.null(num.group.sup)) {
    dist2.group.sup <- dist2.group[num.group.sup]
    dist2.group <- dist2.group[-num.group.sup]
  }
  RV <- sweep(Lg, 2, sqrt(diag(Lg)), "/")
  RV <- sweep(RV, 1, sqrt(diag(Lg)), "/")
  rownames(Lg) <- colnames(Lg) <- rownames(RV) <- colnames(RV) <- c(name.group, 
                                                                    "MFA")
  data.partiel <- vector(mode = "list", length = nbre.group)
  names(data.partiel) <- name.group
  ind.col <- 0
  for (g in 1:nbre.group) {
    if (g %in% group.actif) {
      data.partiel[[g]] <- as.data.frame(matrix(res.globale$call$centre, 
                                                nrow(data), ncol(data), byrow = TRUE, dimnames = dimnames(data)))
      data.partiel[[g]][, (ind.col + 1):(ind.col + group.mod[g])] <- data[, 
                                                                          (ind.col + 1):(ind.col + group.mod[g])]
      ind.col <- ind.col + group.mod[g]
    }
  }
  res.ind.partiel <- vector(mode = "list", length = nbre.group)
  names(res.ind.partiel) <- name.group
  for (g in group.actif) {
    Xis <- t(t(as.matrix(data.partiel[[g]])) - res.globale$call$centre)
    Xis <- t(t(Xis)/res.globale$call$ecart.type)
    coord.ind.sup <- length(group.actif) * as.matrix(Xis)
    coord.ind.sup <- t(t(coord.ind.sup) * res.globale$call$col.w)
    coord.ind.sup <- crossprod(t(coord.ind.sup), res.globale$svd$V)
    res.ind.partiel[[g]]$coord.sup <- coord.ind.sup
  }
  cor.grpe.fact <- as.matrix(matrix(NA, length(group.actif), 
                                    ncp))
  colnames(cor.grpe.fact) <- paste("Dim", c(1:ncp), sep = ".")
  rownames(cor.grpe.fact) <- name.group[group.actif]
  for (f in 1:ncp) {
    for (g in 1:length(group.actif)) cor.grpe.fact[g, f] <- cov.wt(cbind.data.frame(res.ind.partiel[[group.actif[g]]]$coord.sup[1:nb.actif, 
                                                                                                                                f], res.globale$ind$coord[, f]), wt = row.w/sum(row.w), 
                                                                   method = "ML", cor = TRUE)$cor[1, 2]
  }
  It <- vector(length = ncp)
  for (g in group.actif) It <- It + apply(res.ind.partiel[[g]]$coord.sup[1:nb.actif, 
                                                                         ]^2 * row.w, 2, sum)
  rap.inertie <- apply(res.globale$ind$coord^2 * row.w, 2, 
                       sum) * length(group.actif)/It
  res.groupes <- list(Lg = Lg, RV = RV, coord = coord.group[, 
                                                            1:ncp], contrib = contrib.group[, 1:ncp] * 100, cos2 = cos2.group[, 
                                                                                                                              1:ncp], dist2 = dist2.group[-length(dist2.group)], correlation = cor.grpe.fact[, 
                                                                                                                                                                                                             1:ncp])
  if (!is.null(num.group.sup)) {
    res.groupes$coord.sup <- coord.group.sup[, 1:ncp, drop = FALSE]
    res.groupes$cos2.sup <- coord.group.sup[, 1:ncp, drop = FALSE]^2/dist2.group.sup
    res.groupes$dist2.sup <- dist2.group.sup
  }
  inertie.intra.ind.partiel <- matrix(NA, (nb.actif * length(group.actif)), 
                                      ncp)
  nom.ligne <- NULL
  for (i in 1:nb.actif) nom.ligne <- c(nom.ligne, paste(rownames(base)[i], 
                                                        name.group[group.actif], sep = "."))
  tmp <- array(0, dim = c(nrow(res.globale$ind$coord), ncp, 
                          length(group.actif)))
  for (g in 1:length(group.actif)) tmp[, , g] <- (res.ind.partiel[[group.actif[g]]]$coord.sup[1:nb.actif, 
                                                                                              1:ncp, drop = FALSE] - res.globale$ind$coord[, 1:ncp, 
                                                                                                                                           drop = FALSE])^2/length(group.actif)
  tmp <- tmp * row.w
  variab.auxil <- apply(tmp, 2, sum)
  tmp <- sweep(tmp, 2, variab.auxil, FUN = "/") * 100
  inertie.intra.ind <- apply(tmp, c(1, 2), sum)
  for (i in 1:nb.actif) inertie.intra.ind.partiel[((i - 1) * 
                                                     length(group.actif) + 1):(i * length(group.actif)), 
                                                  ] <- t(tmp[i, 1:ncp, ])
  rownames(inertie.intra.ind) <- rownames(res.globale$ind$coord)
  rownames(inertie.intra.ind.partiel) <- nom.ligne
  colnames(inertie.intra.ind) <- colnames(inertie.intra.ind.partiel) <- paste("Dim", 
                                                                              c(1:ncp), sep = ".")
  tab.partial.axes <- matrix(NA, nb.actif, ncp * nbre.group)
  rownames(tab.partial.axes) <- rownames(data)[1:nb.actif]
  nom.axes <- paste("Dim", c(1:ncp), sep = "")
  nom.col <- NULL
  debut <- 0
  for (g in 1:nbre.group) {
    nom.col <- c(nom.col, paste(nom.axes, name.group[g], 
                                sep = "."))
    nbcol <- min(ncp, ncol(res.separe[[g]]$ind$coord))
    tab.partial.axes[, (debut + 1):(debut + nbcol)] <- res.separe[[g]]$ind$coord[, 
                                                                                 1:nbcol]
    debut <- debut + ncp
  }
  colnames(tab.partial.axes) <- nom.col
  indice.col.NA <- which(!is.na(tab.partial.axes[1, ]))
  tab.partial.axes <- tab.partial.axes[, indice.col.NA]
  centre <- apply(tab.partial.axes, 2, moy.p, res.globale$call$row.w)
  tab.partial.axes <- t(t(tab.partial.axes) - centre)
  ecart.type <- apply(tab.partial.axes, 2, ec, res.globale$call$row.w)
  ecart.type[ecart.type <= 1e-08] <- 1
  tab.partial.axes <- t(t(tab.partial.axes)/ecart.type)
  coord.res.partial.axes <- t(tab.partial.axes * res.globale$call$row.w)
  coord.res.partial.axes <- crossprod(t(coord.res.partial.axes), 
                                      res.globale$svd$U[, 1:ncp])
  contrib.res.partial.axes <- coord.res.partial.axes * 0
  debut <- 0
  for (g in 1:nbre.group) {
    nbcol <- min(ncp, ncol(res.separe[[g]]$ind$coord))
    if (g %in% group.actif) 
      contrib.res.partial.axes[(debut + 1):(debut + nbcol), 
                               ] <- coord.res.partial.axes[(debut + 1):(debut + 
                                                                          nbcol), ]^2 * res.separe[[g]]$eig[1:nbcol, 1]/res.separe[[g]]$eig[1, 
                                                                                                                                            1]
    debut <- debut + nbcol
  }
  contrib.res.partial.axes <- t(t(contrib.res.partial.axes)/apply(contrib.res.partial.axes, 
                                                                  2, sum)) * 100
  sigma <- apply(tab.partial.axes, 2, ec, res.globale$call$row.w)
  cor.res.partial.axes <- coord.res.partial.axes/sigma
  colnames(coord.res.partial.axes) <- paste("Dim", c(1:ncol(coord.res.partial.axes)), 
                                            sep = ".")
  dimnames(contrib.res.partial.axes) <- dimnames(cor.res.partial.axes) <- dimnames(coord.res.partial.axes)
  summary.n <- as.data.frame(matrix(NA, 0, 4))
  colnames(summary.n) <- c("group", "variable", "modalite", 
                           "effectif")
  summary.c <- as.data.frame(matrix(NA, 0, 6))
  colnames(summary.c) <- c("group", "variable", "moyenne", 
                           "ecart.type", "minimum", "maximum")
  for (g in 1:nbre.group) {
    if ((type[g] == "c") || (type[g] == "f")) {
      statg <- as.data.frame(matrix(NA, ncol(res.separe[[g]]$call$X), 
                                    6))
      colnames(statg) <- c("group", "variable", "moyenne", 
                           "ecart.type", "minimum", "maximum")
      statg[, "group"] <- rep(g, nrow(statg))
      statg[, "variable"] <- colnames(res.separe[[g]]$call$X)
      statg[, "moyenne"] <- res.separe[[g]]$call$centre
      if (!is.null(res.separe[[g]]$call$ecart.type)) 
        statg[, "ecart.type"] <- res.separe[[g]]$call$ecart.type
      statg[, "minimum"] <- apply(res.separe[[g]]$call$X, 
                                  2, min)
      statg[, "maximum"] <- apply(res.separe[[g]]$call$X, 
                                  2, max)
      if (!is.null(res.separe[[g]]$call$ecart.type)) 
        statg[, -c(1, 2)] <- round(statg[, -c(1, 2)], 
                                   digits = 2)
      else statg[, -c(1, 2, 4)] <- round(statg[, -c(1, 
                                                    2, 4)], digits = 2)
      summary.c <- rbind(summary.c, statg)
    }
    else {
      if (is.null(excl[[g]])) 
        statg <- as.data.frame(matrix(NA, length(res.separe[[g]]$call$marge.col), 
                                      4))
      else statg <- as.data.frame(matrix(NA, length(res.separe[[g]]$call$marge.col[-excl[[g]]]), 
                                         4))
      colnames(statg) <- c("group", "variable", "modalite", 
                           "effectif")
      statg[, "group"] <- rep(g, nrow(statg))
      res.separe[[g]]$call$X <- as.data.frame(res.separe[[g]]$call$X)
      nb.var <- ncol(res.separe[[g]]$call$X)
      nb.mod <- NULL
      nom.mod <- NULL
      nb.mod.orig <- NULL
      nb.mod.high <- 0
      for (v in 1:nb.var) {
        if (is.null(excl[[g]])) {
          nb.mod <- c(nb.mod, nlevels(res.separe[[g]]$call$X[, 
                                                             v]))
        }
        else {
          nb.mod.orig <- c(nb.mod.orig, nlevels(res.separe[[g]]$call$X[, 
                                                                       v]))
          nb.mod.low <- nb.mod.high
          nb.mod.high <- nb.mod.high + nlevels(res.separe[[g]]$call$X[, 
                                                                      v])
          nb.mod.rm.select <- (nb.mod.low < res.separe[[g]]$call$excl & 
                                 res.separe[[g]]$call$excl <= nb.mod.high)
          nb.mod.rm <- res.separe[[g]]$call$excl[nb.mod.rm.select]
          nb.mod <- c(nb.mod, (nb.mod.orig[v] - length(nb.mod.rm)))
        }
        nom.mod <- c(nom.mod, levels(res.separe[[g]]$call$X[, 
                                                            v]))
      }
      if (!is.null(excl[[g]])) 
        nom.mod <- nom.mod[-excl[[g]]]
      statg[, "variable"] <- rep(colnames(res.separe[[g]]$call$X), 
                                 nb.mod)
      statg[, "modalite"] <- nom.mod
      if (is.null(excl[[g]])) 
        statg[, "effectif"] <- res.separe[[g]]$call$marge.col * 
        nbre.ind * nb.var
      else statg[, "effectif"] <- res.separe[[g]]$call$marge.col[-excl[[g]]] * 
        nbre.ind * nb.var
      summary.n <- rbind(summary.n, statg)
    }
  }
  eig <- res.globale$eig
  nom.ligne <- NULL
  for (i in 1:nbre.ind) {
    ind.tmp <- rownames(base)[i]
    nom.ligne <- c(nom.ligne, paste(ind.tmp, name.group[group.actif], 
                                    sep = "."))
  }
  coord.ind.partiel <- matrix(NA, (nbre.ind * length(group.actif)), 
                              ncp)
  rownames(coord.ind.partiel) <- nom.ligne
  colnames(coord.ind.partiel) <- paste("Dim", c(1:ncp), sep = ".")
  coord.ind <- rbind(res.globale$ind$coord[, 1:ncp, drop = FALSE], 
                     res.globale$ind.sup$coord[, 1:ncp, drop = FALSE])
  cos2.ind <- rbind(res.globale$ind$cos2[, 1:ncp, drop = FALSE], 
                    res.globale$ind.sup$cos2[, 1:ncp, drop = FALSE])
  contrib.ind <- res.globale$ind$contrib[, 1:ncp, drop = FALSE]
  liste.ligne <- seq(1, nbre.ind * length(group.actif), by = length(group.actif))
  for (g in 1:length(group.actif)) coord.ind.partiel[liste.ligne + 
                                                       g - 1, ] <- res.ind.partiel[[group.actif[g]]]$coord.sup[, 
                                                                                                               1:ncp, drop = FALSE]
  if (!is.null(ind.sup)) {
    res.ind.sup <- list(coord = coord.ind[(nb.actif + 1):nrow(coord.ind), 
                                          , drop = FALSE], cos2 = cos2.ind[(nb.actif + 1):nrow(coord.ind), 
                                                                           , drop = FALSE], coord.partiel = coord.ind.partiel[(length(group.actif) * 
                                                                                                                                 nb.actif + 1):nrow(coord.ind.partiel), , drop = FALSE])
    res.ind <- list(coord = coord.ind[1:nb.actif, , drop = FALSE], 
                    contrib = contrib.ind, cos2 = cos2.ind[1:nb.actif, 
                                                           , drop = FALSE], within.inertia = inertie.intra.ind[1:nb.actif, 
                                                                                                               1:ncp, drop = FALSE], coord.partiel = coord.ind.partiel[1:(length(group.actif) * 
                                                                                                                                                                            nb.actif), , drop = FALSE], within.partial.inertia = inertie.intra.ind.partiel[1:(length(group.actif) * 
                                                                                                                                                                                                                                                                nb.actif), 1:ncp, drop = FALSE])
  }
  else res.ind <- list(coord = coord.ind, contrib = contrib.ind, 
                       cos2 = cos2.ind, within.inertia = inertie.intra.ind[, 
                                                                           1:ncp, drop = FALSE], coord.partiel = coord.ind.partiel, 
                       within.partial.inertia = inertie.intra.ind.partiel[, 
                                                                          1:ncp, drop = FALSE])
  res.quali.var <- res.quali.var.sup <- NULL
  bool.act <- FALSE
  bool.sup <- FALSE
  if (!is.null(ind.quali)) {
    coord.quali <- res.globale$quali.sup$coord[, 1:ncp, 
                                               drop = FALSE]
    cos2.quali <- res.globale$quali.sup$cos2[, 1:ncp, drop = FALSE]
    val.test.quali <- res.globale$quali.sup$v.test[, 1:ncp, 
                                                   drop = FALSE]
    contrib.quali <- coord.quali * 0
    commun <- intersect(rownames(res.globale$var$contrib), 
                        rownames(contrib.quali))
    if (!is.null(commun)) 
      contrib.quali[commun, ] <- res.globale$var$contrib[commun, 
                                                         1:ncp, drop = FALSE]
    barycentre <- res.globale$call$quali.sup$barycentre
    coord.quali.partiel <- matrix(NA, (nrow(barycentre) * 
                                         length(group.actif)), ncp)
    nom.ligne.bary <- NULL
    for (q in 1:nrow(barycentre)) {
      ind.tmp <- rownames(barycentre)[q]
      nom.ligne.bary <- c(nom.ligne.bary, paste(ind.tmp, 
                                                name.group[group.actif], sep = "."))
    }
    rownames(coord.quali.partiel) <- nom.ligne.bary
    liste.ligne <- seq(1, (nrow(barycentre) * length(group.actif)), 
                       by = length(group.actif))
    inertie.intra.cg.partiel <- matrix(NA, (nrow(barycentre) * 
                                              length(group.actif)), ncp)
    tmp <- array(0, dim = c(nrow(res.globale$quali.sup$coord), 
                            ncp, length(group.actif)))
    ind.col <- 0
    for (g in 1:length(group.actif)) {
      cg.partiel <- as.data.frame(matrix(res.globale$call$centre, 
                                         nrow(barycentre), ncol(barycentre), byrow = TRUE, 
                                         dimnames = dimnames(barycentre)))
      cg.partiel[, (ind.col + 1):(ind.col + group.mod[group.actif[g]])] <- barycentre[, 
                                                                                      (ind.col + 1):(ind.col + group.mod[group.actif[g]])]
      ind.col <- ind.col + group.mod[group.actif[g]]
      Xis <- t((t(cg.partiel) - res.globale$call$centre)/res.globale$call$ecart.type)
      coord.quali.sup <- length(group.actif) * as.matrix(Xis)
      coord.quali.sup <- t(t(coord.quali.sup) * res.globale$call$col.w)
      coord.quali.sup <- crossprod(t(coord.quali.sup), 
                                   res.globale$svd$V)
      coord.quali.partiel[liste.ligne + g - 1, ] <- coord.quali.sup[, 
                                                                    1:ncp]
      tmp[, , g] <- (coord.quali.sup[, 1:ncp, drop = FALSE] - 
                       res.globale$quali.sup$coord[, 1:ncp, drop = FALSE])^2/length(group.actif)
    }
    colnames(coord.quali.partiel) <- paste("Dim", 1:ncp, 
                                           sep = ".")
    tmp <- sweep(tmp, 2, variab.auxil, FUN = "/") * 100
    tmp <- sweep(tmp, 1, poids.bary * sum(row.w), FUN = "*")
    inertie.intra.cg <- apply(tmp, c(1, 2), sum)
    for (i in 1:nrow(barycentre)) inertie.intra.cg.partiel[((i - 
                                                               1) * length(group.actif) + 1):(i * length(group.actif)), 
                                                           ] <- t(tmp[i, 1:ncp, ])
    rownames(inertie.intra.cg) <- rownames(res.globale$quali.sup$coord)
    rownames(inertie.intra.cg.partiel) <- nom.ligne.bary
    colnames(inertie.intra.cg) <- colnames(inertie.intra.cg.partiel) <- paste("Dim", 
                                                                              c(1:ncp), sep = ".")
    ind.col <- 0
    ind.col.act <- NULL
    ind.col.sup <- NULL
    ind.excl <- NULL
    ind.excl.act <- NULL
    for (g in 1:nbre.group) {
      if (type[g] == "n") {
        if (g %in% num.group.sup) 
          ind.col.sup <- c(ind.col.sup, (ind.col + 1):(ind.col + 
                                                         group.mod[g]))
        else ind.col.act <- c(ind.col.act, (ind.col + 
                                              1):(ind.col + group.mod[g]))
        if (!is.null(excl[[g]])) {
          ind.excl <- ind.col + excl[[g]]
          ind.excl.act <- c(ind.excl.act, ind.excl)
        }
        ind.col = ind.col + group.mod[g]
      }
    }
    if (!is.null(ind.excl.act)) 
      ind.col.act <- ind.col.act[-ind.excl.act]
    if (!is.null(ind.col.sup)) {
      coord.quali.sup <- coord.quali[ind.col.sup, , drop = FALSE]
      cos2.quali.sup <- cos2.quali[ind.col.sup, , drop = FALSE]
      val.test.quali.sup <- val.test.quali[ind.col.sup, 
                                           , drop = FALSE]
      coord.quali.partiel.sup <- coord.quali.partiel[unlist(lapply(ind.col.sup, 
                                                                   function(k) seq(length(group.actif) * (k - 1) + 
                                                                                     1, length = length(group.actif)))), ]
      inertie.intra.cg.sup <- inertie.intra.cg[ind.col.sup, 
                                               1:ncp]
      inertie.intra.cg.partiel.sup <- inertie.intra.cg.partiel[unlist(lapply(ind.col.sup, 
                                                                             function(k) seq(length(group.actif) * (k - 1) + 
                                                                                               1, length = length(group.actif)))), 1:ncp]
      bool.sup <- TRUE
    }
    if (!is.null(ind.col.act)) {
      coord.quali.act <- coord.quali[ind.col.act, , drop = FALSE]
      contrib.quali.act <- contrib.quali[ind.col.act, 
                                         , drop = FALSE]
      val.test.quali.act <- NULL
      if (is.null(tab.comp)) {
        cos2.quali.act <- cos2.quali[ind.col.act, , 
                                     drop = FALSE]
        val.test.quali.act <- val.test.quali[ind.col.act, 
                                             , drop = FALSE]
      }
      else {
        cos2.quali.act <- res.globale$quali.var$cos2
      }
      coord.quali.partiel.act <- coord.quali.partiel[unlist(lapply(ind.col.act, 
                                                                   function(k) seq(length(group.actif) * (k - 1) + 
                                                                                     1, length = length(group.actif)))), ]
      inertie.intra.cg.act <- inertie.intra.cg[ind.col.act, 
                                               1:ncp]
      inertie.intra.cg.partiel.act <- inertie.intra.cg.partiel[unlist(lapply(ind.col.act, 
                                                                             function(k) seq(length(group.actif) * (k - 1) + 
                                                                                               1, length = length(group.actif)))), 1:ncp]
      bool.act <- TRUE
    }
    if (bool.act) 
      res.quali.var <- list(coord = coord.quali.act, contrib = contrib.quali.act, 
                            cos2 = cos2.quali.act, v.test = val.test.quali.act, 
                            coord.partiel = coord.quali.partiel.act, within.inertia = inertie.intra.cg.act, 
                            within.partial.inertia = inertie.intra.cg.partiel.act)
    if (bool.sup) 
      res.quali.var.sup <- list(coord = coord.quali.sup, 
                                cos2 = cos2.quali.sup, v.test = val.test.quali.sup, 
                                coord.partiel = coord.quali.partiel.sup, within.inertia = inertie.intra.cg.sup, 
                                within.partial.inertia = inertie.intra.cg.partiel.sup)
  }
  indice.quanti <- NULL
  indice.freq <- NULL
  num.tmp <- 0
  for (g in group.actif) {
    if (type[g] == "c") 
      indice.quanti <- c(indice.quanti, c((num.tmp + 1):(num.tmp + 
                                                           group.mod[g])))
    if (type[g] == "f") 
      indice.freq <- c(indice.freq, c((num.tmp + 1):(num.tmp + 
                                                       group.mod[g])))
    num.tmp <- num.tmp + group.mod[g]
  }
  res.quanti.var <- NULL
  if (!is.null(indice.quanti)) {
    coord.quanti.var <- res.globale$var$coord[indice.quanti, 
                                              1:ncp, drop = FALSE]
    cos2.quanti.var <- res.globale$var$cos2[indice.quanti, 
                                            1:ncp, drop = FALSE]
    contrib.quanti.var <- res.globale$var$contrib[indice.quanti, 
                                                  1:ncp, drop = FALSE]
    cor.quanti.var <- res.globale$var$cor[indice.quanti, 
                                          1:ncp, drop = FALSE]
    res.quanti.var <- list(coord = coord.quanti.var, contrib = contrib.quanti.var, 
                           cos2 = cos2.quanti.var, cor = cor.quanti.var)
  }
  res.freq <- NULL
  if (!is.null(indice.freq)) {
    coord.freq <- res.globale$var$coord[indice.freq, 1:ncp, 
                                        drop = FALSE]
    cos2.freq <- res.globale$var$cos2[indice.freq, 1:ncp, 
                                      drop = FALSE]
    contrib.freq <- res.globale$var$contrib[indice.freq, 
                                            1:ncp, drop = FALSE]
    res.freq <- list(coord = coord.freq, contrib = contrib.freq, 
                     cos2 = cos2.freq)
  }
  res.quanti.var.sup <- NULL
  res.freq.sup <- NULL
  if (!is.null(num.group.sup)) {
    num.tmp <- 0
    indice.quanti <- NULL
    indice.freq <- NULL
    for (g in num.group.sup) {
      if (type[g] == "c") 
        indice.quanti <- c(indice.quanti, c((num.tmp + 
                                               1):(num.tmp + group.mod[g])))
      if (type[g] == "f") 
        indice.freq <- c(indice.freq, c((num.tmp + 1):(num.tmp + 
                                                         group.mod[g])))
      num.tmp <- num.tmp + group.mod[g]
    }
    if (!is.null(indice.quanti)) {
      coord.quanti.var.sup <- res.globale$quanti.sup$coord[indice.quanti, 
                                                           1:ncp, drop = FALSE]
      cos2.quanti.var.sup <- res.globale$quanti.sup$cos2[indice.quanti, 
                                                         1:ncp, drop = FALSE]
      cor.quanti.var.sup <- res.globale$quanti.sup$cor[indice.quanti, 
                                                       1:ncp, drop = FALSE]
      res.quanti.var.sup <- list(coord = coord.quanti.var.sup, 
                                 cos2 = cos2.quanti.var.sup, cor = cor.quanti.var.sup)
    }
    if (!is.null(indice.freq)) {
      coord.freq.sup <- res.globale$quanti.sup$coord[indice.freq, 
                                                     1:ncp, drop = FALSE]
      cos2.freq.sup <- res.globale$quanti.sup$cos2[indice.freq, 
                                                   1:ncp, drop = FALSE]
      res.freq.sup <- list(coord = coord.freq.sup, cos2 = cos2.freq.sup)
    }
  }
  aux <- res.separe[[1]]$ind$coord
  name.aux <- paste(colnames(res.separe[[1]]$ind$coord), name.group[1], 
                    sep = ".")
  for (g in 2:nbre.group) {
    aux <- cbind(aux, res.separe[[g]]$ind$coord)
    name.aux = c(name.aux, paste(colnames(res.separe[[g]]$ind$coord), 
                                 name.group[g], sep = "."))
  }
  cor.partial.axes <- cov.wt(aux, wt = row.w/sum(row.w), method = "ML", 
                             cor = TRUE)$cor
  dimnames(cor.partial.axes) <- list(name.aux, name.aux)
  res.partial.axes <- list(coord = coord.res.partial.axes[, 
                                                          1:ncp], cor = cor.res.partial.axes[, 1:ncp], contrib = contrib.res.partial.axes[, 
                                                                                                                                          1:ncp], cor.between = cor.partial.axes)
  resultats <- list(separate.analyses = res.separe, eig = eig, 
                    group = res.groupes, inertia.ratio = rap.inertie[1:ncp], 
                    ind = res.ind)
  if (!is.null(ind.sup)) 
    resultats$ind.sup <- res.ind.sup
  if (!is.null(c(res.quanti.var, res.quanti.var.sup))) 
    resultats$summary.quanti = summary.c
  if (!is.null(c(bool.act, bool.sup))) 
    resultats$summary.quali = summary.n
  if (!is.null(res.quanti.var)) 
    resultats$quanti.var = res.quanti.var
  if (!is.null(res.quanti.var.sup)) 
    resultats$quanti.var.sup = res.quanti.var.sup
  if (!is.null(res.freq)) 
    resultats$freq = res.freq
  if (!is.null(res.freq.sup)) 
    resultats$freq.sup = res.freq.sup
  if (bool.act) 
    resultats$quali.var = res.quali.var
  if (bool.sup) 
    resultats$quali.var.sup = res.quali.var.sup
  resultats$partial.axes = res.partial.axes
  resultats$call = call
  resultats$call$call <- match.call()
  resultats$global.pca = res.globale
  class(resultats) <- c("MFA", "list")
  if (graph & (ncp > 1)) {
    if (bool.act | bool.sup) {
      cg.plot.partial <- NULL
      if (!is.null(resultats["quali.var"]$quali.var)) {
        max.inertia <- order(apply(resultats$quali.var$within.inertia[, 
                                                                      1:2], 1, sum))
        cg.plot.partial <- rownames(resultats$quali.var$coord)[max.inertia[1:length(max.inertia)]]
      }
      if (!is.null(resultats$quali.var.sup)) {
        max.inertia <- order(apply(resultats$quali.var.sup$within.inertia[, 
                                                                          1:2], 1, sum))
        cg.plot.partial <- c(cg.plot.partial, rownames(resultats$quali.var.sup$coord)[max.inertia[1:length(max.inertia)]])
      }
      p <- plot.MFA(resultats, choix = "ind", invisible = "ind", 
               partial = cg.plot.partial, habillage = "group", lab.var=TRUE,
               axes = axes, new.plot = TRUE)
      
      #ADDITION
      plotlist[["ind1"]] <- p;
      plotlist[["ind1"]]$partial <- cg.plot.partial;
      plotlist[["ind1"]]$axes <- axes;
      #
      
    }
    max.inertia <- order(apply(resultats$ind$within.inertia[, 
                                                            1:2], 1, sum))
    p <- plot.MFA(resultats, choix = "axes", habillage = "group", 
             axes = axes, new.plot = TRUE, shadowtext = TRUE)
    
    #ADDITION
    plotlist[["axes"]] <- p;
    #
    var_partial <- rownames(resultats$ind$coord)[max.inertia[c(1:2, 
                                                               nrow(resultats$ind$coord) - 1, nrow(resultats$ind$coord))]];
    p <- plot.MFA(resultats, choix = "ind", invisible = "quali", 
             partial = var_partial, 
             habillage = "group", axes = axes, new.plot = TRUE)
    
    #ADDITION
    plotlist[["ind2"]] <- p;
    plotlist[["ind2"]]$partial <- var_partial;
    plotlist[["ind2"]]$axes <- axes;
    #
    
    if ("c" %in% type) {
      p <- plot.MFA(resultats, choix = "var", habillage = "group", 
               axes = axes, new.plot = TRUE, shadowtext = TRUE)
      
      #ADDITION
      plotlist[["var"]] <- p;
      #

    }
    
    if ("f" %in% type) {
      p <- plot.MFA(resultats, choix = "freq", habillage = "group", 
               axes = axes, new.plot = TRUE)
      
      #ADDITION
      plotlist[["freq"]] <- p;
      #
      
    }
    
    
    p <- plot.MFA(resultats, choix = "ind", invisible = "quali", 
             habillage = "ind", axes = axes, new.plot = TRUE, 
             col.hab = 1 + 3 * (1:nbre.ind) %in% ind.sup)
    
    #ADDITION
    plotlist[["ind3"]] <- p;
    #
    
    p <- plot.MFA(resultats, choix = "group", axes = axes, new.plot = TRUE)
    
    #ADDITION
    plotlist[["group"]] <- p;
    #
    
  }
  resultats$plotlist <- plotlist;
  private$chem_data_mfa = resultats;
  return(resultats)
}





		



	)

	
)