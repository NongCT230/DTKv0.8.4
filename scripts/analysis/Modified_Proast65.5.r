dtype <- 1 # continuous data type
modeltype <- 5 # model 5 selected
concCol <- 2 # Dose/Concentration column
respCol <- 3 # Response column
inCES <- NA_real_ # 0.05 for 5%, 0.1 for 10%...etc.
includHill <- 0 # include HILL model


#############################################
## To obtain results
## 
## results stored in variable: proastTemp
##
## HILL BMD: proastTemp$HILL$CED
## HILL BMDL: proastTemp$HILL$conf.int[,1]
## HILL BMDU: proastTemp$HILL$conf.int[,2]
##
## EXP BMD: proastTemp$EXP$CED
## EXP BMDL: proastTemp$EXP$conf.int[,1]
## EXP BMDU: proastTemp$EXP$conf.int[,2]
##############################################


##############
## f.proast ##
##############
f.proast <- function (odt = list(), ans.all = 0, er = FALSE, resize = FALSE, 
                      scale.ans = FALSE, const.var = F, show.warnings = F, BMR = 5) 
{
  inCES <<- BMR/100;
  
  if (exists("track")) 
    print("f.proast")
  if (is.list(odt)) 
    if (length(odt$EXP) > 0 || length(odt$MLE) > 0) {
      cat("you probably forgot the comma in f.proast()\n")
      ans.cont <- menu(c("no", "yes"), title = "Do you want to continue nonetheless?")
      if (ans.cont == 1) 
        return(invisible())
    }
  library(assertive, pos = 5)
  dec.point <- is_period_for_decimal_point()
  dec.comma <- is_comma_for_decimal_point()
  if (dec.comma || !dec.point) {
    cat("\nATTENTION: you need to use decimal point, not comma, \n             make sure you are using decimal points in your data, \n             otherwise change this in local settings on your computer\n\n")
    f.press.key.to.continue()
  }
  odt.name <- deparse(substitute(odt))
  f.overlap()
  quit <- F
  if (mode(ans.all) == "list") {
    ans.all <- f.adjust.saved(ans.all)
    ans.all$const.var <- const.var
    ans.all$quick.ans <- 1
    if (ans.all$cont) 
      quit <- f.con(ans.all, list.logic = T)
    else quit <- f.cat(ans.all, list.logic = T)
    if (quit) {
      cat("\n\n   Goodbye,\n    PROAST \n\n")
      return(invisible())
    }
  }
  if (!er) {
    ans.all <- f.ini(odt)
    ans.all$odt.name <- odt.name
    ans.all$scale.ans <- scale.ans
    ans.all$const.var <- const.var
    ans.all$show.warnings <- show.warnings
    f.graphwin.size(resize)
    #dtype <- menu(c("continuous, individual data", "binary", 
    #                "ordinal", "quantal", "clustered continuous, individual data", 
    #                "clustered quantal", "continuous, summary data", 
    #                "clustered continuous, summary data", "quantal, CxT", 
    #                "other"), title = "\nWhat type of response data do you want to consider?")
    if (dtype == 10) 
      dtype <- eval(parse(prompt = "Give value for dtype > "))
    if (dtype == 7) 
      dtype <- 10
    if (dtype == 8) 
      dtype <- 15
    if (dtype == 9) 
      dtype <- 84
    ans.all$dtype <- dtype
    if (dtype %in% c(1, 5, 15, 25, 26)) 
      ans.all$plt.mns <- 1
    if (dtype %in% c(10, 15, 250, 260)) {
      ans.all$plt.mns <- 3
      ans.all$CI.plt <- T
    }
    f.assign(".Pr.last", ans.all)
    if (dtype %in% c(1, 5, 10, 15)) {
      ans.all$quick.ans <- modeltype
      #ans.all$quick.ans <- menu(c("single model", "select model 3 or 5 from nested (exponential/Hill) family of models", 
      #                            "change settings first", "select model 3 from nested (exponential/Hill) family of models", 
      #                            "select model 5 from nested (exponential/Hill) family of models"), 
      #                          title = "\nDo you want to fit a single model or select a model from a nested family of models?")
      if (ans.all$quick.ans == 3) {
        ans.all$change <- rep(F, ans.all$nrQ + 1)
        ans.all <- f.change.settings(ans.all, choose = T)
        ans.all$quick.ans <- menu(c("single model", "select model 3 or 5 from nested (exponential/Hill) family of models", 
                                    "change settings first", "select model 3 from nested (exponential/Hill) family of models", 
                                    "select model 5 from nested (exponential/Hill) family of models"), 
                                  title = "\nDo you want to fit a single model or select a model from a nested family of models?\n")
      }
      if (ans.all$quick.ans %in% 2:5) 
        quit <- f.con(ans.all)
    }
    if (dtype %in% c(2, 3, 4, 6, 84)) {
      ans.all$quick.ans <- menu(c("single model", "set of models", 
                                  "change settings first"), title = "\nDo you want to fit a set of models, or choose a single model?")
      if (ans.all$quick.ans == 3) {
        ans.all$change <- rep(F, ans.all$nrQ + 1)
        ans.all <- f.change.settings(ans.all, choose = T)
        ans.all$quick.ans <- menu(c("single model", "set of models"), 
                                  title = "\nDo you want to fit a set of models, or choose a single model?")
      }
      if (ans.all$quick.ans == 2) 
        quit <- f.cat(ans.all)
    }
    if (quit) {
      #cat("\n\n   Goodbye,\n    PROAST \n\n")
      return(invisible())
    }
    ans.all$change[1] <- T
    ans.all <- f.change.settings(ans.all)
    f.assign(".Pr.last", ans.all)
    if (ans.all$quick.ans == 1) {
      x <- ans.all$odt[, ans.all$xans]
      x <- x[!is.na(x)]
      if (is.list(x)) 
        for (jj in 1:length(x)) f.check.nonneg.num(x[, 
                                                     jj], gui)
      else f.check.nonneg.num(x, gui)
      if (max(x, na.rm = T) > 100) {
        cat("\n(Nonzero) doses range from ", min(x[x != 
                                                     0], na.rm = T), "to ", max(x, na.rm = T))
        cat("\nit is recommended to scale x to prevent numerical problems\n")
        ans.all$sf.x <- eval(parse(prompt = "\nQ10: Give scaling factor for dose (e.g. 1000) > "))
      }
      ans.all$change[2] <- T
      ans.all <- f.change.settings(ans.all)
      f.assign(".Pr.last", ans.all)
    }
  }
  if (er) {
    ans.all <- .Pr.last
    if (ans.all$xans == 0 || ans.all$yans == 0) {
      cat("\nOption er=T cannot be used, start analysis all over\n\n")
      return(invisible())
    }
    ans.all$quick.ans <- 1
  }
  if (ans.all$quick.ans %in% c(1, 3)) 
    ans.all <- f.execute(ans.all)
  if (ans.all$cont) 
    quit <- f.con(ans.all)
  else quit <- f.cat(ans.all)
  if (quit) {
    cat("\n\n   Goodbye,\n    PROAST \n\n")
    return(invisible())
  }
}

###############
## f.overlap ##
###############
f.overlap <- function () 
{
  if (exists("track")) 
    print("f.overlap")
  attached <- search()
  attached <- substring(attached, 9, 14)
  lst1 <- attached == "proast"
  attached <- substring(search(), 1, 4)
  lst2 <- attached == "file"
  pos.pr <- 0
  if (any(lst1)) {
    pos.pr <- (1:length(lst1))[lst1]
  }
  if (any(lst2)) {
    pos.pr <- (1:length(lst2))[lst2]
  }
  pos.pr <- pos.pr[1]
  if (pos.pr > 0) {
    if (exists("last.fit", where = pos.pr)) {
      rm(last.fit, pos = pos.pr)
      rm(newfit, pos = pos.pr)
    }
    if (length(pos.pr) > 1) 
      cat("\n\nATTENTION: you have attached more than one PROAST package\n\n")
    obj.pr <- objects(pos = pos.pr)
    wdir <- getwd()
    dum <- substring(wdir, 30, 36)
    if (dum != "proast") {
      obj.wd <- objects(pos = 1)
      obj.wd <- obj.wd[obj.wd != "newfit"]
      lst <- obj.wd %in% obj.pr
      overlap <- obj.wd[lst]
      if (sum(lst) > 0) {
        cat("\nProast cannot be run !")
        cat("\npress escape to stop this session\n")
        cat("\nand then remove or rename the following object(s):\n")
        for (ii in 1:sum(lst)) {
          cat(overlap[ii])
          cat("\n")
        }
        f.press.key.to.continue(NA)
        stop()
      }
    }
  }
  if (exists("track")) 
    print("f.overlap:   END")
  return(invisible())
}

###########
## f.ini ##
###########
f.ini <- function (odt = NULL, gui = FALSE, no.seed = FALSE) 
{
  if (exists("track")) 
    print("f.ini")
  if (!no.seed) 
    set.seed(125)
  ans.all <- list()
  ans.all$PRversion <- f.version()
  nrQ <- 25
  ans.all$gui <- gui
  if (is.list(odt) && length(odt$nvar) == 1) {
    ans.all$varnames <- odt$varnames
    ans.all$nvar <- odt$nvar
    ans.all$odt <- odt$data
    ans.all$info <- odt$info
  }
  else {
    ans.all$varnames <- colnames(odt)
    ans.all$nvar <- length(ans.all$varnames)
    ans.all$odt <- odt
    ans.all$info <- " "
  }
  ans.all$data.0 <- ans.all$odt
  ans.all$nrQ <- nrQ
  ans.all$change <- rep(F, nrQ + 1)
  ans.all$xans <- 0
  ans.all$yans <- 0
  ans.all$sf.x <- 1
  ans.all$nans <- 0
  ans.all$CED.model <- F
  ans.all$sans <- 0
  ans.all$sd.se <- 0
  ans.all$detlim <- 0
  ans.all$Vdetlim <- 0
  ans.all$ans.detlim <- 1
  ans.all$detlim.col <- 0
  ans.all$auto.detlim <- F
  ans.all$transf.ans <- 3
  ans.all$cens.up <- NA
  ans.all$nest.no <- 0
  ans.all$fct4.no <- 0
  ans.all$fct5.no <- 0
  ans.all$covar.dd <- 0
  ans.all$regr.resid <- NA
  ans.all$pred.value <- NA
  ans.all$x.mn <- NA
  ans.all$ans.m6.sd <- 1
  ans.all$sign.q <- 1
  ans.all$const.var <- FALSE
  ans.all$tans <- 0
  ans.all$sf.t <- 1
  ans.all$model.type <- 1
  ans.all$P.gof <- 0.05
  ans.all$aic.crit <- 2
  ans.all$aic <- NA
  ans.all$ces.ans <- 3
  ans.all$alfa.length <- 0
  ans.all$alfa.start <- 10
  ans.all$alfa.Pvalue <- NA
  ans.all$nr.gr <- 1
  ans.all$CI.plt <- F
  ans.all$th.lvm <- NA
  ans.all$shift.tmp <- 0
  ans.all$remove.no <- 0
  ans.all$remove.lev <- list()
  ans.all$remove.name <- "none"
  ans.all$select.no <- 0
  ans.all$select.lev <- list()
  ans.all$select.name <- "all"
  ans.all$fct1 <- 1
  ans.all$fct2 <- 1
  ans.all$fct3 <- 1
  ans.all$fct4 <- 1
  ans.all$fct5 <- 1
  ans.all$fct1.no <- 0
  ans.all$fct2.no <- 0
  ans.all$fct3.no <- 0
  ans.all$fct4.no <- 0
  ans.all$fct5.no <- 0
  ans.all$displ.no <- 0
  ans.all$covar.no <- 0
  ans.all$nr.aa <- 0
  ans.all$nr.bb <- 0
  ans.all$nr.var <- 0
  ans.all$nr.cc <- 0
  ans.all$nr.dd <- 0
  ans.all$ref.lev <- NA
  ans.all$model.ans <- 0
  ans.all$MLE <- NA
  ans.all$converged <- 0
  ans.all$lb <- NA
  ans.all$ub <- NA
  ans.all$scale.ans <- FALSE
  ans.all$CES <- 0
  ans.all$CES.0 <- 0
  ans.all$NES.ans <- 1
  ans.all$BMR <- 0
  ans.all$CED <- NA
  ans.all$ED50 <- FALSE
  ans.all$RPF1 <- 1
  ans.all$RPF2 <- 1
  ans.all$x1 <- 0
  ans.all$x2 <- 0
  ans.all$x3 <- 0
  ans.all$CED.ref <- NA
  ans.all$BMD <- NA
  ans.all$CED.matr <- matrix(NA, ncol = 2, nrow = 1)
  ans.all$plot.type <- 2
  ans.all$l.ty <- 1
  ans.all$heading <- ""
  ans.all$x.leg <- ""
  ans.all$y.leg <- ""
  ans.all$t.leg <- ""
  ans.all$low.y <- NA
  ans.all$upp.y <- NA
  ans.all$regr.lower <- NA
  ans.all$regr.upper <- NA
  ans.all$x.lim <- rep(NA, 2)
  ans.all$y.lim <- NA
  ans.all$shift <- 0
  ans.all$min.x <- NA
  ans.all$max.x <- NA
  ans.all$plt.mns <- 1
  ans.all$cex.1 <- NA
  ans.all$cex.2 <- NA
  ans.all$identify <- 0
  ans.all$label <- 0
  ans.all$ans.plt <- 0
  ans.all$down <- F
  ans.all$CI.dum <- 1
  ans.all$add.one <- T
  ans.all$LVM <- T
  ans.all$do.MA <- F
  ans.all$combi.ans <- F
  ans.all$svg.plots <- F
  ans.all$ans.scale <- 0
  ans.all$bounds.ans <- 2
  ans.all$modelname <- ""
  assign(".ypos", 0, immediate = T, pos = 1)
  ans.all$store.name <- 0
  ans.all$list.logic <- F
  ans.all$CES.cat <- 1
  ans.all$def.exc <- NA
  ans.all$ced.lines <- F
  ans.all$categ.ans <- 0
  ans.all$cat.ans <- 0
  ans.all$reverse.ans <- 0
  ans.all$f.profile <- F
  ans.all$trace <- T
  ans.all$group <- 0
  ans.all$group.nr <- 1
  ans.all$nruns <- 0
  ans.all$constr.ans <- 1
  ans.all$no.hill <- F
  ans.all$no.CI <- F
  ans.all$one.sided <- F
  ans.all$interrupt <- TRUE
  if (gui) 
    ans.all$interrupt <- FALSE
  ans.all$rank.low <- 0
  ans.all$cc.fix <- NA
  ans.all$increase <- -1
  ans.all$increase.00 <- -1
  ans.all$increase.tmp <- 0
  ans.all$accepted <- T
  ans.all$vv <- 1
  ans.all$hill <- 0
  ans.all$no.refit <- F
  ans.all$WAPP <- FALSE
  ans.all$efsa.tool <- FALSE
  ans.all$simu <- FALSE
  ans.all$output <- TRUE
  ans.all$name.wapp <- NA
  ans.all$plotprefix <- NA
  color <- 1:100
  color <- color[!color %in% c(7, 15)]
  ans.all$color <- color
  ans.all$mark <- c(2, 4:25, 32:255)
  ans.all$x <- NA
  ans.all$y <- NA
  ans.all$yy <- NA
  ans.all$yans.denom <- 0
  ans.all$sd <- NA
  ans.all$mn.log <- 0
  ans.all$sd2.log <- 0
  ans.all$y.add <- 0
  ans.all$nr.dd <- 1
  ans.all$nr.cc <- 1
  ans.all$kk <- NA
  ans.all$ttt <- 0
  ans.all$ttt.txt <- ""
  ans.all$lump.ans <- 2
  ans.all$scores.orig <- 0
  ans.all$displ.ans <- 0
  ans.all$displ.fact <- 1
  ans.all$factor.name <- character()
  ans.all$covariate <- 1
  ans.all$fct1.txt <- ""
  ans.all$fct2.txt <- ""
  ans.all$fct3.txt <- ""
  ans.all$fct4.txt <- ""
  ans.all$fct5.txt <- ""
  ans.all$gr.txt <- ""
  ans.all$covar.txt <- ""
  ans.all$fct1.name <- ""
  ans.all$fct2.name <- ""
  ans.all$fct3.name <- ""
  ans.all$fct4.name <- ""
  ans.all$fct5.name <- ""
  ans.all$covar.name <- ""
  ans.all$twice <- F
  ans.all$tst.ans <- 2
  ans.all$par.start <- NA
  ans.all$eps <- 1e-12
  ans.all$fit.ans <- 1
  ans.all$fitted <- F
  ans.all$new.model <- F
  ans.all$scale.dum <- 0
  ans.all$text.par <- character()
  ans.all$text.par.0 <- NA
  ans.all$nrp <- 0
  ans.all$npar <- NA
  ans.all$adjust.start <- F
  ans.all$cond.ans <- 3
  ans.all$control.lst <- f.control(3)
  ans.all$loglik <- NA
  ans.all$low.yln <- NA
  ans.all$show <- ""
  ans.all$nl <- "\n"
  ans.all$tb <- "\t"
  ans.all$sp <- " "
  ans.all$sep <- FALSE
  ans.all$par.boot <- 0
  ans.all$boot <- F
  ans.all$xy.lim <- NA
  ans.all$y.lim.CI <- NA
  ans.all$replot <- F
  ans.all$pre.16 <- F
  ans.all$pi.full <- NA
  ans.all$x.full <- NA
  ans.all$fct1.full <- NA
  ans.all$fct2.full <- NA
  ans.all$kk.tot <- 0
  ans.all$nn.tot <- 0
  ans.all$regr.start <- numeric()
  ans.all$full6.done <- FALSE
  ans.all$th.start <- NA
  ans.all$th.0.start <- 0
  ans.all$sig.start <- 0.25
  ans.all$regr.par <- numeric()
  ans.all$regr.par.matr <- matrix()
  ans.all$th.par <- NA
  ans.all$sig.par <- NA
  ans.all$nth <- 1
  ans.all$decr.zz <- T
  ans.all$response.matr <- 0
  conf.int <- matrix(NA, ncol = 2)
  ans.all$conf.int <- conf.int
  ans.all$conf.lev <- 0.9
  ans.all$ab.txt <- "-"
  ans.all$quick.ans <- 1
  if (gui) 
    ans.all$quick.ans <- 3
  ans.all$constr <- -Inf
  ans.all$ab.txt <- ""
  ans.all$npar.full <- NA
  ans.all$npar.null <- NA
  ans.all$loglik.full <- NA
  ans.all$loglik.null <- NA
  ans.all$aic.full <- NA
  ans.all$aic.null <- NA
  ans.all$full.ans <- 1
  if (ans.all$gui) {
    ans.all$data <- data.frame()
    ans.all$data.is.from.file <- FALSE
    ans.all$last.upload.dir <- ""
    ans.all$data.name <- ""
    ans.all$comma.for.decimal <- f.use.comma.for.decimal()
  }
  ans.all$dtype.0 <- 0
  ans.all$summ.lst <- list()
  ans.all$notes <- "ATTENTION"
  ans.all$report.pars <- NA
  ans.all$alert.full <- NA
  ans.all$show.no.trend <- NA
  class(ans.all) <- c("list", "proast.session")
  if (exists("track")) 
    print("f.ini : END")
  return(ans.all)
}

###############
## f.version ##
###############
f.version <- function () 
{
  PRversion <- "65.5"
  PRversion
}

###############
## f.control ##
###############
f.control <- function (lev = 3) 
{
  lst <- list()
  switch(lev, {
    lst$eval.max <- 50
    lst$iter.max <- 40
    lst$rel.tol <- 0.001
    lst$x.tol <- 0.015
    lst$step.min <- 0.00022
  }, {
    lst$eval.max <- 100
    lst$iter.max <- 75
    lst$rel.tol <- 1e-06
    lst$x.tol <- 0.00015
    lst$step.min <- 2.2e-07
  }, {
    lst$eval.max <- 1000
    lst$iter.max <- 750
  })
  return(lst)
}

#####################
## f.graphwin.size ##
#####################
f.graphwin.size <- function (resize = F) 
{
  if (exists("track")) 
    print("f.graphwin.size")
  if (!exists(".gw.size")) {
    size.x <- 5
    size.y <- 5
    f.assign(".gw.size", c(size.x, size.y))
  }
  if (resize) {
    f.create.graphwin(.gw.size[1], .gw.size[2])
    oke <- F
    while (!oke) {
      gw.x <- .gw.size[1]
      gw.y <- .gw.size[2]
      cat("\ncurrent values for horizontal and vertical size of the graphical window are:\n", 
          gw.x, gw.y, "\n")
      gw.x <- eval(parse(prompt = paste("\nGive value for horizontal size > ")))
      gw.y <- eval(parse(prompt = paste("Give value for vertical size > ")))
      f.create.graphwin(gw.x, gw.y)
      oke <- 1 - menu(c("no", "yes"), title = "\nIs this window size okay?")
      lst <- dev.list()
      curr <- dev.cur()
      lst <- lst[lst != curr]
      for (ii in lst) dev.off(ii)
    }
    assign(".gw.size", c(gw.x, gw.y), pos = 1, immediate = T)
  }
}

##############
## f.assign ##
##############
f.assign <- function (name, obj) 
{
  if (exists("track")) 
    print("f.assign")
  obj$date <- date()
  if (!is.character(name)) 
    print("f.assign:  names and object are reversed")
  assign(name, obj, pos = 1, immediate = T)
  if (exists("track")) 
    print("f.assign:   END")
}

###########
## f.con ##
###########
f.con <- function (ans.all, list.logic = F) 
{
  if (exists("track")) 
    print("f.con")
  f.assign(".Pr.last", ans.all)
  if (list.logic) {
    cat("\n You have chosen previous results concerning:   ")
    cat(ans.all$odt.name)
    PRversion <- ans.all$PRversion
    ans.all$PRversion <- ans.all$version.old
    if (ans.all$model.fam > 0) {
      cat("\nfitted model: ")
      cat(ans.all$modelname, "\n")
      cat("loglik:", ans.all$loglik, "\n")
      ans.all$show <- f.show.con(ans.all)
      cat(ans.all$show)
    }
    if (ans.all$boot) {
      cat("\nBootstrap data are available")
      nruns <- length(ans.all$CED.boot[, 1])
      cat("\nnruns:", nruns, "\n")
    }
    if (ans.all$model.ans %in% 41:42 && ans.all$nr.cc > 1) {
      with(ans.all, {
        npar <- length(MLE)
        Vcc <- MLE[(npar - nr.dd + 1):npar]
        Vsd <- sqrt(MLE[1:nr.var])
        f.graph.window(1)
        plot(log10(Vcc), log10(Vsd), xlab = "log.sd", 
             ylab = "log.c")
        f.press.key.to.continue()
      })
    }
    if (ans.all$nr.cc > 1) {
      ans.all <- with(ans.all, {
        npar <- length(MLE)
        Vcc <- MLE[(nr.var + nr.aa + nr.bb + 1):(nr.var + 
                                                   nr.aa + nr.bb + nr.cc)]
        Vsd <- sqrt(MLE[1:nr.var])
        if (length(Vsd) == length(Vcc)) {
          gr.txt <- f.pars(ans.all)$gr.txt
          gr.lab <- letters[1:length(gr.txt)]
          show <- "subgroups"
          for (i in 1:length(gr.txt)) show <- paste(show, 
                                                    nl, gr.lab[i], tb, gr.txt[i])
          f.graph.window(1)
          plot(Vsd, log(Vcc), pch = as.character(gr.lab))
          abline(0, 7)
          loc.txt <- max(log(Vcc))
          mtext(show, 4, 1, adj = 0, padj = 1, las = 1, 
                at = loc.txt, cex = 0.7)
          cat("\n\n\nSee plot for relationship between c and sd\n\n")
          f.press.key.to.continue()
          if (model.ans == 6 && fitted) 
            ans.all <- f.fit.model6(ans.all)
        }
        return(ans.all)
      })
    }
    ans.all$par.start <- ans.all$MLE
    if (ans.all$plot.type > 8) 
      ans.all$plot.type <- 2
    if (ans.all$model.fam == 0) 
      f.plot.gui(ans.all)
    else {
      ans.all$show <- f.show.con(ans.all)
      f.plot.all(ans.all)
    }
    ans.all$PRversion <- PRversion
  }
  main.ans <- 2
  if (list.logic) 
    main.ans <- 14
  if (ans.all$model.ans > 0) {
    if (!ans.all$fitted) 
      ans.all <- f.start.con(ans.all)
    main.ans <- 14
  }
  while (T) {
    switch(main.ans, {
      ans.all$change <- c(rep(F, ans.all$nrQ + 1))
      ans.all <- f.change.settings(ans.all, choose = T)
      if (ans.all$fitted) ans.all$show <- f.show.con(ans.all)
      if (!ans.all$fitted) {
        if (ans.all$model.ans == 0) ans.all <- f.choose.model(ans.all)
        if (ans.all$model.ans %in% 38:39) ans.all <- f.quick.con(ans.all) else {
          ans.all <- f.execute(ans.all)
          if (ans.all$ans.m6.sd != 2) ans.all <- f.start.con(ans.all, 
                                                             adjust = F, fitted = F, tmp.quick = F)
          f.plot.all(ans.all)
        }
      } else f.plot.all(ans.all)
      model.names <- f.expect.con(name = TRUE)
      ans.all$modelname <- model.names[ans.all$model.ans]
    }, {
      ans.all <- f.clear(ans.all)
      if (ans.all$quick.ans %in% c(1, 3)) {
        ans.all <- f.choose.model(ans.all)
        if (ans.all$model.ans %in% c(48, 50)) {
          ans.all$change[24] <- T
          ans.all <- f.change.settings(ans.all)
          ans.all <- f.execute(ans.all)
        }
        if (ans.all$model.ans == 6) ans.all$ans.m6.sd <- 1
        ans.all$new.model <- T
      }
      if (ans.all$quick.ans == 1) {
        f.graph.window(1)
        ans.all <- f.start.con(ans.all)
        if (ans.all$adjust.start) {
          ans.all <- f.start.con(ans.all, adjust = T, 
                                 fitted = T, tmp.quick = F)
        }
      }
      if (ans.all$quick.ans > 1) {
        ans.all <- f.quick.con(ans.all)
        ans.all$quick.ans <- 1
        #print (paste("length(ans.all$Vyans): ", length(ans.all$Vyans), "; length(ans.all$HILL): ", length(ans.all$HILL), "; ans.all$EXP:", ans.all$EXP, sep=""))
        if (length(ans.all$Vyans) == 1 && length(ans.all$HILL) > 
            0) 
          model.fam <- 0     
          #model.fam <- menu(c("exponential", "Hill"),
               #                  title = "\nWhich model do you want to continue with?\n                           (type 0 if both)") else model.fam <- 1
        switch(model.fam, ans.all <- f.move.sublist(ans.all, 
                                                    ans.all$EXP), ans.all <- f.move.sublist(ans.all, 
                                                                                            ans.all$HILL))
        switch(model.fam, ans.all$heading <- "Exponential model", 
               ans.all$heading <- "Hill model")
        ans.all$model.fam <- model.fam
        ans.all$modelname <- ans.all$model.names[ans.all$model.ans]
        if (length(ans.all$Vyans) == 1) #f.plot.gui(ans.all) # commented out to avoid plotting figure
        ans.all$ans.plt <- 0
      }
    }, {
      if (ans.all$fitted & !ans.all$new.model) ans.all$par.start <- ans.all$MLE
      ans.dd <- 2
      if (max(ans.all$fct5) > 1) ans.dd <- menu(c("yes", 
                                                  "no"), title = "\nDo you want to use previous MLE from model with one d as start value?\n")
      if (ans.dd == 1) {
        nr.dd <- max(ans.all$fct5)
        dd.start <- ans.all$MLE[length(ans.all$MLE)]
        ans.all$par.start <- c(ans.all$MLE, rep(dd.start, 
                                                nr.dd - 1))
        ans.all$scale.dum <- c(ans.all$scale.dum, rep(1, 
                                                      nr.dd - 1))
      }
      ans.all <- f.start.con(ans.all, adjust = T, fitted = F, 
                             tmp.quick = F)
    }, {
      if (ans.all$model.ans == 0) cat("\n\nfirst choose a model before fitting\n") else {
        ans.all <- f.mm4.con(ans.all)
      }
    }, {
      ans.all <- f.mm5.con(ans.all)
    }, {
      if (ans.all$fitted == F) cat("\nFirst fit the model using option 4\n") else if (ans.all$model.ans == 
                                                                                      11) cat("\nCED not defined for full model\n") else if (ans.all$model.ans == 
                                                                                                                                             1) cat("\nCED not defined for null model\n") else ans.all <- f.mm6.con(ans.all)
    }, {
      if (is.na(ans.all$CED[1])) cat("\nYou did not calculate CED! First use option 6 from main menu\n") else {
        ans.all$nruns <- eval(parse(prompt = "\ngive the number of bootstrap runs you want to be done > "))
        ans.all$plot.ans <- 1
        ans.all$plot.ans <- menu(c("no", "yes (reduces speed)"), 
                                 title = "\ndo you want plots of each sampling run?")
        ans.all <- f.mm7.con(ans.all)
      }
    }, {
      ans.all <- f.mm8.con(ans.all)
    }, {
      with(ans.all, {
        if ((is.na(CED[1])) == T) cat("\nYou did not calculate CED! First use option 6 from main menu\n") else if (length(CED.boot) == 
                                                                                                                   0) cat("\nno bootstrap data available!\nfirst choose option 7 from main menu\n") else {
                                                                                                                     group.nr <- 1
                                                                                                                     if (nr.bb > 1) group.nr <- eval(parse(prompt = paste("\nWhich group of ", 
                                                                                                                                                                          varnames[ans.all$fct2.no], "\ndo you want to consider? > ")))
                                                                                                                     f.RfD.con(CED.boot, CES, sf.x = sf.x, y.leg = y.leg, 
                                                                                                                               group.nr = group.nr, factor.name = varnames[ans.all$fct2.no], 
                                                                                                                               CED.show = show)
                                                                                                                   }
      })
    }, {
      if (length(ans.all$CED.boot) == 0) {
        cat("\n no bootstrap data available! \nfirst choose option 7 from main menu\nif you want probabilistic analysis\n")
        with(ans.all, {
          cat("\nGive the exposure level you want to consider\nin the same units as used in the animal study")
          act.exp <- eval(parse(prompt = "\n(or type 0 to return to the Main Menu) > "))
          if (act.exp != 0) {
            group.nr <- 1
            if (nr.bb > 1) {
              group.nr <- eval(parse(prompt = paste("\nWhich group of ", 
                                                    varnames[ans.all$fct2.no], "\ndo you want to consider? > ")))
              regr./par <- regr.par.matr[group.nr, ]
            }
            ES.anim <- f.es.con(model.ans, act.exp/sf.x, 
                                regr.par, CES = CES)
            ES.anim <- signif(ES.anim, 6)
            cat("the estimated effect size in the animal\nat exposure level ", 
                act.exp, "is: \n", ES.anim, "\n\n")
          }
        })
      } else ans.all <- f.mm10.con(ans.all)
    }, {
      with(ans.all, {
        if (model.ans != 2) {
          cat("\nonly model 2 implemented\n")
          return(invisible())
        }
        regr.par.matr <- f.pars(ans.all)$regr.par.matr
        for (ii in 1:length(regr.par.matr[, 1])) {
          aa <- logb(regr.par.matr[ii, 1])
          bb <- regr.par.matr[ii, 2]
          if (max(fct3) == 1) SD <- sqrt(MLE[1]) else SD <- sqrt(MLE[ii])
          ADI <- eval(parse(prompt = paste("give value for ADI (or MRL)", 
                                           "  > ")))
          kkk <- eval(parse(prompt = paste("give value for k, in tol.lim = x + k*sd", 
                                           "  > ")))
          WP <- (logb(ADI) - kkk * SD - aa)/bb
          cat("\n\nThe withdrawal period is:", WP, "\n")
        }
      })
    }, {
      with(ans.all, {
        if (model.fam == 0 && !is.null(ans.all$HILL)) {
          cat("\nYou first need to choose between the exponential and Hill model\n")
          cat("Close the PROAST session, and restart it using the object with the saved results\n\n")
          f.press.key.to.continue()
        } else {
          cat(paste(1:npar, ": ", text.par[1:npar], "\n"))
          cat("Give number(s) of the (consecutive) parameter(s)\n")
          ans.all$group <- eval(parse(prompt = paste(" -------- > ")))
          MLE <- ans.all$MLE
          if (any(MLE[ans.all$group] < 0)) ans.all$nolog <- T else ans.all$nolog <- F
          conf.int.0 <- ans.all$conf.int
          ans.all$trace <- T
          ans.all <- f.CI(ans.all)
          name.ci <- text.par[ans.all$group]
          confint.spec <- ans.all$conf.int/sf.x
          dimnames(confint.spec) <- list(name.ci, NULL)
          ans.all$conf.int <- conf.int.0
          ans.all$confint.spec <- confint.spec
          cat("\n\n the confidence intervals are:\n")
          print(confint.spec)
          f.store.results(ans.all)
        }
      })
    }, {
      main.ans <- 13
      return(T)
    })
    #cat("\n\n--------\nMAIN MENU :\n--------\n\n")
    main.ans <- 13
    #main.ans <- menu(c("Change settings", "Choose (another) model", 
    #                   "Choose other startvalues", "Fit model", "Plot results", 
    #                   "Calculate CED: point estimate or confidence interval", 
    #                   "Generate bootstrap runs based on fitted model", 
    #                   "Calculate CED distribution for the animal", "", 
    #                   "", "Calculate withdrawal period for veterinary medicines", 
    #                   "Calculate confidence interval for model parameter(s)", 
    #                   "End session"), title = "What do you want to do ? \n")
  }
}

#############
## f.clear ##
#############
f.clear <- function (ans.all) 
{
  if (exists("track")) 
    print("f.clear")
  ans.all$show <- ""
  ans.all$MLE <- NA
  ans.all$regr.par <- NA
  ans.all$regr.par.matr <- NA
  ans.all$CED.matr <- matrix(NA)
  ans.all$regr.par.matr <- matrix(NA)
  ans.all$conf.int <- matrix(NA, ncol = 2)
  ans.all$fitted <- F
  ans.all$ans.plt <- 3
  ans.all$par.start <- NA
  ans.all$scale.dum <- NA
  ans.all$lb <- NA
  ans.all$ub <- NA
  ans.all$lower <- NA
  ans.all$upper <- NA
  ans.all$loglik <- NA
  if (ans.all$quick.ans > 1) {
    ans.all$Vbmdl <- NA
    ans.all$Vbmdu <- NA
  }
  if (exists("track")) 
    print("f.clear:   END")
  return(ans.all)
}

#################
## f.quick.con ##
#################
f.quick.con <- function (ans.all) 
{
  if (exists("track")) 
    print("f.quick.con")
  if (ans.all$WAPP) {
    ans.all$gui <- TRUE
    ans.all$interrupt <- FALSE
  }
  ans.all <- with(ans.all, {
    if (gui) {
      ans.all$quick.ans <- 2
      if (!WAPP && !simu) {
        cat("\nclick in this window to see current output\n")
        f.press.key.to.continue()
      }
      ans.all$data.0 <- odt
      f.overlap()
      Vyans <- yans
      if (dtype %in% c(10, 15)) {
        Vsans <- sans
        Vnans <- nans
      }
      ans.all$const.var <- FALSE
      PRversion <- f.version()
      if (length(ans.all$PRversion) == 0) 
        version.old <- ans.all$version
      else version.old <- ans.all$PRversion
      if (version.old[1] != PRversion) 
        ans.all <- f.adjust.saved(ans.all)
      if (quick.ans == 3) 
        ans.all$quick.ans <- 2
      if (ans.all$WAPP && ans.all$nruns == 0 && dtype %in% 
          c(5, 15)) 
        ans.all$nruns <- 200
    }
    if (!gui) {
      varnames <- ans.all$varnames
      nvar <- ans.all$nvar
      dtype <- ans.all$dtype
      odt <- ans.all$odt
      if (ans.all$xans == 0) {
        ans.all$change[1] <- T
        ans.all <- f.change.settings(ans.all)
      }
      #cat("\n")
      #cat(paste(1:nvar, varnames[1:nvar], "\n"))
      #cat("Give number(s) of the response(s) you want to analyse\n")
      Vyans <- respCol
      #Vyans <- eval(parse(prompt = paste(" -------- > ")))
      zero.ind <- F
      if (length(Vyans) > 1) 
        zero.ind <- any(data.0[, Vyans] == 0, na.rm = T)
      if (dtype %in% c(10, 15, 250)) {
        cat("\n")
        cat(paste(1:nvar, varnames[1:nvar], "\n"))
        cat("Give number(s) of the variation statistic associated to these endpoints\n")
        Vsans <- eval(parse(prompt = paste(" -------- > ")))
        if (length(Vsans) == 1 && length(Vyans) > 1) 
          Vsans <- rep(Vsans, length(Vyans))
        sd.se <- menu(c("standard deviations", "standard errors"), 
                      title = "\nDo you have standard deviations or standard errors associated \n                          with the means?")
        cat("\n")
        cat(paste(1:nvar, varnames[1:nvar], "\n"))
        cat("Give number(s) of the groups size associated to these endpoints\n")
        Vnans <- eval(parse(prompt = paste(" -------- > ")))
        if (length(Vnans) == 1 && length(Vyans) > 1) 
          Vnans <- rep(Vnans, length(Vyans))
        ans.all$Vsans <- Vsans
        ans.all$sd.se <- sd.se
        ans.all$Vnans <- Vnans
      }
      ans.all$covar.dd <- 0
      if (dtype == 5) 
        ans.all$nest.no <- menu(varnames[1:nvar], title = "\nGive number of the nested factor \n (type 0 to skip nested analysis)")
      if (length(ans.all$xans) == 1) 
        covar.no <- 0
        #covar.no <- menu(varnames[1:nvar], title = "\nGive number of factor serving as potential covariate (e.g.sex)\n                    type 0 if none --- \n  ")
      else covar.no <- 0
      ans.all$covar.no <- covar.no
      if (0) 
        if (ans.all$quick.ans == 4 && covar.no > 0) {
          covar.dd <- menu(c("no", "yes"), title = "\nDo you want parameter d to be covariate dependent? \n  ")
          covar.dd <- covar.dd - 1
          ans.all$covar.dd <- covar.dd
          ans.all$CES <- -1
          fix.ans <- 2
          cc.fix <- NA
          if (fix.ans == 1) 
            ans.all$cc.fix <- eval(parse(prompt = "\n\ngive fixed value for parameter c >   "))
        }
      if (!ans.all$covar.dd) {
        if (length(Vyans) > 1) 
          ans.all$NES.ans <- menu(c("no", "yes"), title = "\nDo you want to adjust CES to within group SD?\n")
        if (ans.all$NES.ans == 1) 
          ans.all$CES <- inCES
          #ans.all$CES <- eval(parse(prompt = paste("\nGive value for CES (always positive)\n              type 0 to avoid calculation of CIs  > ")))
        else ans.all$CES <- 0.05
      }
      ans.all$no.hill <- 0 # need to use 0 to get HILL model
      #ans.all$no.hill <- menu(c("yes", "no"), title = "\nDo you want to include Hill model family?") - 
        #1
      if (length(Vyans) == 1) 
        interrupt <- TRUE
      if (length(Vyans) > 1) {
        int.ans <- menu(c("yes (opportunity to save results per endpoint)", 
                          "no"), title = "\nDo you want to interrupt calculations after each endpoint? \n")
        if (int.ans == 1) 
          interrupt <- TRUE
        else interrupt <- FALSE
      }
      ans.all$interrupt <- interrupt
      ans.all$auto.detlim <- F
      if (zero.ind) 
        if (length(Vyans) > 1) {
          auto.ans <- menu(c("yes", "apply automatic detection limit"), 
                           title = "\nDo you want to provide detection limit for each separate response, or automically use smallest nonzero observation?\n")
          switch(auto.ans, ans.all$auto.detlim <- F, 
                 ans.all$auto.detlim <- T)
        }
    }
    ans.all$quick <- T
    ans.all$no.CI <- F
    if (ans.all$CES == 0) {
      ans.all$CES <- 0.05
      ans.all$no.CI <- T
    }
    ans.all$CES.0 <- ans.all$CES
    if (dtype %in% c(5, 15) && !ans.all$no.CI && !ans.all$WAPP) 
      ans.all$nruns <- eval(parse(prompt = "\n\ngive the number of bootstrap runs you want to be done for calculating CED CIs (e.g.200) > "))
    if (covar.no > 0) {
      covariate <- odt[, covar.no]
      covariate <- covariate[!is.na(covariate)]
      covar.txt <- levels(as.factor(covariate))
    }
    else covar.txt <- character(1)
    CED.matr <- matrix(NA, ncol = (2 * length(covar.txt)))
    if (ans.all$NES.ans == 2) {
      Mces <- matrix(ncol = 3)
      Msd <- matrix(ncol = 3)
    }
    else {
      Mces <- NA
      Msd <- NA
    }
    ans.all$CI.plt <- FALSE
    if (dtype %in% c(4, 10)) 
      ans.all$CI.plt <- FALSE
    ans.all$cont <- f.cont(ans.all)
    if (ans.all$cont) {
      ans.all$plot.type <- 4
      ans.all$plt.mns <- 3
    }
    else ans.all$plot.type <- 2
    CED.all <- list(CED.matr = CED.matr, endpoints = character(), 
                    covar.txt = covar.txt, CES = ans.all$CES, Mces = Mces, 
                    Msd = Msd, NES.ans = ans.all$NES.ans, no.hill = ans.all$no.hill, 
                    quick.ans = quick.ans, PRversion = PRversion, date = date())
    ans.all$model.type <- 1
    ans.all$model.list <- c("EXP", "HILL")
    data.0.tmp <- ans.all$data.0
    first.loop <- T
    ii <- 1
    for (yans.tmp in Vyans) {
      ans.all$yans <- yans.tmp
      if (dtype %in% c(10, 15, 250, 260)) {
        ans.all$sans <- Vsans[ii]
        ans.all$nans <- Vnans[ii]
      }
      ans.all$data.0 <- f.remove.NAs(xans = ans.all$xans, 
                                     yans = ans.all$yans, sans = ans.all$sans, nans = ans.all$nans, 
                                     covar.no = 0, dfr = data.0.tmp)
      response <- ans.all$data.0[, yans.tmp]
      skip <- f.check.nonneg.num(response, gui, dtype, 
                                 quick.ans = 2)
      if (skip) 
        cat("\n\nATTENTION: analysis not possible for response", 
            ans.all$varnames[yans.tmp], "\n")
      date.tmp <- date()
      if (!skip) {
        if (length(Vyans) > 1) {
          ans.all$detlim <- 0
          ans.all$Vdetlim <- 0
          ans.all$detlim.col <- 0
          ans.all$covar.no <- covar.no
        }
        ans.all <- f.execute(ans.all, no.plot = T)
        ans.all <- f.clear(ans.all)
        ans.all$twice <- T
        if (first.loop) 
          ans.all$full.ans <- f.full.ans(ans.all, gui = gui)
        first.loop <- F
        if (!gui && yans.tmp == Vyans[1]) 
          if (interrupt) {
            ans.all$lower.dd <- 0.25
            ans.all$upper.dd <- 4
            #cat(paste("\n\nATTENTION: the constraints on parameter d are set at", 
            #          ans.all$lower.dd, "and", ans.all$upper.dd, 
            #          "\n"))
            f.press.key.to.continue()
          }
        if (output) {} 
          #cat(paste("\n\nresponse: ", ans.all$y.leg))
        ans.all$hill <- 0
        ans.all$increase <- 0
        if (ans.all$NES.ans == 2) 
          ans.all$CES <- 0.05
        ans.all <- f.select.m5.con(ans.all)
        if (!ans.all$no.hill) {
          if (output) {}
            #cat(paste("\n\nresponse: ", ans.all$y.leg))
          ans.all$hill <- 10
          if (ans.all$NES.ans == 2) 
            ans.all$CES <- 0.05
          ans.all <- f.select.m5.con(ans.all)
        }
        ans.all$fitted <- T
        bmdl.lowest <- numeric()
        bmdu.highest <- numeric()
        if (!ans.all$no.CI && length(ans.all$HILL) > 
            0) 
          if (!is.na(ans.all$EXP$conf.int[1, 1]) && !is.na(ans.all$HILL$conf.int[1, 
                                                                                 1])) 
            if (identical(ans.all$EXP$fct2, ans.all$HILL$fct2)) {
              nr.gr <- max(ans.all$EXP$fct2)
              if (length(ans.all$xans) > 1) 
                nr.gr <- length(ans.all$xans)
              for (qq in (1:nr.gr)) {
                bmdl.lowest[qq] <- min(ans.all$EXP$conf.int[qq, 
                                                            1], ans.all$HILL$conf.int[qq, 1])
                bmdu.highest[qq] <- max(ans.all$EXP$conf.int[qq, 
                                                             2], ans.all$HILL$conf.int[qq, 2])
              }
              if (output) {}
                #cat("\n\n -----------  CES = ", ans.all$CES.0, 
                #    " --------------------------------- ")
              if (ans.all$NES.ans == 2) 
                if (output) {}
                  #cat("\n-----------  endpoint-specific NES = ", 
                  #    round(ans.all$CES, 4), " --------------------------------- ")
              ans.all$bmdl.lowest <- bmdl.lowest
              ans.all$bmdu.highest <- bmdu.highest
              if (length(bmdl.lowest) == 1) 
                covar.txt <- "all"
              bmdCI.overall <- rbind(bmdl.lowest, bmdu.highest)
              bmdCI.overall <- t(bmdCI.overall)
              if (length(ans.all$xans) > 1) 
                covar.txt <- dimnames(ans.all$EXP$conf.int)[[1]]
              bmdCI.overall <- cbind(covar.txt, bmdCI.overall)
              bmdCI.overall <- as.data.frame(bmdCI.overall)
              ans.all$bmdCI.overall <- bmdCI.overall
            }
        if (output) {
          #cat("\nThe lowest BMDL and highest BMDU from exponential and Hill models are:\n")
          #print(ans.all$bmdCI.overall)
          #cat(" --------------------------------------------------- ")
        }
        if (ans.all$model.ans != 1 && ans.all$no.CI != 
            T && length(Vyans) > 1) {
          CED.all$covar.txt <- ans.all$covar.txt
          CED.all <- f.CED.all(ans.all$EXP, ans.all$HILL, 
                               CED.all, ans.all$y.leg)
        }
        assign("last.fit", ans.all, immed = T, pos = 1)
        if (ans.all$fitted) {
          if (ans.all$dtype %in% c(1, 5)) 
            ans.all$plt.mns <- 1
          #f.plot.gui(ans.all) # commented out to avoid plotting
        }
        if (length(ans.all$xans) > 1) {
          ans.all$conf.int <- f.confint.doseaddition(ans.all)
          cat("\n the confidence intervals (per model) are:\n")
          print(ans.all$conf.int)
        }
        if (!gui) {}
          #cat(paste("end of analysis for response: ", 
          #          ans.all$y.leg))
        if (gui) 
          ans.all$date <- date()
        if (!gui && interrupt) {
          if (length(ans.all$HILL) == 0) 
            ans.all <- f.move.sublist(ans.all, ans.all$EXP)
          ans.all$CED.all <- CED.all
          store.name <- f.store.results(ans.all)
          if (!gui && length(ans.all$notes) > 1) {
            print("f.quick.con")
            cat(ans.all$notes)
            f.press.key.to.continue()
          }
        }
      }
      ii <- ii + 1
    }
    if (ans.all$no.CI == F && length(Vyans) == 1 && length(ans.all$EXP$conf.int[, 
                                                                                1]) > 1) {
      f.plot.CI(ans.all)
    }
    if (ans.all$no.CI == F & length(Vyans) > 1) {
      if (sum(!is.na(CED.all$CED.matr[, 1]) > 0)) {
        dimn.col <- c("BMDL", "BMDU")
        dimn.col <- rep(dimn.col, length(ans.all$covar.txt))
        dimnames(CED.all$CED.matr)[[2]] <- dimn.col
        CED.all$CED.matr <- CED.all$CED.matr[-1, ]
        dimnames(CED.all$CED.matr)[[1]] <- CED.all$endpoints
        dimnames(CED.all$Mces)[[1]] <- NULL
        CED.all$Mces <- as.data.frame(CED.all$Mces[-1, 
                                                   ])
        dimnames(CED.all$Msd)[[1]] <- NULL
        CED.all$Msd <- as.data.frame(CED.all$Msd[-1, 
                                                 ])
        CED.all$max.x <- max(ans.all$x)
        lst.name <- 0
        if (!gui) 
          lst.name <- eval(as.character(parse(prompt = "\n\ngive name for file to store CED table as R object (or 0 to avoid) > ")))
        else lst.name <- "gui.CEDtable"
        if (lst.name != 0) 
          f.assign(eval(lst.name), CED.all)
        f.plot.CED(CED.all)
      }
      else cat("\n\nNone of the endpoints showed a significant trend\n")
    }
    ans.all$Vyans <- Vyans
    ans.all$NES.ans <- 1
    if (WAPP) 
      graphics.off()
    if (output) {
      #cat("\ncalculations started at:", date.tmp)
      #cat("\nand ended at:", date())
      #cat("\n")
    }
    if (exists("track")) 
      print("f.quick.con:  END")
    return(ans.all)
  })
}

#######################
## f.change.settings ##
#######################
f.change.settings <- function (ans.all, choose = F) 
{
  if (exists("track")) 
    print("f.change.settings")
  ans.all$choose <- choose
  with(ans.all, {
    cont <- f.cont(ans.all)
    if (quick.ans == 1 && change[2]) {
      if (dtype %in% c(10, 15, 250, 260)) 
        if (sans[1] == 0) 
          change[3] <- T
        if (dtype %in% c(4, 6, 84) && change[2]) 
          if (nans[1] == 0) 
            change[8] <- T
    }
    while (choose || sum(change[1:(nrQ - 1)]) != 0) {
      if (choose) {
        settings <- f.define.settings(quick.ans, dtype)
        set.ans <- menu(settings, title = "\nWhich setting(s) do you want to change?")
        if (!(set.ans %in% c(9, 18, 23, nrQ))) 
          ans.all$fitted <- F
        change[set.ans] <- T
      }
      if (change[1]) {
        #cat("\nQ1: Which variable do you want to consider as independent variable?\n    (e.g. dose, age)\n")
        #cat(" ")
        #cat(paste(1:nvar, ":", varnames, "\n"))
        #xans <- eval(parse(prompt = paste(" -------- > ")))
        #ans.all$xans <- xans
        ans.all$xans <- concCol
        ans.all$change <- change
        change[1] <- F
        if (length(xans) > 1) {
          ans.all$displ.no <- menu(varnames[1:nvar], 
                                   title = "\nQ9: Which column defines the single doses and mixtures?\n    (type 0 if none)")
        }
      }
      if (change[2]) {
        yans <- menu(varnames[1:nvar], title = "\nQ2: Which response variable do you want to analyse ? \n")
        ans.all$yans <- yans
        ans.all$change <- change
        change[2] <- F
        if (dtype %in% c(4, 6, 84)) 
          change[8] <- T
        if (dtype %in% c(10, 15)) 
          change[3] <- T
        if (dtype %in% c(5)) {
          ans.all$nest.no <- menu(varnames[1:nvar], title = "\nGive column number of the nested factor \n                         (type 0 to skip nested analysis)")
        }
        ans.all$detlim <- 0
        ans.all$Vdetlim <- 0
        ans.all$detlim.col <- 0
        ans.all$reverse.ans <- 0
      }
      if (change[3]) {
        if (dtype == 3) {
          cat.ans <- menu(1:nth, title = paste("\nQ3:  give severity category for dichotimization \n"))
          change[3] <- F
          ans.all$cat.ans <- cat.ans
        }
        else {
          sans <- menu(varnames[1:nvar], title = paste("\nQ3a:  number of the variation statistic associated to", 
                                                       varnames[yans], " \n"))
          sd.se <- menu(c("standard deviations", "standard errors"), 
                        title = "\nQ3b: Do you have standard deviations or standard errors associated with the means?")
          nans <- menu(varnames[1:nvar], title = "\nQ3c: Give the associated sample size\n")
          change[3] <- F
          ans.all$sans <- sans
          ans.all$sd.se <- sd.se
          ans.all$nans <- nans
        }
      }
      if (change[4]) {
        go.on <- 1
        count <- 0
        while (go.on > 0) {
          remove.no <- menu(varnames[1:nvar], title = paste("\nQ4: If you want to remove outliers, type number of factor to indicate them", 
                                                            "\n--- type 0 if you do not want to remove any (more) data --- \n"))
          if (remove.no == 0) {
            go.on <- 0
            tmp.lev <- "none"
          }
          else {
            count <- count + 1
            remove.fact <- data.0[, remove.no]
            if (!is.numeric(remove.fact)) {
              levels.all <- sort(levels(factor(remove.fact)))
              nr.lev <- length(levels.all)
              cat("\nQ4b: Give number(s) of the level(s) you want to remove\n")
              cat(" ")
              cat(paste(1:nr.lev, ":", levels.all, "\n"))
              tmp.no <- eval(parse(prompt = paste(" -------- > ")))
              tmp.lev <- levels.all[tmp.no]
            }
            if (is.numeric(remove.fact)) {
              levels.all <- unique(remove.fact)
              cat("\n the values of this variable are: \n")
              print(sort(levels.all))
              cat("\nQ4b: Give values(s) you want to remove\n")
              tmp.lev <- eval((parse(prompt = paste("      non-numeric levels between quotes! > "))))
              tmp.lev <- levels.all[levels.all %in% tmp.lev]
            }
            ans.all$remove.lev <- sort(tmp.lev)
            ans.all$remove.no <- remove.no
            ans.all$data.0 <- data.0
            data.0 <- f.subset(ans.all, remove = T)
            ans.all$data.0 <- data.0
            remove.tmp <- paste(paste(varnames[remove.no], 
                                      ":", paste(tmp.lev, collapse = " ")), sep = "\n")
            if (remove.name == "none") 
              remove.name <- character()
            remove.name <- paste(remove.name, "\n  ", 
                                 remove.tmp)
          }
          if (count == 0) {
            cat("\nATTENTION: complete dataset will be restored\n")
            data.0 <- odt
            data.0 <- f.remove.NAs(xans, yans, sans, 
                                   dfr = data.0)
            if (length(unique(Vdetlim) == 1)) {
              Vdetlim <- rep(Vdetlim[1], length(data.0[, 
                                                       1]))
              data.0 <- cbind(data.0, Vdetlim)
            }
            ans.all$data.0 <- data.0
            remove.name <- "none"
          }
          ans.all$remove.name <- remove.name
        }
        change[4] <- F
      }
      if (change[5]) {
        yans.denom <- menu(varnames[1:nvar], title = "\nWhich variable do you want as denominator ? (0 = NONE)\n")
        change[5] <- F
        ans.all$yans.denom <- yans.denom
      }
      if (change[6]) {
        if (yans != 0) {
          y <- data.0[, yans]
          #print(sort(y))
          #cat("\nthe lowest nonzero value is: ", min(y[y > 
          #                                               0]), "\n")
          if (min(y) > 0) 
            cat("\n\nNOTE: lowest observation is positive \n\n")
        }
        detlim <- eval(parse(prompt = paste("\nQ6: Give global value of detection limit, \nor type 0 for individual values > ")))
        if (detlim == 0) 
          detlim.col <- menu(varnames[1:nvar], title = "\nGive column with individual detection limits\n")
        if (detlim.col == 0) {
          ans.all$Vdetlim <- rep(detlim, length(x))
          ans.all$data.0 <- cbind(data.0, Vdetlim)
          ans.all$detlim.col <- nvar + 1
          ans.all$nvar <- nvar + 1
          ans.all$varnames[ans.all$nvar] <- "detlim"
        }
        change[6] <- F
      }
      if (change[7]) {
        y <- data.0[, yans]
        y <- y[!is.na(y)]
        cat("\nthe lowest nonzero observation is: ", 
            min(y[y > 0]), "\n")
        y.add <- eval(parse(prompt = "\n\nQ7: Give small value to be added to the data > "))
        change[7] <- F
        ans.all$y.add <- y.add
      }
      if (change[8]) {
        ans.all$nans <- menu(varnames[1:nvar], title = "\nQ3b: Give the associated sample sizes")
        if (dtype == 84) 
          ans.all$tans <- menu(varnames[1:nvar], title = "\nQ3c: Give number of time variable")
        change[8] <- F
      }
      if (change[9]) {
        displ.no <- menu(varnames[1:nvar], title = "\nQ9: Which variable do you want to consider for distinct plotting?\n    (type 0 if none)")
        change[9] <- F
      }
      if (change[10]) {
        if (ans.all$xans != 0) {
          x <- data.0[, xans]
          cat("\n(Nonzero) doses range from ", min(x[x != 
                                                       0], na.rm = T), "to ", max(x))
        }
        ans.all$sf.x <- eval(parse(prompt = "\nQ10: Give scaling factor for dose > "))
        if (tans > 0) {
          ttt <- data.0[, tans]
          cat("\n\n\n(Nonzero) times range from ", min(ttt[ttt != 
                                                             0]), "to ", max(ttt))
          sf.t <- eval(parse(prompt = "\nQ10: Give scaling factor for time > "))
          ans.all$sf.t <- sf.t
        }
        change[10] <- F
        if (model.ans == 46) {
          cat("\nscaling factor on dose not implemented for model 46, it will be ignored\n")
          f.press.key.to.continue()
        }
      }
      if (change[11]) {
        go.on <- 1
        count <- 0
        while (go.on > 0) {
          select.no <- menu(varnames[1:nvar], title = "\nQ11a: Give number of (another) factor for which you want to select data\n                       --- type 0 if you want to analyze all (remaining) data --- \n")
          if (select.no == 0) {
            go.on <- 0
            tmp.lev <- "none"
          }
          else {
            count <- count + 1
            select.fact <- data.0[, select.no]
            if (!is.numeric(select.fact)) {
              levels.all <- sort(levels(factor(select.fact)))
            }
            if (is.numeric(select.fact)) {
              levels.all <- sort(unique(select.fact))
            }
            nr.lev <- length(levels.all)
            cat("\nQ11b: Give number(s) associated with the level(s) you want to select\n")
            cat(" ")
            cat(paste(1:nr.lev, ":", levels.all, "\n"))
            tmp.no <- eval(parse(prompt = paste(" -------- > ")))
            tmp.lev <- levels.all[tmp.no]
            zz.1 <- data.0[, select.no]
            zz.2 <- zz.1[!is.na(zz.1)]
            if (quick.ans == 1 && length(zz.1) > length(zz.2)) 
              cat("\n\nATTENTION: there are missing values in this factor\n")
            check <- F
            for (ii in tmp.lev) check <- c(check, zz.2 == 
                                             ii)
            if (sum(check) == 0) {
              cat("\n\nfactor level(s)", tmp.lev, "not found, re-start f.proast with argument er=T\n\n")
              return(invisible())
            }
            ans.all$select.no <- select.no
            ans.all$select.lev <- tmp.lev
            ans.all$data.0 <- data.0
            data.0 <- f.subset(ans.all, select = T)
            ans.all$data.0 <- data.0
            name.tmp <- paste(varnames[select.no], paste(tmp.lev, 
                                                         collapse = " "))
            if (select.name == "all") 
              select.name <- character()
            select.name <- paste(select.name, "\n  ", 
                                 name.tmp)
          }
        }
        if (count == 0) {
          cat("\nATTENTION: complete dataset will be restored\n")
          data.0 <- odt
          data.0 <- f.remove.NAs(xans, yans, sans, dfr = data.0)
          if (length(unique(Vdetlim) == 1)) {
            Vdetlim <- rep(Vdetlim[1], length(data.0[, 
                                                     1]))
            data.0 <- cbind(data.0, Vdetlim)
          }
          ans.all$data.0 <- data.0
          select.name <- "all"
        }
        change[11] <- F
        ans.all$select.lev <- tmp.lev
        ans.all$select.name <- select.name
      }
      if (change[12]) {
        if (dtype == 2) 
          lump.ans <- menu(c("yes", "no"), title = "\nQ12: Do you want to consider observed fractions per dose group\n \n                  i.e. you have at least 5 observations per dose group?")
        else lump.ans <- 2
        change[12] <- F
        ans.all$lump.ans <- lump.ans
        ans.all$yans <- 0
      }
      if (change[13]) {
        fct1.no <- menu(varnames[1:nvar], title = "\nQ13: Give number of factor serving as covariate with respect to parameter a \n --- type 0 if none --- \n  ")
        change[13] <- F
        ans.all$fct1.no <- fct1.no
        ans.all$fct1.name <- varnames[fct1.no]
      }
      if (change[14]) {
        fct2.no <- menu(varnames[1:nvar], title = "\nQ14: Give number of factor serving as covariate with respect to parameter b \n --- type 0 if none --- \n  ")
        change[14] <- F
        ans.all$fct2.no <- fct2.no
        ans.all$fct2.name <- varnames[fct2.no]
      }
      if (change[15]) {
        if (dtype == 4) 
          fct3.no <- menu(varnames[1:nvar], title = "\nQ15: Give number of factor serving as covariate with respect to parameter theta \n --- type 0 if none --- \n")
        else fct3.no <- menu(varnames[1:nvar], title = "\nQ15: Give number of factor serving as covariate with respect to    \n                    parameter var \n --- type 0 if none --- \n")
        change[15] <- F
        ans.all$fct3.no <- fct3.no
        ans.all$fct3.name <- varnames[fct3.no]
      }
      if (change[16]) {
        if (model.ans %in% c(2, 4, 7, 9, 12, 14, 17, 
                             19, 22, 24)) {
          cat("\nATTENTION: parameter d does not occur in this model\n\n")
          f.press.key.to.continue()
        }
        fct5.no <- menu(varnames[1:nvar], title = "\nQ16: Give number of factor serving as covariate with respect to\n    parameter d \n --- type 0 if none --- \n  ")
        change[16] <- F
        ans.all$fct5.no <- fct5.no
        ans.all$fct5.name <- varnames[fct5.no]
      }
      if (change[17]) {
        if (model.ans %in% c(2, 3, 7, 8, 12, 13, 17, 
                             18, 22, 23)) {
          cat("\nATTENTION: parameter c does not occur in this model\n\n")
          f.press.key.to.continue()
        }
        fct4.no <- menu(varnames[1:nvar], title = "\nQ16: Give number of factor serving as covariate with respect to\n    parameter c \n --- type 0 if none --- \n  ")
        change[17] <- F
        ans.all$fct4.no <- fct4.no
        ans.all$fct4.name <- varnames[fct4.no]
      }
      if (change[18]) {
        if (is.na(x[1])) 
          ans.all <- f.execute(ans.all)
        if (dtype %in% c(1, 5, 10, 15, 25, 26)) {
          res.MS <- with(ans.all, res.MS <- f.con.tst(x, 
                                                      y, yy, dtype, fct1, fct2, fct3, ans.all$sf.x, 
                                                      x.leg, y.leg, info, factor.name, odt, Vdetlim = Vdetlim, 
                                                      select.name = select.name, version = PRversion, 
                                                      sd = sd, nn.gr = nn, ans.all = ans.all, color = color, 
                                                      nrp = 0), return(res.MS))
          change[18] <- F
          ans.all$res.MS <- res.MS
        }
        else if (dtype == 4) {
          f.cat.tst(ans.all)
          change[18] <- F
        }
        else cat("\n\nNOAEL approach not implemented for dtype", 
                 dtype, "\n\n")
      }
      if (change[19] && cont) 
        if (dtype %in% c(1, 5, 10, 15, 25, 26, 250, 260)) {
          transf.ans <- menu(c("no transformation", "sqrt-transformation", 
                               "log-transformation"), title = "\nQ19: Which transformation do you want?")
          ans.all$transf.ans <- transf.ans
          if (dtype %in% c(10, 15, 250, 260)) 
            switch(transf.ans, dtype <- 250, dtype <- 260, 
                   dtype <- 10)
          if (dtype %in% c(1, 25, 26)) 
            switch(transf.ans, dtype <- 25, dtype <- 26, 
                   dtype <- 1)
          if (dtype == 5) 
            switch(transf.ans, dtype <- 25, dtype <- 26, 
                   dtype <- 5)
          if (dtype %in% c(1, 10, 15)) 
            cat("\nData will be analysed on log-scale\n")
          if (dtype == 25 | dtype == 250) 
            cat("\nData will be analysed on original scale\n")
          if (dtype == 26) 
            cat("\nData will be analysed on square root-scale\n")
          ans.all$dtype <- dtype
          if (ans.all$yans > 0) 
            ans.all <- f.execute(ans.all)
          change[19] <- F
        }
      else cat("\nATTENTION : data cannot be transformed\n")
      if (change[19] && !cont) {
        ans.all$constr.ans <- menu(c("no", "yes (not recommended!)"), 
                                   title = "\nQ19: Do you want to constrain the models to have finite slope at zero? ")
        change[19] <- F
      }
      if (change[20]) {
        ans.all$aic.crit <- eval((parse(prompt = paste("Q20: Give value for \n      critical difference in AICs among models (default = 2) > "))))
        change[20] <- F
      }
      if (change[21]) {
        if (cont) {
          if (yans != 0) {
            y <- data.0[, yans]
            print(sort(y))
          }
          cens.up <- eval(parse(prompt = paste("\n\nQ21: Give value for right censoring > ")))
          change[21] <- F
          ans.all$cens.up <- cens.up
        }
        if (!cont) {
          if (decr.zz) {
            cat("DR is increasing, and will be changed to decreasing")
            ans.all$decr.zz <- F
            f.press.key.to.continue()
          }
          if (!decr.zz) {
            cat("DR is decreasing, and will be changed to increasing")
            ans.all$decr.zz <- T
            f.press.key.to.continue()
          }
          change[21] <- F
        }
      }
      if (change[22]) {
        ans.all$fit.ans <- menu(c("error in response only (default)", 
                                  "error in both response and dose (minimize sum of products)"), 
                                title = "\nQ22: What type of regression do you want?")
        change[22] <- F
      }
      if (change[23]) {
        ans.all$cond.ans <- menu(c("mild conditions", 
                                   "moderate conditions", "strict conditions (R default)"), 
                                 title = "\nQ23: How strict do you want the fitting conditions to be?")
        ans.all$control.lst <- f.control(ans.all$cond.ans)
        change[23] <- F
        change[nrQ] <- T
      }
      if (change[24]) {
        ans.all$tans <- menu(varnames[1:nvar], title = "\nQ24: Give column number of time variable ? \n")
        change[24] <- F
      }
      if (change[25]) {
        choose <- F
      }
      ans.all$show <- ""
      ans.all$change <- change
      ans.all$full6.done <- FALSE
      f.assign(".Pr.last", ans.all)
    }
    if (exists("track")) 
      print("f.change.settings:  END")
    return(ans.all)
  })
}

############
## f.cont ##
############
f.cont <- function (ans.all) 
{
  if (exists("track")) 
    print("f.cont")
  if (ans.all$dtype %in% c(1, 5, 15, 10, 25, 26, 250, 260)) 
    cont <- T
  if (ans.all$dtype %in% c(2, 3, 4, 6, 84)) 
    cont <- F
  return(cont)
}

##################
## f.remove.NAs ##
##################
f.remove.NAs <- function (xans = 0, tans = 0, yans = 0, sans = 0, nans = 0, covar.no = 0, 
                          fct1.no = 0, fct2.no = 0, fct3.no = 0, fct4.no = 0, fct5.no = 0, 
                          dfr) 
{
  if (exists("track")) 
    print("f.remove.NAs")
  lst.na.x <- rep(F, length(dfr[, 1]))
  lst.na.x1 <- rep(F, length(dfr[, 1]))
  lst.na.x2 <- rep(F, length(dfr[, 1]))
  lst.na.x3 <- rep(F, length(dfr[, 1]))
  lst.na.t <- rep(F, length(dfr[, 1]))
  lst.na.y <- rep(F, length(dfr[, 1]))
  lst.na.sd <- rep(F, length(dfr[, 1]))
  lst.na.nn <- rep(F, length(dfr[, 1]))
  lst.na.covar <- rep(F, length(dfr[, 1]))
  lst.na.fct1 <- rep(F, length(dfr[, 1]))
  lst.na.fct2 <- rep(F, length(dfr[, 1]))
  lst.na.fct3 <- rep(F, length(dfr[, 1]))
  lst.na.fct4 <- rep(F, length(dfr[, 1]))
  lst.na.fct5 <- rep(F, length(dfr[, 1]))
  if (xans[1] != 0) {
    if (length(xans) == 1) 
      lst.na.x <- is.na(dfr[, xans])
    else {
      lst.na.x1 <- is.na(dfr[, xans[1]])
      lst.na.x2 <- is.na(dfr[, xans[2]])
    }
    if (length(xans) == 3) 
      lst.na.x3 <- is.na(dfr[, xans[3]])
  }
  if (tans != 0) 
    lst.na.t <- is.na(dfr[, tans])
  if (yans != 0) 
    lst.na.y <- is.na(dfr[, yans])
  if (sans != 0) 
    lst.na.sd <- is.na(dfr[, sans])
  if (nans != 0) 
    lst.na.nn <- is.na(dfr[, nans])
  if (covar.no != 0) 
    lst.na.covar <- is.na(dfr[, covar.no])
  if (fct1.no != 0) 
    lst.na.fct1 <- is.na(dfr[, fct1.no])
  if (fct2.no != 0) 
    lst.na.fct2 <- is.na(dfr[, fct2.no])
  if (fct3.no != 0) 
    lst.na.fct3 <- is.na(dfr[, fct3.no])
  if (fct4.no != 0) 
    lst.na.fct4 <- is.na(dfr[, fct4.no])
  if (fct5.no != 0) 
    lst.na.fct5 <- is.na(dfr[, fct5.no])
  if (sum(lst.na.x) > 0) 
    cat("\n\nThere are ", sum(lst.na.x), "missing values in the x column\n")
  if (sum(lst.na.t) > 0) 
    cat("\n\nThere are ", sum(lst.na.t), "missing values in the time column\n")
  if (sum(lst.na.x1) > 0) 
    cat("\n\nThere are ", sum(lst.na.x1), "missing values in the x1 column\n")
  if (sum(lst.na.x2) > 0) 
    cat("\n\nThere are ", sum(lst.na.x2), "missing values in the x2 column\n")
  if (sum(lst.na.y) > 0) 
    cat("\n\nThere are ", sum(lst.na.y), "missing values in the y column\n")
  if (sum(lst.na.sd) > 0) 
    cat("\n\nThere are ", sum(lst.na.sd), "missing values in the sd column\n")
  if (sum(lst.na.nn) > 0) 
    cat("\n\nThere are ", sum(lst.na.nn), "missing values in the group size column\n")
  if (sum(lst.na.covar) > 0) 
    cat("\n\nThere are ", sum(lst.na.covar), "missing values in the covariate column\n")
  if (sum(lst.na.fct1) > 0) 
    cat("\n\nThere are ", sum(lst.na.fct1), "missing values in the covariate column on a\n")
  if (sum(lst.na.fct2) > 0) 
    cat("\n\nThere are ", sum(lst.na.fct2), "missing values in the covariate column in b\n")
  if (sum(lst.na.fct3) > 0) 
    cat("\n\nThere are ", sum(lst.na.fct3), "missing values in the covariate column on var\n")
  if (sum(lst.na.fct4) > 0) 
    cat("\n\nThere are ", sum(lst.na.fct4), "missing values in the covariate column on c\n")
  if (sum(lst.na.fct5) > 0) 
    cat("\n\nThere are ", sum(lst.na.fct5), "missing values in the covariate column on d\n")
  lst.na <- (lst.na.x | lst.na.x1 | lst.na.x2 | lst.na.t | 
               lst.na.y | lst.na.sd | lst.na.nn | lst.na.covar | lst.na.fct1 | 
               lst.na.fct2 | lst.na.fct3 | lst.na.fct4 | lst.na.fct5)
  if (sum(lst.na) > 0) 
    cat("\n these rows are removed from the dataset\n")
  dfr <- dfr[!lst.na, ]
  if (exists("track")) 
    print("f.remove.NAs:  END")
  return(dfr)
}

########################
## f.check.nonneg.num ##
########################
f.check.nonneg.num <- function (vec, gui = FALSE, dtype = NA, quick.ans = 1) 
{
  if (exists("track")) 
    print("f.check.nonneg.num")
  check.out <- F
  if (quick.ans == 1) {
    if (!is.numeric(vec)) {
      print(" values are:")
      print(vec)
      stop("\nATTENTION: not all values are numerical\n")
    }
    if (any(vec < 0)) {
      if (gui) 
        stop("\nvariable contain negative values\n", 
             call. = F)
      if (!gui) {
        cat("\n\nATTENTION: variable contains negative values\n")
        if (!(dtype %in% c(25, 250))) {
          ans.cont <- menu(c("no, this is an error, stop analysis", 
                             "yes, but omit log-transformation"), title = "\nDo you want to continue nonetheless?")
          switch(ans.cont, stop(message = "PROAST interrupted", 
                                call. = F), check.out <- T)
        }
      }
    }
    if (!any(vec > 0)) {
      if (gui) 
        stop("\nvariable contains zeros only\n", call. = F)
      if (!gui) {
        cat("\n\nATTENTION: variable contains zeros only\n")
        ans.cont <- 0
        switch(ans.cont, stop(message = "PROAST interrupted", 
                              call. = F), check.out <- T)
      }
    }
  }
  if (quick.ans > 1) {
    if (!is.numeric(vec) || any(vec < 0) || !any(vec > 0)) {
      cat("\nATTENTION: response data contain negative values, or only zeros,\n              or non-numerical values\nno analysis possible")
      print("the response data are:")
      print(vec)
      check.out <- T
    }
  }
  if (exists("track")) 
    print("f.check.nonneg.num:   END")
  invisible(check.out)
}

###############
## f.execute ##
###############
f.execute <- function (ans.all, no.plot = F) 
{
  if (exists("track")) 
    print("f.execute")
  ans.all$no.plot <- no.plot
  with(ans.all, {
    cont <- f.cont(ans.all)
    if (quick.ans == 1) 
      data.0 <- f.remove.NAs(xans, tans, yans, sans, nans, 
                             covar.no, fct1.no, fct2.no, fct3.no, fct4.no, 
                             fct5.no, dfr = data.0)
    if (dtype %in% c(4, 6, 84, 10, 15, 250, 260)) 
      data.0 <- f.remove.zeros(dtype, sans = sans, nans = nans, 
                               dfr = data.0)
    if (length(xans) > 1) {
      x.leg <- "dose after addition"
      if (sf.x != 1) 
        x.leg <- paste("dose/", sf.x, " after addition", 
                       sep = "")
      x1 <- as.numeric(data.0[, xans[1]])/sf.x
      x2 <- as.numeric(data.0[, xans[2]])/sf.x
      x <- x1 + x2
      if (ans.all$aic.crit == 2) 
        ans.all$aic.crit <- 5
      if (length(xans) == 3) {
        x3 <- as.numeric(data.0[, xans[3]])/sf.x
        x <- x1 + x2 + x3
      }
    }
    else {
      x <- as.numeric(data.0[, xans])/sf.x
      f.check.nonneg.num(x, gui)
      x.leg <- varnames[xans]
      if (sf.x != 1) 
        x.leg <- paste(x.leg, "/", sf.x, sep = "")
    }
    if (length(xans) == 3) {
      data.0 <- data.0[order(x1), ]
      x1 <- as.numeric(data.0[, xans[1]])/sf.x
      x2 <- as.numeric(data.0[, xans[2]])/sf.x
      x3 <- as.numeric(data.0[, xans[3]])/sf.x
      x <- x1 + x2 + x3
    }
    else if (length(xans) == 2) {
      data.0 <- data.0[order(x1), ]
      x1 <- as.numeric(data.0[, xans[1]])/sf.x
      x2 <- as.numeric(data.0[, xans[2]])/sf.x
      x <- x1 + x2
    }
    else {
      data.0 <- data.0[order(x), ]
      x <- as.numeric(data.0[, xans])/sf.x
    }
    if (tans > 0) {
      ttt <- as.numeric(data.0[, tans])/sf.t
      f.check.nonneg.num(ttt, gui)
      t.leg <- varnames[tans]
      ttt.txt <- levels(factor(ttt))
    }
    if (cont) {
      if (yans != 0) {
        y <- as.numeric(data.0[, yans])
        yy <- y
        y.leg <- varnames[yans]
        if (dtype %in% c(1, 5, 10, 15)) {
          data.problem <- f.check.nonneg.num(y, gui, 
                                             dtype)
          if (data.problem) {
            cat("\nATTENTION: analyses will be done without log-transformation, due to negative observations\n\n")
            if (dtype == 1) 
              dtype <- 25
            if (dtype %in% c(10)) 
              dtype <- 250
          }
        }
      }
      if (auto.detlim) {
        if (min(y) == 0) 
          detlim <- min(y[y > 0])
        detlim.col <- nvar + 1
        data.0[, detlim.col] <- detlim
        Vdetlim <- data.0[, detlim.col]
      }
      else if (detlim.col > 0) 
        Vdetlim <- data.0[, detlim.col]
      else if (detlim > 0) {
        detlim.col <- nvar + 1
        data.0[, detlim.col] <- detlim
        Vdetlim <- data.0[, detlim.col]
      }
      if (dtype %in% c(1, 5) && min(y) == 0 && min(Vdetlim) == 
          0) {
        #cat("\n\nATTENTION:  minimum value of observed response data is zero\n")
        if (dtype %in% c(10, 15, 250, 260)) 
          cat("(note that detection limit does not really apply to summary data)\n")
        ans.detlim <- 1
        #ans.detlim <- menu(c("Give detection limit", 
        #                     "Add nonzero value to observations"), title = "Q6-7: What do you want to do?")
        if (ans.detlim == 1) {
          #cat("\nthe lowest nonzero value is: ", min(y[y > 
          #                                               0]), "\n")
          detlim <- min(y[y > 0])
          #detlim <- eval(parse(prompt = paste("\n\nQ6: Give global value of detection limit,\n                       or type 0 if you have individual values > ")))
          if (detlim == 0) 
            Vdetlim <- menu(varnames[1:nvar], title = "\nQ6a: Give column with individual detection limits\n")
          if (detlim != 0 && detlim.col > 0) {
            Vdetlim <- data.0[, detlim.col]
            if (min(Vdetlim <= 0)) {
              cat("\n ATTENTION: there are non-positive detection limits !\n")
              f.press.key.to.continue()
            }
            if (any(is.na(Vdetlim))) {
              cat("\nATTENTION:  column with detection limits contains NAs; \n                       replace NAs with (large) numbers\n\n")
              f.press.key.to.continue()
            }
          }
          if (detlim.col == 0) {
            Vdetlim <- rep(detlim, length(y))
            data.0 <- cbind(data.0, Vdetlim)
            detlim.col <- nvar + 1
            ans.all$nvar <- nvar + 1
          }
          if (ans.detlim == 2) {
            cat("\nThe lowest nonzero observation is: ", 
                min(y[y > 0]), "\n")
            y.add <- eval(parse(prompt = "\n\nQ7: Give small value to be added to the data > "))
          }
        }
      }
      if (min(Vdetlim) > 0) {
        y[y < Vdetlim] <- -1000
        yy <- yy * (y != -1000) + (0.5 * Vdetlim) * (y == 
                                                       -1000)
      }
      if (!is.na(cens.up)) {
        y[y >= cens.up] <- -2000
        yy <- yy * (y != -2000) + (1 * cens.up) * (y == 
                                                     -2000)
      }
      if (max(Vdetlim) > 0) 
        ans.all$low.y <- min(Vdetlim)/2
      else ans.all$low.y <- 0.98 * min(yy)
      ans.all$upp.y <- 1.02 * max(yy)
      if (ans.detlim == 2) {
        y <- y + y.add
        yy <- y
      }
      if (dtype %in% c(1, 5)) {
        y.denom <- 1
        if (yans.denom != 0) {
          y.denom <- data.0[, yans.denom]
          if (min(y.denom) == 0) 
            cat("\nATTENTION: zero values in denominator ! \n\n")
          y <- y/y.denom
          yy <- yy/y.denom
          y.leg <- paste(y.leg, "/", varnames[yans.denom])
        }
      }
      nn <- 0
      if (dtype %in% c(10, 15, 250, 260)) {
        if (sans == 0) 
          stop(message = "no data for variation statistics provided\n")
        if (nans == 0) 
          stop(message = "no data for group sizes provided\n")
        sd <- as.numeric(data.0[, sans])
        nn <- as.numeric(data.0[, nans])
        if (any(y == 0)) {
          cat("\nATTENTION: there are zero (mean) observations in the data\n")
          cat("\nThis is not allowed for mean values\n")
          cat("the analysis will be interrupted\n")
          stop()
        }
        if (any(sd == 0)) {
          cat("\nATTENTION: there are zero SDs in the data\n")
          cat("\nThis is not allowed for SDs\n")
          cat("the analysis will be interrupted\n")
          stop()
        }
        if (sd.se == 2) 
          sd <- sd * sqrt(nn)
        if (0) {
          f.graph.window(2)
          plot(log10(y), log10(sd), xlab = "mean (log10-scale)", 
               ylab = "sd (log10-scale)")
          title(main = "correlation between mean of y and sd of y")
          y.mn <- y
          CV <- sd/y.mn
          mn.log <- logb(y.mn/sqrt(1 + CV^2))
          mn.log <- mn.log/logb(10)
          sd2.log <- logb(CV^2 + 1)
          sd2.log <- sd2.log/logb(10)^2
          plot(mn.log, sd2.log, xlab = "mean of log-y", 
               ylab = "sd of log-y")
          title(main = "correlation between mean of log-y and sd of log-y")
          eval(parse(prompt = "\n\nCheck plot and type c to continue .................. >   "))
        }
        if (dtype %in% c(10, 15)) {
          y.mn <- y
          CV <- sd/y.mn
          mn.log <- logb(y.mn/sqrt(1 + CV^2))
          yy <- exp(mn.log)
          ans.all$mn.log <- mn.log
          ans.all$sd2.log <- logb(CV^2 + 1)
        }
        if (dtype == 250) {
          y.mn <- y
          ans.all$mn.log <- y.mn
          ans.all$sd2.log <- sd^2
          ans.all$yy <- y.mn
        }
      }
    }
    if (dtype %in% c(4, 6, 84)) {
      nn <- data.0[, nans]
      if (yans != 0) {
        kk <- data.0[, yans]
        nn <- data.0[, nans]
      }
      if (sum(kk > nn) > 0) {
        cat("\nNumber of responses is larger than sample size ! \n\n")
        return(ans.all)
      }
      y <- kk/nn
      yy <- y
      if (max(y, na.rm = T) < 0.5) 
        ans.all$y.lim.plt <- c(0, 1.2 * max(y))
      else ans.all$y.lim.plt <- c(0, 1)
    }
    if (!cont) {
      cens <- 0
      if (sum(varnames == "cens") > 0) {
        zzz <- 1:nvar
        cens.col <- zzz[varnames == "cens"]
        cens <- data.0[, cens.col]
      }
    }
    if (dtype %in% 2:3) {
      if (yans != 0) {
        y <- as.numeric(data.0[, yans])
        yy <- y
        y.leg <- varnames[yans]
      }
      nn <- 1
      y.original <- y
      scores.orig <- sort(unique(y))
      ymn <- tapply(y, x, mean)
      xmn <- tapply(x, x, mean)
      ymn <- as.numeric(ymn)
      if (reverse.ans == 0) {
        if (length(xmn) > 1) 
          if (var(xmn != 0) && var(ymn != 0)) 
            if (cor(xmn, ymn) < 0) {
              cat("\nProast assumes zero to represent normal, and higher scores abnormal")
              cat("\nYour data seem to have the opposite direction\n")
              ans.all$reverse.ans <- menu(c("yes", "no"), 
                                          title = "Do your scores have opposite direction?\n")
            }
      }
      if (ans.all$reverse.ans == 1) 
        scores.orig <- rev(scores.orig)
      score <- 0
      n.lev <- length(scores.orig)
      for (ii in 1:n.lev) {
        y[y.original == scores.orig[ii]] <- score
        score <- score + 1
      }
      scores.mtr <- cbind(scores.orig, levels(factor(y)))
      dimnames(scores.mtr) <- list(NULL, c("orig.scores", 
                                           "temp.scores"))
      dum.ord <- sum(scores.mtr[, 1] == scores.mtr[, 2]) != 
        length(scores.mtr[, 1])
      if (dum.ord) {
        cat("\n\nATTENTION: the original scores have been transformed for analysis as follows:\n")
        print(scores.mtr)
      }
      ans.all$y.transf <- y
    }
    fct1.txt <- ""
    fct2.txt <- ""
    fct3.txt <- ""
    fct4.txt <- ""
    fct5.txt <- ""
    fct1 <- rep(1, length(y))
    fct2 <- rep(1, length(y))
    if (cont) 
      fct3 <- rep(1, length(y))
    if (!cont) 
      fct3 <- 1
    fct4 <- rep(1, length(y))
    fct5 <- rep(1, length(y))
    if (quick.ans == 1) {
      if (fct1.no > 0) {
        fct1 <- as.numeric(factor(data.0[, fct1.no]))
        fct1.txt <- levels(factor(data.0[, fct1.no]))
        if (interrupt && max(fct1) == 1) {
          cat("\nAttention 1: The factor you chose as covariate on par a has only one level\n")
          cat("you might have selected a subgroup for this factor\n")
          f.press.key.to.continue()
        }
      }
      if (fct2.no > 0) {
        fct2 <- as.numeric(factor(data.0[, fct2.no]))
        fct2.txt <- levels(factor(data.0[, fct2.no]))
        if (max(fct2) == 1) {
          cat("\nAttention: the factor you chose as covariate on par b has only one level\n")
          cat("you might have selected a subgroup for this factor\n")
          f.press.key.to.continue()
        }
      }
      if (fct3.no > 0) {
        fct3 <- data.0[, fct3.no]
        fct3.txt <- levels(factor(data.0[, fct3.no]))
        if (cont) {
          fct3 <- as.numeric(factor(data.0[, fct3.no]))
          if (max(fct3) == 1) {
            cat("\nAttention: the factor you chose as covariate on var has only one level\n")
            cat("you might have selected a subgroup for this factor\n")
            f.press.key.to.continue()
          }
          ans.all$nr.var <- max(fct3)
        }
        if (!cont) {
          print(levels(factor(data.0[, fct3.no])))
          ans.all$fct3.ref <- min(data.0[, fct3.no])
          ans.all$nr.var <- 0
        }
      }
      if (fct4.no > 0) {
        fct4 <- as.numeric(factor(data.0[, fct4.no]))
        fct4.txt <- levels(factor(data.0[, fct4.no]))
        if (max(fct4) == 1) {
          cat("\nAttention: the factor you chose as covariate on par c has only one level\n")
          cat("you might have selected a subgroup for this factor\n")
          f.press.key.to.continue()
        }
      }
      if (fct5.no > 0) {
        fct5 <- as.numeric(factor(data.0[, fct5.no]))
        fct5.txt <- levels(factor(data.0[, fct5.no]))
        if (max(fct5) == 1) {
          cat("\nAttention: the factor you chose as covariate on par d has only one level\n")
          cat("you might have selected a subgroup for this factor\n")
          f.press.key.to.continue()
        }
      }
      twice <- F
      if (fct1.no != 0) 
        twice <- fct1.no == fct2.no
      if (!twice && fct1.no > 1) 
        twice <- fct1.no == fct5.no
      if (!twice && fct1.no > 1) 
        twice <- fct1.no == fct4.no
    }
    if (quick.ans > 1) {
      twice <- T
      covariate <- rep(1, length(y))
      covar.txt <- ""
      if (covar.no > 0) {
        covariate <- as.numeric(factor(data.0[, covar.no]))
        if (max(covariate, na.rm = T) == 1) {
          cat("\nAttention: the factor you chose as covariate has only one level\n")
          ans.all$covar.no <- 0
        }
        covar.txt <- levels(factor(data.0[, covar.no]))
        covar.name <- varnames[covar.no]
        if (sum(is.na(covariate)) > 0) {
          cat("\nAttention: the factor you chose as covariate contains NAs\n")
          cat("\nThe associated observations will be removed\n")
          f.press.key.to.continue()
          covariate <- data.0[, covar.no]
        }
        fct1.txt <- covar.txt
        fct2.txt <- covar.txt
        if (cont) {
          fct3.txt <- covar.txt
        }
      }
      ans.all$covariate <- covariate
      ans.all$covar.txt <- covar.txt
      ans.all$covar.name <- covar.name
    }
    ans.all$cont <- cont
    if (quick.ans == 1) 
      ans.all$gr.txt <- f.pars(ans.all)$gr.txt
    if (quick.ans > 1) 
      ans.all$gr.txt <- covar.txt
    if (displ.no > 0) {
      ans.all$displ.fact <- as.numeric(factor(data.0[, 
                                                     displ.no]))
      ans.all$displ.txt <- levels(factor(data.0[, displ.no]))
    }
    if (displ.no == 0) {
      ans.all$displ.fact <- 1
      ans.all$displ.txt <- NA
    }
    ans.all$displ.no <- displ.no
    ans.all$nr.aa <- max(fct1)
    ans.all$nr.bb <- max(fct2)
    ans.all$nr.var <- max(fct3)
    ans.all$nr.cc <- max(fct4)
    ans.all$nr.dd <- max(fct5)
    ans.all$fct1 <- fct1
    ans.all$fct2 <- fct2
    ans.all$fct3 <- fct3
    ans.all$fct4 <- fct4
    ans.all$fct5 <- fct5
    ans.all$fct1.txt <- fct1.txt
    ans.all$fct2.txt <- fct2.txt
    if (length(xans) > 1) 
      ans.all$fct2.txt <- substring(varnames[xans[1:length(xans)]], 
                                    1, 2)
    ans.all$fct3.txt <- fct3.txt
    ans.all$fct4.txt <- fct4.txt
    ans.all$fct5.txt <- fct5.txt
    if (quick.ans == 1) 
      trace.plt <- T
    else trace.plt <- F
    if (quick.ans == 1) 
      trace <- T
    else trace <- F
    if (min(x) == 0) 
      dum.contr <- min(x[x > 0])/4
    else dum.contr <- min(x)
    if (min(y) == 0) 
      dum.zero.resp <- min(yy[yy > 0])/2
    else dum.zero.resp <- min(yy)
    xy.lim[1] <- dum.contr
    xy.lim[2] <- min(x)
    xy.lim[3] <- max(x)
    xy.lim[4] <- min(yy)
    xy.lim[5] <- max(yy)
    if (!cont) {
      xy.lim[4] <- 0
      xy.lim[5] <- 1
      if (max(y) < 0.5) 
        xy.lim[5] <- (max(y) + 0.1) * 1.2
    }
    xy.lim[6] <- dum.zero.resp
    if (0) {
      if (sum(yy[yy > 0]) == 0) 
        low.ylog <- -log10(0.3)
      else {
        if (min(yy) == 0) 
          low.ylog <- log10(min(yy[yy > 0])/2/max(yy))
        else low.ylog <- log10(min(yy))
      }
    }
    if (cont) {
      if (plt.mns == 2) 
        cex.1 <- 1
      else cex.1 <- 0.75
      if (plot.type == 9) 
        cex.1 <- 0.6
      if (dtype %in% c(10, 15, 110, 250, 260)) 
        cex.1 <- 1.6
      cex.2 <- 1.5
      if (dtype == 15) 
        cex.2 <- 0.75
    }
    if (!cont) {
      if (dtype == 6) 
        cex.1 <- 0.6
      else cex.1 <- 1.2
      cex.2 <- 1.5
    }
    ans.all$cex.1 <- cex.1
    ans.all$cex.2 <- cex.2
    ans.all$quick.ans <- quick.ans
    if (quick.ans > 1) 
      ans.all$quick <- T
    ans.all$data.0 <- data.0
    ans.all$x <- x
    if (length(xans) > 1) {
      ans.all$x1 <- x1
      ans.all$x2 <- x2
    }
    if (length(xans) > 2) 
      ans.all$x3 <- x3
    ans.all$y <- y
    ans.all$yy <- yy
    ans.all$dtype <- dtype
    ans.all$cont <- cont
    if (cont) 
      ans.all$model.names <- f.expect.con(name = T)
    ans.all$x.leg <- x.leg
    if (yans > 0) 
      ans.all$y.leg <- y.leg
    ans.all$sd <- sd
    ans.all$nn <- nn
    ans.all$detlim <- detlim
    ans.all$Vdetlim <- Vdetlim
    ans.all$detlim.col <- detlim.col
    ans.all$ans.detlim <- ans.detlim
    ans.all$y.add <- y.add
    if (dtype %in% c(4, 6, 84) && yans != 0) {
      ans.all$kk <- kk
      ans.all$y <- kk/nn
      ans.all$cens <- cens
    }
    ans.all$ttt <- ttt
    ans.all$ttt.txt <- ttt.txt
    ans.all$t.leg <- t.leg
    t
    ans.all$xy.lim <- xy.lim
    ans.all$twice <- twice
    ans.all$factor.name <- ""
    ans.all$factor.name <- f.factorname(ans.all)
    ans.all$low.y <- low.y
    ans.all$upp.y <- upp.y
    if (dtype == 11) 
      ans.all$fit.ans <- 2
    if (ans.all$fit.ans == 2) {
      if (!cont) 
        print("fit.ans has value 2 !!")
      ans.all$plot.type <- 4
    }
    ans.all$trace <- trace
    ans.all$trace.plt <- trace.plt
    ans.all$created <- date()
    if (dtype == 3) {
      ans.all$scores.orig <- scores.orig
      ans.all$scores.mtr <- scores.mtr
      ans.all$nth <- max(y)
      if (ans.all$cat.ans > 0) {
        y[y < cat.ans] <- 0
        y[y >= cat.ans] <- 1
        ans.all$dtype <- 2
        ans.all$y <- y
        ans.all$nth <- 1
        ans.all$plot.type <- 2
      }
    }
    if (dtype %in% c(2, 3, 6, 15)) 
      ans.all$CI.plt <- FALSE
    if (dtype %in% c(4, 10, 250)) 
      ans.all$CI.plt <- TRUE
    if (!no.plot) {
      f.graph.window(1)
      if (!fitted) 
        ans.all$heading <- "data"
      if (cont) 
        f.plot.con(ans.all)
      if (dtype %in% c(4, 6, 84)) {
        if (max(fct2) == 1) 
          kk.dum <- kk
        else kk.dum <- 1
        if (dtype == 6) 
          kk.dum <- 1
        if (length(y) > 200) 
          kk.dum <- 1
        if (ans.all$plot.type == 0) 
          ans.all$plot.type <- 1
        f.plot.frq(ans.all)
      }
      if (dtype %in% 2:3) 
        f.plot.cat(plot.type = 1, x, y, x.leg, y.leg, 
                   fct1 = fct1, fct2 = fct2, shift = 0, heading = "", 
                   factor.name = ans.all$factor.name, color = color)
    }
    if (dtype != 6 && !do.MA) {
      ans.all$pi.full <- NA
      ans.all$x.full <- NA
      ans.all$fct1.full <- NA
      ans.all$fct2.full <- NA
      ans.all$loglik.full <- NA
      ans.all$npar.full <- NA
    }
    if (dtype == 6) {
      ans.all$alfa.length <- 1
    }
    if (0) 
      if (!gui) {
        dtype.tmp <- odt$dtype[ans.all$yans]
        if (dtype.tmp != ans.all$dtype) {
          cat("\nThe type of data you indicated does not match the one indicated in the dataset\n")
          cat("indicated in dataset:", dtype.tmp)
          cat("\nintended:", ans.all$dtype, "\n")
          f.press.key.to.continue()
        }
      }
    if (cont && length(xans) == 1) 
      ans.all$x.mn <- mean(x)
    if (cont) {
      dum.nn <- f.nr.replicates(ans.all)
      if (mean(dum.nn) < 3) 
        ans.all$plt.mns <- 2
      if (length(xans) > 1) {
        ans.all$plt.mns <- 3
        ans.all$CI.plt <- TRUE
      }
    }
    ans.all$dtype.0 <- dtype
    if (ans.all$dtype == 3) 
      ans.all$plot.type <- 6
    if (lump.ans == 1) 
      ans.all <- f.lump.bin(ans.all)
    if (quick.ans == 1) 
      ans.all$model.fam <- 1
    else ans.all$model.fam <- 0
    f.assign(".Pr.last", ans.all)
    if (exists("track")) 
      print("f.execute:  END")
    return(ans.all)
  })
}

##################
## f.expect.con ##
##################
f.expect.con <- function (model.ans, x, regr.par = 0, fct1 = 1, fct2 = 1, fct3 = 1, 
                          fct4 = 1, fct5 = 1, name = F, CES = NA, twice = T, ttt = 0, 
                          yy = 0, trace.expect = F, increase, x.mn = NA, ref.lev, ans.m6.sd = 1, 
                          par.start = NA, sign.q = 0) 
{
  if (exists("track2")) 
    print("f.expect.con")
  if (!name) {
    if (model.ans == -6) {
    }
  }
  if (name) {
    model.names <- c("E1: y = a", "E2: y = a*exp(bx)", "E3: y = a*exp(bx^d)", 
                     "E4: y = a*[c-(c-1)exp(-bx)]", "E5: y = a*[c-(c-1)exp(-bx^d)]", 
                     "E5b: y = a*(c^(1 - exp(-(x/b)^d)))", "E2 y = a*exp(x/b)", 
                     "E3 y = a*exp((x/b)^d)", "E4 y = a * [c - (c-1)exp(-(x/b))]", 
                     "E5 y = a * [c - (c-1)exp(-(x/b)^d)]", "Full model: y = group mean", 
                     "E2-CED: y = a*exp(bx)", "E3-CED: y = a*exp(bx^d)", 
                     "E4-CED: y = a * [c-(c-1)exp(-bx)]", "E5-CED: y = a * [c-(c-1)exp(-bx^d)]", 
                     "BMDratio model (M10)", "H2: a * (1 - x/(b+x))", 
                     "H3: a * (1 - x^d/(b^d+x^d))", "H4: a * (1 + (c-1)x/(b+x))", 
                     "H5: a * (1 + (c-1)x^d/(b^d+x^d))", "y = a {[c1-(c1-1)exp(-(x/b1)^d)] * [c2 - (c2-1)exp(-(x/b2)^d)]}", 
                     "H2-CED: a * (1 - x/(b+x))", "H3-CED: a * (1 - x^d/(b^d+x^d))", 
                     "H4-CED: a * (1 + (c-1)x/(b+x))", "H5-CED: a * (1 + (c-1)x^d/(b^d+x^d))", 
                     "y = c + a*x^b", "y = c + exp((x/b)^a)", "y = c + b(x-a)^d*(x>a)", 
                     "y = a + c * (1 - exp( - (x/b)^d))", "y = a + b*x^c", 
                     "y = c + a*exp((x/b)^d)", "y = a {[cc-(cc-1)exp(-(x/b1)^d1)] * [cc - (cc-1)exp(-(x/b2)^d2)]}", 
                     "y = a {[c1-(c1-1)exp(-(x/b1)^d1)] * [c2 - (c2-1)exp(-(x/b2)^d2)]}", 
                     "y = a + bx + cx^2 + dx^3", "y = a/b [1 - exp(-bx)]", 
                     "y = a * (1 + (c-1)(x^d/(b^d+x^d)))", "y =  b*exp(c(x-a)^d) * (x>a) + b * (x<a)", 
                     "select model 5 from nested family", "select model 5 or model 3 from nested family", 
                     "select model 3 from nested family", "y = ax^b", 
                     "", "hill-additive: y = a + cx^d/(b^d+x^d)", "E2 with independent a and b: y = exp(a + b(x - mean.x))", 
                     "y = a {[c1-(c1-1)exp(-(x/b)^d)] * [c2 - (c2-1)exp(-(x/b)^d)]}", 
                     "E15 in terms of RPFs", "E6 in terms of CED", "CxT model E5: y = a * [c - (c-1)exp[-(x/b)^d - (t/g)^h)]", 
                     "E5b-ED50: y = a*(c^(1 - exp(-ln(2)(x/ED50)^d)))", 
                     "CxT model double exponential: y = a * [c - (c-1)exp[-(x/b.time)^d]] with b.time = b * exp( - (t/g)^h)")
    return(model.names)
  }
  nr.aa <- max(fct1)
  nr.bb <- max(fct2)
  nr.var <- max(fct3)
  nr.cc <- max(fct4)
  nr.dd <- max(fct5)
  if (model.ans == 47 || (model.ans == 6 && ans.m6.sd == 2)) {
    sd.tmp <- sqrt(regr.par[1:nr.var])
    sd0 <- rep(0, length(x))
    for (ii in (1:nr.var)) sd0 <- sd0 + sd.tmp[ii] * (fct3 == 
                                                        ii)
    regr.par <- regr.par[-(1:nr.var)]
    qq <- regr.par[nr.aa + nr.bb + 1]
    sign.q0 <- rep(0, length(x))
    for (ii in (1:nr.var)) sign.q0 <- sign.q0 + sign.q[ii] * 
      (fct3 == ii)
  }
  if (model.ans != 11) {
    nrp <- length(regr.par)
    aa0 <- rep(0, length(x))
    aa.tmp <- regr.par[1:nr.aa]
    if (nr.aa == 1) 
      aa0 <- aa0 + aa.tmp[1]
    else for (ii in (1:nr.aa)) aa0 <- aa0 + aa.tmp[ii] * 
      (fct1 == ii)
    bb0 <- rep(0, length(x))
    bb.tmp <- regr.par[(nr.aa + 1):(nr.aa + nr.bb)]
    if (nr.bb == 1) 
      bb0 <- bb0 + bb.tmp[1]
    else for (jj in (1:nr.bb)) bb0 <- bb0 + bb.tmp[jj] * 
      (fct2 == jj)
    par3 <- regr.par[nr.aa + nr.bb + 1]
    if (length(par3) == 0 || is.na(par3)) 
      par3 <- 0
    par4 <- regr.par[nr.aa + nr.bb + nr.cc + 1]
    cc0 <- par3
    dd0 <- par4
    if (model.ans %in% c(4, 5, 6, 9, 10, 14, 15, 19, 20, 
                         24, 25, 41, 42, 46)) {
      cc0 <- par3
      if (max(fct4) > 1) {
        cc0 <- rep(0, length(x))
        cc.tmp <- regr.par[(nr.aa + nr.bb + 1):(nr.aa + 
                                                  nr.bb + nr.cc)]
        for (kk in (1:nr.cc)) cc0 <- cc0 + cc.tmp[kk] * 
          (fct4 == kk)
      }
    }
    if (model.ans %in% c(3, 8, 13, 18, 23)) {
      dd0 <- par3
      if (max(fct5) > 1) {
        dd0 <- rep(0, length(x))
        dd.tmp <- regr.par[(nr.aa + nr.bb + 1):length(regr.par)]
        for (kk in (1:nr.dd)) dd0 <- dd0 + dd.tmp[kk] * 
          (fct5 == kk)
      }
    }
  }
  if (model.ans %in% c(5, 6, 10, 15, 20, 25, 41, 42, 46)) {
    dd0 <- par4
    if (max(fct5) > 1) {
      dd0 <- rep(0, length(x))
      dd.tmp <- regr.par[(nr.aa + nr.bb + nr.cc + 1):length(regr.par)]
      for (kk in (1:nr.dd)) dd0 <- dd0 + dd.tmp[kk] * (fct5 == 
                                                         kk)
    }
  }
  if (model.ans %in% c(48, 50)) {
    gg <- regr.par[nr.aa + nr.bb + nr.cc + nr.dd + 1]
    hh <- regr.par[nr.aa + nr.bb + nr.cc + nr.dd + 2]
  }
  switch(model.ans, y.expect <- aa0, y.expect <- aa0 * exp(bb0 * 
                                                             x + cc0 * ttt), y.expect <- aa0 * exp(bb0 * (x^dd0)), 
         y.expect <- aa0 * (cc0 - (cc0 - 1) * exp(-bb0 * x)), 
         {
           y.expect <- aa0 * (cc0 - (cc0 - 1) * exp(-bb0 * (x^dd0)))
         }, {
           if (ans.m6.sd == 1) y.expect <- exp(log(aa0) + log(cc0) * 
                                                 (1 - exp(-(x/bb0)^dd0)))
           if (ans.m6.sd == 2) y.expect <- exp(log(aa0) + sign.q0 * 
                                                 qq * sd0 * (1 - exp(-(x/bb0)^dd0)))
         }, y.expect <- aa0 * exp(x/bb0), {
           if (bb0[1] > 0) y.expect <- aa0 * exp((x/bb0)^dd0) else y.expect <- aa0 * 
               exp(-(x/-bb0)^dd0)
         }, y.expect <- aa0 * (cc0 - (cc0 - 1) * exp(-x/bb0)), 
         y.expect <- aa0 * (cc0 - (cc0 - 1) * exp(-(x/bb0)^dd0)), 
         {
           x.gr <- levels(factor(x))
           x.fact <- factor(x)
           y.tmp <- rep(NA, length(x))
           y.expect <- rep(0, length(x))
           if (nr.var > 1 & nr.aa == 1 & nr.bb == 1) {
             for (jj in (1:nr.var)) for (ii in (1:length(x.gr))) {
               y.tmp <- yy[x.fact == x.gr[ii] & fct3 == jj]
               if (length(y.tmp) > 0) {
                 y.mn <- exp(mean(logb(y.tmp)))
                 y.expect <- y.expect + y.mn * (x.fact == 
                                                  x.gr[ii]) * (fct3 == jj)
               }
               y.mn <- y.tmp
             }
           } else if (twice) for (jj in (1:nr.bb)) for (ii in (1:length(x.gr))) {
             y.tmp <- yy[x.fact == x.gr[ii] & fct2 == jj]
             if (length(y.tmp) > 0) {
               y.mn <- exp(mean(logb(y.tmp)))
               y.expect <- y.expect + y.mn * (x.fact == x.gr[ii]) * 
                 (fct2 == jj)
             }
           } else if (!twice) for (jj in (1:nr.aa)) for (kk in (1:nr.bb)) for (ii in (1:length(x.gr))) {
             y.tmp <- yy[x.fact == x.gr[ii] & fct1 == jj & 
                           fct2 == kk]
             if (length(y.tmp) > 0) {
               y.mn <- exp(mean(logb(y.tmp)))
               y.expect <- y.expect + y.mn * (x.fact == x.gr[ii]) * 
                 (fct1 == jj) * (fct2 == kk)
             }
           }
         }, {
           y.expect <- aa0 * (CES + 1)^(x/bb0)
         }, y.expect <- aa0 * (CES + 1)^((x/bb0)^dd0), y.expect <- aa0 * 
           (cc0 - (cc0 - 1) * ((CES + 1 - cc0)/(1 - cc0))^(x/bb0)), 
         {
           y.expect <- aa0 * (cc0 - (cc0 - 1) * ((CES + 1 - 
                                                    cc0)/(1 - cc0))^((x/bb0)^(dd0)))
         }, {
           if (!exists("CES.16")) {
             cat("\nATTENTION:  CES.16 not defined; CES2 is set to 0.10 and CES1 to 0.05\n")
             CES2 <- 0.1
             CES1 <- 0.05
             f.press.key.to.continue()
           } else {
             CES2 <- CES.16$CES2
             CES1 <- CES.16$CES1
           }
           if (par4 <= 1) print(c("f.expect.con, value of BMD-ratio:", 
                                  par4))
           dd <- f.uniroot.BMDratio(ratio = par4, cc = cc0, 
                                    CES1 = CES1, CES2 = CES2)
           y.expect <- aa0 * (cc0 - (cc0 - 1) * exp(-(x/bb0)^dd))
           if (is.na(dd)) y.expect <- rep(0, length(x))
         }, {
           y.expect <- aa0 * (1 - x/(bb0 + x))
         }, {
           if (increase == 1) y.expect <- aa0 * (1 - x^dd0/(-(-bb0)^dd0 + 
                                                              x^dd0))
           if (increase == -1) y.expect <- aa0 * (1 - x^dd0/((bb0)^dd0 + 
                                                               x^dd0))
         }, {
           y.expect <- aa0 * (1 + ((cc0 - 1) * x)/(bb0 + x))
         }, {
           y.expect <- aa0 * (1 + ((cc0 - 1) * x^dd0)/(sign(bb0) * 
                                                         (abs(bb0)^dd0 + x^dd0)))
         }, {
           b2 <- regr.par[nr.aa + nr.bb + 1]
           c1 <- regr.par[nr.aa + nr.bb + 2]
           c2 <- regr.par[nr.aa + nr.bb + 3]
           dd <- regr.par[nr.aa + nr.bb + 4]
           y.expect.1 <- c1 - (c1 - 1) * exp(-(x/bb0)^dd)
           y.expect.2 <- c2 - (c2 - 1) * exp(-(x/b2)^dd)
           y.expect <- aa0 * (y.expect.1 * y.expect.2)
         }, {
           dum <- f.bb.con(model.ans, cc = NA, dd = NA, CED = bb0, 
                           CES)
           y.expect <- aa0 * (1 - x/(sign(dum) * (abs(dum) + 
                                                    x)))
           y.expect <- aa0 * (1 - x/(dum + x))
         }, {
           dum <- f.bb.con(model.ans, cc = NA, dd = dd0, CED = bb0, 
                           CES)
           if (increase == 1) y.expect <- aa0 * (1 - x^dd0/(-(-dum)^dd0 + 
                                                              x^dd0))
           if (increase == -1) y.expect <- aa0 * (1 - x^dd0/((dum)^dd0 + 
                                                               x^dd0))
         }, {
           dum <- f.bb.con(model.ans, cc = cc0, dd = NA, CED = bb0, 
                           CES)
           y.expect <- aa0 * (1 + ((cc0 - 1) * x)/(dum + x))
         }, {
           dum <- f.bb.con(model.ans, cc = cc0, dd = dd0, CED = bb0, 
                           CES)
           y.expect <- aa0 * (1 + ((cc0 - 1) * x^dd0)/(dum^dd0 + 
                                                         x^dd0))
         }, {
           y.expect <- cc0 + aa0 * x^bb0
         }, y.expect <- exp(x/bb0)^aa0, y.expect <- cc0 + (bb0 * 
                                                             abs((x - aa0))^dd0) * (x >= aa0), y.expect <- aa0 + 
           cc0 * (1 - exp(-(x/bb0)^dd0)), {
             y.expect <- aa0 + bb0 * x^cc0
           }, y.expect <- aa0 * exp((x/bb0)^dd0) + cc0, {
             b2 <- regr.par[nr.aa + nr.bb + 1]
             cc <- regr.par[nr.aa + nr.bb + 2]
             d1 <- regr.par[nr.aa + nr.bb + 3]
             d2 <- regr.par[nr.aa + nr.bb + 4]
             y.expect.1 <- cc - (cc - 1) * exp(-(x/bb0)^d1)
             y.expect.2 <- cc - (cc - 1) * exp(-(x/b2)^d2)
             y.expect <- aa0 * (y.expect.1 * y.expect.2)
           }, {
             bb0 <- rep(0, length(x))
             bb.tmp <- regr.par[(nr.aa + 1):(nr.aa + nr.bb)]
             for (jj in (1:nr.bb)) bb0 <- bb0 + bb.tmp[jj] * (fct2 == 
                                                                jj)
             b2 <- regr.par[nr.aa + nr.bb + 1]
             c1 <- regr.par[nr.aa + nr.bb + 2]
             c2 <- regr.par[nr.aa + nr.bb + 3]
             d1 <- regr.par[nr.aa + nr.bb + 4]
             d2 <- regr.par[nr.aa + nr.bb + 5]
             y.expect.1 <- c1 - (c1 - 1) * exp(-(x/bb0)^d1)
             y.expect.2 <- c2 - (c2 - 1) * exp(-(x/b2)^d2)
             y.expect <- aa0 * (y.expect.1 * y.expect.2)
           }, {
             y.expect <- aa0 + bb0 * x + cc0 * x^2 + dd0 * x^3
           }, y.expect <- (aa0/bb0) * (1 - exp(-bb0 * x)), y.expect <- aa0 * 
           (1 + (cc0 - 1) * (x^dd0/(bb0^dd0 + x^dd0))), y.expect <- bb0 * 
           exp(cc0 * (x - aa0)^dd0) * (x >= aa0) + bb0 * (x < 
                                                            aa0), cat(""), cat(""), cat(""), y.expect <- aa0 * 
           x^bb0, {
           }, {
             y.expect <- aa0 + cc0 * x^dd0/(bb0^dd0 + x^dd0)
           }, {
             y.expect <- exp(aa0 + bb0 * (x - x.mn))
           }, {
             c1 <- regr.par[nr.aa + nr.bb + 1]
             c2 <- regr.par[nr.aa + nr.bb + 2]
             dd <- regr.par[nr.aa + nr.bb + 3]
             y.expect.1 <- c1 - (c1 - 1) * exp(-(x/bb0)^dd)
             y.expect.2 <- c2 - (c2 - 1) * exp(-(x/bb0)^dd)
             y.expect <- aa0 * (y.expect.1 * y.expect.2)
           }, {
             CED.pars <- regr.par[(nr.aa + 1):(nr.aa + nr.bb)]
             CED.ref <- CED.pars[ref.lev]
             bb0 <- CED.ref/bb0
             bb0[fct2 == ref.lev] <- CED.ref
             y.expect <- aa0 * (cc0 - (cc0 - 1) * ((CES + 1 - 
                                                      cc0)/(1 - cc0))^((x/bb0)^(dd0)))
             if (0) {
               CED.ref <- regr.par[nr.aa + ref.lev]
               dum <- f.bb.con(model.ans, cc = cc0, dd = dd0, 
                               CED = bb0, CES, ref.lev = ref.lev, CED.ref = CED.ref)
               y.expect <- aa0 * (dd0 - (dd0 - 1) * exp(-dum * 
                                                          (x^cc0)))
             }
           }, {
             if (0) {
               print("f.expect.con  pars")
               print(model.ans)
               print(CES)
               print(sd0[1:5])
               print(aa0[1:5])
               print(bb0[1:5])
               print(qq)
               print(dd0[1:5])
             }
             y.expect <- exp(log(aa0) + sign.q0 * qq * sd0 * (1 - 
                                                                exp(x^dd0 * log(1 - CES)/bb0^dd0)))
           }, {
             y.expect = aa0 * (cc0 - (cc0 - 1) * exp(-(x/bb0)^dd0 - 
                                                       (ttt/gg)^hh))
           }, {
             y.expect <- exp(log(aa0) + log(cc0) * (1 - exp(-log(2) * 
                                                              (x/bb0)^dd0)))
           }, {
             bb.time <- bb0 * exp(-(ttt/gg)^hh)
             y.expect <- aa0[1] * (cc0[1] - (cc0[1] - 1) * exp(-(x/bb.time)^dd0[1]))
           }, y.expect <- aa0 * ((1 - cc0) - cc0 * exp(-bb0 * x^dd0)))
  if (exists("track2")) 
    print("f.expect.con:  END")
  return(y.expect)
}

##################
## f.factorname ##
##################
f.factorname <- function (ans.all) 
{
  if (exists("track")) 
    print("f.factorname")
  with(ans.all, {
    factor.name <- ""
    fact1.name <- ""
    fact2.name <- ""
    fact3.name <- ""
    fact4.name <- ""
    fact5.name <- ""
    if (quick.ans > 1) {
      fct1.no <- 0
      fct2.no <- 0
      fct3.no <- 0
      if (nr.aa > 1) 
        fct1.no <- covar.no
      if (nr.bb > 1) 
        fct2.no <- covar.no
      if (nr.var > 1) 
        fct3.no <- covar.no
    }
    if (fct1.no > 0) {
      fact1.name <- paste("\nfact-a: ", varnames[fct1.no])
      factor.name <- paste(factor.name, fact1.name)
    }
    if (fct2.no > 0) {
      fact2.name <- paste("\nfact-b: ", varnames[fct2.no])
      factor.name <- paste(factor.name, fact2.name)
    }
    if (fct3.no > 0) {
      fact3.name <- paste("\nfact-var: ", varnames[fct3.no])
      factor.name <- paste(factor.name, fact3.name)
    }
    if (fct4.no > 0) {
      fact4.name <- paste("\nfact-d: ", varnames[fct4.no])
      factor.name <- paste(factor.name, fact4.name)
    }
    if (fct5.no > 0) {
      fact5.name <- paste("\nfact-c: ", varnames[fct5.no])
      factor.name <- paste(factor.name, fact5.name)
    }
    if (exists("track")) 
      print("f.factorname  END")
    return(factor.name)
  })
}

#####################
## f.nr.replicates ##
#####################
f.nr.replicates <- function (ans.all) 
{
  if (exists("track")) 
    print("f.nr.replicates")
  with(ans.all, {
    if (length(xans) > 1) {
      x <- data.0[, xans[1]]
    }
    if (quick.ans == 1) 
      if (length(fct3) > 1) 
        dum.nn <- table(x, fct1, fct2, fct3, fct5)
      else dum.nn <- table(x, fct1, fct2, fct5)
      if (quick.ans > 1) 
        dum.nn <- table(x, covariate)
      if (dtype %in% c(4, 6, 10, 15, 250, 260)) 
        dum.nn <- nn
      dum.nn <- as.numeric(dum.nn)
      dum.nn <- dum.nn[dum.nn != 0]
      return(dum.nn)
  })
}

################
## f.full.ans ##
################
f.full.ans <- function (ans.all, gui) 
{
  if (exists("track")) 
    print("f.full.ans")
  if (gui == T) 
    ans.all$gui <- T
  with(ans.all, {
    if (dtype %in% c(2, 3)) {
      full.ans <- 2
      return(full.ans)
    }
    if (dtype %in% c(6)) {
      full.ans <- 1
      return(full.ans)
    }
    dum.nn <- f.nr.replicates(ans.all)
    if (cont) {
      ans.all$model.ans <- 11
      regr.par <- f.regr.par.full(ans.all)
      nr.dosegr <- length(regr.par)
    }
    if (!cont) {
      ans.all$model.type <- 1
      ans.all$model.ans <- 14
      ans.all <- f.start.bin(ans.all, tmp.quick = T)
      nr.dosegr <- length(ans.all$regr.par)
    }
    full.ans <- ans.all$full.ans
    if (!gui && dtype != 15) {
      if (mean(dum.nn) < 4) {
        #cat("\n\naverage group size is only", round(mean(dum.nn), 
        #                                            1))
        full.ans <- 2
        #full.ans <- menu(c("yes", "no"), title = "\nDo you nonetheless want to fit full model?\n")
      }
      else if (nr.dosegr > 100 || length(dum.nn) > 100) {
        cat("\ndataset is large (total number of dose groups: ", 
            nr.dosegr, "), fitting full model may take long")
        full.ans <- menu(c("yes", "no"), title = "\nDo you nonetheless want to fit full model?\n")
      }
    }
    if (gui) {
      full.ans <- 1
      if (mean(dum.nn) < 4) {
        #cat("\n\naverage group size is only", round(mean(dum.nn), 
        #                                            1))
        #cat("\ntherefore full model will not be fitted\n")
        full.ans <- 2
      }
      else if (nr.dosegr > 100 || length(dum.nn) > 100) {
        cat("\ndataset is large, therefore full model will not be fitted\n")
        full.ans <- 2
      }
    }
    if (exists("track")) 
      print("f.full.ans:   END")
    return(full.ans)
  })
}

#####################
## f.regr.par.full ##
#####################
f.regr.par.full <- function (ans.all) 
{
  if (exists("track")) 
    print("f.regr.par.full")
  with(ans.all, {
    if (dtype %in% c(10, 15, 250, 260)) {
      if (dtype %in% c(10, 15)) 
        regr.par <- exp(mn.log)
      if (dtype == 250) 
        regr.par <- (mn.log)
      if (dtype == 260) 
        regr.par <- (mn.log)^2
    }
    else {
      regr.par <- numeric()
      nn <- numeric()
      if (length(xans) > 1) 
        x <- data.0[, xans[1]]
      if (nr.var == 1) {
        if (twice) 
          for (jj in levels(factor(fct2))) {
            x.tmp <- x[fct2 == jj]
            y.tmp <- yy[fct2 == jj]
            regr.tmp <- as.numeric(exp(tapply(logb(y.tmp), 
                                              x.tmp, mean)))
            regr.par <- c(regr.par, regr.tmp)
            nn.tmp <- tapply(y.tmp, x.tmp, length)
            nn <- c(nn, nn.tmp)
          }
        if (!twice) 
          for (jj in levels(factor(fct2))) for (ii in levels(factor(fct1))) {
            x.tmp <- x[fct1 == ii & fct2 == jj]
            y.tmp <- yy[fct1 == ii & fct2 == jj]
            regr.tmp <- as.numeric(exp(tapply(logb(y.tmp), 
                                              x.tmp, mean)))
            regr.par <- c(regr.par, regr.tmp)
            nn.tmp <- tapply(y.tmp, x.tmp, length)
            nn <- c(nn, nn.tmp)
          }
      }
      if (nr.var > 1) {
        regr.par <- numeric()
        for (jj in levels(factor(fct3))) {
          x.tmp <- x[fct3 == jj]
          y.tmp <- yy[fct3 == jj]
          regr.tmp <- as.numeric(exp(tapply(logb(y.tmp), 
                                            x.tmp, mean)))
          regr.par <- c(regr.par, regr.tmp)
          nn.tmp <- tapply(y.tmp, x.tmp, length)
          nn <- c(nn, nn.tmp)
        }
      }
      ans.all$nn <- nn
    }
    return(regr.par)
  })
}

#############################
## f.press.key.to.continue ##
#############################
f.press.key.to.continue <- function (prompt = NA) 
{
  #if (!is.na(prompt)) 
  #  prompt <- paste(prompt, "\nPress enter key to continue ... > ")
  #else prompt <- "\nPress enter key to continue ... > "
  #cat(readline(prompt = prompt))
  #return(invisible())
}

#####################
## f.select.m5.con ##
#####################
f.select.m5.con <- function (ans.all, output = TRUE) 
{
  if (exists("track")) 
    print("f.select.m5.con")
  gui <- ans.all$gui
  hill <- ans.all$hill
  cont <- ans.all$cont
  quick.ans <- ans.all$quick.ans
  x <- ans.all$x
  y <- ans.all$y
  yy <- ans.all$yy
  nn <- ans.all$nn
  if (ans.all$dtype.0 == 5) 
    ans.all$dtype <- 5
  if (ans.all$dtype.0 == 15) 
    ans.all$dtype <- 15
  dtype <- ans.all$dtype
  if (dtype == 5) 
    ans.all$dtype.0 <- 5
  nth <- ans.all$nth
  covariate <- ans.all$covariate
  covar.no <- ans.all$covar.no
  covar.txt <- ans.all$covar.txt
  CES <- abs(ans.all$CES)
  ces.ans <- ans.all$ces.ans
  if (dtype == 3) 
    ans.all$CES.cat <- 1
  CED <- NA
  cens.up <- ans.all$cens.up
  select.name <- ans.all$select.name
  factor.name <- ans.all$factor.name
  if (cont) 
    plot.type <- ans.all$plot.type
  else plot.type <- 0
  color <- ans.all$color
  th.start <- ans.all$th.start
  sig.start <- ans.all$sig.start
  alfa.length <- ans.all$alfa.length
  alfa.start <- ans.all$alfa.start
  full.ans <- ans.all$full.ans
  interrupt <- ans.all$interrupt
  nr.covar <- max(covariate)
  alert.full <- ans.all$alert.full
  output <- ans.all$output
  if (cont) 
    nr.var <- 1
  vv <- 1
  if (cont && hill == 10) {
    vv <- ans.all$vv
    nr.var <- length(ans.all$vv)
  }
  if (!cont) 
    nr.var <- 0
  if (hill == 0) {
    skip.b <- F
    skip.0 <- F
  }
  if (nr.covar > 10) 
    skip.0 <- T
  dd.start <- 1
  if (hill == 10) {
    var.start <- ans.all$var.start
    Vaa.start <- ans.all$Vaa.start
    aa.start <- ans.all$aa.start
    Vced.start <- ans.all$Vced.start
    ced.start <- ans.all$ced.start
    cc.start <- ans.all$cc.start
    ths.start <- ans.all$ths.start
    if (vv == 2) 
      fct3 <- covariate
    else fct3 <- rep(1, length(x))
    if (cont) 
      nr.var <- max(fct3)
    if (dtype == 6) 
      alfa.start <- ans.all$alfa.start
    aic.m1 <- ans.all$aic.m1
    row.sel.m1 <- ans.all$row.sel.m1
    skip.0 <- ans.all$skip.0
    skip.b <- ans.all$skip.b
  }
  ii <- 0
  if (hill == 0) {
    #title.tmp <- "nested exponential models"
  }
  if (hill == 10) {
    #title.tmp <- "nested Hill models"
  }
  name.tmp <- ""
  if (!ans.all$WAPP && output && interrupt && cont) 
    #f.graph.window(12, WAPP = ans.all$WAPP, title = title.tmp, 
    #               name.wapp = name.tmp, plotprefix = ans.all$plotprefix)
  tb <- "\t"
  model.txt <- character()
  Vconverged <- numeric()
  Vloglik <- numeric()
  Vaic <- numeric()
  Vnpar <- integer()
  if (cont) 
    MLE.all <- matrix(NA, 18, 3 * nr.covar + 4)
  if (dtype %in% 2:3) 
    MLE.all <- matrix(NA, 18, 2 * nr.covar + 2 + nth + 1)
  if (dtype == 4) 
    MLE.all <- matrix(NA, 18, 2 * nr.covar + 2 + 1 + 1)
  if (dtype == 6) 
    MLE.all <- matrix(NA, 18, 1 + 2 * nr.covar + 2 + 1 + 
                        1)
  if (length(ans.all$xans) > 1) 
    MLE.all <- cbind(MLE.all, rep(NA, 18))
  if (length(ans.all$xans) > 2) 
    MLE.all <- cbind(MLE.all, rep(NA, 18))
  alert.full <- NA
  ans.all.fit <- ans.all
  ans.all.fit$par.start <- NA
  if (ans.all$displ.ans == 0) 
    ans.all.fit$displ.fact <- covariate
  nr.full <- 0
  if (output) {
    if (hill == 0) {}
      #cat("\nANALYSIS WITH EXPONENTIAL MODELS")
    else {} #cat("\nANALYSIS WITH HILL MODELS")
    #cat("\nmodel", tb, "converged", tb, "npar", tb, "loglik", 
    #    tb, "aic")
  }
  if (hill == 0) {
    fct3 <- rep(1, length(x))
    if (full.ans == 1) {
      ii <- ii + 1
      model.txt[ii] <- "full   "
      ans.all.fit$model.ans <- 11
      if (dtype %in% c(4, 6)) {
        ans.all.fit$model.type <- 1
        ans.all.fit$model.ans <- 14
      }
      ans.all.fit$modelname <- model.txt[ii]
      ans.all.fit$fct1 <- covariate
      ans.all.fit$fct2 <- covariate
      ans.all.fit$fct3 <- fct3
      if (dtype != 6) 
        ans.all.fit <- f.qfit(ans.all.fit, plot.type, 
                              output = output)
      else {
        ans.all.fit <- f.dtype6.mn(ans.all.fit)
        ans.all$full6.done <- TRUE
        ans.all$alfa <- ans.all.fit$alfa
        ans.all$Pvalue.alfa <- ans.all.fit$Pvalue.alfa
        ans.all$x.full <- ans.all.fit$x.full
        ans.all$pi.full <- ans.all.fit$pi.full
        ans.all$fct1.full <- ans.all.fit$fct1.full
        ans.all$fct2.full <- ans.all.fit$fct2.full
      }
      ans.all.fit$kk.tot <- ans.all$kk.tot
      ans.all.fit$nn.tot <- ans.all$nn.tot
      Vnpar[ii] <- ans.all.fit$npar.aic
      Vloglik[ii] <- round(ans.all.fit$loglik, 2)
      Vaic[ii] <- ans.all.fit$aic
      Vconverged[ii] <- ans.all.fit$converged
      MLE <- ans.all.fit$MLE
      MLE.all[ii, ] <- signif(MLE[1:length(MLE.all[1, ])], 
                              7)
      if (output) {}
        #cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
        #          tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
      if (cont) {
        var.start <- MLE[1]
        ans.all$var.start <- var.start
      }
      nr.full <- 1
      results.full <- f.model.specific.results(ans.all.fit)
      if (ans.all$const.var == TRUE) 
        vv <- 1
      else if (cont && covar.no > 0) {
        ii <- ii + 1
        model.txt[ii] <- paste("full-", "v", sep = "")
        ans.all.fit$model.ans <- 11
        ans.all.fit$modelname <- model.txt[ii]
        fct3 <- covariate
        ans.all.fit$fct3 <- fct3
        if (cont) 
          nr.var <- max(fct3)
        ans.all.fit$par.start <- c(rep(var.start, nr.var), 
                                   ans.all.fit$par.start[-1])
        ans.all.fit <- f.qfit(ans.all.fit, plot.type = 0, 
                              output = output)
        Vnpar[ii] <- ans.all.fit$npar.aic
        Vloglik[ii] <- round(ans.all.fit$loglik, 2)
        Vaic[ii] <- ans.all.fit$aic
        Vconverged[ii] <- ans.all.fit$converged
        MLE <- ans.all$MLE
        if (dtype == 6) 
          alfa <- MLE[1]
        if (Vaic[ii] > Vaic[ii - 1]) {
          vv <- 1
          fct3 <- rep(1, length(x))
          nr.var <- 1
        }
        else {
          vv <- 2
          var.start <- ans.all.fit$MLE[1:nr.var]
          ans.all$var.start <- var.start
          results.full <- f.model.specific.results(ans.all.fit)
        }
        MLE <- ans.all.fit$MLE
        MLE.all[ii, ] <- signif(MLE[1:length(MLE.all[1, 
                                                     ])], 7)
        if (output) {}
          #cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
          #          tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
        nr.full <- 2
      }
      ans.all$FULL.model <- results.full
    }
    if (full.ans == 2) {
      ii <- 1
      model.txt[1] <- "full"
      if (cont) 
        nr.var <- 1
      if (output) {
        if (dtype == 6) {
          Vnpar[ii] <- ans.all$npar
          Vloglik[ii] <- round(ans.all$loglik, 2)
          Vaic[ii] <- ans.all$aic
          Vconverged[ii] <- ans.all$converged
          #cat(paste("\n", model.txt[ii], tb, tb, Vconverged[ii], 
          #          tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
          ans.all$loglik.full <- Vloglik[ii]
          ans.all$aic.full <- Vaic[ii]
          ans.all$npar.aic.full <- Vnpar[ii]
        }
        else {
          Vloglik[1] <- NA
          Vaic[1] <- NA
          Vnpar[1] <- NA
          MLE.all[1, 1:4] <- 0
          #cat(paste("\n", model.txt[ii], tb, "not evaluated"))
        }
      }
    }
    if (!cont) 
      ans.all.fit$model.type <- 2
    ii <- ii + 1
    row.sel.m1 <- ii
    model.txt[ii] <- "null-"
    if (max(fct3) > 1) 
      model.txt[ii] <- paste(model.txt[ii], "v", sep = "")
    model.ans <- 1
    ans.all.fit$model.ans <- model.ans
    ans.all.fit$modelname <- model.txt[ii]
    fct1 <- rep(1, length(x))
    fct2 <- rep(1, length(x))
    ans.all.fit$fct1 <- fct1
    ans.all.fit$fct2 <- fct2
    ans.all.fit$fct3 <- fct3
    ans.all.fit$par.start <- NA
    if (cont && full.ans == 1) 
      ans.all.fit$par.start <- c(var.start, exp(mean(log(y))))
    ans.all.fit <- f.qfit(ans.all.fit, plot.type, output = output)
    if (is.na(ans.all.fit$loglik)) {
      if (output) 
        cat("\nATTENTION:  model 1 did not result in finite log-likelihood value \n")
      ans.all.fit$loglik <- NA
    }
    Vnpar[ii] <- ans.all.fit$npar.aic
    Vloglik[ii] <- round(ans.all.fit$loglik, 2)
    Vaic[ii] <- ans.all.fit$aic
    Vconverged[ii] <- ans.all.fit$converged
    MLE <- ans.all.fit$MLE
    MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
    if (output) {}
      #cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
      #          tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
    aa.start <- MLE[nr.var + 1]
    if (cont && full.ans == 2) 
      var.start <- MLE[1]
    if (!cont) {
      npar.tmp <- length(MLE)
      ths.start <- MLE[(alfa.length + 1 + 1):(npar.tmp)]
    }
    aic.m1 <- ans.all.fit$aic
    results.m1 <- f.model.specific.results(ans.all.fit)
    if (covar.no > 0) {
      ii <- ii + 1
      model.txt[ii] <- "null-a"
      if (max(fct3) > 1) 
        model.txt[ii] <- paste(model.txt[ii], "v", sep = "")
      model.ans <- 1
      ans.all.fit$model.ans <- model.ans
      ans.all.fit$modelname <- model.txt[ii]
      fct1 <- covariate
      fct2 <- rep(1, length(x))
      nr.aa <- max(fct1)
      ans.all.fit$fct1 <- fct1
      ans.all.fit$fct2 <- fct2
      ans.all.fit$fct3 <- fct3
      if (cont) 
        par.start <- c(var.start, rep(aa.start, nr.aa))
      if (!cont) 
        par.start <- c(rep(aa.start, nr.aa), ths.start)
      if (dtype == 6) 
        par.start <- c(alfa.start, par.start)
      ans.all.fit$par.start <- par.start
      ans.all.fit <- f.qfit(ans.all.fit, plot.type, output = output)
      Vnpar[ii] <- ans.all.fit$npar.aic
      Vloglik[ii] <- round(ans.all.fit$loglik, 2)
      Vaic[ii] <- ans.all.fit$aic
      Vconverged[ii] <- ans.all.fit$converged
      MLE <- ans.all.fit$MLE
      MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
      if (is.na(Vloglik[ii])) {
        if (output) 
          cat("\nATTENTION:  model 1 did not result in finite log-likelihood value \n")
        Vloglik[ii] <- NA
      }
      if (output) {
        #cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
        #          tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
      }
      if (cont && full.ans == 2) 
        var.start <- MLE[1:nr.var]
      Vaa.start <- MLE[(nr.var + 1):(nr.var + nr.aa)]
      if (Vaic[ii] < Vaic[ii - 1]) {
        results.m1 <- f.model.specific.results(ans.all.fit)
        row.sel.m1 <- ii
        aic.m1 <- ans.all.fit$aic
      }
    }
    ans.all$NULL.model <- results.m1
    ans.all.fit$model.ans <- 2
    ans.all.fit$modelname <- "m2"
    if (max(covariate) > 1) 
      ans.all.fit$modelname <- "m2-ab"
    if (max(fct3) > 1) 
      ans.all.fit$modelname <- paste(ans.all.fit$modelname, 
                                     "v", sep = "")
    fct1 <- covariate
    fct2 <- covariate
    nr.aa <- max(fct1)
    nr.bb <- max(fct2)
    if (cont) 
      nr.var <- max(fct3)
    ans.all.fit$fct1 <- fct1
    ans.all.fit$fct2 <- fct2
    if (vv == 1) 
      ans.all.fit$fct3 <- rep(1, length(x))
    if (vv == 2) 
      ans.all.fit$fct3 <- covariate
    ans.all.fit$nr.aa <- nr.aa
    ans.all.fit$nr.bb <- nr.bb
    if (cont) 
      ans.all.fit$nr.var <- max(ans.all.fit$fct3)
    if (cont) 
      par.start.tmp <- f.start.con(ans.all.fit, tmp.quick = T)$par.start
    if (!cont) 
      par.start.tmp <- f.start.cat(ans.all.fit, tmp.quick = T)$par.start
    if (covar.no == 0) 
      aa.start.tmp <- aa.start
    if (covar.no > 0) 
      aa.start.tmp <- Vaa.start
    bb.start <- par.start.tmp[(alfa.length + nr.var + nr.aa + 
                                 1):(alfa.length + nr.var + nr.aa + nr.bb)]
    if (cont) {
      ans.all.fit$par.start <- c(var.start, aa.start.tmp, 
                                 bb.start)
    }
    else ans.all.fit$par.start <- NA
    ans.all.fit <- f.qfit(ans.all.fit, plot.type, output = output)
    MLE <- ans.all.fit$MLE
    if (0) {
      if (vv == 2) {}
        #cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
        #          tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
      else {}#cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
                     #tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
    }
    Vaa.start <- MLE[(alfa.length + nr.var + 1):(alfa.length + 
                                                   nr.var + nr.aa)]
    aa.start <- exp(mean(log(Vaa.start)))
    ans.all$Vaa.start <- Vaa.start
    ans.all$aa.start <- aa.start
    Vbb.start <- MLE[(alfa.length + nr.var + max(fct1) + 
                        1):(alfa.length + nr.var + max(fct1) + max(fct2))]
    if (!cont) {
      npar.tmp <- length(MLE)
      ths.start <- MLE[(alfa.length + nr.aa + nr.bb + 1):npar.tmp]
    }
    if (output) 
      if (dtype != 6 & alfa.length != 0) {
        print("f.select.m5.con")
        cat("\nalfa.length not equal to zero\n")
        f.press.key.to.continue()
      }
    if (cont) {
      if (abs(sum(sign(Vbb.start))) != length(Vbb.start)) 
        cat("\nSlopes have opposite signs")
      if (mean(Vbb.start) > 0) {
        ans.all.fit$increase <- 1
        ans.all$increase <- 1
        ans.all.fit$CES <- CES
        ans.all$CES <- CES
      }
      else {
        ans.all.fit$increase <- -1
        ans.all.fit$CES <- -CES
        ans.all$CES <- -CES
        ans.all$increase <- -1
      }
    }
    if (cont) {
      regr.par.matr <- ans.all.fit$regr.par.matr
      bb.tmp <- regr.par.matr[, 2]
      if (ans.all$increase == 1) {
        bb.tmp[bb.tmp < 0] <- min(bb.tmp[bb.tmp > 0])
        regr.par.matr[, 2] <- bb.tmp
      }
      if (ans.all$increase == -1) {
        bb.tmp[bb.tmp > 0] <- max(bb.tmp[bb.tmp < 0])
        regr.par.matr[, 2] <- bb.tmp
      }
      ced.lst <- f.ced.con(model.ans = 2, regr.par.matr, 
                           CES = ans.all.fit$CES, nr.aa = nr.aa, nr.bb = nr.bb, 
                           nr.var = 1)
      Vced.start <- ced.lst$CED
      ced.start <- exp(mean(log(Vced.start)))
    }
    if (!cont) {
      ced.lst <- f.ced.cat(ans.all.fit)
      Vced.start <- ced.lst$CED.matr[, 1]
      ced.start <- exp(mean(log(Vced.start)))
      Vced.start[Vced.start == Inf] <- Vced.start
    }
    ans.all$Vced.start <- Vced.start
    ans.all$ced.start <- ced.start
    if (ans.all$increase == -1) 
      cc.start <- 0.1
    if (ans.all$increase == 1) 
      cc.start <- 10
    if (!cont) 
      cc.start <- 0.1
    ans.all$cc.start <- cc.start
  }
  ans.all.fit$increase <- ans.all$increase
  if (hill == 0) 
    nr.first.models <- ii
  else nr.first.models <- 0
  ii.m3 <- 0
  m3x <- numeric()
  Vaic.m3 <- numeric()
  if (!skip.0) {
    model.ans <- 13 + hill
    ii <- ii + 1
    ii.m3 <- ii.m3 + 1
    if (hill == 0) 
      model.txt[ii] <- "Expon. m3-"
    if (hill == 10) 
      model.txt[ii] <- "Hill m3-"
    fct1 <- rep(1, length(x))
    fct2 <- rep(1, length(x))
    if (vv == 2) 
      fct3 <- covariate
    nr.aa <- max(fct1)
    nr.bb <- max(fct2)
    if (cont) 
      nr.var <- max(fct3)
    if (max(fct3) > 1) 
      model.txt[ii] <- paste(model.txt[ii], "v", sep = "")
    ans.all.fit$model.ans <- model.ans
    ans.all.fit$modelname <- model.txt[ii]
    ans.all.fit$fct1 <- fct1
    ans.all.fit$fct2 <- fct2
    ans.all.fit$fct3 <- fct3
    ans.all.fit$nr.aa <- nr.aa
    ans.all.fit$nr.bb <- nr.bb
    ans.all.fit$nr.var <- nr.var
    if (cont) 
      par.start <- c(var.start, aa.start, ced.start, 1)
    if (!cont) 
      par.start <- c(aa.start, ced.start, 1, ths.start)
    if (dtype == 6) 
      par.start <- c(alfa.start, par.start)
    ans.all.fit$par.start <- par.start
    ans.all.fit <- f.qfit(ans.all.fit, plot.type, output = output)
    Vnpar[ii] <- ans.all.fit$npar.aic
    Vloglik[ii] <- round(ans.all.fit$loglik, 2)
    Vaic[ii] <- ans.all.fit$aic
    Vconverged[ii] <- ans.all.fit$converged
    MLE <- ans.all.fit$MLE
    MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
    Vaic.m3[ii.m3] <- Vaic[ii]
    m3x[ii.m3] <- "0"
    if (output) {
      #cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
      #          tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
    }
    aa.start <- MLE[alfa.length + nr.var + 1]
    if (full.ans == 2) 
      var.start <- MLE[1:nr.var]
    results.m3.0 <- f.model.specific.results(ans.all.fit)
  }
  if (covar.no > 0) {
    model.ans <- 13 + hill
    ii <- ii + 1
    ii.m3 <- ii.m3 + 1
    if (hill == 0) 
      model.txt[ii] <- "Expon. m3-a"
    if (hill == 10) 
      model.txt[ii] <- "Hill m3-a"
    if (max(fct3) > 1) 
      model.txt[ii] <- paste(model.txt[ii], "v", sep = "")
    ans.all.fit$model.ans <- model.ans
    ans.all.fit$modelname <- model.txt[ii]
    fct1 <- covariate
    ans.all.fit$fct1 <- fct1
    ans.all.fit$fct2 <- rep(1, length(x))
    ans.all.fit$fct3 <- fct3
    nr.aa <- max(fct1)
    nr.bb <- 1
    ans.all.fit$nr.aa <- nr.aa
    ans.all.fit$nr.bb <- nr.bb
    ans.all.fit$nr.var <- nr.var
    if (cont) 
      par.start <- c(var.start, Vaa.start, ced.start, 1)
    if (!cont) 
      par.start <- c(aa.start, Vced.start, 1, ths.start)
    if (dtype == 6) 
      par.start <- c(alfa.start, par.start)
    ans.all.fit$par.start <- par.start
    ans.all.fit <- f.qfit(ans.all.fit, plot.type, output = output)
    Vnpar[ii] <- ans.all.fit$npar.aic
    Vloglik[ii] <- round(ans.all.fit$loglik, 2)
    Vaic[ii] <- ans.all.fit$aic
    Vconverged[ii] <- ans.all.fit$converged
    MLE <- ans.all.fit$MLE
    MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
    Vaic.m3[ii.m3] <- Vaic[ii]
    m3x[ii.m3] <- "a"
    if (output) {}
      #cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
      #          tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
    Vaa.start <- MLE[(alfa.length + nr.var + 1):(alfa.length + 
                                                   nr.var + max(fct1))]
    results.m3.a <- f.model.specific.results(ans.all.fit)
    if (hill == 0) 
      if (Vaic[ii - 1] - Vaic[ii] > 10) {
        skip.b <- T
        skip.0 <- T
      }
    else {
      skip.b <- F
    }
    if (!skip.b) {
      ii <- ii + 1
      ii.m3 <- ii.m3 + 1
      if (hill == 0) 
        model.txt[ii] <- "Expon. m3-b"
      if (hill == 10) 
        model.txt[ii] <- "Hill m3-b"
      if (max(fct3) > 1) 
        model.txt[ii] <- paste(model.txt[ii], "v", sep = "")
      ans.all.fit$modelname <- model.txt[ii]
      fct1 <- rep(1, length(x))
      fct2 <- covariate
      ans.all.fit$fct1 <- fct1
      ans.all.fit$fct2 <- fct2
      ans.all.fit$fct3 <- fct3
      nr.aa <- 1
      nr.bb <- max(fct2)
      ans.all.fit$nr.aa <- nr.aa
      ans.all.fit$nr.bb <- nr.bb
      ans.all.fit$nr.var <- nr.var
      if (cont) 
        par.start <- c(var.start, aa.start, Vced.start, 
                       1)
      if (!cont) 
        par.start <- c(aa.start, Vced.start, 1, ths.start)
      if (dtype == 6) 
        par.start <- c(alfa.start, par.start)
      ans.all.fit$par.start <- par.start
      ans.all.fit <- f.qfit(ans.all.fit, plot.type, output = output)
      Vnpar[ii] <- ans.all.fit$npar.aic
      Vloglik[ii] <- round(ans.all.fit$loglik, 2)
      Vaic[ii] <- ans.all.fit$aic
      Vconverged[ii] <- ans.all.fit$converged
      MLE <- ans.all.fit$MLE
      MLE.2.b <- MLE
      MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
      Vaic.m3[ii.m3] <- Vaic[ii]
      m3x[ii.m3] <- "b"
      if (output) {}
        #cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
        #          tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
      aa.start.mb <- MLE[alfa.length + nr.var + 1]
      results.m3.b <- f.model.specific.results(ans.all.fit)
    }
    ii <- ii + 1
    ii.m3 <- ii.m3 + 1
    if (hill == 0) 
      model.txt[ii] <- "Expon. m3-ab"
    if (hill == 10) 
      model.txt[ii] <- "Hill m3-ab"
    if (max(fct3) > 1) 
      model.txt[ii] <- paste(model.txt[ii], "v", sep = "")
    ans.all.fit$modelname <- model.txt[ii]
    fct1 <- covariate
    fct2 <- covariate
    ans.all.fit$fct1 <- fct1
    ans.all.fit$fct2 <- fct2
    ans.all.fit$fct3 <- fct3
    nr.aa <- max(fct1)
    nr.bb <- max(fct2)
    if (cont) 
      nr.var <- max(fct3)
    ans.all.fit$nr.aa <- nr.aa
    ans.all.fit$nr.bb <- nr.bb
    ans.all.fit$nr.var <- nr.var
    if (cont) 
      par.start <- c(var.start, rep(aa.start, nr.aa), Vced.start, 
                     1)
    if (!cont) 
      par.start <- c(rep(aa.start, nr.aa), Vced.start, 
                     1, ths.start)
    if (dtype == 6) 
      par.start <- c(alfa.start, par.start)
    ans.all.fit$par.start <- par.start
    ans.all.fit <- f.qfit(ans.all.fit, plot.type, output = output)
    Vnpar[ii] <- ans.all.fit$npar.aic
    Vloglik[ii] <- round(ans.all.fit$loglik, 2)
    Vaic[ii] <- ans.all.fit$aic
    Vconverged[ii] <- ans.all.fit$converged
    MLE <- ans.all.fit$MLE
    MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
    Vaic.m3[ii.m3] <- Vaic[ii]
    m3x[ii.m3] <- "ab"
    if (output) {
      if (vv == 2) {}
        #cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
        #          tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
      else {}#cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
              #       tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
    }
    Vaa.start <- MLE[(alfa.length + nr.var + 1):(alfa.length + 
                                                   nr.var + max(fct1))]
    Vced.start <- MLE[(alfa.length + nr.var + max(fct1) + 
                         1):(alfa.length + nr.var + max(fct1) + max(fct2))]
    if (!cont) {
      npar.tmp <- length(MLE)
      ths.start <- MLE[(alfa.length + nr.aa + nr.bb + 2):(npar.tmp)]
    }
    if (dtype != 6 & alfa.length != 0) {
      print("f.select.m5.con")
      cat("\nalfa.length not equal to zero\n")
      f.press.key.to.continue()
    }
    results.m3.ab <- f.model.specific.results(ans.all.fit)
  }
  rank.m3 <- order(Vaic.m3)[1]
  row.sel.m3 <- nr.first.models + rank.m3
  if (quick.ans == 4) 
    if (covar.no > 0) {
      if (output) {
        cat("\n---------------------------------------------------\n")
        cat("Best model 3 with covariates is:", model.txt[nr.first.models + 
                                                            rank.m3], "\n")
      }
      if (ii.m3 > 1) {
        rank.m3.2 <- order(Vaic.m3)[2]
        if (output && Vaic.m3[rank.m3.2] - Vaic.m3[rank.m3] <= 
            2) 
          cat("However", model.txt[nr.first.models + 
                                     rank.m3.2], "is a reasonable model as well\n")
      }
      if (ii.m3 > 2) {
        rank.m3.3 <- order(Vaic.m3)[3]
        if (output && Vaic.m3[rank.m3.3] - Vaic.m3[rank.m3] <= 
            2) 
          cat("and so is", model.txt[nr.first.models + 
                                       rank.m3.3], "\n")
      }
      if (ii.m3 > 3) {
        rank.m3.4 <- order(Vaic.m3)[3]
        if (output && Vaic.m3[rank.m3.4] - Vaic.m3[rank.m3] <= 
            2) 
          cat("and so is", model.txt[nr.first.models + 
                                       rank.m3.4], "\n")
      }
      if (output) 
        cat("---------------------------------------------------\n")
    }
  if (quick.ans != 4) {
    ii.m5 <- 0
    m5x <- numeric()
    Vaic.m5 <- numeric()
    nr.first.models.plus3 <- ii
    if (!skip.0) {
      model.ans <- 15 + hill
      ii <- ii + 1
      ii.m5 <- ii.m5 + 1
      if (hill == 0) 
        model.txt[ii] <- "Expon. m5-"
      if (hill == 10) 
        model.txt[ii] <- "Hill m5-"
      fct1 <- rep(1, length(x))
      fct2 <- rep(1, length(x))
      if (max(fct3) > 1) 
        model.txt[ii] <- paste(model.txt[ii], "v", sep = "")
      ans.all.fit$model.ans <- model.ans
      ans.all.fit$modelname <- model.txt[ii]
      ans.all.fit$fct1 <- fct1
      ans.all.fit$fct2 <- fct2
      ans.all.fit$fct3 <- fct3
      ans.all.fit$nr.aa <- max(fct1)
      ans.all.fit$nr.bb <- max(fct2)
      ans.all.fit$nr.var <- max(fct3)
      if (cont) 
        nr.var <- max(fct3)
      ans.all.fit$nr.var <- nr.var
      if (cont) 
        par.start <- c(var.start, aa.start, ced.start, 
                       cc.start, dd.start)
      if (!cont) 
        par.start <- c(aa.start, ced.start, cc.start, 
                       dd.start, ths.start)
      if (dtype == 6) 
        par.start <- c(alfa.start, par.start)
      ans.all.fit$par.start <- par.start
      ans.all.fit <- f.qfit(ans.all.fit, plot.type, output = output)
      Vnpar[ii] <- ans.all.fit$npar.aic
      Vloglik[ii] <- round(ans.all.fit$loglik, 2)
      Vaic[ii] <- ans.all.fit$aic
      Vconverged[ii] <- ans.all.fit$converged
      MLE <- ans.all.fit$MLE
      MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
      Vaic.m5[ii.m5] <- Vaic[ii]
      m5x[ii.m5] <- "0"
      if (output) {}
        #cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
        #          tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
      results.m5.0 <- f.model.specific.results(ans.all.fit)
    }
    if (covar.no > 0) {
      ii <- ii + 1
      ii.m5 <- ii.m5 + 1
      model.ans <- 15 + hill
      if (hill == 0) 
        model.txt[ii] <- "Expon. m5-a"
      if (hill == 10) 
        model.txt[ii] <- "Hill m5-a"
      if (max(fct3) > 1) 
        model.txt[ii] <- paste(model.txt[ii], "v", sep = "")
      ans.all.fit$model.ans <- model.ans
      ans.all.fit$modelname <- model.txt[ii]
      fct1 <- covariate
      fct2 <- rep(1, length(x))
      ans.all.fit$fct1 <- fct1
      ans.all.fit$fct2 <- fct2
      ans.all.fit$fct3 <- fct3
      nr.aa <- max(fct1)
      nr.bb <- max(fct2)
      ans.all.fit$nr.aa <- nr.aa
      ans.all.fit$nr.bb <- nr.bb
      ans.all.fit$nr.var <- nr.var
      if (cont) 
        par.start <- c(var.start, Vaa.start, ced.start, 
                       cc.start, dd.start)
      if (!cont) 
        par.start <- c(Vaa.start, ced.start, cc.start, 
                       dd.start, ths.start)
      if (dtype == 6) 
        par.start <- c(alfa.start, par.start)
      ans.all.fit$par.start <- par.start
      ans.all.fit <- f.qfit(ans.all.fit, plot.type, output = output)
      Vnpar[ii] <- ans.all.fit$npar.aic
      Vloglik[ii] <- round(ans.all.fit$loglik, 2)
      Vaic[ii] <- ans.all.fit$aic
      Vconverged[ii] <- ans.all.fit$converged
      MLE <- ans.all.fit$MLE
      MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
      Vaic.m5[ii.m5] <- Vaic[ii]
      m5x[ii.m5] <- "a"
      if (output) {}
        #cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
        #          tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
      Vaa.start <- MLE[(alfa.length + nr.var + 1):(alfa.length + 
                                                     nr.var + max(fct1))]
      results.m5.a <- f.model.specific.results(ans.all.fit)
      if (!skip.b) {
        ii <- ii + 1
        ii.m5 <- ii.m5 + 1
        model.ans <- 15 + hill
        if (hill == 0) 
          model.txt[ii] <- "Expon. m5-b"
        if (hill == 10) 
          model.txt[ii] <- "Hill m5-b"
        if (max(fct3) > 1) 
          model.txt[ii] <- paste(model.txt[ii], "v", 
                                 sep = "")
        ans.all.fit$modelname <- model.txt[ii]
        fct1 <- rep(1, length(x))
        fct2 <- covariate
        ans.all.fit$fct1 <- fct1
        ans.all.fit$fct2 <- fct2
        ans.all.fit$fct3 <- fct3
        nr.aa <- 1
        nr.bb <- max(fct2)
        ans.all.fit$nr.aa <- nr.aa
        ans.all.fit$nr.bb <- nr.bb
        ans.all.fit$nr.var <- nr.var
        if (cont) 
          par.start <- c(var.start, aa.start.mb, Vced.start, 
                         cc.start, dd.start)
        if (!cont) 
          par.start <- c(aa.start.mb, Vced.start, cc.start, 
                         dd.start, ths.start)
        if (dtype == 6) 
          par.start <- c(alfa.start, par.start)
        ans.all.fit$par.start <- par.start
        ans.all.fit <- f.qfit(ans.all.fit, plot.type, 
                              output = output)
        Vnpar[ii] <- ans.all.fit$npar.aic
        Vloglik[ii] <- round(ans.all.fit$loglik, 2)
        Vaic[ii] <- ans.all.fit$aic
        Vconverged[ii] <- ans.all.fit$converged
        MLE <- ans.all.fit$MLE
        MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
        Vaic.m5[ii.m5] <- Vaic[ii]
        m5x[ii.m5] <- "b"
        if (output) {}
          #cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
          #          tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
        results.m5.b <- f.model.specific.results(ans.all.fit)
      }
      ii <- ii + 1
      ii.m5 <- ii.m5 + 1
      model.ans <- 15 + hill
      if (hill == 0) 
        model.txt[ii] <- "Expon. m5-ab"
      if (hill == 10) 
        model.txt[ii] <- "Hill m5-ab"
      if (max(fct3) > 1) 
        model.txt[ii] <- paste(model.txt[ii], "v", sep = "")
      ans.all.fit$modelname <- model.txt[ii]
      fct1 <- covariate
      fct2 <- covariate
      ans.all.fit$fct1 <- fct1
      ans.all.fit$fct2 <- fct2
      ans.all.fit$fct3 <- fct3
      nr.aa <- max(fct1)
      nr.bb <- max(fct2)
      ans.all.fit$nr.aa <- nr.aa
      ans.all.fit$nr.bb <- nr.bb
      ans.all.fit$nr.var <- nr.var
      if (cont) 
        nr.var <- max(fct3)
      ans.all.fit$nr.var <- nr.var
      if (cont) 
        par.start <- c(var.start, Vaa.start, Vced.start, 
                       cc.start, dd.start)
      if (!cont) 
        par.start <- c(Vaa.start, Vced.start, cc.start, 
                       dd.start, ths.start)
      if (dtype == 6) 
        par.start <- c(alfa.start, par.start)
      ans.all.fit$par.start <- par.start
      ans.all.fit <- f.qfit(ans.all.fit, plot.type, output = output)
      Vnpar[ii] <- ans.all.fit$npar.aic
      Vloglik[ii] <- round(ans.all.fit$loglik, 2)
      Vaic[ii] <- ans.all.fit$aic
      Vconverged[ii] <- ans.all.fit$converged
      MLE <- ans.all.fit$MLE
      MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
      Vaic.m5[ii.m5] <- Vaic[ii]
      m5x[ii.m5] <- "ab"
      if (output) {
        if (vv == 2) {}
          #cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
          #          tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
        else {}#cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
               #        tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
      }
      results.m5.ab <- f.model.specific.results(ans.all.fit)
    }
    rank.m5 <- order(Vaic.m5)[1]
    row.sel.m5 <- nr.first.models.plus3 + rank.m5
    if (quick.ans == 5) 
      if (covar.no > 0) {
        if (output) {
          cat("\n---------------------------------------------------\n")
          cat("Best model 5 with covariates is:", model.txt[nr.first.models + 
                                                              rank.m5], "\n")
        }
        if (ii.m5 > 1) {
          rank.m5.2 <- order(Vaic.m5)[2]
          if (output && Vaic.m5[rank.m5.2] - Vaic.m5[rank.m5] <= 
              2) 
            cat("However", model.txt[nr.first.models + 
                                       rank.m5.2], "is a reasonable model as well\n")
        }
        if (ii.m5 > 2) {
          rank.m5.3 <- order(Vaic.m5)[3]
          if (output && Vaic.m5[rank.m5.3] - Vaic.m5[rank.m5] <= 
              2) 
            cat("and so is", model.txt[nr.first.models + 
                                         rank.m5.3], "\n")
        }
        if (ii.m5 > 3) {
          rank.m5.4 <- order(Vaic.m5)[3]
          if (output && Vaic.m5[rank.m5.4] - Vaic.m5[rank.m5] <= 
              2) 
            cat("and so is", model.txt[nr.first.models. + 
                                         rank.m5.4], "\n")
        }
        if (output) 
          cat("---------------------------------------------------\n")
      }
  }
  if (quick.ans == 4) 
    Vaic.models <- Vaic.m3
  if (quick.ans == 5) 
    Vaic.models <- Vaic.m5
  if (quick.ans == 2) 
    Vaic.models <- c(Vaic.m3, Vaic.m5)
  model.sel <- 1000
  # if (all(Vaic.models - aic.m1 > -2)) {
  #   cat("\nNo significant trend according to AIC criterion for model", 
  #       model.ans, "\n\n")
  #   model.sel <- 1
  #   row.sel <- row.sel.m1
  # }
  if (quick.ans == 2 && covar.no > 0 && model.sel > 1) {
    rank.first <- order(Vaic.models)[1]
    row.sel <- rank.first + nr.first.models
    if (covar.no > 0) {
      if (output) {
        cat("\n---------------------------------------------------\n")
        cat("Best model with covariates is:", model.txt[nr.first.models + 
                                                          rank.first], "\n")
      }
      rank.second <- order(Vaic.models)[2]
      if (output && Vaic.models[rank.second] <= Vaic.models[rank.first] + 
          2) 
        cat("However", model.txt[nr.first.models + rank.second], 
            "is a reasonable model as well\n")
      rank.third <- order(Vaic.models)[3]
      if (output && Vaic.models[rank.third] <= Vaic.models[rank.first] + 
          2) 
        cat("and so is", model.txt[nr.first.models + 
                                     rank.third], "\n")
      if (length(Vaic.models) > 3) {
        rank.fourth <- order(Vaic.models)[4]
        if (output && Vaic.models[rank.fourth] <= Vaic.models[rank.first] + 
            2) 
          cat("and so is", model.txt[nr.first.models + 
                                       rank.fourth], "\n")
      }
      if (length(Vaic.models) > 4) {
        rank.fifth <- order(Vaic.models)[5]
        if (output && Vaic.models[rank.fifth] <= Vaic.models[rank.first] + 
            2) 
          cat("and so is", model.txt[nr.first.models + 
                                       rank.fifth], "\n")
      }
      if (output) 
        cat("---------------------------------------------------\n")
    }
  }
  if (quick.ans == 4) {
    m3x.sel <- m3x[rank.m3]
    row.sel <- row.sel.m3
    if (m3x.sel == "0") 
      ans.all.spec <- results.m3.0
    if (m3x.sel == "a") 
      ans.all.spec <- results.m3.a
    if (m3x.sel == "b") 
      ans.all.spec <- results.m3.b
    if (m3x.sel == "ab") 
      ans.all.spec <- results.m3.ab
  }
  if (quick.ans == 5) {
    m5x.sel <- m5x[rank.m5]
    row.sel <- row.sel.m5
    if (m5x.sel == "0") 
      ans.all.spec <- results.m5.0
    if (m5x.sel == "a") 
      ans.all.spec <- results.m5.a
    if (m5x.sel == "b") 
      ans.all.spec <- results.m5.b
    if (m5x.sel == "ab") 
      ans.all.spec <- results.m5.ab
  }
  aic.tmp <- 0
  if (quick.ans == 2) {
    if (!ans.all$cont) 
      aic.tmp <- 5
    m5x.sel <- m5x[rank.m5]
    m3x.sel <- m3x[rank.m3]
    aic.m5.sel <- Vaic.m5[rank.m5]
    aic.m3.sel <- Vaic.m3[rank.m3]
    if (aic.m3.sel <= aic.m5.sel + aic.tmp) {
      row.sel <- row.sel.m3
      if (m3x.sel == "0") 
        ans.all.spec <- results.m3.0
      if (m3x.sel == "a") 
        ans.all.spec <- results.m3.a
      if (m3x.sel == "b") 
        ans.all.spec <- results.m3.b
      if (m3x.sel == "ab") 
        ans.all.spec <- results.m3.ab
    }
    if (aic.m5.sel + aic.tmp < aic.m3.sel) {
      row.sel <- row.sel.m5
      if (m5x.sel == "0") 
        ans.all.spec <- results.m5.0
      if (m5x.sel == "a") 
        ans.all.spec <- results.m5.a
      if (m5x.sel == "b") 
        ans.all.spec <- results.m5.b
      if (m5x.sel == "ab") 
        ans.all.spec <- results.m5.ab
    }
  }
  if (model.sel == 1) 
    ans.all.spec$trend <- FALSE
  else ans.all.spec$trend <- TRUE
  ans.all.tmp <- f.move.sublist(ans.all, ans.all.spec)
  if (!ans.all$const.var) 
    if (cont && full.ans == 2 && max(fct3) == 1 && covar.no > 
        0) {
      ii <- ii + 1
      model.txt[ii] <- paste(model.txt[row.sel], "v", sep = "")
      ans.all.tmp$fct3 <- covariate
      nr.var <- max(covariate)
      ans.all.tmp$nr.var <- nr.var
      par.start <- c(rep(var.start, nr.var), ans.all.tmp$MLE[-1])
      ans.all.tmp$par.start <- par.start
      if (length(Vloglik) > 9) 
        plot.type.tmp <- 0
      else plot.type.tmp <- plot.type
      ans.all.tmp <- f.qfit(ans.all.tmp, plot.type.tmp, 
                            output = output)
      Vnpar[ii] <- ans.all.tmp$npar
      Vloglik[ii] <- round(ans.all.tmp$loglik, 2)
      Vaic[ii] <- ans.all.tmp$aic
      Vconverged[ii] <- ans.all.tmp$converged
      MLE <- ans.all.tmp$MLE
      MLE.all[ii, 1:length(MLE)] <- signif(MLE, 7)
      if (output) {}
        #cat(paste("\n", model.txt[ii], tb, Vconverged[ii], 
        #          tb, Vnpar[ii], tb, Vloglik[ii]), tb, Vaic[ii])
      if (Vaic[ii] >= Vaic[row.sel]) {
        fct3 <- rep(1, length(x))
        nr.var <- 1
      }
      else {
        ans.all.spec <- f.model.specific.results(ans.all.tmp)
        row.sel <- ii
        vv <- 2
        var.start <- MLE[1:nr.var]
        if (model.sel == 1) {
          cat("ATTENTION: Without covariate-dependent variances no trend was assessed (based on AIC).\n                     You might consider to repeat the analysis including the full model to check if this changes\n\n")
        }
      }
    }
  if (output && model.sel != 1) {
    #cat("\n------------------------------------------------------\n")
    #cat("selected model:", model.txt[row.sel], "\n")
    if (aic.tmp != 0) 
      cat("using crit diff in AIC of", aic.tmp, "\n")
  }
  if (full.ans == 1) {
    if (dtype %in% c(5, 15)) 
      cat("\n\nnote that the AIC of a fitted model for continuous data with litter effects \n  cannot be compare to the full model\n\n")
    else if (!ans.all$do.MA && full.ans == 1 && model.sel > 
             1) {
      aic.full <- ans.all$FULL.model$aic
      npar.aic.full <- ans.all$FULL.model$npar.aic
      if (npar.aic.full >= Vnpar[row.sel]) {
        if (output && (Vaic[row.sel] > aic.full + 2)) {
          alert.full <- f.alert.full()
        }
      }
    }
  }
  ans.all.tmp <- f.move.sublist(ans.all, ans.all.spec)
  if (output) {
    #cat(ans.all.tmp$report.pars)
    #cat("------------------------------------------------------\n")
  }
  #if (!is.na(alert.full)) 
  #  cat(alert.full)
  if (dtype %in% c(1, 5)) {
    ans.all.tmp <- f.resid.con(ans.all.tmp)
    regr.resid <- ans.all.tmp$regr.resid
    sig2 <- mean(MLE[1:nr.var])
    mu <- mean(regr.resid)
    length.old <- length(regr.resid)
    limit.tmp <- f.grubb(length.old)
    limit.low <- -limit.tmp * sqrt(sig2) + mu
    limit.up <- limit.tmp * sqrt(sig2) + mu
    regr.resid.min <- regr.resid[(regr.resid > limit.low) & 
                                   (regr.resid < limit.up)]
    length.new <- length(regr.resid.min)
    if (length.old > length.new) {}
      #cat("\nATTENTION: ", length.old - length.new, "outliers detected\n")
  }
  if (dtype %in% c(5, 15)) {
    if (max(ans.all.tmp$Vdetlim, na.rm = T) != 0) 
      cat("\nlitter effects not implemented for nonzero detection limit\n")
    if (dtype == 5 && ans.all.tmp$nest.no == 0) 
      dtype <- 1
    else {
      ans.all.tmp <- f.nested.con(ans.all.tmp)
      ans.all.tmp$inter.var <- ans.all.tmp$inter.var
      ans.all.tmp$intra.var <- ans.all.tmp$intra.var
    }
  }
  if (ans.all.tmp$NES.ans == 2 && model.sel == 1) {
    ans.all.tmp <- f.nes(ans.all.tmp)
  }
  ans.all.tmp$factor.name <- f.factorname(ans.all.tmp)
  ans.all.tmp$model.txt <- model.txt
  ans.all.tmp$Vloglik <- Vloglik
  ans.all.tmp$Vnpar <- Vnpar
  ans.all.tmp$Vaic <- Vaic
  ans.all.tmp$Vconverged <- Vconverged
  ans.all.tmp$MLE.all <- MLE.all
  ans.all.tmp$regr.par.matr <- f.pars(ans.all.tmp)$regr.par.matr
  ans.all.tmp$row.sel <- row.sel
  ans.all.tmp$conf.int <- matrix(rep(NA, 2), ncol = 2)
  if (model.sel != 1) {
    if (!ans.all$no.CI && dtype == 3) 
      ans.all.tmp <- f.CI.sel.ord(ans.all.tmp)
    else if (ans.all$NES.ans == 2) {
      ans.all.tmp <- f.refit.nes(ans.all.tmp)
      if (output) {
        #cat(ans.all.tmp$report.pars)
        #cat("------------------------------------------------------\n")
      }
    }
    else if (!ans.all$no.CI && cont) {
      ans.all.tmp <- f.CI.sel(ans.all.tmp)
    }
    else if (!cont) 
      ans.all.tmp <- f.CI.sel(ans.all.tmp)
  }
  if (model.sel == 1) {
    ans.all.tmp$CED.matr <- matrix(rep(NA, 2), ncol = 2)
    ans.all.tmp$rank.low <- 1
  }
  loglik.summary <- data.frame(model = ans.all.tmp$model.txt, 
                               converged = ans.all.tmp$Vconverged, loglik = ans.all.tmp$Vloglik, 
                               npar = ans.all.tmp$Vnpar, AIC = ans.all.tmp$Vaic)
  ans.all.tmp$loglik.summary <- loglik.summary
  if (hill == 0) {
    ans.all[["EXP"]] <- f.model.specific.results(ans.all.tmp)
  }
  if (hill == 10) {
    ans.all[["HILL"]] <- f.model.specific.results(ans.all.tmp)
  }
  if (hill == 0) 
    ans.all$logliks <- ans.all$EXP$loglik.summary
  if (hill == 10) 
    ans.all$logliks <- rbind(ans.all$EXP$loglik.summary, 
                             ans.all$HILL$loglik.summary)
  ans.all$npar.null <- ans.all$NULL.model$npar
  ans.all$npar.aic.null <- ans.all$NULL.model$npar.aic
  ans.all$loglik.null <- ans.all$NULL.model$loglik
  ans.all$aic.null <- ans.all$NULL.model$aic
  if (length(ans.all$FULL.model) > 0) {
    ans.all$npar.full <- ans.all$FULL.model$npar
    ans.all$npar.aic.full <- ans.all$FULL.model$npar.aic
    ans.all$loglik.full <- ans.all$FULL.model$loglik
    ans.all$aic.full <- ans.all$FULL.model$aic
  }
  ans.all$fitted <- T
  if (hill == 0) {
    ans.all$row.sel.m1 <- row.sel.m1
    ans.all$aic.m1 <- aic.m1
    if (covar.no > 0) {
      ans.all$vv <- vv
    }
    if (cont) 
      ans.all$var.start <- var.start
    if (!cont) 
      ans.all$ths.start <- ths.start
    ans.all$skip.b <- skip.b
    ans.all$skip.0 <- skip.0
    ans.all$nr.first.models <- nr.first.models
    ans.all$model.type <- 2
    if (cont && !(model.sel == 1)) 
      if (max(ans.all$fct3) != ans.all$nr.var) {
        print("error in f.select.m5.con")
        print(ans.all$fct3)
        print(ans.all$nr.var)
        f.press.key.to.continue()
      }
  }
  if (0) 
    if (!(model.sel == 1 & hill == 10)) {
      if (max(fct1) > 1) 
        ans.all$fct1.no <- covar.no
      if (max(fct2) > 1) 
        ans.all$fct2.no <- covar.no
      if (max(fct3) > 1) 
        ans.all$fct3.no <- covar.no
    }
  if (0) 
    if (model.sel != 1 && hill == 10 && !ans.all$do.MA) {
      ans.all$loglik <- NULL
      ans.all$aic <- NULL
      ans.all$MLE <- NULL
      ans.all$MLE.all <- NULL
      ans.all$regr.par <- NULL
      ans.all$Vlogliks <- NULL
      ans.all$Vaic <- NULL
      ans.all$fct1 <- NULL
      ans.all$fct2 <- NULL
      if (!(dtype %in% c(2, 4, 6))) {
        ans.all$fct3 <- NULL
        ans.all$nr.var <- NULL
      }
      else {
        ans.all$fct3 <- 0
        ans.all$nr.var <- 0
      }
      ans.all$nr.aa <- NULL
      ans.all$nr.bb <- NULL
      ans.all$nrp <- NULL
      ans.all$Vnpar <- NULL
      ans.all$npar <- NULL
    }
  if (ans.all$no.hill) 
    ans.all$model.fam <- 1
  else ans.all$model.fam <- 0
  ans.all$CED.model <- TRUE
  if (dtype %in% c(5, 15)) {
    ans.all$inter.var <- ans.all.tmp$inter.var
    ans.all$intra.var <- ans.all.tmp$intra.var
  }
  ans.all$alert.full <- alert.full
  if (exists("track")) 
    print("f.select.m5.con:   END")
  return(ans.all)
}

####################
## f.graph.window ##
####################
f.graph.window <- function (nr.pl = 1, nr.gr = 1, WAPP = FALSE, title = "", name.wapp = NA, 
                            plotprefix = NA, quick.ans, svg.plots = FALSE, output = TRUE) 
{
  if (exists("track")) 
    print("f.graph.window")
  if (nr.pl > 64) 
    cat("\nATTENTION:  nr of plots too large to plot together\n")
  if (!exists(".gw.size")) 
    .gw.size <- c(5, 5)
  if (WAPP && !svg.plots) 
    .gw.size <- c(400, 400)
  else if (WAPP && svg.plots) 
    .gw.size <- c(6, 6)
  f.create.graphwin(.gw.size[1], .gw.size[2], nr.gr = nr.gr, 
                    WAPP = WAPP, title = title, name.wapp = name.wapp, plotprefix = plotprefix, 
                    svg.plots = svg.plots, output = output)
  if (nr.pl == 1) {
    par(mfcol = c(1, 1))
    par(mex = 0.7)
    par(mar = c(6, 5, 6, 13))
  }
  if (nr.pl %in% 2) {
    par(mfcol = c(2, 1))
    par(mar = c(5.5, 4.8, 4, 12))
    par(cex.main = 1.2)
    par(cex.sub = 1)
    par(cex.lab = 1)
    par(cex = 0.6)
  }
  if (nr.pl %in% 3:4) {
    par(mfcol = c(2, 2))
    par(mar = c(3.2, 3.7, 2, 1.5))
    par(cex.main = 1.2)
    par(cex.sub = 0.9)
    par(cex.lab = 0.9)
    par(cex = 0.9)
  }
  if (nr.pl %in% 5:6) {
    par(mfcol = c(3, 2))
    par(mar = c(5, 3.5, 2.5, 2.5))
    par(cex.main = 1.2)
    par(cex.sub = 0.6)
    par(cex.lab = 1.1)
    par(cex = 0.6)
  }
  if (nr.pl %in% 7:9) {
    par(mfcol = c(3, 3))
    par(mar = c(2.5, 3, 2.5, 0.5))
    par(cex.main = 1.2)
    par(cex.sub = 0.6)
    par(cex.lab = 1.3)
    par(cex = 0.6)
  }
  if (nr.pl %in% 10:12) {
    par(mfcol = c(4, 3))
    par(mar = c(2.5, 2.5, 1.5, 1.5))
    par(cex.main = 1.2)
    par(cex.sub = 0.6)
    par(cex.lab = 0.6)
    par(cex = 0.6)
  }
  if (nr.pl %in% 13:16) {
    par(mfcol = c(4, 4))
    par(mar = c(2, 2, 1.2, 0.5))
  }
  if (nr.pl %in% 17:20) {
    par(mfcol = c(4, 5))
    par(mar = c(3.3, 2.5, 2, 1.5))
  }
  if (nr.pl %in% 21:25) {
    par(mfcol = c(5, 5))
    par(mar = c(3, 2, 1.8, 1.2))
  }
  if (nr.pl %in% 26:30) {
    par(mfcol = c(6, 5))
    par(mar = c(2, 1.5, 1.5, 1))
  }
  if (nr.pl %in% 31:36) {
    par(mfcol = c(6, 6))
    par(mar = c(2, 1.5, 1.5, 1))
  }
  if (nr.pl %in% 37:42) {
    par(mfcol = c(7, 6))
    par(mar = c(2, 1.5, 1.5, 1))
  }
  if (nr.pl %in% 43:49) {
    par(mfcol = c(7, 7))
    par(mar = c(2, 1.5, 1.5, 1))
  }
  if (nr.pl %in% 50:56) {
    par(mfcol = c(8, 7))
    par(mar = c(2, 1.5, 1.5, 1))
  }
  if (nr.pl %in% 57:64) {
    par(mfcol = c(8, 8))
    par(mar = c(2, 1.5, 1.5, 1))
  }
  if (exists("track")) 
    print("f.graph.window   END")
  return(invisible())
}

#######################
## f.create.graphwin ##
#######################
f.create.graphwin <- function (aa, bb, title = "", name.wapp = NA, WAPP = FALSE, plotprefix = NA, 
                               nr.gr = 1, svg.plots = FALSE, output = TRUE) 
{
  if (exists("track")) 
    print("f.create.graphwin")
  lst <- dev.list()
  if (length(lst) > 60) {
    if (nr.gr > 25) 
      graphics.off()
    else {
      cat("\n ATTENTION:  too many graphical windows are opened")
      cat("\nthey will be deleted, this is the last chance to save any of them")
      if (output) 
        f.press.key.to.continue("\npress any key if they can be deleted ... > ")
      graphics.off()
      assign(".ypos", 0)
    }
  }
  if (WAPP) {
    if (svg.plots) 
      svg(paste(plotprefix, name.wapp, ".svg", sep = ""), 
          aa, bb)
    else png(paste(plotprefix, name.wapp, ".png", sep = ""), 
             aa, bb)
  }
  else {
    title <- paste(dev.cur(), "PROAST:", title)
    y.pos <- .ypos + 5
    if (y.pos > 100) 
      y.pos <- 0
    assign(".ypos", y.pos, immediate = T, pos = 1)
    if (.Platform$OS.type == "windows") 
      windows(aa, bb, 10, ypos = y.pos, restoreConsole = T, 
              title = title)
    else X11(, aa, bb, 10, title = title, ypos = y.pos)
  }
  if (exists("track")) 
    print("f.create.graphwin:   END")
  return(invisible())
}

############
## f.qfit ##
############
f.qfit <- function (ans.all.fit, plot.type, plt.mns = 3, output = TRUE) 
{
  if (exists("track")) 
    print("f.qfit")
  plot.type.save <- ans.all.fit$plot.type
  ans.all.fit$plot.type <- plot.type
  ans.all.fit$plt.mns <- plt.mns
  with(ans.all.fit, {
    if (cont) {
      increase.00 <- ans.all.fit$increase
      ans.all.fit$text.par <- f.text.par(ans.all.fit)
      par.start.tmp <- ans.all.fit$par.start
      if (length(xans) > 1) 
        par.start.tmp <- NA
      if (any(is.na(par.start.tmp))) {
        ans.all.fit <- f.start.con(ans.all.fit, adjust = F, 
                                   fitted = F, tmp.quick = T)
        if (ans.all.fit$adjust.start) {
          cat("\nadjust start values")
          ans.all.fit <- f.start.con(ans.all.fit, adjust = T, 
                                     fitted = F, tmp.quick = T)
          par.start <- ans.all.fit$par.start
        }
        if (ans.all.fit$par.start[1] == 0) {
          ans.all.fit <- f.start.con(ans.all.fit, adjust = T)
        }
      }
      else {
        ans.all.fit$regr.par <- par.start.tmp[-(1:max(fct3))]
        ans.all.fit$par.start <- par.start.tmp
        ans.all.fit$npar <- length(par.start.tmp)
      }
      ans.all.fit <- f.constr.con(ans.all.fit, tmp.quick = T)
      ans.all.fit$increase <- increase.00
      ans.all.fit <- f.nlminb(ans.all.fit, tmp.quick = T)
      if (model.ans == -25 || model.ans == -13) {
        print("f.qfit")
        print(ans.all.fit$par.start)
        print(ans.all.fit$MLE)
        print(ans.all.fit$loglik)
        print(ans.all.fit$CES)
        print(ans.all.fit$sf.x)
      }
      if (ans.all.fit$MLE[1] == 0) 
        return(ans.all.fit)
      ans.all.fit <- f.pars(ans.all.fit)
      if (output) 
        ans.all.fit$report.pars <- f.report.pars(ans.all.fit)
      MLE <- ans.all.fit$MLE
      ans.all.fit$regr.par <- MLE[(max(fct3) + 1):length(MLE)]
      if (interrupt && plot.type > 0 && !WAPP) {
        ans.all.fit$heading <- ""
        ans.all.fit$xy.lim[2:3] <- c(min(x), max(x))
        ignore <- f.plot.con(ans.all.fit)
        modelname <- ans.all.fit$modelname
        title(main = paste(modelname))
        ans.all.fit$l.ty = 1
        if (model.ans != 11) 
          f.lines.con(ans.all.fit)
      }
    }
    if (!cont) {
      if (model.type == 1) {
        if (model.ans != 14) {
          cat("\n error in f.qfit: model.type and model.ans do not match\n")
          stop()
        }
        ans.all.fit$nr.aa <- max(fct1)
        ans.all.fit$nr.bb <- max(fct2)
        ans.all.fit <- f.start.bin(ans.all.fit, tmp.quick = T)
        ans.all.fit <- f.nlminb(ans.all.fit, tmp.quick = T)
        if (dtype == 6) {
          ans.all.fit$pi.full <- ans.all.fit$MLE[-1]
          ans.all.fit$alfa.start <- ans.all.fit$MLE[1]
        }
        if (interrupt && plot.type > 0 && hill == 0) {
          ans.all.fit$heading <- ""
          ans.all.fit$xy.lim[2:3] <- c(min(x), max(x))
          ans.all.fit$shift.tmp <- f.plot.frq(ans.all.fit)$shift.tmp
          if (is.R()) 
            title(main = paste(modelname))
          else title(main = paste(modelname), cex = 0.5)
        }
        ans.all.fit$npar <- length(ans.all.fit$MLE)
      }
      if (model.type == 2) {
        par.start.tmp <- ans.all.fit$par.start
        if (any(is.na(par.start.tmp)) || length(xans) > 
            1) {
          ans.all.fit <- f.start.cat(ans.all.fit, tmp.quick = T)
          nrp <- ans.all.fit$nrp
        }
        else {
          nrp <- length(par.start.tmp) - nth - 1
          if (dtype == 6) 
            nrp <- nrp - 1
          loglik.start <- -f.lik.cat(par.start.tmp, x, 
                                     y, kk, nn, dtype, fct1, fct2, nrp, nth, nr.aa, 
                                     nr.bb, model.ans, model.type, CES.cat = CES.cat, 
                                     ttt = ttt, twice = twice, CES = CES, decr.zz = decr.zz, 
                                     alfa.length = alfa.length, cens.up = cens.up, 
                                     fct3 = fct3, ces.ans = ces.ans, kk.tot = kk.tot, 
                                     nn.tot = nn.tot, quick.ans = quick.ans)
          if (is.na(loglik.start) | loglik.start == -1e+12) {
            ans.all.fit <- f.start.cat(ans.all.fit, tmp.quick = T)
          }
          else {
            ans.all.fit$par.start <- par.start.tmp
            ans.all.fit$npar <- length(par.start.tmp)
          }
        }
        ans.all.fit$nrp <- nrp
        ans.all.fit <- f.constr.con(ans.all.fit, tmp.quick = T)
        ans.all.fit <- f.nlminb(ans.all.fit, tmp.quick = T)
        par.lst <- f.split.par(ans.all.fit$MLE, ans.all.fit$nrp, 
                               ans.all.fit$nth, dtype)
        regr.par <- par.lst$regr.par
        th.par <- par.lst$th.par
        sig.par <- par.lst$sig.par
        ans.all.fit <- f.pars(ans.all.fit)
        ans.all.fit$report.pars <- f.report.pars(ans.all.fit)
        if (interrupt && plot.type > 0 && hill == 0) {
          print("f.qfit 000000")
          print("Attention to programmer: this part of the code is used, which is unexpected")
          f.press.key.to.continue()
          ans.all.fit$heading <- ""
          ans.all.fit$xy.lim[2:3] <- c(min(x), max(x))
          ans.all.fit$regr.par <- regr.par
          ans.all.fit$th.par <- th.par
          ans.all.fit$sig.par <- sig.par
          if (dtype %in% 2:3) {
            ans.all.fit$plot.type <- plot.type.save
            f.lines.cat(ans.all.fit)
            ans.all.fit$regr.par <- regr.par
          }
          else {
            kk.dum <- 1
            f.plot.frq(ans.all.fit)
            f.lines.frq(ans.all.fit)
            ans.all.fit$regr.par <- regr.par
          }
          title(main = paste(modelname))
        }
      }
    }
    if (!is.finite(ans.all.fit$loglik)) 
      ans.all.fit$loglik <- -1e+12
    ans.all.fit$plot.type <- plot.type.save
    ans.all.fit$fitted <- TRUE
    if (exists("track")) 
      print("f.qfit:  END")
    return(ans.all.fit)
  })
}

################
## f.text.par ##
################
f.text.par <- function (ans.all, brief = F) 
{
  if (exists("track")) 
    print("f.text.par")
  with(ans.all, {
    if (brief && model.ans != 46) {
      nr.var <- 1
      nr.aa <- 1
      nr.bb <- 1
      nr.cc <- 1
      nr.dd <- 1
    }
    else {
      nr.var <- max(fct3)
      nr.aa <- max(fct1)
      nr.bb <- max(fct2)
      nr.cc <- max(fct4)
      nr.dd <- max(fct5)
    }
    if (!cont) 
      nr.var <- 0
    if (nr.aa == 1) {
      fct1.txt <- NULL
      sep.a <- NULL
    }
    if (nr.bb == 1 && length(xans) == 1 && model.ans != 46) {
      fct2.txt <- NULL
    }
    if (nr.var == 1) {
      fct3.txt <- NULL
    }
    if (nr.cc == 1) {
      fct4.txt <- NULL
    }
    if (nr.dd == 1) {
      fct5.txt <- NULL
    }
    if (!cont && model.type == 1) {
      switch(model.ans, text.par <- paste("a", fct1.txt[1:nr.aa], 
                                          sep = "-"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                          sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-")), 
             text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                           paste("b", fct2.txt[1:nr.bb], sep = "-"), "c"), 
             text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                           paste("b", fct2.txt[1:nr.bb], sep = "-"), "c", 
                           "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                     sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                                               "c"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                         sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                                                                   "c"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                             sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                                                                                       "c", "dd"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                       sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                                                                                                                 "c"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                           sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                                                                                                                                     "c"), text.par <- c(rep("aa", nr.aa), rep("bb", 
                                                                                                                                                                               nr.bb)), text.par <- c(rep("aa", nr.aa), rep("bb", 
                                                                                                                                                                                                                            nr.bb), "cc"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                               sep = "-"), paste("bb", fct2.txt[1:nr.bb], 
                                                                                                                                                                                                                                                                                 sep = "-")), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                                                                                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-")), 
             {
               if (length(ans.all$x.full) == 0) x.full <- x
               if (dtype %in% c(4, 6, 84)) text.par <- paste("group", 
                                                             1:length(x.full))
             }, text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                    sep = "-"), paste("BMD", fct2.txt[1:nr.bb], 
                                                      sep = "-")), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                       sep = "-"), paste("BMD", fct2.txt[1:nr.bb], 
                                                                                                         sep = "-"), "c"), text.par <- c(paste("a", 
                                                                                                                                               fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                                                                                                                                                                                    fct2.txt[1:nr.bb], sep = "-"), "c", "d"), text.par <- c(paste("a", 
                                                                                                                                                                                                                                                  fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                                                                                                                                                                                                                                                                                       fct2.txt[1:nr.bb], sep = "-"), "c"), text.par <- c(paste("a", 
                                                                                                                                                                                                                                                                                                                                                fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                                                                                                                                                                                                                                                                                                                                                                                     fct2.txt[1:nr.bb], sep = "-"), "c"), text.par <- c(paste("a", 
                                                                                                                                                                                                                                                                                                                                                                                                                                              fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   fct2.txt[1:nr.bb], sep = "-"), "c", "dd"), 
             text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                           paste("BMD", fct2.txt[1:nr.bb], sep = "-"), 
                           "c"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                     sep = "-"), paste("BMD", fct2.txt[1:nr.bb], 
                                                                       sep = "-"), "c"), text.par <- c(paste("a", 
                                                                                                             fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                                                                                                                                                  fct2.txt[1:nr.bb], sep = "-")), text.par <- c(paste("a", 
                                                                                                                                                                                                      fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                                                                                                                                                                                                                                           fct2.txt[1:nr.bb], sep = "-"), "cc"), {
                                                                                                                                                                                                                                             if (ces.ans == 1) text.par <- c(paste("a", 
                                                                                                                                                                                                                                                                                   fct1.txt[1:nr.aa], sep = "-"), paste("ED50", 
                                                                                                                                                                                                                                                                                                                        fct2.txt[1:nr.bb], sep = "-"))
                                                                                                                                                                                                                                             if (ces.ans %in% 2:3) text.par <- c(paste("a", 
                                                                                                                                                                                                                                                                                       fct1.txt[1:nr.aa], sep = "-"), paste("BMD", 
                                                                                                                                                                                                                                                                                                                            fct2.txt[1:nr.bb], sep = "-"))
                                                                                                                                                                                                                                           }, text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                                  sep = "-"), paste("BMD", fct2.txt[1:nr.bb], 
                                                                                                                                                                                                                                                                                    sep = "-")), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                                                                                     sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                                                                                                                                                                                                                                                                                                               "c", "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                                                                                                              sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                                                                                                                                                                                                                                                                                                                                        "c", "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                                                                                                                                       sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                                                                                                                                                                                                                                                                                                                                                                 "c", "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                                                                                                                                                                sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                                                                                                                                                                                                                                                                                                                                                                                          "BMDratio"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                                                                                                                                                                                           sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                                                                                                                                                                                                                                                                                                                                                                                                                     "c", "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                                                                                                                                                                                                                    sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                                                                                                                                                                                                                                                                                                                                                                                                                                              "c"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                  "c"))
      if (!(model.ans == 14)) 
        if (length(xans) > 1) {
          if (length(xans) == 2) 
            text.par <- c(text.par, paste("RPF", fct2.txt[2], 
                                          sep = "-"))
          if (length(xans) == 3) 
            text.par <- c(text.par, paste("RPF", fct2.txt[2:3], 
                                          sep = "-"))
        }
    }
    if (cont || model.type == 2) {
      switch(model.ans, text.par <- paste("a", fct1.txt[1:nr.aa], 
                                          sep = "-"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                          sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-")), 
             text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                           paste("b", fct2.txt[1:nr.bb], sep = "-"), paste("d", 
                                                                           fct5.txt[1:nr.dd], sep = "-")), text.par <- c(paste("a", 
                                                                                                                               fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                                                                                                                                                                    sep = "-"), paste("c", fct4.txt[1:nr.cc], sep = "-")), 
             text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                           paste("b", fct2.txt[1:nr.bb], sep = "-"), paste("c", 
                                                                           fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                                                                                                                fct5.txt[1:nr.dd], sep = "-")), text.par <- c(paste("a", 
                                                                                                                                                                    fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                                                                                                                                                                                                         sep = "-"), paste("c", fct4.txt[1:nr.cc], sep = "-"), 
                                                                                                                                                              paste("d", fct5.txt[1:nr.dd], sep = "-")), 
             text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                           paste("b", fct2.txt[1:nr.bb], sep = "-")), 
             text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                           paste("b", fct2.txt[1:nr.bb], sep = "-"), paste("d", 
                                                                           fct5.txt[1:nr.dd], sep = "-")), text.par <- c(paste("a", 
                                                                                                                               fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                                                                                                                                                                    sep = "-"), paste("c", fct4.txt[1:nr.cc], sep = "-")), 
             text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                           paste("b", fct2.txt[1:nr.bb], sep = "-"), paste("c", 
                                                                           fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                                                                                                                fct5.txt[1:nr.dd], sep = "-")), text.par <- rep("GM", 
                                                                                                                                                                length(regr.par)), text.par <- c(paste("a", 
                                                                                                                                                                                                       fct1.txt[1:nr.aa], sep = "-"), paste("CED", 
                                                                                                                                                                                                                                            fct2.txt[1:nr.bb], sep = "-")), text.par <- c(paste("a", 
                                                                                                                                                                                                                                                                                                fct1.txt[1:nr.aa], sep = "-"), paste("CED", 
                                                                                                                                                                                                                                                                                                                                     fct2.txt[1:nr.bb], sep = "-"), paste("d", fct5.txt[1:nr.dd], 
                                                                                                                                                                                                                                                                                                                                                                          sep = "-")), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                                                                                                                                                                           sep = "-"), paste("CED", fct2.txt[1:nr.bb], 
                                                                                                                                                                                                                                                                                                                                                                                                                             sep = "-"), paste("c", fct4.txt[1:nr.cc], sep = "-")), 
             text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                           paste("CED", fct2.txt[1:nr.bb], sep = "-"), 
                           paste("c", fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                                                                           fct5.txt[1:nr.dd], sep = "-")), text.par <- c(paste("a", 
                                                                                                                               fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                                                                                                                                                                    sep = "-"), "c", "BMDratio"), text.par <- c(paste("a", 
                                                                                                                                                                                                                      fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                                                                                                                                                                                                                                                           sep = "-")), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                                                            sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                                                                                                                                                                                                                                                                                      paste("d", fct5.txt[1:nr.dd], sep = "-")), 
             text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                           paste("b", fct2.txt[1:nr.bb], sep = "-"), paste("c", 
                                                                           fct4.txt[1:nr.cc], sep = "-")), text.par <- c(paste("a", 
                                                                                                                               fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                                                                                                                                                                    sep = "-"), paste("c", fct4.txt[1:nr.cc], sep = "-"), 
                                                                                                                         paste("d", fct5.txt[1:nr.dd], sep = "-")), 
             text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                           paste("bb1", fct2.txt[1:nr.bb], sep = "-"), 
                           "bb2", "cc1", "cc2", "dd"), text.par <- c(paste("a", 
                                                                           fct1.txt[1:nr.aa], sep = "-"), paste("CED", 
                                                                                                                fct2.txt[1:nr.bb], sep = "-")), text.par <- c(paste("a", 
                                                                                                                                                                    fct1.txt[1:nr.aa], sep = "-"), paste("CED", 
                                                                                                                                                                                                         fct2.txt[1:nr.bb], sep = "-"), paste("d", fct5.txt[1:nr.dd], 
                                                                                                                                                                                                                                              sep = "-")), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                                               sep = "-"), paste("CED", fct2.txt[1:nr.bb], 
                                                                                                                                                                                                                                                                                                 sep = "-"), paste("c", fct4.txt[1:nr.cc], sep = "-")), 
             text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                           paste("CED", fct2.txt[1:nr.bb], sep = "-"), 
                           paste("c", fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                                                                           fct5.txt[1:nr.dd], sep = "-")), text.par <- c(paste("a", 
                                                                                                                               fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                                                                                                                                                                    sep = "-"), "c"), text.par <- c(paste("a", 
                                                                                                                                                                                                          fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                                                                                                                                                                                                                                               sep = "-")), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                                                sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                                                                                                                                                                                                                                                                          "c", "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                                                                         sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                                                                                                                                                                                                                                                                                                   "c", "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                                                                                                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                                                                                                                                                                                                                                                                                                                            "c"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                                                                                                                      sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                                                                                                                                                                                                                                                                                                                                                "c", "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                                                                                                                                               sep = "-"), paste("bb1", fct2.txt[1:nr.bb], 
                                                                                                                                                                                                                                                                                                                                                                                                 sep = "-"), "bb2", "cc", "dd1", "dd2"), text.par <- c(paste("aa", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                             fct1.txt[1:nr.aa], sep = "-"), paste("bb1", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  fct2.txt[1:nr.bb], sep = "-"), "bb2", "cc1", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                       "cc2", "dd1", "dd2"), text.par <- c(paste("a", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      sep = "-"), "cc", "dd"), text.par <- c(paste("a", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   fct1.txt[1:nr.aa], sep = "-"), paste("b", fct2.txt[1:nr.bb], 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        sep = "-")), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                   "c", "d"), text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            "c", "d"), , , , text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-")), 
             , text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                   sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                             "c", paste("d", fct5.txt[1:nr.dd], sep = "-")), 
             text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                           paste("b", fct2.txt[1:nr.bb], sep = "-")), 
             text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                           paste("b", fct2.txt[1:nr.bb], sep = "-"), "c1", 
                           "c2", "d"), {
                             nr.bb <- length(fct2.txt)
                             ced.txt <- paste("RPF", fct2.txt[1:nr.bb], 
                                              sep = "-")
                             ced.txt[ref.lev] <- paste("CED", fct2.txt[ref.lev], 
                                                       sep = "-")
                             text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                 sep = "-"), ced.txt, paste("c", fct4.txt[1:nr.cc]), 
                                           paste("d", fct5.txt[1:nr.dd], sep = "-"))
                           }, text.par <- c(paste("a", fct1.txt[1:nr.aa], 
                                                  sep = "-"), paste("b", fct2.txt[1:nr.bb], sep = "-"), 
                                            "c", paste("d", fct5.txt[1:nr.dd], sep = "-")))
      if (model.ans == 6 && ans.m6.sd == 2) 
        text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                      "b", "q", "d ")
      if (model.ans == 47) 
        text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                      "GCED", "q", "d ")
      if (model.ans %in% c(48, 50)) 
        text.par <- c(paste("a", fct1.txt[1:nr.aa], sep = "-"), 
                      paste("b", fct2.txt[1:nr.bb], sep = "-"), paste("c", 
                                                                      fct4.txt[1:nr.cc], sep = "-"), paste("d", 
                                                                                                           fct5.txt[1:nr.dd], sep = "-"), "gg", "hh")
      if (cont && fit.ans == 1) 
        text.par <- c(paste("var", fct3.txt[1:nr.var], 
                            sep = "-"), text.par)
      if (model.ans != 11) 
        if (length(xans) > 1) {
          if (length(xans) == 2) 
            text.par <- c(text.par, paste("RPF", fct2.txt[2], 
                                          sep = "-"))
          if (length(xans) == 3) 
            text.par <- c(text.par, paste("RPF", fct2.txt[2:3], 
                                          sep = "-"))
        }
      if (!cont && ces.ans != 4) {
        if (model.ans %in% c(12:15, 22:25) && nr.aa > 
            1 && nr.bb == 1) 
          for (ii in 2:(nr.aa + 1)) text.par[ii] <- paste("CED", 
                                                          ii - 1, sep = "")
      }
      if (!cont) {
        if (dtype == 4) 
          text.th <- "th"
        else text.th <- paste("th-", 1:nth, sep = "")
        text.par <- c(text.par, text.th, "sigma")
      }
    }
    if (dtype == 6) 
      text.par <- c("alfa", text.par)
    if (exists("track")) 
      print("f.text.par:  END")
    return(text.par)
  })
}

#################
## f.start.con ##
#################
f.start.con <- function (ans.all, adjust = F, fitted = F, tmp.quick = F) 
{
  if (exists("track")) 
    print("f.start.con")
  if (ans.all$model.ans %in% 38:39) {
    cat("\nmodel.ans is 38 or 39 in f.start.con\n")
    f.press.key.to.continue
    return(ans.all)
  }
  if (ans.all$model.ans == 0) {
    cat("choose model first")
    f.press.key.to.continue()
    return(ans.all)
  }
  if (ans.all$cont && ans.all$model.ans == 10 & ans.all$list.logic) {
    ans.dum <- menu(c("yes", "no"), title = "\ndo you want to use previous fitted model 5 for start values\n?")
    if (ans.dum == 1) {
      model.cur <- as.numeric(readline(prompt = "\ngive model number of previous fit >"))
      MLE <- ans.all$MLE
      nr.var <- ans.all$nr.var
      nr.aa <- ans.all$nr.aa
      nr.bb <- ans.all$nr.bb
      par.start <- MLE
      cc <- MLE[length(MLE) - 1]
      dd <- MLE[length(MLE)]
      tmp.bb <- MLE[(nr.var + nr.aa + 1):(nr.var + nr.aa + 
                                            nr.bb)]
      if (model.cur == 15) {
        tmp.bb <- f.bb.con(15, cc, dd, tmp.bb, ans.all$CES)
      }
      bb.new <- (1/tmp.bb)^(1/dd)
      par.start[(nr.var + nr.aa + 1):(nr.var + nr.aa + 
                                        nr.bb)] <- bb.new
      ans.all$par.start <- par.start
      ans.all$regr.par <- ans.all$par.start[-(1:nr.var)]
      ignore <- f.plot.con(ans.all)
      ans.all$nr.gr <- f.lines.con(ans.all)
      with(ans.all, {
        loglik.first <- -f.lik.con(par.start, x, y, dtype, 
                                   fct1, fct2, fct3, model.ans, mn.log, sd2.log, 
                                   nn, Vdetlim = Vdetlim, CES = CES, twice = twice, 
                                   ttt = ttt, trace = F, fct4 = fct4, fct5 = fct5, 
                                   cens.up = cens.up, par.tmp = NA, increase = increase, 
                                   ref.lev = ref.lev, sign.q = sign.q, ans.m6.sd = ans.m6.sd)
        cat("\n\nlog-likelihood at start values: ", round(loglik.first, 
                                                          2), "\n\n")
      })
      ans.all <- f.constr.con(ans.all, tmp.quick = T)
      return(ans.all)
    }
  }
  if (!adjust) 
    if (ans.all$model.ans != 16 && !tmp.quick && max(ans.all$fct5) == 
        1) 
      ans.all <- f.clear(ans.all)
    ans.all$adjust <- adjust
    ans.all$fitted <- fitted
    ans.all$tmp.quick <- tmp.quick
    ans.all$nr.aa <- max(ans.all$fct1)
    ans.all$nr.bb <- max(ans.all$fct2)
    ans.all$nr.var <- max(as.numeric(ans.all$fct3))
    ans.all$nr.cc <- max(ans.all$fct4)
    ans.all$nr.dd <- max(ans.all$fct5)
    ans.all$heading <- "start values"
    with(ans.all, {
      if (fit.ans == 2) 
        nr.var <- 0
      cc.upp <- 10000
      dd.low <- eps
      dd.upp <- 10
      if (dtype == 5) 
        dtype <- 1
      trace.tmp <- F
      if (quick.ans == 1) 
        trace.tmp <- T
      y.sort <- yy[order(x)]
      x.sort <- sort(x)
      if (adjust == F) {
        if (dtype %in% c(10, 15, 250, 260) & length(x) < 
            10) {
          xtmp <- x.sort
          ytmp <- y.sort
        }
        else {
          sub <- round(length(x.sort)/5)
          if (sub == 0) 
            sub <- 1
          xtmp <- numeric()
          ytmp <- numeric()
          regr.par <- numeric()
          e1 <- 0
          e2 <- 0
          for (k in 1:5) {
            e2 <- e2 + sub
            qq <- x.sort[e1:e2]
            xtmp[k] <- mean(qq[is.na(qq) == F])
            qq <- y.sort[e1:e2]
            if (dtype %in% c(25, 250)) 
              ytmp[k] <- mean(qq[is.na(qq) == F])
            if (dtype %in% c(26, 260)) 
              ytmp[k] <- mean(sqrt(qq[is.na(qq) == F]))^2
            else if (!(dtype %in% c(25, 26, 250, 260))) 
              ytmp[k] <- exp(mean(logb(qq[is.na(qq) == 
                                            F])))
            e1 <- e1 + sub
          }
        }
        if (length(ytmp[!is.na(ytmp)]) == 3) {
          ytmp[4:5] <- ytmp[3]
          xtmp[4:5] <- xtmp[3]
        }
        if (length(ytmp[!is.na(ytmp)]) == 4) {
          ytmp[5] <- ytmp[4]
          xtmp[5] <- xtmp[4]
        }
        top <- length(ytmp)
        aa <- NA
        bb <- NA
        cc <- NA
        dd <- NA
        if (model.ans != 1) {
          lm.res <- f.start.lm.con(xtmp, ytmp, model.ans, 
                                   dtype)
          aa <- abs(lm.res[1])
          if (model.ans %in% c(2, 3, 7, 8, 17, 18)) 
            bb <- lm.res[2]
          else bb <- abs(lm.res[2])
          if (quick.ans == 1 && model.type == 1) 
            increase <- lm.res[3]
          if (model.ans %in% c(26, 30, 46)) 
            increase <- -increase
          if (model.ans == 6 && max(fct4) > 1) 
            increase <- 0
          if (cont && model.ans %in% c(c(12:15, 22:25, 
                                         46)) && quick.ans == 1 && increase == -1 && 
              CES > 0) {
            incr.ans <- menu(c("yes", "no"), title = "\n you entered a positive value for CES, \n                          is the function increasing?\n")
            if (incr.ans == 1) 
              increase <- 1
            if (incr.ans == 2) 
              CES <- -CES
          }
          if (cont && model.ans %in% c(c(12:15, 22:25, 
                                         46)) && quick.ans == 1 && increase == 1 && 
              CES < 0) {
            ces.ans <- menu(c("yes", "no"), title = "\n you entered a negative value for CES, is the function decreasing?\n")
            if (ces.ans == 1) 
              increase <- -1
            if (ces.ans == 2) 
              CES <- -CES
          }
          ans.all$CES <- CES
        }
        switch(model.ans, {
          if (nr.bb > 1) {
            cat("\nThis model is not defined for different b parameters!")
            f.press.key.to.continue()
            ans.all$fct2 <- rep(1, length(x))
          }
          if (nr.cc > 1) {
            cat("\nThis model is not defined for different c parameters!")
            return(ans.all)
          }
          if (nr.dd > 1) {
            cat("\nThis model is not defined for different d parameters!")
            return(ans.all)
          }
          aa <- mean(ytmp)
          regr.par <- rep(aa, nr.aa)
        }, {
          if (fitted == T) {
            aa <- par.start[nr.var + 1]
            bb <- par.start[nr.var + nr.aa + 1]
          }
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb))
        }, {
          if (fitted) {
            aa <- par.start[nr.var + 1]
            bb <- par.start[nr.var + nr.aa + 1]
          }
          dd <- 1
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(dd, nr.dd))
        }, {
          if (fitted == F) {
            cc <- ytmp[top]/aa
            if (!is.finite(bb)) bb <- xtmp[5]/100
          } else {
            aa <- par.start[nr.var + 1]
            bb <- par.start[nr.var + nr.aa + 1]
            cc <- ytmp[top]/aa
          }
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(cc, nr.cc))
          if (cc < 1) increase <- -1 else increase <- 1
          if (increase.tmp != 0) increase <- increase.tmp
          incr <- increase + 2
        }, {
          dd <- 1
          if (fitted == F) {
            cc <- ytmp[top]/aa
            if (!is.finite(bb)) bb <- xtmp[5]/100
          } else {
            aa <- par.start[nr.var + 1]
            bb <- par.start[nr.var + nr.aa + 1]
            cc <- ytmp[top]/aa
          }
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(cc, nr.cc), rep(dd, nr.dd))
          if (cc < 1) increase <- -1 else increase <- 1
          if (increase.tmp != 0) increase <- increase.tmp
        }, {
          bb <- 1/bb
          dd <- 1
          if (fitted == F) {
            cc <- 1.01 * ytmp[top]/aa
            if (!is.finite(bb)) bb <- xtmp[5]/100
          } else {
            aa <- par.start[nr.var + 1]
            bb <- par.start[nr.var + nr.aa + 1]
            cc <- ytmp[top]/aa
          }
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(cc, nr.cc), rep(dd, nr.dd))
        }, {
          if (fitted == T) {
            aa <- par.start[nr.var + 1]
            bb <- par.start[nr.var + nr.aa + 1]
          } else bb <- 1/bb
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb))
        }, {
          if (fitted == F) bb <- 1/bb else {
            aa <- par.start[nr.var + 1]
            bb <- par.start[nr.var + nr.aa + 1]
          }
          dd <- 1
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(dd, nr.dd))
        }, {
          if (!fitted) {
            bb <- abs(1/bb)
            cc <- ytmp[top]/aa
          } else {
            aa <- par.start[nr.var + 1]
            bb <- par.start[nr.var + nr.aa + 1]
            cc <- ytmp[top]/aa
          }
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(cc, nr.cc))
        }, {
          cc <- ytmp[top]/aa
          bb <- 1/bb
          dd <- 1
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(cc, nr.cc), rep(dd, nr.dd))
        }, {
          regr.par <- f.regr.par.full(ans.all)
          if (0) {
            if (dtype %in% c(10, 15, 250, 260)) {
              if (dtype %in% c(10, 15)) regr.par <- exp(mn.log)
              if (dtype == 250) regr.par <- (mn.log)
              if (dtype == 260) regr.par <- (mn.log)^2
            } else {
              regr.par <- numeric()
              nn <- numeric()
              if (nr.var == 1) {
                if (twice) for (jj in levels(factor(fct2))) {
                  x.tmp <- x[fct2 == jj]
                  y.tmp <- yy[fct2 == jj]
                  regr.tmp <- as.numeric(exp(tapply(logb(y.tmp), 
                                                    x.tmp, mean)))
                  regr.par <- c(regr.par, regr.tmp)
                  nn.tmp <- tapply(y.tmp, x.tmp, length)
                  nn <- c(nn, nn.tmp)
                }
                if (!twice) for (jj in levels(factor(fct2))) for (ii in levels(factor(fct1))) {
                  x.tmp <- x[fct1 == ii & fct2 == jj]
                  y.tmp <- yy[fct1 == ii & fct2 == jj]
                  regr.tmp <- as.numeric(exp(tapply(logb(y.tmp), 
                                                    x.tmp, mean)))
                  regr.par <- c(regr.par, regr.tmp)
                  nn.tmp <- tapply(y.tmp, x.tmp, length)
                  nn <- c(nn, nn.tmp)
                }
              }
              if (nr.var > 1) {
                regr.par <- numeric()
                for (jj in levels(factor(fct3))) {
                  x.tmp <- x[fct3 == jj]
                  y.tmp <- yy[fct3 == jj]
                  regr.tmp <- as.numeric(exp(tapply(logb(y.tmp), 
                                                    x.tmp, mean)))
                  regr.par <- c(regr.par, regr.tmp)
                  nn.tmp <- tapply(y.tmp, x.tmp, length)
                  nn <- c(nn, nn.tmp)
                }
              }
              ans.all$nn <- nn
            }
          }
          ans.all$regr.par <- regr.par
        }, {
          if (fitted == T) {
            aa <- par.start[nr.var + 1]
            bb <- par.start[nr.var + nr.aa + 1]
          }
          if (fitted == F || is.na(bb) || bb == 0) {
            if (is.na(bb) || bb == 0) bb <- mean(xtmp)
            CED <- (1/bb) * logb(CES + 1)
            CED <- abs(CED)
            if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                                                        min(x))/3
          }
          regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb))
        }, {
          if (fitted == T) {
            aa <- par.start[nr.var + 1]
            bb <- par.start[nr.var + nr.aa + 1]
          }
          if (fitted == F || is.na(bb) || bb == 0) {
            if (is.na(bb) || bb == 0) bb <- mean(xtmp)
            CED <- (1/bb) * logb(CES + 1)
            CED <- abs(CED)
            if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                                                        min(x))/3
          }
          dd <- 1
          regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                        rep(dd, nr.dd))
        }, {
          if (fitted == T) {
            aa <- par.start[nr.var + 1]
            bb <- par.start[nr.var + nr.aa + 1]
            cc <- 0
          }
          if (fitted == F || is.na(bb) || bb == 0) {
            cc <- ytmp[top]/aa
            if (cont) {
              if (CES < 0 & cc > 1 + CES) {
                cc <- 1 + CES - (1 + CES)/100
              }
              if (CES > 0 & cc < 1 + CES) {
                cc <- 1 + CES + CES/100
              }
            }
            if (is.na(bb) || bb == 0) bb <- mean(xtmp)
            if (cc == 0) cc <- 0.1
            if (cc == 1) cc <- 0.9
            CED <- -(1/bb) * logb((CES + 1 - cc)/(1 - cc))
            CED <- abs(CED)
            if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                                                        min(x))/3
          }
          regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                        rep(cc, nr.cc))
        }, {
          dd <- 1
          cc <- ytmp[top]/aa
          if (cont) {
            if ((CES < 0) && (cc > 1 + CES)) {
              cc <- 1 + CES - (1 + CES)/100
            }
            if ((CES > 0) & (cc < 1 + CES)) {
              cc <- 1 + CES + CES/100
            }
          }
          if (cc == 0) cc <- 0.1
          if (cc == 1) cc <- 0.9
          if (!cont) cc <- 0.1
          CED <- (-(1/bb) * logb((-CES + 1 - cc)/(1 - cc)))^(1/dd)
          CED <- abs(CED)
          if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                                                      min(x))/3
          if (is.na(CED) || CED == 0) CED <- mean(xtmp)
          regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                        rep(cc, nr.cc), rep(dd, nr.dd))
        }, {
          regr.par <- MLE[-(1:nr.var)]
          if (!exists("CES.16")) {
            cat("\nCES.16 does not exist and is created (0.10 and 0.05)\n")
            CES.16 <- list(CES2 = 0.1, CES1 = 0.05)
            f.assign("CES.16", CES.16)
          }
          CES2 <- CES.16$CES2
          CES1 <- CES.16$CES1
          aa.vec <- regr.par[1:nr.aa]
          bb.vec <- regr.par[(nr.aa + 1):(nr.aa + nr.bb)]
          cc <- regr.par[nr.aa + nr.bb + 1]
          dd <- regr.par[length(regr.par)]
          if (cc < 1) {
            CES2 <- -CES2
            CES1 <- -CES1
            if (cc > 1 + CES2) {
              cat("\nATTENTION: value of c does not allow BMD ratio")
              f.press.key.to.continue()
              return(ans.all)
            }
          }
          if (cc > 1) {
            if (cc < 1 + CES2) {
              cat("\nATTENTION: value of c does not allow BMD ratio")
              return(ans.all)
            }
          }
          BMDratio <- (log((CES2 + 1 - cc)/(1 - cc))/log((CES1 + 
                                                            1 - cc)/(1 - cc)))^(1/dd)
          regr.par <- c(aa.vec, bb.vec, rep(cc, nr.cc), 
                        BMDratio)
        }, {
          if (bb == Inf) bb <- 1000
          if (model.type == 2) increase <- -1
          if (increase == 1) bb <- -abs(bb)
          if (increase == -1) bb <- abs(bb)
          if (increase == 1 & bb > max(x)) bb <- -max(x)
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb))
          if (cont) {
            par.dum <- c(rep(0.1, nr.var), regr.par)
            loglik.first <- -f.lik.con(par.dum, x, y, dtype, 
                                       fct1, fct2, fct3, model.ans, mn.log, sd2.log, 
                                       nn, Vdetlim = Vdetlim, CES = CES, twice = twice, 
                                       ttt = ttt, trace = trace.tmp, fct4 = fct4, 
                                       fct5 = fct5, cens.up = cens.up, par.tmp = NA, 
                                       increase = increase, ref.lev = ref.lev, sign.q = sign.q, 
                                       ans.m6.sd = ans.m6.sd)
            count <- 0
            while (!is.finite(loglik.first) & count < 10) {
              bb <- bb * 1.3
              regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb))
              par.dum <- c(rep(0.1, nr.var), regr.par)
              loglik.first <- -f.lik.con(par.dum, x, y, 
                                         dtype, fct1, fct2, fct3, model.ans, mn.log, 
                                         sd2.log, nn, Vdetlim = Vdetlim, CES = CES, 
                                         twice = twice, ttt = ttt, trace = trace.tmp, 
                                         fct4 = fct4, fct5 = fct5, cens.up = cens.up, 
                                         par.tmp = NA, increase = increase, ref.lev = ref.lev, 
                                         sign.q = sign.q, ans.m6.sd = ans.m6.sd)
              count <- count + 1
            }
          }
          if (cont) {
            y.expect <- f.expect.con(model.ans, x, regr.par, 
                                     CES = CES, increase = increase, fct1 = fct1, 
                                     fct2 = fct2)
            count <- 0
            while (any(y.expect < 0) && count < 20) {
              bb <- 1.5 * bb
              regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb))
              y.expect <- f.expect.con(model.ans, x, regr.par, 
                                       CES = CES, increase = increase, fct1 = fct1, 
                                       fct2 = fct2)
              count <- count + 1
            }
          }
        }, {
          if (bb == Inf) bb <- 1000
          if (model.type == 2) increase <- -1
          if (quick.ans == 1 && cont) if (bb < 0) increase <- 1 else increase <- -1
          if (increase == 1) bb <- -abs(bb)
          if (increase == -1) bb <- abs(bb)
          if (increase == 1 && bb > max(x)) bb <- -max(x)
          dd <- 1
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(dd, nr.dd))
          if (0) if (cont) {
            par.dum <- c(rep(0.1, nr.var), regr.par)
            loglik.first <- -f.lik.con(par.dum, x, y, dtype, 
                                       fct1, fct2, fct3, model.ans, mn.log, sd2.log, 
                                       nn, Vdetlim = Vdetlim, CES = CES, twice = twice, 
                                       ttt = ttt, trace = trace.tmp, fct4 = fct4, 
                                       fct5 = fct5, cens.up = cens.up, par.tmp = NA, 
                                       increase = increase, ref.lev = ref.lev)
            count <- 0
            while (!is.finite(loglik.first) & count < 10) {
              bb <- bb * 1.3
              regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                            rep(dd, nr.dd))
              par.dum <- c(rep(0.1, nr.var), regr.par)
              loglik.first <- -f.lik.con(par.dum, x, y, 
                                         dtype, fct1, fct2, fct3, model.ans, mn.log, 
                                         sd2.log, nn, Vdetlim = Vdetlim, CES = CES, 
                                         twice = twice, ttt = ttt, trace = trace.tmp, 
                                         fct4 = fct4, fct5 = fct5, cens.up = cens.up, 
                                         par.tmp = NA, increase = increase, ref.lev = ref.lev)
              count <- count + 1
            }
          }
          if (cont) {
            y.expect <- f.expect.con(model.ans, x, regr.par, 
                                     CES = CES, increase = increase, fct1 = fct1, 
                                     fct2 = fct2)
            y.expect <- y.expect[!is.na(y.expect)]
            count <- 0
            while (any(y.expect < 0) && count < 20) {
              bb <- 1.5 * bb
              regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                            rep(dd, nr.dd))
              y.expect <- f.expect.con(model.ans, x, regr.par, 
                                       CES = CES, increase = increase, fct1 = fct1, 
                                       fct2 = fct2)
              count <- count + 1
            }
          }
        }, {
          cc <- ytmp[top]/aa
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(cc, nr.cc))
        }, {
          cc <- ytmp[top]/aa
          dd <- 1
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(cc, nr.cc), rep(dd, nr.dd))
        }, {
          nr.doses <- length(ytmp)
          aa <- ytmp[1]
          bb1 <- mean(xtmp[1:trunc(nr.doses/2)])
          bb2 <- mean(xtmp[trunc(nr.doses/2):nr.doses])
          dd <- 1
          cc1 <- ytmp[nr.doses]/ytmp[1]
          cc2 <- ytmp[nr.doses]/ytmp[1]
          regr.par <- c(rep(aa, nr.aa), rep(bb1, nr.bb), 
                        bb2, cc1, cc2, dd)
          ans.all$lower <- c(eps, eps, eps, eps, eps, dd.low)
          ans.all$upper <- c(Inf, Inf, Inf, Inf, Inf, 100)
        }, {
          if (CES > 0) bb <- -bb
          CED <- f.inv.con(model.ans = 17, c(aa, bb), CES)
          if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                                                      min(x))/3
          regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb))
          if (cont) {
            y.expect <- f.expect.con(model.ans, x, regr.par, 
                                     CES = CES, increase = increase, fct1 = fct1, 
                                     fct2 = fct2)
            count <- 0
            while (any(y.expect < 0) && count < 20) {
              CED <- 1.5 * CED
              regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb))
              y.expect <- f.expect.con(model.ans, x, regr.par, 
                                       CES = CES, increase = increase, fct1 = fct1, 
                                       fct2 = fct2)
              count <- count + 1
            }
          }
        }, {
          if (CES > 0) bb <- -bb
          if (CES < 0) increase <- -1 else increase <- 1
          dd <- 1
          CED <- f.inv.con(model.ans = 18, c(aa, bb, dd), 
                           CES)
          regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                        rep(dd, nr.dd))
          if (cont) {
            y.expect <- f.expect.con(model.ans, x, regr.par, 
                                     CES = CES, increase = increase, fct1 = fct1, 
                                     fct2 = fct2)
            count <- 0
            while (any(y.expect < 0) && count < 20) {
              CED <- 1.5 * CED
              regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                            rep(dd, nr.dd))
              y.expect <- f.expect.con(model.ans, x, regr.par, 
                                       CES = CES, increase = increase, fct1 = fct1, 
                                       fct2 = fct2)
              count <- count + 1
            }
          }
        }, {
          cc <- ytmp[top]/aa
          if (cont) {
            if (CES < 0 & cc > 1 + CES) {
              cc <- 1 + CES - (1 + CES)/100
            }
            if (CES > 0 & cc < 1 + CES) {
              cc <- 1 + CES + CES/100
            }
          }
          bb <- (aa * (cc - 1) * mean(xtmp))/(mean(ytmp) - 
                                                aa) - mean(xtmp)
          if (is.na(bb) || bb == 0) bb <- mean(xtmp)
          if (cc == 0) cc <- 0.1
          if (cc == 1) cc <- 0.9
          CED <- f.inv.con(model.ans = 19, c(aa, bb, cc), 
                           CES)
          if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                                                      min(x))/3
          regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                        rep(cc, nr.cc))
        }, {
          cc <- ytmp[top]/aa
          dd <- 1
          if (cont) {
            if (CES < 0 & cc > 1 + CES) {
              cc <- 1 + CES - (1 + CES)/100
            }
            if (CES > 0 & cc < 1 + CES) {
              cc <- 1 + CES + CES/100
            }
          }
          if (is.na(bb) || bb == 0) bb <- mean(xtmp)
          if (cc == 0) cc <- 0.1
          if (cc == 1) cc <- 0.9
          CED <- f.inv.con(model.ans = 20, c(aa, bb, cc, 
                                             dd), CES)
          if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                                                      min(x))/3
          regr.par <- c(rep(aa, nr.aa), rep(CED, nr.bb), 
                        rep(cc, nr.cc), rep(dd, nr.dd))
        }, {
          dt.fr <- data.frame(yyy = ytmp, dose.tmp = xtmp)
          res.lm <- lm(yyy ~ dose.tmp, dt.fr)
          cc <- res.lm[[1]][1]
          aa <- res.lm[[1]][2]
          bb <- 1
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(cc, nr.cc))
          if (increase == 1) {
            ans.all$lower <- c(eps, eps, eps)
            ans.all$upper <- c(Inf, Inf, Inf)
          }
          if (increase == -1) {
            ans.all$lower <- c(-Inf, eps, eps)
            ans.all$upper <- c(eps, Inf, Inf)
          }
        }, {
          aa <- 1
          bb <- -(logb(aa) - logb(ytmp[top]))/(xtmp[top] - 
                                                 xtmp[1])
          bb <- 1/bb
          cc <- ytmp[1]
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb))
          ans.all$lower <- c(nr.aa, Inf, nr.bb)
          ans.all$upper <- c(Inf, Inf)
        }, {
          cc <- ytmp[1]
          aa <- xtmp[2]
          bb <- (ytmp[3] - cc)/(xtmp[3] - aa)
          dd <- 1
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(cc, nr.cc), dd)
          ans.all$lower <- c(eps, -Inf, eps, eps)
          ans.all$upper <- c(Inf, Inf, Inf, 100)
        }, {
          aa <- ytmp[1]
          cc <- ytmp[top]
          dd <- 1
          bb <- 1
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(cc, nr.cc), dd)
          ans.all$lower <- c(eps, eps, eps, eps)
          ans.all$upper <- c(Inf, Inf, Inf, 100)
        }, {
          dt.fr <- data.frame(yyy = ytmp, dose.tmp = xtmp)
          res.lm <- lm(yyy ~ dose.tmp, dt.fr)
          aa <- res.lm[[1]][1]
          bb <- res.lm[[1]][2]
          cc <- 1
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(cc, nr.cc))
          if (increase == 1) {
            ans.all$lower <- c(eps, eps, eps)
            ans.all$upper <- c(Inf, Inf, Inf)
          }
          if (increase == -1) {
            ans.all$lower <- c(eps, -Inf, eps)
            ans.all$upper <- c(Inf, eps, Inf)
          }
        }, {
          bb <- -(logb(aa) - logb(ytmp[top]))/(xtmp[top] - 
                                                 xtmp[1])
          bb <- 1/bb
          cc <- 0
          dd <- 1
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(cc, nr.cc), dd)
          ans.all$lower <- c(eps, eps, eps, eps)
          ans.all$upper <- c(Inf, Inf, Inf, 100)
        }, {
          nr.doses <- length(ytmp)
          aa <- ytmp[1]
          bb1 <- mean(xtmp[1:trunc(nr.doses/2)])
          bb2 <- mean(xtmp[trunc(nr.doses/2):nr.doses])
          dd1 <- 1
          dd2 <- 1
          cc <- ytmp[nr.doses]/ytmp[1]
          regr.par <- c(rep(aa, nr.aa), rep(bb1, nr.bb), 
                        bb2, rep(cc, nr.cc), dd1, dd2)
          ans.all$lower <- c(eps, eps, eps, eps, eps, eps)
          ans.all$upper <- c(Inf, Inf, Inf, Inf, 100, 100)
        }, {
          nr.doses <- length(ytmp)
          aa <- ytmp[1]
          bb1 <- mean(xtmp[1:trunc(nr.doses/2)])
          bb2 <- mean(xtmp[trunc(nr.doses/2):nr.doses])
          dd1 <- 1
          dd2 <- 1
          cc1 <- ytmp[nr.doses]/ytmp[1]
          cc2 <- ytmp[nr.doses]/ytmp[1]
          regr.par <- c(rep(aa, nr.aa), rep(bb1, nr.bb), 
                        bb2, cc1, cc2, dd1, dd2)
          ans.all$lower <- c(eps, eps, eps, eps, eps, eps, 
                             eps)
          ans.all$upper <- c(Inf, Inf, Inf, Inf, Inf, 100, 
                             100)
        }, {
          aa <- ytmp[1]
          bb <- (ytmp[3] - aa)/(xtmp[3])
          cc <- 0
          dd <- 0
          if (!increase) bb <- -bb
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(cc, nr.cc), dd)
          ans.all$lower <- c(-Inf, -Inf, -Inf, -Inf)
          ans.all$upper <- c(Inf, Inf, Inf, Inf)
        }, {
          aa <- eval(parse(prompt = paste("give value for parameter a", 
                                          "  > ")))
          bb <- aa/ytmp[top]
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb))
          ans.all$lower <- c(-Inf, -Inf)
          ans.all$upper <- c(Inf, Inf)
        }, {
          aa <- ytmp[1] + 0.001
          cc <- ytmp[top]/aa
          dd <- 1
          bb <- mean(x)
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(cc, nr.cc), dd)
          ans.all$lower <- c(eps, -Inf, eps, eps)
          ans.all$upper <- c(Inf, Inf, Inf, Inf)
        }, {
          bb <- ytmp[1] + 0.001
          dd <- 1
          aa <- 0
          cc <- 1
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(cc, nr.cc), dd)
          ans.all$lower <- c(eps, -Inf, eps, eps)
          ans.all$upper <- c(Inf, Inf, Inf, Inf)
        }, cat(""), cat(""), cat(""), {
          bb <- 1
          aa <- (ytmp[5] - aa)/(xtmp[5] - xtmp[1])
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb))
          ans.all$lower <- c(eps, eps)
          ans.all$upper <- c(Inf, Inf)
        }, {
        }, {
          dd <- 1
          cc <- ytmp[5]
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(cc, nr.cc), rep(dd, nr.dd))
          ans.all$lower <- c(eps, eps, rep(eps, nr.cc), 
                             eps)
          ans.all$upper <- c(Inf, Inf, rep(cc.upp, nr.cc), 
                             100)
        }, {
          if (fitted == T) {
            aa <- par.start[nr.var + 1]
            bb <- par.start[nr.var + nr.aa + 1]
          }
          if (fitted == F || is.na(bb) || bb == 0) {
            aa <- mean(logb(yy))
            bb <- logb(ytmp[5]/ytmp[1])/(xtmp[top] - xtmp[1])
          }
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb))
          ans.all$lower <- c(-Inf, -Inf)
          ans.all$upper <- c(Inf, Inf)
        }, {
          nr.doses <- length(ytmp)
          aa <- ytmp[1]
          bb <- mean(xtmp[1:trunc(nr.doses/2)])
          dd <- 1
          cc1 <- ytmp[nr.doses]/ytmp[1]
          cc2 <- ytmp[nr.doses]/ytmp[1]
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        cc1, cc2, dd)
          ans.all$lower <- c(eps, eps, eps, eps, dd.low)
          ans.all$upper <- c(Inf, Inf, Inf, Inf, 100)
        }, {
          dd <- 1
          cc <- ytmp[top]/aa
          if (cont) {
            if ((CES < 0) && (cc > 1 + CES)) {
              cc <- 1 + CES - (1 + CES)/100
            }
            if ((CES > 0) & (cc < 1 + CES)) {
              cc <- 1 + CES + CES/100
            }
          }
          if (cc == 0) cc <- 0.1
          if (cc == 1) cc <- 0.9
          if (!cont) cc <- 0.1
          CED <- (-(1/bb) * logb((-CES + 1 - cc)/(1 - cc)))^(1/dd)
          CED <- abs(CED)
          if (CED < (max(x) - min(x))/1000) CED <- (max(x) - 
                                                      min(x))/3
          if (is.na(CED) || CED == 0) CED <- mean(xtmp)
          CED.vec <- rep(1, nr.bb)
          CED.vec[ref.lev] <- CED
          regr.par <- c(rep(aa, nr.aa), CED.vec, rep(cc, 
                                                     nr.cc), rep(dd, nr.dd))
          if (0) {
            CED.vec <- regr.par[(nr.aa + 1):(nr.aa + nr.bb)]
            CED.ref <- CED.vec[ref.lev]
            RPF.vec <- CED.vec/CED.ref
            RPF.vec[ref.lev] <- CED.ref
            regr.par[(nr.aa + 1):(nr.aa + nr.bb)] <- RPF.vec
            ans.all$CED.ref <- CED.ref
          }
        }, {
          dd <- 1
          qq <- 6
          GCED <- mean(x)
          var <- 0.05
          par.start <- c(rep(var, nr.var), rep(aa, nr.aa), 
                         GCED, qq, dd)
          regr.par <- par.start
        }, {
          dd <- 1
          cc <- ytmp[top]/aa
          if (!is.finite(bb)) bb <- xtmp[5]/100
          gg <- bb
          hh <- 1
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(cc, nr.cc), rep(dd, nr.dd), gg, hh)
        }, {
          bb <- 1/bb
          dd <- 1
          bb <- bb/(log(2))^(1/dd)
          if (fitted == F) {
            cc <- 1.01 * ytmp[top]/aa
            if (!is.finite(bb)) bb <- xtmp[5]/100
          } else {
            aa <- par.start[nr.var + 1]
            bb <- par.start[nr.var + nr.aa + 1]
            cc <- ytmp[top]/aa
          }
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(cc, nr.cc), rep(dd, nr.dd))
        }, {
          bb <- 1
          dd <- 1
          cc <- ytmp[top]/aa
          gg <- mean(ttt)
          hh <- 1
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(cc, nr.cc), rep(dd, nr.dd), gg, hh)
        }, {
          dd <- 1
          if (fitted == F) {
            cc <- ytmp[top]/aa
            if (!is.finite(bb)) bb <- xtmp[5]/100
          } else {
            aa <- par.start[nr.var + 1]
            bb <- par.start[nr.var + nr.aa + 1]
            cc <- ytmp[top]/aa
          }
          regr.par <- c(rep(aa, nr.aa), rep(bb, nr.bb), 
                        rep(cc, nr.cc), rep(dd, nr.dd))
          if (cc < 1) increase <- -1 else increase <- 1
          incr <- increase + 2
          switch(incr, {
            dum.1 <- eps
            dum.2 <- 1
          }, {
            dum.1 <- eps
            dum.2 <- cc.upp
          }, {
            dum.1 <- 1
            dum.2 <- cc.upp
          })
          ans.all$lower <- c(eps, eps, dum.1, dd.low)
          ans.all$upper <- c(Inf, Inf, dum.2, dd.upp)
        })
        text.par <- f.text.par(ans.all)
        if (!(model.ans == 11 && model.type == 1)) {
          if (length(xans) == 2) {
            ans.all$RPF1 <- 1
            regr.par <- c(regr.par, RPF1)
            x <- x1 + RPF1 * x2
            ans.all$x <- x
            ans.all$xy.lim[2:3] <- c(min(ans.all$x), max(ans.all$x))
          }
          if (length(xans) == 3) {
            ans.all$RPF1 <- 1
            ans.all$RPF2 <- 1
            regr.par <- c(regr.par, RPF1, RPF2)
            x <- x1 + RPF1 * x2 + RPF2 * x3
            ans.all$x <- x
            ans.all$xy.lim[2:3] <- c(min(ans.all$x), max(ans.all$x))
          }
        }
        ans.all$regr.par <- regr.par
        if (!cont) {
          if (exists("track")) 
            print("f.start.con:  END.sub")
          ans.all$regr.par <- regr.par
          return(ans.all)
        }
        ans.all$nrp <- length(text.par) - nr.var
        if (!(dtype %in% c(10, 15, 250, 260))) {
          expect.trans <- f.expect.con(model.ans, x, regr.par, 
                                       fct1 = fct1, fct2 = fct2, fct3 = fct3, CES = CES, 
                                       twice = twice, ttt = ttt, y = yy, increase = increase, 
                                       x.mn = x.mn, ref.lev = ref.lev, trace.expect = F, 
                                       par.start = par.start, sign.q = sign.q, ans.m6.sd = ans.m6.sd)
          yy.trans <- yy
          # if (sum(expect.trans <= 0, na.rm = TRUE) > 0 && !(dtype %in% 
          #                                     c(25, 250))) {
          #   cat("\nATTENTION:  model predicts negative values \n")
          #   f.press.key.to.continue()
          # }
          if (dtype == 26) {
            expect.trans <- sqrt(expect.trans)
            yy.trans <- sqrt(yy)
          }
          if (dtype %in% c(1, 5, 15)) {
            expect.trans <- logb(expect.trans)
            yy.trans <- logb(yy)
          }
          if (sum(is.na(expect.trans)) > 0) {
            ans.all$par.start <- c(rep(NA, nr.var), regr.par)
            if (quick.ans > 1) {
              cat("\n there may be a fitting problem, try single model\n")
              ans.all$par.start[1] <- 0
              f.press.key.to.continue()
            }
            else {
              cat("\n adjust start values\n")
              ans.all$text.par <- text.par
              ans.all$adjust.start <- T
            }
            return(ans.all)
          }
          var.start <- var(yy.trans - expect.trans)
        }
        else var.start <- mean(sd2.log)
        if (is.na(var.start)) 
          var.start <- eval(parse(prompt = paste("\ngive value for start value of var (e.g. 0.5)\n")))
        if (model.ans != 47) 
          par.start <- c(rep(var.start, nr.var), regr.par)
        if (!model.ans %in% c(1:10, 12:20, 22:25)) {
          ans.all$lower <- c(1e-06, ans.all$lower)
          ans.all$upper <- c(1, ans.all$upper)
        }
        if (ans.all$model.ans == 16) {
          par.start <- MLE
          par.start[length(MLE)] <- BMDratio
          ans.all$regr.par <- par.start[-(1:nr.var)]
          ans.all$nrp <- length(MLE) - nr.var
        }
        ans.all$par.start <- par.start
        ans.all$text.par <- text.par
        ans.all$text.par.0 <- text.par
        if (plot.type > 0) 
          ans.all$xy.lim <- xy.lim
        if (plot.type > 8) 
          plot.type <- 1
        if (!tmp.quick && plot.type > 0) {
          ans.all$nr.gr <- f.plot.con(ans.all)
          ans.all$l.ty <- 2
          ans.all$increase <- increase
          if (ans.all$model.ans != 11) {
            ans.all$nr.gr <- f.lines.con(ans.all)
          }
        }
        switch(fit.ans, loglik.first <- -f.lik.con(ans.all$par.start, 
                                                   x, y, dtype, fct1, fct2, fct3, model.ans, mn.log, 
                                                   sd2.log, nn, Vdetlim = Vdetlim, CES = CES, twice = twice, 
                                                   ttt = ttt, trace = trace.tmp, fct4 = fct4, fct5 = fct5, 
                                                   cens.up = cens.up, par.tmp = NA, increase = increase, 
                                                   x.mn = x.mn, ref.lev = ref.lev, sign.q = sign.q, 
                                                   ans.m6.sd = ans.m6.sd), loglik.first <- f.lik.biv(par.start, 
                                                                                                     x, yy, fct1, fct2, model.ans, trace = F, twice = twice))
        if (is.na(loglik.first) | loglik.first < -1e+20) {
          cat("\nATTENTION: You need to adjust start values before fitting the model...\n")
          ans.all$adjust.start <- T
          f.press.key.to.continue()
          return(ans.all)
        }
        if (!tmp.quick) {
          if (model.ans == 11) {
            cat("\n\napproximate value of log-likelihood for saturated model: ", 
                round(loglik.first, 2), "\n")
            cat("\nfit the model for a more precise value of the log-likelihood\n\n")
            ans.all$loglik.first <- loglik.first
          }
          else cat("\n\nlog-likelihood at start values: ", 
                   round(loglik.first, 2), "\n\n")
        }
      }
      # if (adjust == T) {
      #   ans.all.plt <- ans.all
      #   switch(fit.ans, ans.all.plt$regr.par <- par.start[-(1:ans.all$nr.var)], 
      #          ans.all.plt$regr.par <- par.start)
      #   nrp <- length(ans.all.plt$regr.par)
      #   ans <- 1
      #   loglik.old <- 0
      #   if (sum(is.na(par.start)) == 0) {
      #     if (ans.plt == 1) {
      #       ans.all.plt$heading <- ""
      #       f.plot.sep(ans.all.plt)
      #     }
      #     else {
      #       f.plot.con(ans.all.plt)
      #     }
      #   }
      #   while (ans[1] != (nrp + (nr.var + 1))) {
      #     cat("\n=======================================================\n")
      #     cat(paste("model: ", model.names[model.ans]))
      #     cat("\n\nold log-likelihood: ", round(loglik.old, 
      #                                           3))
      #     loglik.new <- NA
      #     if (sum(is.finite(par.start)) == length(par.start)) 
      #       switch(fit.ans, loglik.new <- -f.lik.con(par.start, 
      #                                                x, y, dtype, fct1, fct2, fct3, model.ans, 
      #                                                mn.log, sd2.log, nn, Vdetlim = Vdetlim, CES = CES, 
      #                                                twice = twice, ttt = ttt, trace = trace.tmp, 
      #                                                fct4 = fct4, fct5 = fct5, cens.up = cens.up, 
      #                                                par.tmp = NA, increase = increase, x.mn = x.mn, 
      #                                                ref.lev = ref.lev, sign.q = sign.q, ans.m6.sd = ans.m6.sd), 
      #              loglik.new <- f.lik.biv(par.start, x, yy, 
      #                                      fct1, fct2, model.ans, trace = F, twice = twice))
      #     cat("\nnew log-likelihood: ", round(loglik.new, 
      #                                         3), "\n")
      #     cat("=======================================================\n\n")
      #     space <- rep("    current value:      ", (nrp + 
      #                                                 nr.var))
      #     if (model.ans == 11) 
      #       cat("\nATTENTION: log-likelihood not implemented for full model\n")
      #     par.list <- "start values"
      #     for (ii in 1:length(text.par)) {
      #       par.list <- paste(par.list, "\n", ii, "  ", 
      #                         text.par[ii], " ", signif(par.start[ii], 
      #                                                   3))
      #     }
      #     par.list <- paste(par.list, "\n", length(text.par) + 
      #                         1, "  ", "Continue\n")
      #     cat(par.list)
      #     ans <- eval(parse(prompt = paste("Which parameter(s) do you want to change?  > ")))
      #     if (ans[1] != (nrp + nr.var + 1)) {
      #       start.tmp <- eval(parse(prompt = paste("Please, give new value for ", 
      #                                              text.par[ans], "  >")))
      #       par.start[ans] <- start.tmp
      #       if (model.ans == 46) {
      #         CED.vec <- par.start[(nr.var + nr.aa + 1):(nr.var + 
      #                                                      nr.aa + nr.bb)]
      #         CED.ref <- CED.vec[ref.lev]
      #         ans.all$CED.ref <- CED.ref
      #         ans.all.plt$CED.ref <- CED.ref
      #       }
      #     }
      #     ans.all.plt$par.start <- par.start
      #     if (length(xans) == 2) {
      #       RPF1 <- par.start[length(par.start)]
      #       ans.all.plt$x <- x1 + RPF1 * x2
      #       ans.all.plt$xy.lim[2:3] <- c(min(ans.all.plt$x), 
      #                                    max(ans.all.plt$x))
      #     }
      #     if (length(xans) == 3) {
      #       RPF1 <- par.start[length(par.start) - 1]
      #       RPF2 <- par.start[length(par.start)]
      #       ans.all.plt$x <- x1 + RPF1 * x2 + RPF2 * x3
      #       ans.all.plt$xy.lim[2:3] <- c(min(ans.all.plt$x), 
      #                                    max(ans.all.plt$x))
      #     }
      #     switch(fit.ans, ans.all.plt$regr.par <- par.start[-(1:nr.var)], 
      #            ans.all.plt$regr.par <- par.start)
      #     if (ans.all.plt$plot.type == 0) 
      #       ans.all.plt$plot.type <- 1
      #     if (ans.plt == 1) {
      #       ans.all.plt$heading <- ""
      #       f.plot.sep(ans.all.plt)
      #     }
      #     else {
      #       f.plot.con(ans.all.plt)
      #       f.lines.con(ans.all.plt)
      #     }
      #     loglik.old <- loglik.new
      #   }
      #   ans.all$par.start <- par.start
      #   ans.all$adjust.start <- F
      #   ans.all$loglik <- loglik.new
      # }
      ans.all$nr.var <- nr.var
      ans.all$npar <- length(ans.all$par.start)
      ans.all$CED <- NA
      ans.all$increase <- increase
      ans.all <- f.constr.con(ans.all, tmp.quick = T)
      ans.all$l.ty <- 1
      if (exists("track")) 
        print("f.start.con:  END")
      return(ans.all)
    })
}

###############
## f.lik.con ##
###############
f.lik.con <- function (theta, x, y, dtype, fct1, fct2, fct3, model.ans, mn.log, 
                       sd2.log, nn, Vdetlim, CES, twice = T, ttt = 0, trace = F, 
                       fct4 = 1, fct5 = 1, cens.up = NA, lb = -Inf, ub = Inf, par.tmp, 
                       increase = increase, x.mn = NA, ref.lev = ref.lev, ans.m6.sd = 1, 
                       sign.q, x1 = 0, x2 = 0, x3 = 0) 
{
  if (model.ans == -13) {
    print("f.lik.con:  begin")
    print(theta)
  }
  if (sum(is.na(theta)) > 0) {
    if (trace) {
      cat("\nnote: problem in f.lik.con: NAs in theta\n")
      print(theta)
    }
    return(NA)
  }
  if (0) 
    if ((length(fct3) > 1) & (length(fct3) != length(x))) {
      print("f.lik.con")
      cat("\nATTENTION: fct3 incorrect\n\n")
      print(length(fct3))
      stop()
    }
  variance <- 0
  for (jj in 1:max(fct3)) variance <- variance + theta[jj] * 
    (fct3 == jj)
  if (!model.ans %in% c(6, 47) || ans.m6.sd != 2) 
    regr.par <- theta[(max(fct3) + 1):length(theta)]
  else regr.par <- theta
  if (model.ans == 11) {
    yy <- y
    yy[y == -1000] <- Vdetlim[y == -1000]
    yy[y == -2000] <- cens.up
  }
  if (sum(!is.finite(theta)) > 0) {
    theta[!is.finite(theta)] <- par.tmp[!is.finite(theta)]
  }
  if (any(x3 > 0)) {
    RPF1 <- regr.par[length(regr.par) - 1]
    RPF2 <- regr.par[length(regr.par)]
    x <- x1 + RPF1 * x2 + RPF2 * x3
  }
  else if (any(x2 > 0)) {
    RPF1 <- regr.par[length(regr.par)]
    x <- x1 + RPF1 * x2
  }
  theta <- theta[!is.na(theta)]
  if (sum(theta <= lb) > 0) {
    theta[theta < lb] <- 1.1 * par.tmp[theta < lb]
  }
  if (sum(theta >= ub) > 0) {
    theta[theta > ub] <- 0.9 * par.tmp[theta > ub]
  }
  expect <- f.expect.con(model.ans, x, regr.par, fct1 = fct1, 
                         fct2 = fct2, fct3 = fct3, fct4 = fct4, fct5 = fct5, CES = CES, 
                         twice = twice, ttt = ttt, y = yy, trace.expect = F, increase = increase, 
                         x.mn = x.mn, ref.lev = ref.lev, sign.q = sign.q, ans.m6.sd = ans.m6.sd)
  if (dtype == 1 | dtype == 5) {
    if (trace & sum(is.na(expect)) > 0) {
      cat("\nnote: NAs in predicted response\nat parameter values:\n")
      cat(signif(theta, 4), "\n")
    }
    expect <- logb(expect)
    y.log <- logb(y)
    y.log[y < 0] <- 0
    score1 <- (y > 0) * (-0.5 * logb(2 * pi * variance) - 
                           ((y.log - expect)^2)/(2 * variance))
    score.detlim <- logb(pnorm((logb(Vdetlim) - expect)/sqrt(variance)))
    score.detlim[!is.finite(score.detlim)] <- 0
    score.censup <- logb(1 - pnorm((logb(cens.up) - expect)/sqrt(variance)))
    score.censup[!is.finite(score.censup)] <- 0
    score2 <- (y == -1000) * score.detlim + (y == -2000) * 
      score.censup
    score <- score1 + score2
  }
  if (dtype == 25) {
    score1 <- (y > 0) * (-0.5 * logb(2 * pi * variance) - 
                           ((y - expect)^2)/(2 * variance))
    score.detlim <- logb(pnorm(((Vdetlim) - expect)/sqrt(variance)))
    score.detlim[!is.finite(score.detlim)] <- 0
    score.censup <- logb(1 - pnorm(((cens.up) - expect)/sqrt(variance)))
    score.censup[!is.finite(score.censup)] <- 0
    score2 <- (y == -1000) * score.detlim + (y == -2000) * 
      score.censup
    score <- score1 + score2
  }
  if (dtype == 26) {
    expect <- sqrt(expect)
    y.sqrt <- sqrt(yy)
    score1 <- (y > 0) * (-0.5 * logb(2 * pi * variance) - 
                           ((y.sqrt - expect)^2)/(2 * variance))
    score.detlim <- logb(pnorm((sqrt(Vdetlim) - expect)/sqrt(variance)))
    score.detlim[!is.finite(score.detlim)] <- 0
    score.censup <- logb(1 - pnorm((sqrt(cens.up) - expect)/sqrt(variance)))
    score.censup[!is.finite(score.censup)] <- 0
    score2 <- (y == -1000) * score.detlim + (y == -2000) * 
      score.censup
    score <- score1 + score2
  }
  if (dtype %in% c(10, 15)) {
    if (model.ans != 11) {
      expect <- logb(f.expect.con(model.ans, x, regr.par, 
                                  fct1 = fct1, fct2 = fct2, fct3 = fct3, fct4 = fct4, 
                                  fct5 = fct5, CES = CES, twice = twice, ttt = ttt, 
                                  y = y, increase = increase, x.mn = x.mn, ref.lev = ref.lev, 
                                  sign.q = sign.q, ans.m6.sd = ans.m6.sd))
      dum <- nn * (mn.log - expect)^2 + (nn - 1) * sd2.log
    }
    else dum <- (nn - 1) * sd2.log
    score <- -(nn * logb(sqrt(2 * pi * variance)) + dum/(2 * 
                                                           variance))
  }
  if (dtype == 250) {
    if (model.ans != 11) {
      expect <- f.expect.con(model.ans, x, regr.par, fct1 = fct1, 
                             fct2 = fct2, fct3 = fct3, fct4 = fct4, fct5 = fct5, 
                             CES = CES, twice = twice, ttt = ttt, y = y, increase = increase, 
                             x.mn = x.mn, ref.lev = ref.lev, sign.q = sign.q, 
                             ans.m6.sd = ans.m6.sd)
      dum <- nn * (mn.log - expect)^2 + (nn - 1) * sd2.log
    }
    else dum <- (nn - 1) * sd2.log
    score <- -(nn * logb(sqrt(2 * pi * variance)) + dum/(2 * 
                                                           variance))
  }
  if (dtype == 260) {
    if (model.ans != 11) {
      expect <- sqrt(f.expect.con(model.ans, x, regr.par, 
                                  fct1 = fct1, fct2 = fct2, fct3 = fct3, fct4 = fct4, 
                                  fct5 = fct5, CES = CES, twice = twice, ttt = ttt, 
                                  y = y, increase = increase, x.mn = x.mn, ref.lev = ref.lev, 
                                  sign.q = sign.q, ans.m6.sd = ans.m6.sd))
      dum <- nn * (mn.log - expect)^2 + (nn - 1) * sd2.log
    }
    else dum <- (nn - 1) * sd2.log
    score <- -(nn * logb(sqrt(2 * pi * variance)) + dum/(2 * 
                                                           variance))
  }
  if (model.ans == -13) {
    print("f.lik.con END")
    print(sum(score))
  }
  return(-sum(score))
}

##################
## f.constr.con ##
##################
f.constr.con <- function (ans.all, tmp.quick = F) 
{
  if (exists("track")) 
    print("f.constr.con")
  ans.all$tmp.quick <- tmp.quick
  with(ans.all, {
    if (!cont && model.type == 1) {
      cat("\n f.constr.con used for model.type = 1\n")
      f.press.key.to.continue()
    }
    nr.aa <- max(fct1)
    nr.bb <- max(fct2)
    nr.cc <- max(fct4)
    nr.dd <- max(fct5)
    nr.var <- max(fct3)
    lower.var <- 1e-06
    upper.var <- 10
    if (dtype %in% c(25, 250)) 
      upper.var <- 1e+10
    if (dtype %in% c(26, 260)) 
      upper.var <- 10000
    if (model.ans %in% c(1:10, 12:20, 22:25, 46:48, 50)) {
      if (cont) {
        lower.aa <- min(yy/100)
        upper.aa <- max(yy * 100)
      }
      if (!cont) {
        if (model.type == 2) {
          lower.aa <- 1
          if (model.ans == 1) 
            lower.aa <- 0.01
          if (nth > 4) 
            upper.aa <- exp(2 * sig.start)
          else switch(nth, upper.aa <- exp(6 * sig.start), 
                      upper.aa <- exp(5 * sig.start), upper.aa <- exp(4 * 
                                                                        sig.start), upper.aa <- exp(3 * sig.start))
        }
      }
      lower.bb <- -Inf
      upper.bb <- Inf
      if (model.ans %in% c(4:6, 9:10, 13:15, 16, 19:20, 
                           23:25, 46:48, 50)) 
        lower.bb <- 1e-06
      if (model.ans %in% 17:18) {
        if (increase == 1) {
          lower.bb <- -Inf
          upper.bb <- -1.01 * max(x)
        }
        else {
          lower.bb <- eps
          upper.bb <- Inf
        }
      }
      if (0) 
        if (!cont && dtype != 3 && quick.ans == 1 && 
            model.ans %in% c(4, 14, 19, 24)) {
          incr.ans <- 0
          if (increase == 1 && par.start[npar] <= 1) 
            incr.ans <- menu(c("increasing", "decreasing"), 
                             title = "\nis the dose-response decreasing or increasing?\n")
          if (increase == -1 && par.start[npar] > 1) 
            incr.ans <- menu(c("increasing", "decreasing"), 
                             title = "\nis the dose-response decreasing or increasing?\n")
          if (incr.ans == 1) 
            increase <- 1
          if (incr.ans == 2) 
            increase <- -1
          ans.all$increase <- increase
        }
      if (model.ans %in% 3) {
        if (increase == -1) {
          lower.bb <- -Inf
          upper.bb <- 0
        }
        if (increase == 1) {
          lower.bb <- 0
          upper.bb <- Inf
        }
      }
      if (model.ans %in% 18) {
        if (increase == 1) {
          lower.bb <- -Inf
          upper.bb <- 0
        }
        if (increase == -1) {
          lower.bb <- 0
          upper.bb <- Inf
        }
      }
      if (increase == 1) {
        lower.cc <- 1.1
        upper.cc <- 1e+06
      }
      if (increase == -1) {
        lower.cc <- 1e-06
        upper.cc <- 0.9
      }
      if (increase == 0) {
        lower.cc <- 1e-06
        upper.cc <- 1000
      }
      if (model.ans == 6 && ans.m6.sd == 2) {
        lower.cc <- 1e-06
        upper.cc <- 1000
      }
      if (cont && model.ans %in% c(4, 5, 14, 15, 19, 20, 
                                   24, 25, 46)) {
        if (increase == 1) {
          lower.cc <- (1 + abs(CES)) * 1.01
          if (lower.cc > 1.1) {
            #cat("\nATTENTION: lower constraint on parameter cc has been increased to", 
            #    lower.cc, "\ndue to value of CES\n")
            #cat("This might result in suboptimal fit of model", 
            #    model.ans, "\n")
          }
        }
        if (increase == -1) {
          upper.cc <- (1 - abs(CES))/1.01
          if (upper.cc < 0.9) {
            #cat("\nATTENTION: upper constraint on parameter cc has been decreased to", 
            #    upper.cc, "\ndue to value of CES\n")
            #cat("This might result in suboptimal fit of model", 
            #    model.ans, "\n")
          }
        }
      }
      lower.dd <- 0.25
      upper.dd <- 4
      lower.RPF.1 <- 1e-12
      upper.RPF.1 <- 1e+12
      lower.RPF.2 <- 1e-12
      upper.RPF.2 <- 1e+12
      if (model.ans == 1) {
        lower <- c(lower.aa)
        upper <- c(upper.aa)
      }
      if (model.ans %in% c(2, 7, 12, 17, 22)) {
        lower <- c(lower.aa, lower.bb)
        upper <- c(upper.aa, upper.bb)
      }
      if (model.ans %in% c(3, 8, 13, 18, 23)) {
        lower <- c(lower.aa, lower.bb, lower.dd)
        upper <- c(upper.aa, upper.bb, upper.dd)
      }
      if (model.ans %in% c(4, 9, 14, 19, 24)) {
        lower <- c(lower.aa, lower.bb, lower.cc)
        upper <- c(upper.aa, upper.bb, upper.cc)
      }
      if (model.ans %in% c(5, 6, 10, 15, 16, 20, 25, 46)) {
        lower <- c(lower.aa, lower.bb, lower.cc, lower.dd)
        upper <- c(upper.aa, upper.bb, upper.cc, upper.dd)
      }
      if (model.ans == 47) {
        lower.qq <- 5
        upper.qq <- 10
        lower <- c(lower.aa, lower.bb, lower.qq, lower.dd)
        upper <- c(upper.aa, upper.bb, upper.qq, upper.dd)
      }
      if (model.ans %in% c(48, 50)) {
        lower.gg <- lower.bb
        upper.gg <- upper.bb
        lower.hh <- lower.dd
        upper.hh <- upper.dd
        lower <- c(lower.aa, lower.bb, lower.cc, lower.dd, 
                   lower.gg, lower.hh)
        upper <- c(upper.aa, upper.bb, upper.cc, upper.dd, 
                   upper.gg, upper.hh)
      }
      if (cont && fit.ans == 1) {
        lower <- c(lower.var, lower)
        upper <- c(upper.var, upper)
      }
    }
    else if (model.ans == 11) {
      lower <- c(rep(lower.var, nr.var), regr.par)
      upper <- c(rep(upper.var, nr.var), regr.par)
    }
    else {
      if (cont) 
        dum <- 0
      else dum <- 1
      lower.var <- lower[1 - dum]
      lower.aa <- lower[2 - dum]
      lower.bb <- lower[3 - dum]
      lower.cc <- lower[4 - dum]
      lower.dd <- lower[5 - dum]
      lower.ee <- lower[6 - dum]
      lower.ff <- lower[7 - dum]
      lower.gg <- lower[8 - dum]
      lower.hh <- lower[9 - dum]
      upper.var <- upper[1 - dum]
      upper.aa <- upper[2 - dum]
      upper.bb <- upper[3 - dum]
      upper.cc <- upper[4 - dum]
      upper.dd <- upper[5 - dum]
      upper.ee <- upper[6 - dum]
      upper.ff <- upper[7 - dum]
      upper.gg <- upper[8 - dum]
      upper.hh <- upper[9 - dum]
    }
    if (tmp.quick == F) {
      if (model.ans == 46) {
        cat("\n parameter constraints cannot be canged \n               (not implemented for this model)\n")
        f.press.key.to.continue()
      }
      else {
        vabcd <- f.text.par(ans.all, brief = T)
        if (length(xans) > 1) {
          lower <- c(lower, lower.RPF.1)
          upper <- c(upper, upper.RPF.1)
        }
        if (length(xans) > 2) {
          lower <- c(lower, lower.RPF.2)
          upper <- c(upper, upper.RPF.2)
        }
        ans <- 1
        while (ans <= length(vabcd)) {
          ans <- menu(c(paste(vabcd, ":     ", lower, 
                              " --- ", upper), "continue"), title = "\n Which parameter do you want to constraint? ")
          if (ans <= length(vabcd)) {
            cat("\n give infinite constraints as Inf or -Inf\n\n")
            lower[ans] <- eval(parse(prompt = paste("give lower bound for ", 
                                                    vabcd[ans], "  ", "  > ")))
            upper[ans] <- eval(parse(prompt = paste("give upper bound for ", 
                                                    vabcd[ans], "  > ")))
          }
        }
        lower.var <- lower[1]
        lower.aa <- lower[2]
        lower.bb <- lower[3]
        upper.var <- upper[1]
        upper.aa <- upper[2]
        upper.bb <- upper[3]
        if (model.ans %in% c(4, 5, 6, 9, 10, 14, 15, 
                             19, 20, 24, 25, 26, 30)) {
          lower.cc <- lower[4]
          upper.cc <- upper[4]
          lower.dd <- lower[5]
          upper.dd <- upper[5]
          if (length(xans) > 1) {
            lower.RPF.1 <- lower[6]
            upper.RPF.1 <- upper[6]
          }
          if (length(xans) > 2) {
            lower.RPF.2 <- lower[7]
            upper.RPF.2 <- upper[7]
          }
        }
        if (model.ans %in% c(3, 8, 13, 18, 23)) {
          lower.dd <- lower[4]
          upper.dd <- upper[4]
          if (length(xans) > 1) {
            lower.RPF.1 <- lower[5]
            upper.RPF.1 <- upper[5]
          }
          if (length(xans) > 2) {
            lower.RPF.2 <- lower[6]
            upper.RPF.2 <- upper[6]
          }
        }
        if (model.ans %in% c(2, 7, 12, 17, 22)) {
          if (length(xans) > 1) {
            lower.RPF.1 <- lower[4]
            upper.RPF.1 <- upper[4]
          }
          if (length(xans) > 2) {
            lower.RPF.2 <- lower[5]
            upper.RPF.2 <- upper[5]
          }
        }
        if (model.ans == 47) {
          lower.qq <- lower[4]
          upper.qq <- upper[4]
          lower.dd <- lower[5]
          upper.dd <- upper[5]
        }
      }
    }
    if (model.ans == 1) {
      lb <- c(rep(lower.aa, nr.aa))
      ub <- c(rep(upper.aa, nr.aa))
    }
    if (model.ans %in% c(2, 7, 12, 17, 22, 35, 41, 44)) {
      lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb))
      ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb))
    }
    if (model.ans %in% c(3, 8, 13, 18, 23)) {
      lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
              rep(lower.dd, nr.dd))
      ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
              rep(upper.dd, nr.dd))
    }
    if (model.ans %in% c(4, 9, 14, 19, 24, 26, 27, 30)) {
      lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
              rep(lower.cc, nr.cc))
      ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
              rep(upper.cc, nr.cc))
    }
    if (model.ans %in% c(5, 6, 10, 15, 16, 20, 25, 28, 29, 
                         31, 34, 36, 37, 41, 42, 43)) {
      lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
              rep(lower.cc, nr.cc), rep(lower.dd, nr.dd))
      ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
              rep(upper.cc, nr.cc), rep(upper.dd, nr.dd))
    }
    if (model.ans == 11) {
      lb <- regr.par
      ub <- regr.par
    }
    if (model.ans == 21) {
      lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
              lower.cc, lower.dd, lower.ee, rep(lower.ff, nr.dd))
      ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
              upper.cc, upper.dd, upper.ee, rep(upper.ff, nr.dd))
    }
    if (model.ans == 32) {
      lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
              lower.cc, lower.dd, lower.ee, lower.ff)
      ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
              upper.cc, upper.dd, upper.ee, upper.ff)
    }
    if (model.ans == 33) {
      lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
              lower.cc, lower.dd, lower.ee, lower.ff, lower.gg)
      ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
              upper.cc, upper.dd, upper.ee, upper.ff, upper.gg)
    }
    if (model.ans == 45) {
      lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
              lower.cc, lower.dd, lower.ee)
      ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
              upper.cc, upper.dd, upper.ee)
    }
    if (model.ans == 46) {
      lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
              lower.cc, lower.dd)
      ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
              upper.cc, upper.dd)
    }
    if (model.ans == 47) {
      lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
              lower.qq, lower.dd)
      ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
              upper.qq, upper.dd)
    }
    if (model.ans %in% c(48, 50)) {
      lb <- c(rep(lower.aa, nr.aa), rep(lower.bb, nr.bb), 
              rep(lower.cc, nr.cc), rep(lower.dd, nr.dd), lower.gg, 
              lower.hh)
      ub <- c(rep(upper.aa, nr.aa), rep(upper.bb, nr.bb), 
              rep(upper.cc, nr.cc), rep(upper.dd, nr.dd), upper.gg, 
              upper.hh)
    }
    if (!(model.ans == 11 && model.type == 1)) {
      if (length(xans) > 1) {
        lb <- c(lb, lower.RPF.1)
        ub <- c(ub, upper.RPF.1)
      }
      if (length(xans) > 2) {
        lb <- c(lb, lower.RPF.2)
        ub <- c(ub, upper.RPF.2)
      }
    }
    if (cont && fit.ans == 1) {
      lb <- c(rep(lower.var, nr.var), lb)
      ub <- c(rep(upper.var, nr.var), ub)
    }
    if (!cont) {
      if (model.ans %in% c(2:3, 7:8)) 
        ub[(nr.aa + 1):(nr.aa + nr.bb)] <- 0
      if (model.ans %in% c(4:5, 9:10)) 
        lb[(nr.aa + 1):(nr.aa + nr.bb)] <- 0
      if (model.ans %in% c(17:20)) 
        lb[(nr.aa + 1):(nr.aa + nr.bb)] <- 0
      if (model.ans %in% c(4, 5, 9, 10, 14, 15, 19, 20, 
                           24, 25, 46)) {
        lb[nr.aa + nr.bb + 1] <- 1e-06
        ub[nr.aa + nr.bb + 1] <- 0.999999
      }
      if (model.ans %in% c(12:15, 22:25)) {
        if (nr.aa > 1 && nr.bb == 1) {
          lb[2:(nr.aa + 1)] <- eps
          ub[2:(nr.aa + 1)] <- Inf
        }
        else {
          lb[(nr.aa + 1):(nr.aa + nr.bb)] <- eps
          ub[(nr.aa + 1):(nr.aa + nr.bb)] <- Inf
        }
      }
      lb[nrp + 1] <- th.0.start
      ub[nrp + 1] <- th.0.start
      if (dtype == 6) 
        npar <- npar - 1
      if (nth > 1) {
        lb[(nrp + 2):npar] <- -Inf
        ub[(nrp + 2):npar] <- 1e-10
      }
      if (model.type == 1 && dtype == 6) 
        npar <- npar + 1
      lb[npar] <- sig.start
      ub[npar] <- sig.start
      if (max(as.numeric(fct3)) > 1) {
        lb[npar] <- 0
        ub[npar] <- Inf
      }
      if (model.ans %in% c(3, 5, 18, 20, 13, 15, 23, 25)) 
        if (constr != -Inf) 
          lb[length(par.start) - nth - 1] <- constr
      if (model.ans %in% c(4, 5, 9, 10, 14, 15, 24, 25)) {
        lb[nr.aa + nr.bb + 1] <- 0
        ub[nr.aa + nr.bb + 1] <- 1
      }
      if (dtype == 6) {
        lb <- c(1e-10, lb)
        ub <- c(Inf, ub)
      }
    }
    ans.all$lb <- lb
    ans.all$ub <- ub
    ans.all$lower <- lower
    ans.all$upper <- upper
    text.par <- f.text.par(ans.all)
    for (ii in (1:length(lower))) {
      if (lower[ii] == upper[ii]) 
        text.par[ii] <- paste(text.par[ii], "(fixed)", 
                              sep = "")
    }
    if (exists("track")) 
      print("f.constr.con:  END")
    ans.all$text.par <- text.par
    return(ans.all)
  })
}

##############
## f.nlminb ##
##############
f.nlminb <- function (ans.all, tmp.quick = F) 
{
  if (exists("track")) 
    print("f.nlminb")
  ans.all$tmp.quick <- tmp.quick
  with(ans.all, {
    if (0) {
      print("f.nlminb")
      print(par.start)
    }
    scale.dum <- abs(1/par.start)
    scale.dum[scale.dum == Inf] <- 1000
    if (dtype == 6) 
      scale.dum[1] <- 1
    if (scale.ans && quick.ans == 1 && model.ans != 1 && 
        !((model.ans == 14 & model.type == 1) || (model.ans == 
                                                  11 & model.type == 2))) {
      ans <- 1
      while (ans <= npar) {
        ans <- menu(c(paste(text.par, ":      ", signif(scale.dum, 
                                                        3)), "continue"), title = "\n For what parameter do you want to change scaling used in nlminb? ")
        if (ans <= npar) 
          scale.dum[ans] <- eval(parse(prompt = paste("Give value larger than one ", 
                                                      text.par[ans], ": ")))
      }
    }
    ans.all$scale.dum <- scale.dum
    trace.tmp <- F
    if (quick.ans == 1) 
      trace.tmp <- T
    if (cont) {
      if (model.ans %in% c(4, 5, 9, 10, 14, 15, 19, 20, 
                           24, 25)) {
        cc.start <- par.start[(nr.var + nr.aa + nr.bb + 
                                 1):(nr.var + nr.aa + nr.bb + nr.cc)]
        if (increase == -1 && cc.start >= (1 - abs(CES))) 
          cc.start <- 0.9 * (1 - abs(CES))
        if (increase == 1 && cc.start <= (1 + abs(CES))) 
          cc.start <- 1.1 * (1 + abs(CES))
        par.start[(nr.var + nr.aa + nr.bb + 1):(nr.var + 
                                                  nr.aa + nr.bb + nr.cc)] <- cc.start
      }
      loglik.check <- f.lik.con(par.start, x, y, dtype, 
                                fct1, fct2, fct3, model.ans, mn.log, sd2.log, 
                                nn, Vdetlim = Vdetlim, CES = CES, twice = twice, 
                                ttt = ttt, trace = trace.tmp, fct4 = fct4, fct5 = fct5, 
                                cens.up = cens.up, par.tmp = NA, increase = increase, 
                                x.mn = x.mn, ref.lev = ref.lev, sign.q = sign.q, 
                                ans.m6.sd = ans.m6.sd, x1 = x1, x2 = x2, x3 = x3)
      count <- 0
      while (!is.finite(loglik.check) && model.ans %in% 
             c(17, 18, 22, 23) && count < 30) {
        ced.start <- par.start[(nr.var + nr.aa + 1):(nr.var + 
                                                       nr.aa + nr.bb)]
        ced.start <- 1.5 * ced.start
        par.start[(nr.var + nr.aa + 1):(nr.var + nr.aa + 
                                          nr.bb)] <- ced.start
        loglik.check <- f.lik.con(par.start, x, y, dtype, 
                                  fct1, fct2, fct3, model.ans, mn.log, sd2.log, 
                                  nn, Vdetlim = Vdetlim, CES = CES, twice = twice, 
                                  ttt = ttt, trace = trace.tmp, fct4 = fct4, 
                                  fct5 = fct5, cens.up = cens.up, par.tmp = NA, 
                                  increase = increase, x.mn = x.mn, ref.lev = ref.lev, 
                                  sign.q = sign.q, ans.m6.sd = ans.m6.sd, x1 = x1, 
                                  x2 = x2, x3 = x3)
        count <- count + 1
      }
      if (model.ans == -25) {
        cat("\n")
        print("f.nlminb start")
        print(model.ans)
        print(CES)
        print(increase)
        print(cbind(text.par, par.start, lb, ub, scale.dum))
        print(loglik.check)
        f.press.key.to.continue()
      }
      fit.res <- nlminb(par.start, f.lik.con, scale = scale.dum, 
                        lower = lb, upper = ub, control = control.lst, 
                        x = x, y = y, dtype = dtype, fct1 = fct1, fct2 = fct2, 
                        fct3 = fct3, fct4 = fct4, fct5 = fct5, model.ans = model.ans, 
                        mn.log = mn.log, sd2.log = sd2.log, nn = nn, 
                        Vdetlim = Vdetlim, CES = CES, twice = twice, 
                        trace = F, cens.up = cens.up, lb = lb, ub = ub, 
                        par.tmp = par.start, increase = increase, x.mn = x.mn, 
                        ref.lev = ref.lev, sign.q = sign.q, ans.m6.sd = ans.m6.sd, 
                        x1 = x1, x2 = x2, x3 = x3)
      if (is.na(fit.res$obj)) 
        fit.res$obj <- 1e-12
      count <- 0
      while (is.na(fit.res$par[1]) || (fit.res$obj == 0) || 
             abs(fit.res$obj) == 1e-12 || !is.finite(fit.res$obj) && 
             count <= 20) {
        cat("\n      model is refitted with other scaling parameter")
        scale.dum <- rep(0.5, length(par.start))
        scale.dum <- 2 * scale.dum
        fit.res <- nlminb(par.start, f.lik.con, scale = scale.dum, 
                          lower = lb, upper = ub, control = control.lst, 
                          x = x, y = y, dtype = dtype, fct1 = fct1, fct2 = fct2, 
                          fct3 = fct3, fct4 = fct4, fct5 = fct5, model.ans = model.ans, 
                          mn.log = mn.log, sd2.log = sd2.log, nn = nn, 
                          Vdetlim = Vdetlim, CES = CES, twice = T, cens.up = cens.up, 
                          lb = lb, ub = ub, par.tmp = fit.res$par, increase = increase, 
                          x.mn = x.mn, ref.lev = ref.lev, sign.q = sign.q, 
                          ans.m6.sd = ans.m6.sd, x1 = x1, x2 = x2, x3 = x3)
        if (is.na(fit.res$obj)) 
          fit.res$obj <- 1e-12
        count <- count + 1
        if (count > 20) {
          fit.res$obj <- -1e+10
          fit.res$par <- 0
        }
      }
      ans.all$MLE <- fit.res$par
      if (model.ans == 47 || (model.ans == 6 && ans.m6.sd == 
                              2)) 
        ans.all$regr.par <- ans.all$MLE
      else ans.all$regr.par <- ans.all$MLE[-(1:max(fct3))]
      if (0) {
        loglik.tmp <- -f.lik.con(MLE, x, y, dtype, fct1, 
                                 fct2, fct3, model.ans, mn.log, sd2.log, nn, 
                                 Vdetlim = Vdetlim, CES = CES, twice = twice, 
                                 ttt = ttt, trace = T, fct4 = fct4, fct5 = fct5, 
                                 cens.up = cens.up, par.tmp = NA, increase = increase, 
                                 x.mn = x.mn, ref.lev = ref.lev, sign.q = sign.q, 
                                 ans.m6.sd = ans.m6.sd, x1 = x1, x2 = x2, x3 = x3)
        print(loglik.tmp)
        f.press.key.to.continue()
      }
    }
    if (!cont) {
      date.tmp <- date()
      loglik.check <- f.lik.cat(par.start, x = x, y = y, 
                                kk = kk, nn = nn, dtype = dtype, fct1 = fct1, 
                                fct2 = fct2, nrp = nrp, nth = nth, nr.aa = nr.aa, 
                                nr.bb = nr.bb, model.ans = model.ans, model.type = model.type, 
                                ans.nobg = 0, CES = CES, CES.cat = CES.cat, ttt = ttt, 
                                twice = twice, cens.up = cens.up, fct3 = fct3, 
                                fct4 = fct4, fct5 = fct5, x.full = x.full, fct1.full = fct1.full, 
                                fct2.full = fct2.full, alfa.length = alfa.length, 
                                trace = F, ces.ans = ces.ans, decr.zz = decr.zz, 
                                CES1 = CES1, CES2 = CES2, nn.tot = nn.tot, kk.tot = kk.tot, 
                                fct3.ref = fct3.ref, ref.lev = ref.lev, x1 = x1, 
                                x2 = x2, x3 = x3, output = output)
      if (model.ans == -13 || model.ans == -15) {
        cat("\n")
        print("f.nlminb start")
        print(model.ans)
        print(CES)
        print(increase)
        print(cbind(text.par, par.start, lb, ub, scale.dum))
        print(kk[1:10])
        print(loglik.check)
        f.press.key.to.continue()
      }
      retry <- FALSE
      if (!is.finite(loglik.check) || abs(loglik.check) == 
          2e+10) 
        retry <- TRUE
      count <- 0
      while (retry && (model.type == 2 && model.ans %in% 
                       c(13, 18, 23)) && count < 15) {
        ced.start <- par.start[(alfa.length + nr.aa + 
                                  1):(alfa.length + nr.aa + nr.bb)]
        ced.start <- 1.5 * ced.start
        par.start[(alfa.length + nr.aa + 1):(alfa.length + 
                                               nr.aa + nr.bb)] <- ced.start
        if (count < 3) 
          cat("\nModel re-fitted with new start values\n")
        loglik.check <- f.lik.cat(par.start, x = x, y = y, 
                                  kk = kk, nn = nn, dtype = dtype, fct1 = fct1, 
                                  fct2 = fct2, nrp = nrp, nth = nth, nr.aa = nr.aa, 
                                  nr.bb = nr.bb, model.ans = model.ans, model.type = model.type, 
                                  ans.nobg = 0, CES = CES, CES.cat = CES.cat, 
                                  ttt = ttt, twice = twice, cens.up = cens.up, 
                                  fct3 = fct3, fct4 = fct4, fct5 = fct5, x.full = x.full, 
                                  fct1.full = fct1.full, fct2.full = fct2.full, 
                                  alfa.length = alfa.length, trace = F, ces.ans = ces.ans, 
                                  decr.zz = decr.zz, CES1 = CES1, CES2 = CES2, 
                                  nn.tot = nn.tot, kk.tot = kk.tot, fct3.ref = fct3.ref, 
                                  ref.lev = ref.lev, x1 = x1, x2 = x2, x3 = x3, 
                                  output = output)
        if (is.finite(loglik.check) && abs(loglik.check) < 
            2e+10) 
          retry <- FALSE
        count <- count + 1
      }
      fit.res <- nlminb(par.start, f.lik.cat, scale = scale.dum, 
                        lower = lb, upper = ub, control = control.lst, 
                        x = x, y = y, kk = kk, nn = nn, dtype = dtype, 
                        fct1 = fct1, fct2 = fct2, nrp = nrp, nth = nth, 
                        nr.aa = nr.aa, nr.bb = nr.bb, model.ans = model.ans, 
                        model.type = model.type, ans.nobg = 0, CES = CES, 
                        CES.cat = CES.cat, ttt = ttt, twice = twice, 
                        cens.up = cens.up, fct3 = fct3, fct4 = fct4, 
                        fct5 = fct5, x.full = x.full, fct1.full = fct1.full, 
                        fct2.full = fct2.full, alfa.length = alfa.length, 
                        trace = F, ces.ans = ces.ans, decr.zz = decr.zz, 
                        CES1 = CES1, CES2 = CES2, nn.tot = nn.tot, kk.tot = kk.tot, 
                        fct3.ref = fct3.ref, ref.lev = ref.lev, x1 = x1, 
                        x2 = x2, x3 = x3, output = output)
      if (model.ans == -15) {
        print(" f.nlminb 3333")
        print(model.ans)
        print(CES)
        print(CES.cat)
        print(cbind(text.par, par.start, lb, ub, scale.dum))
        print(fit.res$objective)
        print(fit.res$par)
      }
      while (is.na(fit.res$par[1]) | (fit.res$obj == 0) | 
             (fit.res$obj == 1e-12) | is.na(fit.res$obj)) {
        cat("\n     model", model.ans, "refitted with other scaling parameter")
        if (model.ans == 14 && model.type == 1 && dtype == 
            6) 
          par.start[1] <- eval(parse(prompt = paste("give start value for alfa", 
                                                    "  > ")))
        scale.dum <- rep(0.5, length(par.start))
        scale.dum <- 2 * scale.dum
        count <- count + 1
        fit.res <- nlminb(par.start, f.lik.cat, scale = scale.dum, 
                          lower = lb, upper = ub, control = control.lst, 
                          x = x, y = y, kk = kk, nn = nn, dtype = dtype, 
                          fct1 = fct1, fct2 = fct2, nrp = nrp, nth = nth, 
                          nr.aa = nr.aa, nr.bb = nr.bb, model.ans = model.ans, 
                          model.type = model.type, ans.nobg = 0, CES = CES, 
                          CES.cat = CES.cat, ttt = ttt, twice = twice, 
                          cens.up = cens.up, fct3 = fct3, fct5 = fct5, 
                          x.full = x.full, fct1.full = fct1.full, fct2.full = fct2.full, 
                          alfa.length = alfa.length, trace = F, ces.ans = ces.ans, 
                          decr.zz = decr.zz, CES1 = CES1, CES2 = CES2, 
                          fct3.ref = fct3.ref, ref.lev = ref.lev, x1 = x1, 
                          x2 = x2, x3 = x3, output = output)
        if (count > 50) {
          fit.res$obj <- -1e+10
          fit.res$par <- 0
        }
      }
      ans.all$MLE <- fit.res$par
      if (model.type == 2) {
        par.lst <- f.split.par(ans.all$MLE, nrp, nth, 
                               dtype)
        ans.all$regr.par <- par.lst$regr.par
        ans.all$th.par <- par.lst$th.par
        ans.all$sig.par <- par.lst$sig.par
      }
      else ans.all$regr.par <- ans.all$MLE
    }
    ans.all$loglik <- round(-fit.res$objective, 2)
    if ((model.type == 1 && model.ans == 14) || (cont && 
                                                 model.ans == 11)) 
      npar.aic <- npar
    else npar.aic <- length(ans.all$MLE) - sum(lb == ub)
    ans.all$aic <- 2 * npar.aic - 2 * ans.all$loglik
    ans.all$npar.aic <- npar.aic
    if (!is.finite(ans.all$loglik)) 
      ans.all$loglik <- -1e+12
    mess <- fit.res$message
    ans.all$converged <- f.converged(mess, fit.res$conv, 
                                     tmp.quick = tmp.quick)
    if (dtype == 6 && (model.ans == 14 & model.type == 1)) {
      dum <- list()
      dum$loglik <- ans.all$loglik
      dum$MLE <- ans.all$MLE
      ans.all$full.model <- dum
      alfa.start <- ans.all$MLE[1]
      nn.dum <- table(x, x)
      nn.dum <- nn.dum[nn.dum != 0]
      if (length(nn.dum) > 15 && alfa.start > 20) 
        alfa.start <- 5
      ans.all$alfa.start <- alfa.start
    }
    ans.all$fitted <- T
    ans.all$l.ty <- 1
    ans.all$new.model <- F
    if (0) {
      ans.all$fit.res <- fit.res
      f.vcov(ans.all)
    }
    if (exists("track")) 
      print("f.nlminb : END")
    return(ans.all)
  })
}

#################
## f.converged ##
#################
f.converged <- function (mess, conv.out, tmp.quick = F) 
{
  if (exists("track2")) 
    print("f.converged")
  converged <- F
  if (mode(mess) == "NULL") {
    converged <- F
    if (!tmp.quick) 
      cat("\nATTENTION: convergence message is NULL\n")
  }
  else converged <- 1 - conv.out
  if (!tmp.quick) 
    cat("\n CONVERGION TYPE: ", mess, "\n")
  if (exists("track2")) 
    print("f.converged:  END")
  return(converged)
}

############
## f.pars ##
############
f.pars <- function (ans.all) 
{
  if (exists("track2")) 
    print("f.pars")
  if (ans.all$fit.ans == 1) 
    if (ans.all$model.ans == 47 || (ans.all$model.ans == 
                                    6 && ans.all$ans.m6.sd == 2)) 
      return(f.pars.m6(ans.all))
  with(ans.all, {
    nrp <- length(regr.par)
    regr.par.matr <- numeric()
    nr.aa <- max(fct1)
    nr.bb <- max(fct2)
    if (length(ans.all$nr.cc) == 0) 
      nr.cc <- 1
    if (length(ans.all$nr.dd) == 0) 
      nr.dd <- 1
    if (length(ans.all$nr.var) == 0) 
      nr.var <- 1
    if (length(fct4) == 1) 
      fct4 <- rep(1, length(x))
    nr.subgr <- max(nr.aa, nr.bb, nr.cc, nr.dd)
    if (nr.aa == 1) 
      fct1.txt <- rep("", nr.subgr)
    if (nr.bb == 1) 
      fct2.txt <- rep("", nr.subgr)
    if (nr.cc == 1) 
      fct4.txt <- rep("", nr.subgr)
    if (nr.dd == 1) 
      fct5.txt <- rep("", nr.subgr)
    if (identical(fct1, fct2)) 
      fct2.txt <- rep("", nr.subgr)
    if (identical(fct1, fct4)) 
      fct4.txt <- rep("", nr.subgr)
    if (identical(fct4, fct2)) 
      fct4.txt <- rep("", nr.subgr)
    if (nr.aa == 1 && nr.bb == 1 && nr.cc == 1 && nr.dd == 
        1 && nr.var < 2) {
      ans.all$regr.par.matr <- matrix(regr.par, nrow = 1)
      ans.all$nr.gr <- 1
      ans.all$gr.txt <- ""
      if (exists("track2")) 
        print("end of f.pars")
      return(ans.all)
    }
    if (cont || model.type == 2) {
      if (model.ans %in% c(3, 8, 13, 18, 23)) 
        nr.cc <- 0
      if (model.ans %in% c(4, 9, 14, 19, 24)) 
        nr.dd <- 0
      if (model.ans %in% c(2, 12, 17, 22)) {
        nr.cc <- 0
        nr.dd <- 0
      }
      cc <- regr.par[nr.aa + nr.bb + 1]
      dd <- regr.par[nr.aa + nr.bb + nr.cc + 1]
    }
    if (!cont && model.type == 1) {
      nr.dd <- 0
    }
    gr.txt <- character(0)
    if (!cont && model.type == 2 && (model.ans %in% c(12:15, 
                                                      22:25)) && ces.ans %in% c(1:3, 5)) {
      if (nrp > (nr.aa + nr.bb)) 
        par.rest <- regr.par[(nr.aa + nr.bb + 1):nrp]
      else par.rest <- numeric()
      if (nr.aa > 1 && nr.bb == 1) 
        for (ii in 1:nr.aa) {
          par.tmp <- c(regr.par[1], regr.par[ii + 1], 
                       par.rest)
          regr.par.matr <- rbind(regr.par.matr, par.tmp)
        }
      if (nr.aa == 1 & nr.bb > 1) 
        for (jj in 1:nr.bb) {
          par.tmp <- c(regr.par[1], regr.par[jj + 1], 
                       par.rest)
          regr.par.matr <- rbind(regr.par.matr, par.tmp)
        }
      if (nr.aa > 1 & nr.bb > 1) 
        for (jj in 1:nr.bb) {
          par.tmp <- c(regr.par[jj], regr.par[jj + nr.aa], 
                       par.rest)
          regr.par.matr <- rbind(regr.par.matr, par.tmp)
        }
      if (nr.aa == 1 & nr.bb == 1) 
        regr.par.matr <- matrix(regr.par, nrow = 1)
      ans.all$regr.par.matr <- regr.par.matr
      if (quick.ans > 1 && length(ans.all$covar.txt) > 
          0) 
        fct1.txt <- covar.txt
      kk <- 1
      for (jj in 1:nr.bb) {
        f1 <- fct1[fct2 == jj]
        f1.lev <- levels(factor(f1))
        for (ii.index in f1.lev) {
          ii <- as.numeric(ii.index)
          if (!twice) 
            gr.txt[kk] <- paste(fct1.txt[ii], fct2.txt[jj], 
                                sep = "-")
          if (twice) {
            gr.txt[kk] <- paste(fct1.txt[kk], sep = "-")
          }
          kk <- kk + 1
        }
      }
      ans.all$nr.gr <- length(ans.all$regr.par.matr[, 1])
      ans.all$gr.txt <- gr.txt
      if (exists("track2")) 
        print("f.pars LVM: END ")
      return(ans.all)
    }
    else if (nr.dd < 2) {
      if ((dtype %in% c(1, 5)) && (model.ans == 11) && 
          (nr.var > 1)) 
        fct1 <- fct3
      par.tmp <- rep(NA, length(regr.par))
      kk <- 1
      for (jj in 1:nr.bb) {
        f1 <- fct1[fct2 == jj]
        f1.lev <- levels(factor(f1))
        for (ii.index in f1.lev) {
          ii <- as.numeric(ii.index)
          f4 <- fct4[fct1 == ii & fct2 == jj]
          f4.lev <- levels(factor(f4))
          for (mm.index in f4.lev) {
            mm <- as.numeric(mm.index)
            par.tmp <- c(regr.par[ii], regr.par[nr.aa + 
                                                  jj])
            if (nr.cc > 0) 
              par.tmp <- c(par.tmp, regr.par[nr.aa + 
                                               nr.bb + mm])
            if (nr.dd > 0) 
              par.tmp <- c(par.tmp, dd)
            if (length(xans) == 2) 
              par.tmp <- c(par.tmp, regr.par[length(regr.par)])
            if (length(xans) == 3) 
              par.tmp <- c(par.tmp, regr.par[length(regr.par) - 
                                               1], regr.par[length(regr.par)])
            regr.par.matr <- rbind(regr.par.matr, par.tmp)
            gr.txt[kk] <- paste(fct1.txt[ii], fct2.txt[jj], 
                                fct4.txt[mm], sep = "-")
            kk <- kk + 1
          }
        }
      }
      if (!cont && model.ans == 14 && model.type == 1) {
        ans.all$regr.par.matr <- NA
        ans.all$gr.txt <- gr.txt
        return(ans.all)
      }
      else if (cont && model.ans == 11) {
        ans.all$regr.par.matr <- NA
        ans.all$gr.txt <- gr.txt
        return(ans.all)
      }
      else {
        n.col <- length(par.tmp)
        regr.par.matr <- matrix(regr.par.matr[, 1:n.col], 
                                ncol = n.col)
      }
    }
    if (exists("track2")) 
      print("f.pars, continue 2")
    if (nr.cc == 1 && nr.dd > 1) {
      if ((dtype %in% c(1, 5)) && (model.ans == 11) && 
          (nr.var > 1)) 
        fct1 <- fct3
      par.tmp <- rep(NA, length(regr.par))
      gr.txt <- character(0)
      for (jj in 1:nr.bb) {
        f1 <- fct1[fct2 == jj]
        f1.lev <- levels(factor(f1))
        for (ii.index in f1.lev) {
          ii <- as.numeric(ii.index)
          f5 <- fct5[fct1 == ii & fct2 == jj]
          f5.lev <- levels(factor(f5))
          for (mm.index in f5.lev) {
            mm <- as.numeric(mm.index)
            par.tmp <- c(regr.par[ii], regr.par[nr.aa + 
                                                  jj], regr.par[nr.aa + nr.bb + 1])
            if (nr.dd > 0) 
              par.tmp <- c(par.tmp, regr.par[nr.aa + 
                                               nr.bb + 1 + mm])
            if (length(xans) == 2) 
              par.tmp <- c(par.tmp, regr.par[length(regr.par)])
            if (length(xans) == 3) 
              par.tmp <- c(par.tmp, regr.par[length(regr.par) - 
                                               1], regr.par[length(regr.par)])
            regr.par.matr <- rbind(regr.par.matr, par.tmp)
            gr.txt[kk] <- paste(fct1.txt[ii], fct2.txt[jj], 
                                fct5.txt[mm], sep = "-")
            kk <- kk + 1
          }
        }
      }
      n.col <- length(par.tmp)
      if (model.ans == 14 && model.type == 1) {
        ans.all$regr.par.matr <- NA
        ans.all$gr.txt <- gr.txt
        ans.all$nr.gr <- nr.gr
        return(ans.all)
      }
      else regr.par.matr <- matrix(regr.par.matr[, 1:n.col], 
                                   ncol = n.col)
    }
    if (exists("track2")) 
      print("f.pars, continue 3")
    if (0) {
      if (identical(fct1, fct2) && identical(fct2, fct4)) 
        if (nr.aa > 1 && nr.bb > 1 && nr.cc > 1) {
          ii <- 0
          jj <- nr.aa
          if (model.ans %in% c(5, 6, 10, 15, 16, 20, 
                               25, 41, 42)) 
            dd <- regr.par[nr.aa + nr.bb + nr.cc + 1]
          else dd <- NULL
          kk <- nr.aa + nr.bb
          zz <- 1
          for (jj in 1:nr.cc) {
            if (nr.aa > 1) 
              ii <- ii + 1
            else ii <- 1
            if (nr.bb > 1) 
              jj <- jj + nr.aa
            else jj <- nr.aa + 1
            kk <- kk + 1
            par.tmp <- c(regr.par[ii], regr.par[jj], 
                         regr.par[kk], dd)
            regr.par.matr <- rbind(regr.par.matr, par.tmp)
            gr.txt[zz] <- paste(fct4.txt[zz])
            zz <- zz + 1
          }
        }
      if (identical(fct1, fct2) && identical(fct2, fct5)) 
        if (nr.aa > 1 && nr.bb > 1 && nr.dd > 1) {
          ii <- 0
          jj <- nr.aa
          if (model.ans %in% c(3, 8, 13, 18, 23, 25, 
                               41, 42)) {
            cc <- NULL
            kk <- nr.aa + nr.bb
          }
          if (model.ans %in% c(5, 6, 10, 15, 16, 20, 
                               25, 41, 42)) {
            cc <- regr.par[nr.aa + nr.bb + 1]
            kk <- nr.aa + nr.bb + 1
          }
          zz <- 1
          for (jj in 1:nr.dd) {
            if (nr.aa > 1) 
              ii <- ii + 1
            else ii <- 1
            if (nr.bb > 1) 
              jj <- jj + nr.aa
            else jj <- nr.aa + 1
            kk <- kk + 1
            par.tmp <- c(regr.par[ii], regr.par[jj], 
                         cc, regr.par[kk])
            regr.par.matr <- rbind(regr.par.matr, par.tmp)
            gr.txt[zz] <- paste(fct5.txt[zz])
            zz <- zz + 1
          }
        }
    }
    ans.all$regr.par.matr <- regr.par.matr
    ans.all$gr.txt <- gr.txt
    ans.all$nr.gr <- length(regr.par.matr[, 1])
    if (exists("track2")) 
      print("f.pars: END ")
    return(ans.all)
  })
}

###################
## f.report.pars ##
###################
f.report.pars <- function (ans.all) 
{
  if (exists("track")) 
    print("f.report.pars")
  with(ans.all, {
    report.pars <- paste(" estimate for ", text.par, ": ", 
                         signif(MLE, 4), "\n")
    if (quick.ans == 1) {
      report.pars <- c(report.pars, "\nlog-likelihood: ", 
                       round(ans.all$loglik, 3))
      report.pars <- c(report.pars, "\nAIC: ", round(ans.all$aic, 
                                                     3))
      if (cont && CES != 0) 
        report.pars <- c(report.pars, "\nCES: ", ans.all$CES)
      if (!cont && BMR != 0) 
        report.pars <- c(report.pars, "\nBMR: ", ans.all$BMR)
      if (cont && fit.ans == 1) {
        if (max(fct3) == 1) 
          report.pars <- c(report.pars, "\nestimate for residual C.V.", 
                           fct3.txt[1], ":", round(100 * sqrt(exp(MLE[1]) - 
                                                                1), 2))
        else for (jj in 1:max(fct3)) report.pars <- c(report.pars, 
                                                      "\nestimate for residual C.V. in group", fct3.txt[jj], 
                                                      ":", round(sqrt(exp(MLE[jj]) - 1), 2))
      }
      if (model.type == 2 && !cont) {
        if (dtype == 6) 
          nrp <- nrp + 1
        report.pars[nrp + 1] <- " fixed value for th-1 : 0 \n"
        report.pars[nrp + nth + 1] <- " fixed value for sigma : 0.25 \n"
      }
    }
    return(report.pars)
  })
}

################
## f.plot.con ##
################
f.plot.con <- function (ans.all) 
{
  if (exists("track")) 
    print("f.plot.con")
  ans.all <- f.data.plt.con(ans.all)
  with(ans.all, {
    if (plot.type == 0) 
      cat("\nWARNING: plot type = 0 (f.plot.con)\n")
    if (plot.type == 9) {
      qq.out <- qqnorm(regr.resid, plot = FALSE)
      qqnorm(regr.resid, xlab = "theoretical quantiles", 
             ylab = "regression residuals", main = "", cex = cex.1, 
             mkh = 0.03, plot = TRUE)
      xxx <- f.distr(regr.resid, Vdetlim, plot.ans = FALSE, 
                     ans.distr = 3)
      mu <- xxx[1]
      sig2 <- xxx[2]
      abline(mu, sqrt(sig2))
      f.dtype.txt(dtype)
      return(ans.all)
    }
    bbb <- (x.lim.plt[2] - x.lim.plt[1])/70
    dfr <- points.plt.lst[[1]]
    par(mar = rep(1,4))
    plot(dfr$x, dfr$y, xlim = x.lim.plt, ylim = y.lim.plt, 
         xlab = xleg, ylab = yleg, cex = 0.5, col = 1, type = "n")
    title(main = heading)
    nr.gr <- length(points.plt.lst)
    if (plt.mns < 3) 
      for (ii in 1:nr.gr) {
        dfr <- points.plt.lst[[ii]]
        x.tmp <- dfr$x
        y.tmp <- dfr$y
        if (mark[ii] == 1) 
          cex.1 <- cex.1 * 2
        points(x.tmp, y.tmp, col = color[ii], pch = mark[ii], 
               cex = cex.1)
        if (mark[ii] == 1) 
          cex.1 <- cex.1/2
      }
    if (plt.mns %in% c(1, 3)) 
      for (ii in 1:nr.gr) {
        dfr <- means.plt.lst[[ii]]
        x.tmp <- dfr$x
        y.tmp <- dfr$y
        type.dum <- "p"
        l.ty <- 1
        if (model.ans == 11) {
          l.ty <- 2
          type.dum <- "b"
        }
        points(x.tmp, y.tmp, col = color[ii], pch = mark[ii], 
               type = type.dum, cex = cex.2)
        if (CI.plt) {
          L025 <- dfr$CI.low
          L975 <- dfr$CI.upp
          if (!is.na(L025[1])) 
            for (qq in 1:length(x.tmp)) {
              for (qq in 1:length(x.tmp)) {
                lines(rep(x.tmp[qq], 2), c(L025[qq], 
                                           L975[qq]), col = color[ii])
                lines(c(x.tmp[qq] - bbb, x.tmp[qq] + 
                          bbb), rep(L025[qq], 2), col = color[ii])
                lines(c(x.tmp[qq] - bbb, x.tmp[qq] + 
                          bbb), rep(L975[qq], 2), col = color[ii])
              }
            }
        }
      }
    if (plot.type == 10) {
      f.dtype.txt(dtype)
      abline(0, 0, lty = 2)
    }
    if (0) 
      if (plot.type == 11) {
        x.plt <- x
        if (x.leg != "") 
          xleg <- paste("log10 of normalized ", x.leg, 
                        sep = "")
        else xleg <- ""
        y.plt <- yy
        x.plt <- x.plt[!is.na(y.plt)]
        y.plt <- y.plt[!is.na(y.plt)]
        yleg <- "residual"
        yleg <- "residual"
        x.lim.tmp <- log10(x.lim)
        if (!is.finite(x.lim.tmp)[1]) 
          x.lim.tmp[1] <- log10(dum.contr)
        if (is.na(xy.lim[4])) 
          y.lim <- c(min(y.plt, na.rm = T), max(y.plt, 
                                                na.rm = T))
        y.lim.tmp <- c(min(y.plt, na.rm = T), max(y.plt, 
                                                  na.rm = T))
      }
    if (max(Vdetlim, na.rm = T) > 0 & plot.type < 9) {
      abline(detlim.plt.1, 0, lty = 3)
      abline(detlim.plt.2, 0, lty = 3)
    }
    if (!is.na(cens.up) & plot.type < 9) {
      lines(c(min(xy.plt$x.plt, na.rm = T), max(xy.plt$x.plt, 
                                                na.rm = T)), rep(cens.up, 2), lty = 3)
    }
    if (any(label != 0)) {
      if (plot.type == 10) 
        cat("\nthis option is not possible for residuals\n")
      else {
        cat("\n\nUse LEFT mouse button to click on data points that you want to identify\n\n")
        cat("and RIGHT mouse button to stop\n")
        flush.console()
        identify(xy.plt$x.plt, xy.plt$y.plt, labels = label)
      }
    }
    ans.all$nr.gr <- nr.gr
    if (exists("track")) 
      print("f.plot.con:  END")
    return(ans.all)
  })
}

####################
## f.data.plt.con ##
####################
f.data.plt.con <- function (ans.all) 
{
  if (exists("track")) 
    print("f.data.plt.con")
  with(ans.all, {
    if (0) {
      if (is.na(cex.1)) {
        if (plt.mns == 2) 
          cex.1 <- 1
        else cex.1 <- 0.6
        if (plot.type == 9) 
          cex.1 <- 0.6
        if (dtype %in% c(10, 15, 110, 250, 260)) 
          cex.1 <- 1.6
      }
    }
    if (is.na(cex.2)) {
      if (plt.mns == 2) 
        cex.2 <- 1
      else cex.2 <- 1.5
    }
    if (dtype %in% c(10, 15, 110, 250, 260)) 
      plt.mns <- 3
    if (covar.no > 0) 
      nr.bb <- max(covariate)
    if (length(xans) > 1) {
      fct2 <- factor(data.0[, displ.no])
      gr.txt <- levels(fct2)
      fct2 <- as.numeric(fct2)
      nr.bb <- max(fct2)
    }
    else {
      gr.txt <- f.pars(ans.all)$gr.txt
      if (nr.aa == 1 && nr.bb == 1 && nr.var == 1) 
        if (max(displ.fact) > 1 && quick.ans == 1) {
          fct2 <- displ.fact
          nr.bb <- max(fct2)
          cat("\n\nno covariates, but subgroups are plotted distinctly\n\n")
        }
    }
    dum.contr <- xy.lim[1]
    dum.zero.resp <- xy.lim[6]
    x.lim <- xy.lim[2:3]
    y.lim <- xy.lim[4:5]
    if (!is.na(cens.up)) 
      y.lim <- c(xy.lim[4], max(xy.lim[5], 1.5 * cens.up))
    low.x <- x.lim[1]
    upp.x <- x.lim[2]
    yy <- yy[x >= low.x & x <= upp.x]
    nn <- nn[x >= low.x & x <= upp.x]
    sd2.log <- sd2.log[x >= low.x & x <= upp.x]
    x <- x[x >= low.x & x <= upp.x]
    if (length(fct1) > 1) 
      fct1 <- fct1[x >= low.x & x <= upp.x]
    if (length(fct2) > 1) 
      fct2 <- fct2[x >= low.x & x <= upp.x]
    if (length(fct3) > 1) 
      fct3 <- fct3[x >= low.x & x <= upp.x]
    mark <- c(1:2, 4:25, 33:500)
    if (plot.type == 0) 
      cat("\nWARNING: plot type = 0 (f.data.plt.con)\n")
    f.means <- function(xx, yy, dtype, plot.type, CI = F, 
                        sd2.log, nn) {
      if (exists("track")) 
        print("f.means within f.data.plt.con")
      mean.y <- yy
      if (dtype %in% c(1, 5, 101)) 
        mean.y <- exp(tapply(logb(yy), xx, mean))
      if (dtype %in% c(26)) 
        mean.y <- (tapply(sqrt(yy), xx, mean))^2
      if (dtype %in% c(25) || plot.type == 10) 
        mean.y <- tapply(yy, xx, mean)
      if (plot.type %in% c(1, 2, 5, 10, 11, 15)) 
        mean.y.plt <- mean.y
      if (plot.type %in% c(3, 4, 6)) 
        mean.y.plt <- log10(mean.y)
      if (plot.type %in% c(7, 8)) 
        mean.y.plt <- sqrt(mean.y)
      if (!CI) 
        return(mean.y.plt)
      if (CI) {
        f.var <- function(dtype, xx, yy, sd2.log, nn) {
          if (exists("track")) 
            print("f.var within f.data.plt.con")
          if (dtype %in% c(10, 15, 110, 250)) {
            SS <- sum(sd2.log * (nn - 1))
            df <- sum(nn - 1)
            var.within <- SS/df
          }
          else {
            if (dtype %in% c(1, 5)) 
              yy <- log(yy)
            if (dtype %in% c(26)) 
              yy <- sqrt(yy)
            mn <- tapply(yy, xx, mean)
            nn <- tapply(yy, xx, length)
            Vmn <- rep(mn[1], nn[1])
            if (dtype %in% c(1, 5)) 
              for (ii in 2:length(mn)) Vmn <- c(Vmn, 
                                                rep(mn[ii], nn[ii]))
            resid <- yy - Vmn
            df <- sum(nn - 1)
            if (df == 0) 
              cat("\nconf intervals cannot be calculated, df = 0\n")
            var.within <- (sum(nn) - 1) * var(resid)/df
          }
          Vsem <- sqrt(var.within/nn)
          out.lst <- list(Vsem = Vsem, df = df, var.within = var.within)
          return(out.lst)
        }
        var.out <- f.var(dtype, xx, yy, sd2.log, nn)
        Vsem <- var.out$Vsem
        df <- var.out$df
        Vconf <- qt(0.975, df) * Vsem
        if (dtype %in% c(1, 5, 10, 15)) {
          conf.L <- mean.y/exp(Vconf)
          conf.U <- mean.y * exp(Vconf)
        }
        if (dtype %in% c(25, 250)) {
          conf.L <- mean.y - Vconf
          conf.U <- mean.y + Vconf
        }
        if (dtype %in% c(26)) {
          conf.L <- (sqrt(mean.y) - Vconf)^2
          conf.U <- (sqrt(mean.y) + Vconf)^2
        }
        if (plot.type %in% c(3, 4, 6)) {
          conf.L <- log10(conf.L)
          conf.U <- log10(conf.U)
        }
        if (plot.type %in% c(7, 8)) {
          conf.L <- sqrt(conf.L)
          conf.U <- sqrt(conf.U)
        }
        y.lim.plt <- c(min(conf.L), max(conf.U))
        out.lst <- list(conf.L = conf.L, conf.U = conf.U, 
                        y.lim.plt = y.lim.plt, var.within = var.out$var, 
                        within)
        return(out.lst)
      }
    }
    detlim.plt.1 <- NA
    detlim.plt.2 <- NA
    if (max(Vdetlim, na.rm = T) > 0) {
      detlim.plt.1 <- min(Vdetlim[Vdetlim > 0], na.rm = T)
      detlim.plt.2 <- max(Vdetlim[Vdetlim > 0], na.rm = T)
    }
    x.plt <- NA
    y.plt <- NA
    switch(plot.type, {
      x.plt <- x
      y.plt <- yy
      xleg <- x.leg
      yleg <- y.leg
      x.lim.plt <- x.lim
      y.lim.plt <- y.lim
    }, {
      x.plt <- x
      x.plt[x.plt == 0] <- dum.contr
      x.plt <- log10(x.plt)
      y.plt <- yy
      if (x.leg != "") xleg <- paste("log10-", x.leg, sep = "") else xleg <- ""
      yleg <- y.leg
      y.lim.plt <- y.lim
      x.lim.plt <- log10(x.lim)
      if (!is.finite(x.lim.plt)[1]) x.lim.plt[1] <- log10(dum.contr)
    }, {
      x.plt <- x
      y.plt <- log10(yy)
      cens.up <- log10(cens.up)
      xleg <- x.leg
      if (y.leg != "") yleg <- paste("log10-", y.leg, sep = "") else yleg <- ""
      if (max(Vdetlim, na.rm = T) > 0) {
        detlim.plt.1 <- log10(detlim.plt.1)
        detlim.plt.2 <- log10(detlim.plt.2)
      }
      x.lim.plt <- x.lim
      y.lim.plt <- log10(y.lim)
      if (!is.finite(y.lim.plt)[1]) y.lim.plt[1] <- log10(dum.zero.resp)
    }, {
      x.plt <- x
      x.plt[x.plt == 0] <- dum.contr
      x.plt <- log10(x.plt)
      y.plt <- log10(yy)
      if (is.na(y.lim[1])) y.lim <- xy.lim[4:5]
      cens.up <- log10(cens.up)
      if (x.leg != "") xleg <- paste("log10-", x.leg, sep = "") else xleg <- ""
      if (y.leg != "") yleg <- paste("log10-", y.leg, sep = "") else yleg <- ""
      if (max(Vdetlim, na.rm = T) > 0) {
        detlim.plt.1 <- log10(detlim.plt.1)
        detlim.plt.2 <- log10(detlim.plt.2)
      }
      y.lim.plt <- log10(y.lim)
      x.lim.plt <- log10(x.lim)
      if (!is.finite(x.lim.plt)[1]) x.lim.plt[1] <- log10(dum.contr)
      if (!is.finite(y.lim.plt)[1]) y.lim.plt[1] <- log10(dum.zero.resp)
    }, {
      y.plt <- yy
      x.plt <- sqrt(x)
      if (x.leg != "") xleg <- paste("sqrt-", x.leg, sep = "") else xleg <- ""
      yleg <- y.leg
      y.lim.plt <- y.lim
      x.lim.plt <- sqrt(x.lim)
    }, {
      x.plt <- sqrt(x)
      y.plt <- log10(yy)
      cens.up <- log10(cens.up)
      if (is.na(y.lim[1])) y.lim <- xy.lim[4:5]
      if (y.leg != "") yleg <- paste("log10-", y.leg, sep = "") else yleg <- ""
      if (x.leg != "") xleg <- paste("sqrt-", x.leg, sep = "") else xleg <- ""
      if (max(Vdetlim, na.rm = T) > 0) {
        detlim.plt.1 <- log10(detlim.plt.1)
        detlim.plt.2 <- log10(detlim.plt.2)
      }
      x.lim.plt <- sqrt(x.lim)
      y.lim.plt <- log10(y.lim)
      if (!is.finite(y.lim.plt)[1]) y.lim.plt[1] <- log10(dum.zero.resp)
    }, {
      x.plt <- x
      y.plt <- sqrt(yy)
      cens.up <- sqrt(cens.up)
      if (is.na(y.lim[1])) y.lim <- xy.lim[4:5]
      xleg <- x.leg
      if (y.leg != "") yleg <- paste("sqrt-", y.leg, sep = "") else yleg <- ""
      if (max(Vdetlim, na.rm = T) > 0) {
        detlim.plt.1 <- sqrt(detlim.plt.1)
        detlim.plt.2 <- sqrt(detlim.plt.2)
      }
      x.lim.plt <- x.lim
      y.lim.plt <- sqrt(y.lim)
    }, {
      x.plt <- x
      x.plt[x.plt == 0] <- dum.contr
      x.plt <- log10(x.plt)
      y.plt <- sqrt(yy)
      cens.up <- sqrt(cens.up)
      if (is.na(y.lim[1])) y.lim <- xy.lim[4:5]
      if (y.leg != "") yleg <- paste("sqrt-", y.leg, sep = "") else yleg <- ""
      if (x.leg != "") xleg <- paste("log10-", x.leg, sep = "") else xleg <- ""
      if (max(Vdetlim, na.rm = T) > 0) {
        detlim.plt.1 <- sqrt(detlim.plt.1)
        detlim.plt.2 <- sqrt(detlim.plt.2)
      }
      y.lim.plt <- sqrt(y.lim)
      x.lim.plt <- log10(x.lim)
      if (!is.finite(x.lim.plt)[1]) x.lim.plt[1] <- log10(dum.contr)
    }, {
      ans.all$y.plt <- regr.resid
      ans.all$y.lim.plt <- c(min(ans.all$y.plt, na.rm = T), 
                             max(ans.all$y.plt, na.rm = T))
      ans.all$cex.1 <- cex.1
      if (exists("track")) print("f.data.plt.con:   END")
      return(ans.all)
    }, {
      y.plt <- regr.resid
      x.plt <- x
      x.plt[x.plt == 0] <- dum.contr
      x.plt <- log10(x.plt)
      x.plt <- x.plt[!is.na(y.plt)]
      if (x.leg != "") xleg <- paste("log10-", x.leg, sep = "") else xleg <- ""
      yleg <- "residual"
      x.lim.plt <- log10(x.lim)
      if (!is.finite(x.lim.plt)[1]) x.lim.plt[1] <- log10(dum.contr)
      y.lim.plt <- c(min(y.plt, na.rm = T), max(y.plt, 
                                                na.rm = T))
      if (plt.mns == 3) plt.mns <- 2
    }, {
      x.plt <- x
      if (x.leg != "") xleg <- paste("log10 of normalized ", 
                                     x.leg, sep = "") else xleg <- ""
                                     y.plt <- yy
                                     x.plt <- x.plt[!is.na(y.plt)]
                                     y.plt <- y.plt[!is.na(y.plt)]
                                     yleg <- "residual"
                                     yleg <- "residual"
                                     x.lim.plt <- log10(x.lim)
                                     if (!is.finite(x.lim.plt)[1]) x.lim.plt[1] <- log10(dum.contr)
                                     if (is.na(xy.lim[4])) y.lim <- c(min(y.plt, na.rm = T), 
                                                                      max(y.plt, na.rm = T))
                                     y.lim.plt <- c(min(y.plt, na.rm = T), max(y.plt, 
                                                                               na.rm = T))
    })
    shift.plt <- 0
    shift.tmp <- 0
    if (shift > 0) {
      shift.plt <- (max(x.plt, na.rm = T) - min(x.plt, 
                                                na.rm = T))/shift
      x.lim.plt[2] <- x.lim.plt[2] + shift.plt
    }
    ans.all$Vvar.within <- numeric()
    if (CI.plt) 
      for (jj in 1:nr.bb) for (ii in 1:nr.aa) {
        x.part <- x.plt[fct1 == ii & fct2 == jj]
        y.part <- y.plt[fct1 == ii & fct2 == jj]
        y.tmp <- yy[fct1 == ii & fct2 == jj]
        if (plot.type == 11) {
          x.part <- (x.part/max(x.part, na.rm = T))
          x.part[x.part == 0] <- dum.contr
          x.part <- log10(x.part)
        }
        if (length(y.part) > 0) {
          if (dtype %in% c(10, 15, 250)) {
            sd2.log.part <- sd2.log[fct1 == ii & fct2 == 
                                      jj]
            nn.part <- nn[fct1 == ii & fct2 == jj]
          }
          out.lst <- f.means(x.part, y.tmp, dtype, plot.type, 
                             CI = TRUE, sd2.log = sd2.log.part, nn = nn.part)
          ans.all$Vvar.within <- c(ans.all$Vvar.within, 
                                   out.lst$var.within)
          if (!is.na(out.lst$conf.L[1])) {
            y.lim.CI <- out.lst$y.lim.plt
            y.lim.plt[1] <- min(y.lim.CI[1], y.lim.plt[1], 
                                na.rm = T)
            y.lim.plt[2] <- max(y.lim.CI[2], y.lim.plt[2], 
                                na.rm = T)
          }
        }
      }
    if (nr.aa == 1 & nr.bb == 1 & nr.var > 1) {
      nr.bb <- nr.var
      fct2 <- fct3
    }
    if (nr.aa > 1 & nr.bb > 1 && sum(fct1 != fct2) == 0) {
      nr.aa <- 1
      fct1 <- rep(1, length(x))
    }
    points.plt.lst <- list()
    means.plt.lst <- list()
    if (gr.txt[1] == "") 
      element.txt <- "all"
    else element.txt <- gr.txt
    if (quick.ans > 2) {
      fct2 <- covariate
    }
    zz <- 0
    for (jj in 1:nr.bb) for (ii in 1:nr.aa) {
      x.part <- x.plt[fct1 == ii & fct2 == jj]
      y.part <- y.plt[fct1 == ii & fct2 == jj]
      y.tmp <- yy[fct1 == ii & fct2 == jj]
      if (plot.type == 11) {
        x.part <- (x.part/max(x.part, na.rm = T))
        x.part[x.part == 0] <- dum.contr
        x.part <- log10(x.part)
      }
      if (length(y.part) > 0) {
        zz <- zz + 1
        if (zz > 500) 
          cat("\nATTENTION: number of subgroups too large for plotting\n\n")
        x.part <- x.part + shift.tmp
        mean.x <- NA
        if (plt.mns %in% c(1, 3)) {
          if (dtype %in% c(1, 5, 15, 25, 26)) {
            mean.x <- as.numeric(tapply(x.part, x.part, 
                                        mean))
            mean.y <- f.means(x.part, y.tmp, dtype, plot.type)
          }
          if (dtype %in% c(10, 15, 250)) {
            mean.x <- x.part
            mean.y <- y.part
          }
          if (CI.plt) {
            if (dtype %in% c(10, 15, 250)) {
              sd2.log.part <- sd2.log[fct1 == ii & fct2 == 
                                        jj]
              nn.part <- nn[fct1 == ii & fct2 == jj]
            }
            out.lst <- f.means(x.part, y.tmp, dtype, 
                               plot.type, CI = TRUE, sd2.log = sd2.log.part, 
                               nn = nn.part)
            L025 <- out.lst$conf.L
            L975 <- out.lst$conf.U
          }
        }
        if (plt.mns == 2 || !CI.plt) {
          L025 <- rep(NA, length(mean.x))
          L975 <- rep(NA, length(mean.x))
        }
        points.plt.lst[[element.txt[zz]]] <- data.frame(x = x.part, 
                                                        y = y.part)
        if (plt.mns %in% c(1, 3)) 
          means.plt.lst[[element.txt[zz]]] <- data.frame(x = mean.x, 
                                                         y = mean.y, CI.low = L025, CI.upp = L975)
        shift.tmp <- shift.tmp + shift.plt
      }
    }
    ans.all$xy.plt <- data.frame(x.plt = x.plt, y.plt = y.plt)
    ans.all$points.plt.lst <- points.plt.lst
    ans.all$means.plt.lst <- means.plt.lst
    ans.all$cex.1 <- cex.1
    ans.all$cex.2 <- cex.2
    ans.all$xleg <- xleg
    ans.all$yleg <- yleg
    ans.all$x.lim.plt <- x.lim.plt
    ans.all$y.lim.plt <- y.lim.plt
    ans.all$detlim.plt.1 <- detlim.plt.1
    ans.all$detlim.plt.2 <- detlim.plt.2
    ans.all$gr.txt <- gr.txt
    ans.all$xy.lim <- signif(xy.lim, 3)
    if (exists("track")) 
      print("f.data.plt.con:  END")
    return(ans.all)
  })
}

#################
## f.lines.con ##
#################
f.lines.con <- function (ans.all) 
{
  if (exists("track")) 
    print("f.lines.con")
  with(ans.all, {
    lines.plt.lst <- f.lines.plt.con(ans.all)
    nr.gr <- length(lines.plt.lst)
    for (ii in 1:nr.gr) {
      dfr <- lines.plt.lst[[ii]]
      x.tmp <- dfr$xline
      y.tmp <- dfr$expect
      lines(x.tmp[1:2], y.tmp[1:2], col = color[ii], pch = mark[ii], 
            lty = 2)
      lines(x.tmp[-(1:2)], y.tmp[-(1:2)], col = color[ii], 
            pch = mark[ii], lty = l.ty)
    }
    if (exists("track")) 
      print("f.lines.con : END")
    return(lines.plt.lst)
  })
  if (0) {
    if (sum(is.na(regr.par)) > 0) {
      print("f.lines.con")
      cat("\nregr.par contains NAs\n")
      print(regr.par)
      return(invisible())
    }
    if (plot.type == 12) 
      return(invisible())
    count <- 0
    if ((model.ans != 6) && (nr.aa == 1) && (nr.bb == 1) && 
        (nr.var > 1)) 
      for (ii in 1:nr.var) {
        count <- count + 1
        ans.all.tmp$regr.par <- regr.par.matr[1, ]
        ans.all.tmp$col.tmp <- color[count]
        f.lines.tmp(ans.all.tmp)
      }
    else {
      kk.max <- length(regr.par.matr[, 1])
      par.bb.old <- NA
      par.bb.new <- 0
      for (kk in 1:kk.max) {
        ans.all.tmp$regr.par <- regr.par.matr[kk, ]
        if (model.ans > 1) 
          par.bb.new <- ans.all.tmp$regr.par[2]
        if (ans.plt == 1) {
          if (!identical(par.bb.new, par.bb.old) || (model.ans == 
                                                     1)) {
            count <- count + 1
            ans.all.tmp$col.tmp <- color[count]
          }
        }
        else {
          count <- count + 1
          ans.all.tmp$col.tmp <- color[count]
          if (ans.m6.sd == 2) 
            ans.all.tmp$sign.q <- sign.q[kk]
        }
        f.lines.tmp(ans.all.tmp)
        par.bb.old <- par.bb.new
      }
    }
    if (exists("track")) 
      print("f.lines.con : END")
    return(nr.gr)
  }
}

#####################
## f.lines.plt.con ##
#####################
f.lines.plt.con <- function (ans.all) 
{
  if (exists("track")) 
    print("f.lines.plt.con")
  with(ans.all, {
    if (sum(is.na(regr.par)) > 0) {
      print("f.lines.plt.con")
      cat("\nregr.par contains NAs\n")
      print(regr.par)
      return(invisible())
    }
    if (plot.type == 12) 
      return(invisible())
    f.lines.tmp <- function(ans.all.tmp) {
      dfr <- with(ans.all.tmp, {
        nbins <- 1000
        nbins.0 <- nbins/10
        if (plot.type %in% c(2, 4, 8)) 
          xline <- 10^(seq(from = log10(dum.contr), to = log10(max.x), 
                           length = nbins))
        else if (plot.type %in% 5:6) 
          xline <- (seq(from = sqrt(min.x), to = sqrt(max.x), 
                        length = nbins))^2
        else xline <- seq(from = min.x, to = max.x, length = nbins)
        twice <- F
        expect <- f.expect.con(model.ans, xline, regr.par, 
                               fct1 = rep(1, length(xline)), fct2 = rep(1, 
                                                                        length(xline)), fct3 = 1, CES = CES, twice = twice, 
                               trace = F, increase = increase, x.mn = x.mn, 
                               ref.lev = ref.lev, sign.q = sign.q, ans.m6.sd = ans.m6.sd)
        if (plot.type %in% c(2, 4, 8)) {
          if (min.x == 0) {
            xline <- log10(c(dum.contr, xline[nbins.0:nbins]))
            expectdum <- f.expect.con(model.ans, 0, regr.par, 
                                      fct1 = 1, fct2 = 1, fct3 = 1, CES = CES, 
                                      twice = twice, increase = increase, x.mn = x.mn, 
                                      ref.lev = ref.lev, sign.q = sign.q, ans.m6.sd = ans.m6.sd)
            expect <- c(expectdum, expect[nbins.0:nbins])
          }
          else xline <- log10(xline)
        }
        if (plot.type %in% 5:6) 
          xline <- sqrt(xline)
        if (plot.type %in% c(3, 4, 6)) 
          expect <- log10(expect)
        if (plot.type %in% 7:8) 
          expect <- sqrt(expect)
        dfr <- data.frame(xline = xline, expect = expect)
        if (plot.type %in% c(2, 4, 6) && min(x) == 0) 
          attr(dfr, "dummy.index") <- 1
        else attr(dfr, "dummy.index") <- 0
        return(dfr)
      })
    }
    ans.all.tmp <- ans.all
    ans.all.tmp$dum.contr <- xy.lim[1]
    ans.all.tmp$min.x <- xy.lim[2]
    ans.all.tmp$max.x <- xy.lim[3]
    if (1) 
      if (model.ans == 46) {
        ans.all.tmp$model.ans <- 15
        CED <- regr.par[(nr.aa + 1):(nr.aa + nr.bb)]
        CED.ref <- CED[ref.lev]
        CED <- CED.ref/CED
        CED[ref.lev] <- CED.ref
        regr.par[(nr.aa + 1):(nr.aa + nr.bb)] <- CED
        ans.all.tmp$regr.par <- regr.par
      }
    if (ans.plt == 1) {
      regr.par.matr <- matrix(regr.par, nrow = 1)
    }
    else {
      ans.all.tmp <- f.pars(ans.all.tmp)
      regr.par.matr <- ans.all.tmp$regr.par.matr
    }
    lines.plt.lst <- list()
    if ((model.ans != 6) && (nr.aa == 1) && (nr.bb == 1) && 
        (nr.var > 1)) {
      ans.all.tmp$regr.par <- regr.par.matr[1, ]
      for (ii in 1:nr.var) {
        dfr <- f.lines.tmp(ans.all.tmp)
        lines.plt.lst$all <- dfr
      }
    }
    else {
      kk.max <- length(regr.par.matr[, 1])
      for (kk in 1:kk.max) {
        ans.all.tmp$regr.par <- regr.par.matr[kk, ]
        if (ans.m6.sd == 2) 
          ans.all.tmp$sign.q <- sign.q[kk]
        dfr <- f.lines.tmp(ans.all.tmp)
        if (kk.max == 1) 
          lines.plt.lst$all <- dfr
        else lines.plt.lst[[gr.txt[kk]]] <- dfr
      }
    }
    if (exists("track")) 
      print("f.lines.plt.con : END")
    return(lines.plt.lst)
  })
}

##############################
## f.model.specific.results ##
##############################
f.model.specific.results <- function (ans.all) 
{
  if (exists("track")) 
    print("f.model.specific.results")
  sublist <- list()
  sublist$model.ans <- ans.all$model.ans
  sublist$model.type <- ans.all$model.type
  sublist$modelname <- ans.all$modelname
  sublist$MLE <- ans.all$MLE
  sublist$report.pars <- ans.all$report.pars
  if (length(ans.all$xans) > 1) {
    if (length(ans.all$xans) == 2) {
      RPF <- ans.all$regr.par[length(ans.all$regr.par)]
      sublist$x <- ans.all$x1 + RPF * ans.all$x2
    }
    if (length(ans.all$xans) == 3) {
      RPF1 <- ans.all$regr.par[length(ans.all$regr.par) - 
                                 1]
      RPF2 <- ans.all$regr.par[length(ans.all$regr.par)]
      RPF <- c(RPF1, RPF2)
      sublist$x <- ans.all$x1 + RPF[1] * ans.all$x2 + RPF[2] * 
        ans.all$x3
    }
    sublist$RPF <- RPF
    xy.lim <- ans.all$xy.lim
    xy.lim[2:3] <- c(min(sublist$x), max(sublist$x))
    sublist$xy.lim <- xy.lim
  }
  sublist$sd <- ans.all$sd
  sublist$par.start <- ans.all$MLE
  sublist$loglik.old <- ans.all$loglik
  sublist$nrp <- ans.all$nrp
  sublist$npar <- ans.all$npar
  sublist$npar.aic <- ans.all$npar.aic
  sublist$text.par <- ans.all$text.par
  sublist$lb <- ans.all$lb
  sublist$ub <- ans.all$ub
  sublist$loglik <- ans.all$loglik
  sublist$aic <- ans.all$aic
  sublist$converged <- ans.all$converged
  sublist$fct1 <- ans.all$fct1
  sublist$fct2 <- ans.all$fct2
  sublist$fct3 <- ans.all$fct3
  sublist$fct4 <- ans.all$fct4
  sublist$fct5 <- ans.all$fct5
  sublist$nr.aa <- max(ans.all$fct1)
  sublist$nr.bb <- max(ans.all$fct2)
  sublist$nr.var <- max(ans.all$fct3)
  sublist$nr.cc <- max(ans.all$fct4)
  sublist$nr.dd <- max(ans.all$fct5)
  sublist$nr.gr <- ans.all$nr.gr
  sublist$gr.txt <- ans.all$gr.txt
  sublist$factor.name <- ans.all$factor.name
  sublist$show <- ans.all$show
  if (ans.all$dtype == 5) {
    sublist$inter.var <- ans.all$inter.var
    sublist$intra.var <- ans.all$intra.var
  }
  model.ans <- ans.all$model.ans
  model.type <- ans.all$model.type
  if (!(model.ans == 14 & model.type == 1)) {
    sublist$regr.par <- ans.all$regr.par
    sublist$regr.par.matr <- ans.all$regr.par.matr
    sublist$MLE.all <- ans.all$MLE.all
    if (!ans.all$cont) {
      sublist$th.par <- ans.all$th.par
      sublist$sig.par <- ans.all$sig.par
      sublist$nth <- ans.all$nth
      sublist$sens.lev <- ans.all$sens.lev
      sublist$rank.low <- ans.all$rank.low
    }
    sublist$CED <- ans.all$CED
    sublist$CED.matr <- ans.all$CED.matr
    sublist$CED.boot <- ans.all$CED.boot
    sublist$CES <- ans.all$CES
    sublist$CES.0 <- ans.all$CES.0
    sublist$ced.table <- ans.all$ced.table
    sublist$loglik.summary <- ans.all$loglik.summary
    sublist$conf.int <- ans.all$conf.int
    if (ans.all$dtype == 3) 
      sublist$conf.int.matr <- ans.all$conf.int.matr
    if (!ans.all$cont) 
      sublist$pi.bg.matr <- ans.all$pi.bg.matr
    sublist$response.matr <- ans.all$response.matr
    sublist$trend <- ans.all$trend
    sublist$regr.resid <- ans.all$regr.resid
    sublist$lines.plt.lst <- ans.all$lines.plt.lst
    sublist$cedes.plt.lst <- ans.all$cedes.plt.lst
    sublist$report.pars <- ans.all$report.pars
    if (ans.all$dtype %in% c(1, 5)) 
      sublist$outliers <- ans.all$outliers
  }
  if (exists("track")) 
    print("f.model.specific.results:  END")
  return(sublist)
}

####################
## f.start.lm.con ##
####################
f.start.lm.con <- function (xtmp, ytmp, model.ans, dtype) 
{
  if (exists("track")) 
    print("f.start.lm.con")
  if (dtype %in% c(25, 250)) 
    dt.fr <- data.frame(yyy = ytmp, dose.tmp = xtmp)
  else if (dtype %in% c(26, 260)) 
    dt.fr <- data.frame(yyy = sqrt(ytmp), dose.tmp = xtmp)
  else dt.fr <- data.frame(yyy = log(ytmp), dose.tmp = xtmp)
  res.lm <- lm(yyy ~ dose.tmp, dt.fr)
  aa <- as.numeric(res.lm[[1]][1])
  if (dtype %in% c(1, 5, 10, 15)) 
    aa <- exp(aa)
  if (dtype %in% c(26, 260)) 
    aa <- aa^2
  bb <- as.numeric(res.lm[[1]][2])
  if (!(model.ans %in% 2:15)) {
    dt.fr <- data.frame(yyy = 1/ytmp, dose.tmp = xtmp)
    res.lm <- lm(yyy ~ dose.tmp, dt.fr)
    aa <- as.numeric(1/res.lm[[1]][1])
    bb <- as.numeric(1/res.lm[[1]][2]/aa)
  }
  if (bb < 0) 
    increase <- -1
  else increase <- 1
  if (model.ans %in% c(17:20, 22:25)) {
    if (bb < 0) 
      increase <- 1
    else increase <- -1
  }
  if (exists("track")) 
    print("f.start.lm.con:  END")
  return(c(aa, bb, increase))
}

###############
## f.ced.con ##
###############
f.ced.con <- function (model.ans, regr.par.matr, CES, nr.aa, nr.bb, nr.var, 
                       max.x = NA, ES.abs = 1, ans.m6.sd = 1, ED50 = FALSE) 
{
  if (exists("track")) 
    print("f.ced.con")
  if (0) {
    print("f.ced.con")
    print(ED50)
    print(model.ans)
    print(regr.par.matr)
    print(CES)
    print(nr.aa)
    print(nr.bb)
    print(nr.var)
  }
  CED <- numeric()
  if (model.ans == 16) 
    CED <- rep(0, nr.bb)
  else if (model.ans == 6 && ans.m6.sd == 2) 
    CED <- f.inv.con(model.ans, regr.par.matr[1, -1], CES, 
                     ans.m6.sd = 2)
  else if (model.ans == 47) 
    CED <- regr.par.matr[, 3]
  else {
    nr.row <- length(regr.par.matr[, 1])
    for (ii in 1:nr.row) {
      par.tmp <- regr.par.matr[ii, ]
      CED[ii] <- f.inv.con(model.ans, par.tmp, CES, max.x = max.x, 
                           ES.abs = ES.abs, ED50 = ED50)
    }
  }
  CED.lst <- list(CED = CED)
  if (exists("track")) 
    print("f.ced.con: END")
  return(CED.lst)
}

###############
## f.inv.con ##
###############
f.inv.con <- function (model.ans, params, CES, max.x = NA, ES.abs = 1, ans.m6.sd = 1, 
                       dtype = 1, ED50 = FALSE) 
{
  if (exists("track")) 
    print("f.inv.con")
  CED <- NA
  CES.tmp <- CES
  aa <- params[1]
  bb <- params[2]
  par3 <- params[3]
  par4 <- params[4]
  if (ED50) {
    CED <- bb * (log(2)^(1/par4))
    return(CED)
  }
  switch(model.ans, {
    CED <- NA
    cat("\nfor null model CED is not defined\n")
  }, {
    if (bb < 0) CES.tmp <- -abs(CES)
    CED <- (1/bb) * logb(CES.tmp + 1)
  }, {
    if (bb < 0) CES.tmp <- -abs(CES)
    CED <- (1/bb) * logb(CES.tmp + 1)
    CED <- CED^(1/par3)
  }, {
    if ((par3 < 1)) CES.tmp <- -abs(CES)
    dum <- (CES.tmp + 1 - par3)/(1 - par3)
    dum[dum < 0] <- NA
    if (dtype == 3) CED <- -(1/bb) * logb(dum) else {
      if (sum(is.na(dum)) > 0) {
        cat("\nATTTENTION: parameter c does not allow chosen value for CES\n")
        CED <- NA
      } else CED <- -(1/bb) * logb(dum)
    }
  }, {
    if ((par3 < 1)) CES.tmp <- -abs(CES)
    dum <- (CES.tmp + 1 - par3)/(1 - par3)
    dum[dum < 0] <- NA
    if (dtype == 3) CED <- (-(1/bb) * logb(dum))^(1/par4) else {
      if (sum(is.na(dum)) > 0) {
        cat("\nATTTENTION: parameter c does not allow chosen value for CES\n")
        CED <- NA
      } else CED <- (-(1/bb) * logb(dum))^(1/par4)
    }
  }, {
    if (ans.m6.sd == 1) CED <- bb * (-log(1 - log(CES + 1)/log(par3)))^(1/par4)
    if (ans.m6.sd == 2) CED <- bb * (-log(1 - CES))^(1/par4)
  }, {
    if (bb < 0) CES.tmp <- -abs(CES)
    CED <- bb * (logb(CES.tmp + 1))
  }, {
    if (bb > 0) CED <- bb * (logb(CES.tmp + 1))^(1/par3) else CED <- bb * 
        (-logb(-CES.tmp + 1))^(1/par3)
  }, {
    if (par3 < 1) CES.tmp <- -abs(CES)
    dum <- (CES.tmp + 1 - par3)/(1 - par3)
    dum[dum < 0] <- NA
    if (dtype == 3) CED <- bb * ((-logb(dum))) else {
      if (sum(is.na(dum)) > 0) {
        cat("\nATTTENTION: parameter c does not allow chosen value for CES\n")
        CED <- NA
      } else CED <- bb * ((-logb(dum)))
    }
  }, {
    if (par3 < 1) CES.tmp <- -abs(CES)
    dum <- (CES.tmp + 1 - par3)/(1 - par3)
    dum[dum < 0] <- NA
    if (dtype == 3) {
      CED <- (-(1/bb) * logb(dum))^(1/par4)
      CED <- bb * ((-logb(dum))^(1/par4))
    } else {
      if (sum(is.na(dum)) > 0) {
        cat("\nATTTENTION: parameter c does not allow chosen value for CES\n")
        CED <- NA
      } else {
        CED <- (-(1/bb) * logb(dum))^(1/par4)
        CED <- bb * ((-logb(dum))^(1/par4))
      }
    }
  }, , CED <- bb, CED <- bb, CED <- bb, CED <- bb, {
    if (aa < 0) CES.tmp <- -abs(CES)
    CED <- ((CES.tmp * par3)/aa)^(1/bb)
  }, {
    if (bb > 0) CES.tmp <- -abs(CES)
    CED <- -bb * (CES.tmp/(CES.tmp + 1))
  }, {
    if (bb > 0) CES.tmp <- -abs(CES)
    CED <- ((-sign(bb) * CES.tmp)/(1 + CES.tmp))^(1/par3) * 
      (abs(bb))
  }, {
    if ((par3 < 1)) CES.tmp <- -abs(CES)
    CED <- sign(bb) * bb * ((CES.tmp)/(par3 - 1 - CES.tmp))
  }, {
    if ((par3 < 1)) CES.tmp <- -abs(CES)
    if (par3 == 0) CED <- ((-(bb^par4) * CES.tmp)/(1 + CES.tmp))^(1/par4) else CED <- sign(bb) * 
        bb * ((CES.tmp)/(par3 - 1 - CES.tmp))^(1/par4)
  }, {
    f.dum <- function(x, params, CES) {
      b1 <- params[2]
      b2 <- params[3]
      c1 <- params[4]
      c2 <- params[5]
      dd <- params[6]
      y.expect.1 <- c1 - (c1 - 1) * exp(-(x/b1)^dd)
      y.expect.2 <- c2 - (c2 - 1) * exp(-(x/b2)^dd)
      y.expect <- (y.expect.1 * y.expect.2)
      return(y.expect - CES - 1)
    }
    CED <- NA
    uniroot.out <- uniroot(f.dum, interval = c(0, max.x), 
                           params = params, CES = CES)
    CED <- uniroot.out$root
  }, CED <- bb, CED <- bb, CED <- bb, CED <- bb, {
    if (ES.abs == 1) {
      CED <- ((CES.tmp - par3)/aa)^(1/bb)
      CED[CED <= 0] <- 1e-06
    } else {
      if (aa < 0) CES.tmp <- -abs(CES)
      CED <- ((par3 * CES.tmp)/aa)^(1/bb)
    }
  }, , , , {
    if (ES.abs == 1) {
      CED <- ((CES.tmp - aa)/bb)^(1/par3)
      CED[CED <= 0] <- 1e-06
    } else {
      if (bb < 0) CES.tmp <- -abs(CES)
      CED <- ((aa * CES.tmp)/bb)^(1/par3)
    }
  }, cat(""), cat(""), cat(""), cat(""), cat(""), cat(""), 
  cat(""), cat(""), cat(""), cat(""), {
    CED <- (CES.tmp/aa)^(1/bb)
  }, {
    cat("")
  }, cat(""), cat(""), cat(""), {
    CED <- bb
  }, {
  }, {
    CED <- bb
  })
  if (is.na(CED) && bb < 1e-10) 
    CED <- 1e+10
  if (exists("track")) 
    print("f.inv.con : END")
  return(CED)
}

####################
## f.move.sublist ##
####################
f.move.sublist <- function (lst, sub) 
{
  if (exists("track")) 
    print("f.move.sublist")
  lst$model.ans <- sub$model.ans
  if (length(sub$model.type) > 0) 
    lst$model.type <- sub$model.type
  lst$modelname <- sub$modelname
  lst$MLE <- sub$MLE
  lst$MLE.all <- sub$MLE.all
  if (length(sub$CED) > 0) {
    lst$CED <- sub$CED
    lst$CED.matr <- sub$CED.matr
  }
  lst$response.matr <- sub$response.matr
  lst$report.pars <- sub$report.pars
  if (length(lst$xans) > 1) {
    lst$RPF <- sub$RPF
    lst$x <- sub$x
    lst$xy.lim <- sub$xy.lim
  }
  if (lst$NES.ans == 2) 
    lst$CES <- sub$CES
  lst$par.start <- sub$MLE
  lst$loglik.old <- sub$loglik
  lst$regr.par <- sub$regr.par
  lst$regr.par.matr <- sub$regr.par.matr
  if (!lst$cont) {
    lst$th.par <- sub$th.par
    lst$sig.par <- sub$sig.par
    lst$nth <- sub$nth
    lst$sens.lev <- sub$sens.lev
    lst$rank.low <- sub$rank.low
  }
  if (length(sub$nrp) > 0) 
    lst$nrp <- sub$nrp
  lst$npar <- sub$npar
  lst$npar.aic <- sub$npar.aic
  lst$text.par <- sub$text.par
  lst$lb <- sub$lb
  lst$ub <- sub$ub
  lst$loglik <- sub$loglik
  lst$logliks <- sub$logliks
  lst$aic <- sub$aic
  lst$converged <- sub$converged
  lst$conf.int <- sub$conf.int
  if (lst$dtype == 3) 
    lst$conf.int.matr <- sub$conf.int.matr
  lst$fct1 <- sub$fct1
  lst$fct2 <- sub$fct2
  lst$fct3 <- sub$fct3
  lst$fct4 <- sub$fct4
  lst$fct5 <- sub$fct5
  lst$nr.aa <- max(sub$fct1)
  lst$nr.bb <- max(sub$fct2)
  lst$nr.var <- max(sub$fct3)
  lst$nr.cc <- max(sub$fct4)
  lst$nr.dd <- max(sub$fct5)
  lst$nr.gr <- sub$nr.gr
  lst$show <- sub$show
  lst$trend <- sub$trend
  lst$regr.resid <- sub$regr.resid
  if (lst$dtype == 5) {
    lst$inter.var <- sub$inter.var
    lst$intra.var <- sub$intra.var
  }
  lst$gr.txt <- sub$gr.txt
  lst$lines.plt.lst <- sub$lines.plt.lst
  lst$cedes.lst <- sub$cedes.lst
  if (lst$dtype %in% c(1, 5)) 
    lst$outliers <- sub$outliers
  if (exists("track")) 
    print("f.move.sublist:  END")
  return(lst)
}

#################
## f.resid.con ##
#################
f.resid.con <- function (ans.all) 
{
  if (exists("track")) 
    print("f.resid.con")
  with(ans.all, {
    pred.value <- f.expect.con(model.ans, x, regr.par = regr.par, 
                               fct1 = fct1, fct2 = fct2, fct3 = fct3, fct4 = fct4, 
                               fct5 = fct5, name = F, CES = CES, twice = twice, 
                               ttt = 0, y = yy, increase = increase, ref.lev = ref.lev, 
                               x.mn = x.mn, ans.m6.sd = ans.m6.sd, sign.q = sign.q)
    if (dtype %in% c(1, 5, 10, 15)) 
      regr.resid <- logb(yy) - logb(pred.value)
    if (dtype %in% c(25, 250)) 
      regr.resid <- pred.value - yy
    if (dtype %in% c(26, 260)) 
      regr.resid <- sqrt(yy) - sqrt(pred.value)
    ans.all$regr.resid <- regr.resid
    ans.all$pred.value <- pred.value
    if (exists("track")) 
      print("f.resid.con:  END")
    ans.all <- f.outliers(ans.all)
    return(ans.all)
  })
}

################
## f.outliers ##
################
f.outliers <- function (ans.all) 
{
  with(ans.all, {
    if (nr.var > 1) 
      cat("\nATTENTION: mean of estimated variances is used\n")
    if (!fitted && quick.ans == 1) 
      f.press.key.to.continue("\n first fit model\n")
    length.old <- length(regr.resid)
    sig2 <- mean(MLE[1:nr.var])
    mu <- mean(regr.resid)
    limit.tmp <- f.grubb(length.old)
    if (plot.type %in% c(10, 15)) {
      cat("\nATTENTION:  identifying outliers is not valid for summary data\n")
      f.press.key.to.continue()
    }
    #cat("\nThe Grubb outlier for this sample size is:", round(limit.tmp, 
    #                                                          5), "\n")
    limit.low <- -limit.tmp * sqrt(sig2) + mu
    limit.up <- limit.tmp * sqrt(sig2) + mu
    outlier.vec <- regr.resid[regr.resid < limit.low]
    outlier.vec <- c(outlier.vec, regr.resid[regr.resid > 
                                               limit.up])
    regr.resid.min <- regr.resid[(regr.resid > limit.low) & 
                                   (regr.resid < limit.up)]
    length.new <- length(regr.resid.min)
    if (length.old == length.new) {}
      #cat("\nNo outliers detected\n")
    ans.all$out.ans <- 2
    if (length.old > length.new) {
      lst <- (abs(regr.resid) >= min(abs(outlier.vec))) & 
        (y > -1000)
      #cat("\n", sum(lst), "outliers detected, with values for x and y:\n")
      numbering <- 1:length(lst)
      rrr <- numbering[lst]
      for (ii in rrr) cat(paste(x[ii], y[ii]), "\n")
      ans.all$outliers <- data.0[rrr, ]
      if (quick.ans == 1 && length(outlier.vec) > 0 && 
          length(xans) == 1) 
        ans.all$out.ans <- menu(c("yes", "no"), title = "\nDo you want to remove these outliers and refit the model?")
      if (ans.all$out.ans == 2 && ans.all$plot.type %in% 
          c(9:10)) 
        title(main = "outliers included")
      if (ans.all$out.ans == 1) {
        cat("\n\nThe following data are removed as outliers:\n")
        cat(paste(x[lst], y[lst]))
        ans.all$remove.no <- ans.all$yans
        ans.all$remove.lev <- ans.all$y[lst]
        if (ans.all$remove.name == "none") 
          ans.all$remove.name <- "outliers"
        else ans.all$remove.name <- c(ans.all$remove.name, 
                                      "outliers", sep = "\n")
        ans.all$y <- y[!lst]
        ans.all$yy <- yy[!lst]
        ans.all$xy.lim[4] <- min(ans.all$yy)
        ans.all$x <- x[!lst]
        ans.all$fct1 <- fct1[!lst]
        ans.all$fct2 <- fct2[!lst]
        ans.all$fct3 <- fct3[!lst]
        fct4 <- fct4[!lst]
        ans.all$fct4 <- fct4[!is.na(fct4)]
        fct5 <- fct5[!lst]
        ans.all$fct5 <- fct5[!is.na(fct5)]
        ans.all$regr.resid <- regr.resid[!lst]
        f.press.key.to.continue("\n\nNow, refit the model \n")
        return(ans.all)
      }
    }
    return(ans.all)
  })
}

#############
## f.grubb ##
#############
f.grubb <- function (ss = 25, alfa = 0.05) 
{
  alfa <- alfa/2
  t.crit <- qt(alfa/ss, ss - 2)
  G <- (ss - 1)/sqrt(ss)
  G <- G * sqrt(t.crit^2/(ss - 2 + t.crit^2))
  return(G)
}

##############
## f.CI.sel ##
##############
f.CI.sel <- function (ans.all) 
{
  if (exists("track")) 
    print("f.CI.sel")
  ans.all <- with(ans.all, {
    nr.aa <- max(fct1)
    nr.bb <- max(fct2)
    nr.var <- max(fct3)
    if (cont) {
      cc.OK <- f.check.cc(ans.all)
      if (!cc.OK) {
        note.tmp <- paste(y.leg, ": parameter cc too close to CES; no BMD CI calculated")
        ans.all$notes <- paste(ans.all$notes, "\n", note.tmp)
        ans.all$conf.int <- matrix(NA, ncol = 2, nrow = max(fct2))
        return(ans.all)
      }
      CED <- regr.par[(nr.aa + 1):(nr.aa + nr.bb)]
      from <- nr.var + nr.aa + 1
      until <- nr.var + nr.aa + nr.bb
      ans.all$group <- from:until
    }
    if (!cont) {
      CED.lst <- f.ced.cat(ans.all)
      ans.all$CED.matr <- CED.lst$CED.matr
      ans.all$response.matr <- CED.lst$response.matr
      CED <- CED.lst$CED
      rank.low <- order(CED)[1]
      if (max(fct1) > 1 && max(fct2) == 1) 
        group <- rank.low + 1
      else group <- rank.low + max(fct1)
      if (dtype == 6) 
        group <- group + 1
      if (length(xans) == 1) 
        ans.all$group <- group
      ans.all$rank.low <- rank.low
      if (model.ans > 1) 
        ans.all$sens.lev <- covar.txt[rank.low]
      else ans.all$sens.lev <- NA
    }
    if (dtype %in% c(2, 4, 6)) 
      no.CI <- T
    if (!no.CI) {
      if (dtype %in% c(5, 15)) {
        ans.all$plot.ans <- 1
        ans.all$nrp <- length(regr.par)
        ans.all$CED <- CED
        ans.all <- f.mm7.con(ans.all)
        conf.int <- ans.all$conf.int
      }
      else {
        if (ans.all$loglik == 1e-12) 
          cat("\nno finite log-likelihood found, CI will not be calculated\n")
        else {
          ans.all$trace <- F
          ans.all$trace.plt <- F
          if (length(xans) > 1) {
            group <- nr.var + nr.aa + 1
            if (length(xans) == 2) 
              group <- c(group, npar)
            if (length(xans) == 3) 
              group <- c(group, npar - 1, npar)
            ans.all$group <- group
          }
          if (ans.all$output) {}
            #cat("\n calculating confidence intervals ....\n")
          ans.all <- f.CI(ans.all)
          conf.int <- ans.all$conf.int
          if (ans.all$update) {
            if (cont) {
              CED <- ans.all$regr.par[(nr.aa + 1):(nr.aa + 
                                                     nr.bb)]
            }
            if (!cont) {
              CED.lst <- f.ced.cat(ans.all)
              ans.all$CED.matr <- CED.lst$CED.matr
              ans.all$response.matr <- CED.lst$response.matr
              CED <- CED.lst$CED
              if (nr.bb == 1) 
                CED <- rep(CED, nr.aa)
            }
            ans.all$regr.par.matr <- f.pars(ans.all)$regr.par.matr
          }
        }
      }
      if (model.ans > 1) {
        CED <- signif(CED, 2)
        conf.int <- signif(conf.int, 3)
        CED.uniq <- unique(CED)
        #if (cont && ans.all$output) 
        #  if (length(CED.uniq) > 1) 
        #    for (qq in 1:max(fct2)) cat("\nthe CED (in orig. units) and the", 
        #                                100 * conf.lev, "% confidence interval for group", 
        #                                covar.txt[qq], "is: \n", sf.x * CED[qq], 
        #                                "\n", conf.int[qq, 1], "\n", conf.int[qq, 
        #                                                                      2], "\n")
        #else cat("\nthe CED (in orig. units) and the", 
        #         100 * conf.lev, "% confidence interval is: \n", 
        #         sf.x * CED[1], "\n", conf.int[1, 1], "\n", 
        #         conf.int[1, 2], "\n")
      }
      if (length(CED) > 1) 
        ced.table <- data.frame(subgroup = covar.txt, 
                                BMDL = conf.int[, 1], BMDU = conf.int[, 2], 
                                BMD = CED)
      else ced.table <- data.frame(BMDL = conf.int[, 1], 
                                   BMDU = conf.int[, 2], BMD = CED)
      ans.all$ced.table <- ced.table
    }
    if (0) 
      if (!gui && length(ans.all$MLE) > 10) {
        f.store.results(ans.all, "refit.tmp")
        cat("\nNOTE: confidence intervals were stored in refit.tmp\n")
      }
    ans.all$CED <- signif(CED, 2)
    if (exists("track")) 
      print("f.CI.sel:  END")
    return(ans.all)
  })
}

################
## f.check.cc ##
################
f.check.cc <- function (ans.all) 
{
  with(ans.all, {
    if (exists("track")) 
      print("f.check.cc")
    cc.OK <- T
    cc <- MLE[nr.var + nr.aa + nr.bb + 1]
    if (model.ans %in% c(4, 5, 9, 10, 14, 15, 24, 25)) {
      if (increase == 1) 
        if (is.na(cc) | cc < 0.02 + (1 + CES)) {
          #cat("\n\nthe value of parameter c is too close to CES\n                  this indicates (nonrandom) errors in the data,\n                  or too high value of CES")
          cc.OK <- F
          return(cc.OK)
        }
      if (increase == -1) 
        if (is.na(cc) | cc > -0.02 + (1 - abs(CES))) {
          #cat("\n\nthe value of parameter c is too close to CES\n                  this indicates (nonrandom) errors in the data,\n                  or too high value of CES")
          cc.OK <- F
          return(cc.OK)
        }
    }
    return(cc.OK)
  })
}

##########
## f.CI ##
##########
f.CI <- function (ans.all, nolog = F) 
{
  if (exists("track")) 
    print("f.CI")
  with(ans.all, {
    if (quick.ans == 1) 
      ans.all$trace.plt <- T
    if (length(ans.all$MLE) > 25) {
      cat("\nnumber of parameters exceeds 25\n")
      ans.all$change[nrQ] <- T
      ans.all$change[nrQ + 1] <- F
      ans.all <- f.change.settings(ans.all)
    }
    profile.out <- f.profile.all(ans.all, nolog)
    update <- F
    count <- 0
    loglik.tmp <- profile.out$loglik
    while (sum(profile.out$MLE.new) != 0) {
      count <- count + 1
      if (interrupt) {
        cat(" \n.... calculation of CI is re-started\n")
        cat(" improved fit will be stored in R-object called newfit\n")
      }
      ans.all$MLE <- profile.out$MLE.new
      if (cont) 
        ans.all$regr.par <- ans.all$MLE[-(1:max(fct3))]
      ans.all$loglik <- profile.out$loglik
      if (ans.all$loglik < loglik.tmp - 0.02) {
        cat("\nATTENTION: log-likelihood fluctuates, CI is not calculated\n")
        profile.out$conf.int <- matrix(NA, ncol = 2)
        profile.out$MLE.new <- 0
      }
      else {
        loglik.tmp <- ans.all$loglik
        cat("\n new log-likelihood: ", ans.all$loglik, 
            "\n")
        f.store.results(ans.all, "newfit")
        profile.out <- f.profile.all(ans.all, nolog)
        update <- T
        if (count > 50) {
          cat("\nno global optimum found\n")
          profile.out$MLE.new <- 0
          profile.out$conf.int <- matrix(NA, ncol = 2)
        }
      }
    }
    if (update) {
      if (!cont) {
        CED <- ans.all$MLE[(alfa.length + nr.aa):(alfa.length + 
                                                    nr.aa + nr.bb)]
        CED.matr <- matrix(CED, ncol = 1)
        if (model.type == 1) {
          ans.all$regr.par <- ans.all$MLE[1:nrp]
          if (dtype == 6) 
            ans.all$regr.par <- ans.all$regr.par[-1]
        }
        if (model.type == 2) {
          par.lst <- f.split.par(ans.all$MLE, nrp, nth, 
                                 dtype)
          ans.all$regr.par <- par.lst$regr.par
          ans.all$th.par <- par.lst$th.par
          ans.all$sig.par <- par.lst$sig.par
        }
        ans.all$CED <- signif(CED, 2)
        ans.all$CED.matr <- CED.matr
      }
    }
    if (nolog) 
      conf.int <- signif(profile.out$conf.int, 3)
    else conf.int <- signif(10^profile.out$conf.int, 3)
    low <- conf.int[, 1]
    upp <- conf.int[, 2]
    lst1 <- is.na(low)
    lst2 <- is.na(upp)
    low[lst1] <- 0
    upp[lst2] <- Inf
    conf.int[, 1] <- low
    conf.int[, 2] <- upp
    lst <- lst1 & lst2
    conf.int[lst, ] <- rep(NA, 2)
    ans.all$conf.int <- conf.int
    ans.all$update <- update
    ans.all$profile <- profile.out$profile
    if (quick.ans == 1) {
      cat("\n\n\nCheck the log-likelihood plot")
    }
    if (exists("track")) 
      print("f.CI END")
    return(ans.all)
  })
}

###################
## f.profile.all ##
###################
f.profile.all <- function (ans.all, nolog = F) 
{
  if (exists("track")) 
    print("f.profile.all")
  date.start <- date()
  with(ans.all, {
    max.runs <- 200
    count <- 0
    count.2 <- 0
    CI.NA <- FALSE
    count.max <- 100
    crit <- 0.5 * qchisq(conf.lev, 1)
    jump.crit <- crit * 3
    dist <- 1.01
    par.nr <- 0
    loglik.max <- loglik
    lb.orig <- lb
    ub.orig <- ub
    profile.out <- list()
    MLE.new <- 0
    tb <- "\t"
    if (dtype == 3) {
      trace.plt <- F
      if (0) 
        f.graph.window(nr.bb, nr.gr = nr.gr)
    }
    if (cont) {
      sig2 <- mean(MLE[1:nr.var])
      if (group[1] == 0) {
        if (model.ans == 16) {
          from <- nr.var + nr.aa + nr.bb + 2
          until <- from
        }
        else {
          from <- nr.var + nr.aa + 1
          until <- nr.var + nr.aa + nr.bb
        }
        group <- from:until
      }
      displ.fact <- 1
      if (nr.var > 1) 
        displ.fact <- fct3
      if (nr.bb > 1) 
        displ.fact <- fct2
      if (nr.aa > 1) 
        displ.fact <- fct1
      if (dtype == 5) {
        cat("\nATTENTION: confidence intervals are calculated without accounting for \n             nested structure in the data (e.g. intralitter correlations) \n")
        f.press.key.to.continue()
      }
    }
    if (!cont) {
      if (dtype %in% c(4, 6)) 
        y <- kk/nn
      if (group[1] == 0) {
        from <- nr.aa + 1
        until <- nr.aa + nr.bb
        if (model.type == 2) {
          if (nr.aa == 1) {
            until <- until - nr.bb + 1
          }
          if (nr.aa > 1 && nr.bb == 1 && ces.ans != 4) {
            from <- 2
            until <- nr.aa + 1
          }
          if (nr.aa == 1 & nr.bb > 1) {
            from <- 2
            until <- nr.bb + 1
          }
          if (nr.aa > 1 & nr.bb > 1) {
            from <- nr.aa + 1
            until <- nr.aa + nr.bb
          }
        }
        if (model.ans == 30 & model.type == 1) {
          from <- nr.aa + nr.bb + 1
          until <- nr.aa + nr.bb + 1
        }
        if (dtype == 6) {
          from <- from + 1
          until <- until + 1
        }
        group <- from:until
      }
      displ.fact <- 1
      if (nr.bb > 1) 
        displ.fact <- fct2
      if (nr.aa > 1) 
        displ.fact <- fct1
    }
    conf.int <- matrix(NA, nrow = length(group), ncol = 2)
    for (jj in group) {
      par.nr <- par.nr + 1
      if (trace) {
        if (model.type == 2 && nr.aa > 1 && nr.bb == 
            1) 
          cat("\ncalculating C.I.for group", fct1.txt[par.nr], 
              "......\n")
        else if (nr.bb == 1 || quick.ans == 2) 
          cat("\nCalculating C.I. ......\n")
        else if (group[1] == 0) 
          cat("\nCalculating C.I.for group", fct2.txt[par.nr], 
              "......\n")
        else cat("\nCalculating C.I.for paramter", jj, 
                 "......\n")
      }
      CED.upp.inf <- FALSE
      small.step.low <- F
      small.step.upp <- F
      if (trace) {
        cat("\n=========== lower limit =================\n")
        cat(text.par, "loglik", "loglik boundary", sep = tb)
        cat("\n")
        cat(signif(MLE, 3), round(loglik.max, 2), round(loglik.max - 
                                                          crit, 2), sep = tb)
      }
      lb <- lb.orig
      ub <- ub.orig
      lb[jj] <- MLE[jj]
      # if (lb[jj] == lb.orig[jj]) {
      #   cat("\n\nATTENTION:  MLE of this parameter is equal to its lower constraint (", 
      #       lb.orig[jj], ")\nrefit the model with another lower constraint\n\n")
      #   return(ans.all)
      # }
      # if (lb[jj] == ub.orig[jj]) {
      #   cat("\n\nATTENTION:  MLE of this parameter is equal to its upper constraint (", 
      #       ub.orig[jj], ")\nrefit the model with another upper constraint\n\n")
      #   return(ans.all)
      # }
      if (lb[jj] < 0) 
        nolog <- T
      start <- dist * MLE
      loglik.low <- loglik.max
      CED.low <- MLE[jj]
      CED.low.old <- numeric()
      loglik.low.old <- numeric()
      loglik.current <- loglik.max
      MLE.current <- MLE
      if (trace.plt) {
        ans.all.plt <- ans.all
        ans.all.plt$color <- rep(1, length(x))
        ans.all.plt$color[jj] <- 3
        f.graph.window(2, nr.gr = nr.gr)
        if (!cont & dtype != 3) {
          f.plot.frq(ans.all.plt)
          par.lst <- f.split.par(MLE, nrp, nth, dtype)
          ans.all.plt$regr.par <- par.lst$regr.par
          f.lines.frq(ans.all.plt)
        }
        else {
          ans.all.plt$plt.mns <- 3
          ans.all.plt$cex.2 <- 2.5
          ignore <- f.plot.con(ans.all.plt)
        }
      }
      if (cont && dtype %in% c(1, 5, 10)) 
        step.start <- 1.02 + 0.11 * sig2
      else step.start <- 1.08
      step <- step.start
      stop <- F
      run <- 0
      while (stop == F) {
        run <- run + 1
        lb[jj] <- lb[jj]/step
        ub[jj] <- lb[jj]
        start[jj] <- lb[jj]
        step.tmp <- step
        # while (lb[jj] < lb.orig[jj]) {
        #   cat("\nf.profile.all: value of parameter exceeds lower constraint, attempting to avoid this\n")
        #   step.tmp <- sqrt(step.tmp)
        #   lb[jj] <- lb.orig[jj] * step.tmp
        #   ub[jj] <- lb[jj]
        #   start[jj] <- lb[jj]
        #   step <- step.tmp
        # }
        start[start < lb] <- lb[start < lb]
        start[start > ub] <- ub[start > ub]
        ans.all$par.start <- start
        ans.all$lb <- lb
        ans.all$ub <- ub
        if (cont) 
          loglik.try <- -f.lik.con(start, x, y, dtype, 
                                   fct1, fct2, fct3, model.ans, mn.log, sd2.log, 
                                   nn, Vdetlim = Vdetlim, CES = CES, twice = twice, 
                                   ttt = ttt, trace = F, fct4 = fct4, fct5 = fct5, 
                                   cens.up = cens.up, par.tmp = NA, increase = increase, 
                                   x.mn = x.mn, ref.lev = ref.lev, sign.q = sign.q, 
                                   ans.m6.sd = ans.m6.sd, x1 = x1, x2 = x2)
        else loglik.try <- 0
        if (!is.finite(loglik.try)) {
          lb[jj] <- lb[jj] * step
          step <- sqrt(step)
          cat("\ntemporary fitting problem\n")
          start <- start * runif(length(start), 0.9, 
                                 1.1)
        }
        else {
          nlminb.out <- f.nlminb(ans.all, tmp.quick = T)
          MLE.current <- nlminb.out$MLE
          loglik.current <- nlminb.out$loglik
          start <- MLE.current * dist
          if (!is.finite(loglik.current) || loglik.current <= 
              -1e+10) {
            cat("\nf.profile.all:  bad fit, new try with adjusted start values \n")
            start <- dist * nlminb.out$MLE
            loglik.current <- -Inf
          }
          if (loglik.current > loglik.max + 0.03) {
            ans.all$par.start <- MLE.current * dist
            cat("\nlocal optimum found....\n")
            profile.out <- list()
            profile.out$MLE.new <- nlminb.out$MLE
            profile.out$loglik <- loglik.current
            profile.out$conf.int <- conf.int
            return(profile.out)
          }
          if (length(loglik.low) <= 2 || run > max.runs) {
            loglik.low <- c(loglik.current, loglik.low)
            CED.low <- c(lb[jj], CED.low)
          }
          if (length(loglik.low) > 2) {
            if (loglik.current > loglik.low[2]) {
              loglik.low <- c(loglik.current, loglik.low[-(1:2)])
              CED.low <- c(lb[jj], CED.low[-(1:2)])
            }
            else if (loglik.current > loglik.low[1]) {
              loglik.low <- c(loglik.current, loglik.low[-1])
              CED.low <- c(lb[jj], CED.low[-1])
            }
            else if (loglik.current <= loglik.low[1]) {
              loglik.low <- c(loglik.current, loglik.low)
              CED.low <- c(lb[jj], CED.low)
            }
          }
          if (loglik.current == -Inf) {
            loglik.low <- loglik.low[-1]
            CED.low <- CED.low[-1]
          }
          if (length(loglik.low) > 1) {
            if (abs(loglik.low[1] - loglik.low[2]) > 
                0.25 * crit) 
              step <- sqrt(step)
            if (abs(loglik.low[1] - loglik.low[2]) < 
                0.1 * crit) 
              step <- step^1.26
          }
          if (trace) {
            cat("\n")
            cat(signif(MLE.current, 3), round(loglik.current, 
                                              2), round(loglik.max - crit, 2), sep = tb)
          }
          if (trace.plt & (loglik.current > (loglik.max - 
                                             crit))) {
            if (!cont & dtype != 3) {
              par.lst <- f.split.par(MLE.current, nrp, 
                                     nth, dtype)
              ans.all.plt$regr.par <- par.lst$regr.par
              f.lines.frq(ans.all.plt)
            }
            else {
              if (ans.m6.sd != 2) 
                ans.all.plt$regr.par <- MLE.current[-c(1:nr.var)]
              f.lines.con(ans.all.plt)
            }
          }
        }
        if (step < 1.00001) {
          stop <- T
          small.step.low <- T
        }
        # if (abs(CED.low[1]) < 1e-20) {
        #   cat("\ncurrent parameter value is too small:", 
        #       CED.low[1], "\n")
        #   stop <- T
        # }
        if (length(loglik.low) > 30) 
          if (max(abs(diff(loglik.low[1:10]))) < 0.001) {
            stop <- T
          }
        if (loglik.max - loglik.current > crit) {
          if (length(unique(loglik.low)) < 6) {
            if (trace) 
              cat("\n------- not enough points: ", length(unique(loglik.low)), 
                  " ------------------")
            CED.low <- c(CED.low, CED.low.old)
            CED.low.old <- CED.low
            loglik.low <- c(loglik.low, loglik.low.old)
            loglik.low.old <- loglik.low
            power <- 1.1
            step <- step.start^power
            step.start <- step
            lb[jj] <- MLE[jj]
            start <- dist * MLE
            run <- 0
            count.2 <- count.2 + 1
            if (count.2 > 10) {
              cat("\n did not succeed to collect large enough number of points of the profile\n")
              CI.NA <- TRUE
              stop <- T
              loglik.low <- NA
            }
          }
          else if ((loglik.max - loglik.current) > jump.crit) {
            lb[jj] <- lb[jj] * step
            step <- step.start^0.5
            step.start <- step
            count <- count + 1
            if (step.start < 1.01) {
              CED.low[1] <- lb[jj]
              loglik.low[1] <- loglik.current
              count <- count.max
            }
            if (count == count.max) {
              stop <- T
            }
            run <- 0
            if (0) {
              MLE.current[jj] <- lb[jj]
              par.lst <- f.split.par(MLE.current, nrp, 
                                     nth, dtype)
              ans.all.plt$regr.par <- par.lst$regr.par
              f.lines.frq(ans.all.plt)
              f.press.key.to.continue()
            }
          }
          else {
            loglik.diff <- diff(loglik.low[order(CED.low)])
            if (sum(loglik.diff < -0.1) > length(CED.low) - 
                3) {
              step.start <- 1 + (step.start - 1)/2
              step <- step.start
              run <- 0
            }
            else {
              stop <- T
            }
          }
        }
        if (run > max.runs) {
          cat("\nATTENTION: Maximum number of runs reached in establishing profile\n")
          if (abs(loglik.low[1] - loglik.low[2]) < 0.01) 
            stop <- T
        }
        jump.low <- F
        if (count == count.max) {
          #cat("\nATTENTION:  jump in log-likelihood\n")
          jump.low <- TRUE
        }
      }
      if (trace) {
        cat("\n\n=========== upper limit =================\n")
        cat(text.par, "loglik", "loglik boundary", sep = tb)
        cat("\n")
        cat(signif(MLE, 3), round(loglik.max, 2), round(loglik.max - 
                                                          crit, 2), sep = tb)
      }
      if (MLE.new[1] != 0) 
        MLE <- MLE.new
      lb <- lb.orig
      ub <- ub.orig
      lb[jj] <- MLE[jj]
      start <- dist * MLE
      loglik.upp <- loglik.max
      CED.upp <- MLE[jj]
      CED.upp.old <- numeric()
      loglik.upp.old <- numeric()
      if (cont && dtype %in% c(1, 5, 10)) 
        step.start <- 1.02 + 0.11 * sig2
      else step.start <- 1.08
      step <- step.start
      run <- 0
      stop <- F
      count <- 0
      count.2 <- 0
      while (stop == F) {
        run <- run + 1
        lb[jj] <- lb[jj] * step
        ub[jj] <- lb[jj]
        start[jj] <- lb[jj]
        step.tmp <- step
        # while (lb[jj] > ub.orig[jj]) {
        #   cat("f.profile.all: value of parameter exceeds upper constraint, attempting to avoid this\n")
        #   step.tmp <- sqrt(step.tmp)
        #   lb[jj] <- lb[jj]/step.tmp
        #   ub[jj] <- lb[jj]
        #   start[jj] <- lb[jj]
        #   step <- step.tmp
        # }
        start[start < lb] <- lb[start < lb]
        start[start > ub] <- ub[start > ub]
        ans.all$par.start <- start
        ans.all$lb <- lb
        ans.all$ub <- ub
        if (cont) 
          loglik.try <- -f.lik.con(ans.all$par.start, 
                                   x, y, dtype, fct1, fct2, fct3, model.ans, 
                                   mn.log, sd2.log, nn, Vdetlim = Vdetlim, CES = CES, 
                                   twice = twice, ttt = ttt, trace = F, fct4 = fct4, 
                                   fct5 = fct5, cens.up = cens.up, par.tmp = NA, 
                                   increase = increase, x.mn = x.mn, ref.lev = ref.lev, 
                                   sign.q = sign.q, ans.m6.sd = ans.m6.sd, x1 = x1, 
                                   x2 = x2)
        else loglik.try <- 0
        if (!is.finite(loglik.try)) {
          lb[jj] <- lb[jj]/step
          step <- sqrt(step)
          cat("\ntemporary fitting problem\n")
          start <- start * runif(length(start), 0.9, 
                                 1.1)
        }
        else {
          nlminb.out <- f.nlminb(ans.all, tmp.quick = T)
          MLE.current <- nlminb.out$MLE
          loglik.current <- nlminb.out$loglik
          start <- MLE.current * dist
          if (!is.finite(loglik.current) || loglik.current <= 
              -1e+10) {
            start <- dist * nlminb.out$MLE
            cat("\nf.profile.all:  bad fit, new try with adjusted start values \n")
            loglik.current <- -Inf
          }
          if (loglik.current > loglik.max + 0.03) {
            ans.all$par.start <- MLE.current * dist
            cat("\nlocal optimum found....\n")
            profile.out <- list()
            profile.out$MLE.new <- nlminb.out$MLE
            profile.out$loglik <- loglik.current
            profile.out$conf.int <- conf.int
            return(profile.out)
          }
          if (length(loglik.upp) <= 2 || (run > max.runs)) {
            loglik.upp <- c(loglik.current, loglik.upp)
            CED.upp <- c(ub[jj], CED.upp)
          }
          if (length(loglik.upp) > 2) {
            if (loglik.current > loglik.upp[2]) {
              loglik.upp <- c(loglik.current, loglik.upp[-(1:2)])
              CED.upp <- c(ub[jj], CED.upp[-(1:2)])
            }
            else if (loglik.current > loglik.upp[1]) {
              loglik.upp <- c(loglik.current, loglik.upp[-1])
              CED.upp <- c(ub[jj], CED.upp[-1])
            }
            else if (loglik.current <= loglik.upp[1]) {
              loglik.upp <- c(loglik.current, loglik.upp)
              CED.upp <- c(ub[jj], CED.upp)
            }
          }
          if (loglik.current == -Inf) {
            loglik.upp <- loglik.upp[-1]
            CED.upp <- CED.upp[-1]
          }
          if (length(loglik.upp) > 1) {
            if (abs(loglik.upp[1] - loglik.upp[2]) > 
                0.25 * crit) 
              step <- sqrt(step)
            if (abs(loglik.upp[1] - loglik.upp[2]) < 
                0.1 * crit) 
              step <- step^1.26
          }
          if (0) {
            if ((run <= max.runs) & loglik.current <= 
                (loglik.upp[1] - 5)) {
              lb[jj] <- lb[jj]/step
              step <- 1.01
              loglik.current <- loglik.upp[1]
            }
            else {
              if (abs(loglik.upp[1] - loglik.upp[2]) > 
                  0.25 * crit) {
                step <- sqrt(step)
              }
              if (abs(loglik.upp[1] - loglik.upp[2]) < 
                  0.1 * crit) {
                step <- step^1.26
              }
              if (loglik.current == -Inf) {
                loglik.upp <- loglik.upp[-1]
                CED.upp <- CED.upp[-1]
              }
            }
          }
          if (trace) {
            cat("\n")
            cat(signif(MLE.current, 3), round(loglik.current, 
                                              2), round(loglik.max - crit, 2), sep = tb)
          }
          if (trace.plt & (loglik.current > (loglik.max - 
                                             crit))) {
            if (!cont & dtype != 3) {
              par.lst <- f.split.par(MLE.current, nrp, 
                                     nth, dtype)
              ans.all.plt$regr.par <- par.lst$regr.par
              f.lines.frq(ans.all.plt)
            }
            else {
              if (ans.m6.sd != 2) 
                ans.all.plt$regr.par <- MLE.current[-c(1:nr.var)]
              f.lines.con(ans.all.plt)
            }
          }
        }
        if (step < 1.00001) {
          stop <- T
          small.step.upp <- T
          cat("\nstopped due to small step size\n")
        }
        if (length(loglik.upp) > 30) 
          if (max(abs(diff(loglik.upp[1:10]))) < 0.001) 
            stop <- T
        if (CED.upp[1] > 1e+20) 
          stop <- T
        if (ub[jj] > 1e+50) {
          stop <- T
        }
        if (loglik.max - loglik.current > crit) {
          if (length(unique(loglik.upp)) < 6) {
            if (trace) 
              cat("\n------- not enough points: ", length(unique(loglik.upp)), 
                  " ------------------")
            CED.upp <- c(CED.upp, CED.upp.old)
            loglik.upp <- c(loglik.upp, loglik.upp.old)
            CED.upp.old <- CED.upp
            loglik.upp.old <- loglik.upp
            power <- 0.9
            step <- step.start^power
            step.start <- step
            lb[jj] <- MLE[jj]
            start <- dist * MLE
            run <- 0
            count.2 <- count.2 + 1
            if (count.2 > 10) {
              cat("\n did not succeed to collect large enough number of points of the profile\n")
              CI.NA <- TRUE
              stop <- T
              loglik.upp <- NA
              CED.upp <- NA
            }
          }
          else if ((loglik.max - loglik.current) > jump.crit) {
            lb[jj] <- lb[jj]/step
            step <- step.start^0.5
            step.start <- step
            count <- count + 1
            if (step.start < 1.01) {
              CED.upp[1] <- lb[jj]
              loglik.upp[1] <- loglik.current
              count <- count.max
            }
            if (count == count.max) 
              stop <- T
            run <- 0
            if (0) {
              MLE.current[jj] <- lb[jj]
              par.lst <- f.split.par(MLE.current, nrp, 
                                     nth, dtype)
              ans.all.plt$regr.par <- par.lst$regr.par
              f.lines.frq(ans.all.plt)
              f.press.key.to.continue()
            }
          }
          else {
            loglik.diff <- diff(loglik.upp[order(CED.upp)])
            if (sum(loglik.diff < -0.1) > length(CED.upp) - 
                3) {
              step.start <- 1 + (step.start - 1)/2
              step <- step.start
              run <- 0
            }
            else stop <- T
          }
        }
        if (run > max.runs) {
          cat("\nATTENTION: Maximum number of runs reached in establishing profile\n")
          if (abs(loglik.upp[1] - loglik.upp[2]) < 0.01) 
            stop <- T
        }
        if (0) 
          if (model.ans %in% c(14, 15, 24, 25)) {
            cc <- MLE.current[nr.var + nr.aa + nr.bb + 
                                1]
            if (increase == 1) 
              if (cc <= 1 + CES) {
                print("f.profile.all:   c does not reach CES")
                print(MLE.current)
                print(cc)
                stop <- T
                CED.upp.inf <- TRUE
              }
            if (increase == -1) 
              if (cc >= 1 - abs(CES)) {
                print("f.profile.all decrease:   c does not reach CES")
                print(MLE.current)
                print(cc)
                print(1 - abs(CES))
                stop <- T
                CED.upp.inf <- TRUE
              }
          }
        jump.upp <- FALSE
        if (count == count.max) {
          #cat("\nATTENTION:  jump in log-likelihood\n")
          jump.upp <- TRUE
        }
      }
      if (MLE[jj] < 0) {
        loglik.low.neg <- loglik.low
        CED.low.neg <- CED.low
        loglik.low <- loglik.upp
        CED.low <- CED.upp
        loglik.upp <- loglik.low.neg
        CED.upp <- CED.low.neg
      }
      if (!nolog) 
        CED.low <- log10(CED.low)
      if (max(abs(diff(loglik.low))) > 0.01) {
        loglik.low <- loglik.low[order(CED.low)]
        CED.low <- sort(CED.low)
        #stop -Inf
        CED.low[CED.low == -Inf] <- min(CED.low[CED.low != -Inf], na.rm=TRUE) - 10;
        spline.low <- spline(CED.low, loglik.low)
        if (min(spline.low$x) > min(CED.low)) {
          spline.low$x <- c(min(CED.low), spline.low$x)
          spline.low$y <- c(loglik.low[1], spline.low$y)
        }
        if (length(unique(spline.low$y)) == 1) {
          cat("\n\nonly one point in lower spline\n")
          ci.low <- list(y = -Inf)
          conf.int[par.nr, 1] <- ci.low$y
        }
        else {
          ci.low <- approx(spline.low$y, spline.low$x, 
                           xout = loglik.max - crit)
          conf.int[par.nr, 1] <- ci.low$y
        }
        if (small.step.low) 
          conf.int[par.nr, 1] <- (CED.low[1])
        if (loglik.low[1] > loglik.max - crit) 
          conf.int[par.nr, 1] <- -Inf
      }
      else {
        spline.low <- list()
        spline.low$x <- CED.low
        spline.low$y <- rep(loglik.low[1], length(CED.low))
        ci.low <- list()
        ci.low$y <- min(CED.low)
      }
      if (!nolog) 
        CED.upp <- log10(CED.upp)
      if (max(abs(diff(loglik.upp))) > 0.01) {
        loglik.upp <- loglik.upp[order(CED.upp)]
        CED.upp <- sort(CED.upp)
        spline.upp <- spline(CED.upp, loglik.upp)
        if (max(spline.upp$x) < max(CED.upp)) {
          spline.upp$x <- c(spline.upp$x, max(CED.upp))
          spline.upp$y <- c(spline.upp$y, loglik.upp[length(loglik.upp)])
        }
        if (length(unique(spline.upp$y)) == 1) {
          cat("\n\nonly one point in upper spline\n")
          ci.upp <- list(y = Inf)
          conf.int[par.nr, 2] <- ci.upp$y
        }
        else {
          ci.upp <- approx(spline.upp$y, spline.upp$x, 
                           xout = loglik.max - crit)
          conf.int[par.nr, 2] <- ci.upp$y
        }
        if (small.step.upp) 
          conf.int[par.nr, 2] <- log10(CED.upp[1])
        if (loglik.upp[length(loglik.upp)] > loglik.max - 
            crit) 
          conf.int[par.nr, 2] <- Inf
        if (CED.upp.inf) 
          conf.int[par.nr, 2] <- Inf
      }
      else {
        spline.upp <- list()
        spline.upp$x <- CED.upp
        spline.upp$y <- rep(loglik.upp[1], length(CED.upp))
        ci.upp <- list()
        ci.upp$y <- max(CED.upp)
      }
      CED.vec <- c(CED.low, CED.upp)
      loglik.vec <- c(loglik.low, loglik.upp)
      if (!CI.NA) 
        if (trace.plt) {
          y.low <- min(spline.low$y, spline.upp$y)
          y.upp <- max(spline.low$y, spline.upp$y)
          y.lim <- c(min(y.low, loglik.max - crit), y.upp)
          if (nolog) 
            x.lab = text.par[jj]
          else x.lab = paste("log10(", text.par[jj], 
                             ")")
          plot(CED.vec, loglik.vec, xlab = x.lab, ylim = y.lim, 
               ylab = "log(likelihood)", col = 1, pch = 16)
          title("calculation of confidence interval")
          if (nolog) 
            MLE.tmp <- MLE[jj]
          else MLE.tmp <- log10(MLE[jj])
          points(MLE.tmp, loglik.max, col = 8, pch = 16, 
                 cex = 2.5)
          lines(spline.low$x, spline.low$y, col = 6)
          lines(spline.upp$x, spline.upp$y, col = 6)
          segments(min(CED.vec), loglik.max - crit, max(CED.vec), 
                   loglik.max - crit, col = 1, lty = 2)
          lines(rep(ci.low$y, 2), c(min(loglik.vec), 
                                    loglik.max), lty = 2)
          lines(rep(ci.upp$y, 2), c(min(loglik.vec), 
                                    loglik.max), lty = 2)
          y.low.old <- y.low
          y.upp.old <- y.upp
          if (y.low < loglik.max/1000) 
            y.low <- max(loglik.max - 5, min(spline.low$y, 
                                             spline.upp$y))
          if (y.upp > loglik.max * 1000) 
            y.upp <- min(loglik.max + 5, max(spline.low$y, 
                                             spline.upp$y))
          if (y.low != y.low.old | y.upp != y.upp.old) {
            cat("\n\nplot with outliers omitted\n\n")
            f.graph.window(1, nr.gr = nr.gr)
            plot(CED.vec, loglik.vec, xlab = x.lab, ylim = c(y.low, 
                                                             y.upp), ylab = "log(likelihood)", col = 1, 
                 pch = 16)
            points(MLE.tmp, loglik.max, col = 8, pch = 16, 
                   cex = 2.5)
            lines(spline.low$x, spline.low$y, col = 6)
            lines(spline.upp$x, spline.upp$y, col = 6)
            segments(min(CED.vec), loglik.max - crit, 
                     max(CED.vec), loglik.max - crit, col = 1)
            lines(rep(ci.low$y, 2), c(min(loglik.vec), 
                                      loglik.max), lty = 2)
            lines(rep(ci.upp$y, 2), c(min(loglik.vec), 
                                      loglik.max), lty = 2)
          }
        }
      if (jump.low) 
        conf.int[par.nr, 1] <- CED.vec[1]
      if (jump.upp) 
        conf.int[par.nr, 2] <- CED.vec[length(CED.vec)]
    }
    {
      if (max(abs(diff(loglik.low))) < 0.01) 
        conf.int[par.nr, 1] <- -Inf
      if (max(abs(diff(loglik.upp))) < 0.01) 
        conf.int[par.nr, 2] <- Inf
      if (CI.NA) {
        cat("\nCI could not be established, log-likelihood profile did not change\n")
        conf.int[par.nr, 1] <- NA
        conf.int[par.nr, 2] <- NA
      }
      if (length(conf.int[, 1]) > 1) 
        dimnames(conf.int)[[1]] <- text.par[group]
    }
    profile.out <- list(conf.int = conf.int, MLE.new = MLE.new, 
                        loglik = nlminb.out$loglik, profile = cbind(CED.vec, 
                                                                    loglik.vec))
    if (quick.ans == 1) {
      cat("\nstarted at:", date.start)
      cat("\n")
      cat(date())
    }
    if (exists("track")) 
      print("end of profile.all")
    return(profile.out)
  })
}

################
## f.plot.gui ##
################
f.plot.gui <- function (ans.all, HTML = FALSE, model.summ = TRUE) 
{
  if (exists("track")) 
    print("f.plot.gui")
  WAPP <- ans.all$WAPP
  if (length(ans.all$xans) == 1) 
    ans.all$displ.ans <- ans.all$covar.no
  if (ans.all$cont || ans.all$dtype == 3) {
    if (ans.all$plot.type %in% c(9, 10)) {
      if (ans.all$dtype %in% c(10, 15, 250, 260)) {
        cat("\n residual plots not possible for summary data\n")
        return(invisible())
      }
      if (ans.all$model.fam == 0) {
        cat("\n\nATTENTION: first select a single model when plotting residuals\n\n")
        return(invisible())
      }
    }
    if (length(ans.all$xans) > 1) {
      ans.all$plt.mns <- 3
      ans.all$CI.plt <- TRUE
    }
    if (ans.all$model.fam == 0) {
      if (ans.all$dtype != 3) 
        if (!HTML) {
          if (ans.all$no.hill) 
            f.graph.window(1)
          else {
            if (WAPP & !ans.all$svg.plots) 
              .gw.size <- c(800, 400)
            else if (WAPP & ans.all$svg.plots) 
              .gw.size <- c(8, 4)
            else .gw.size <- c(7, 3)
            f.create.graphwin(.gw.size[1], .gw.size[2], 
                              WAPP = WAPP, title = "", name.wapp = "exphill", 
                              plotprefix = ans.all$plotprefix, svg.plots = ans.all$svg.plots)
            par(mfcol = c(1, 2))
            par(mar = c(5.5, 4.8, 4, 12))
            par(cex.main = 1.2)
            par(cex.sub = 1)
            par(cex.lab = 1.5)
            par(cex = 0.6)
          }
          y.pos <- 0
          assign(".ypos", y.pos, immediate = T, pos = 1)
        }
      if (ans.all$dtype == 3) 
        if (ans.all$no.hill) 
          f.graph.window(1, WAPP = WAPP, name.wapp = "exphill", 
                         plotprefix = ans.all$plot.refix, svg.plots = ans.all$svg.plots)
      else f.graph.window(2, WAPP = WAPP, name.wapp = "exphill", 
                          plotprefix = ans.all$plot.refix, svg.plots = ans.all$svg.plots)
      ans.all.plt <- ans.all
      ans.all.plt$heading <- ""
      if (ans.all$dtype == 3) 
        ans.all.plt$show <- f.show.con(ans.all.plt)
      else ans.all.plt$show <- f.show.con(ans.all.plt)
      if (ans.all$EXP$trend == FALSE) 
        ans.all.plt <- f.move.sublist(ans.all.plt, ans.all$NULL.model)
      else ans.all.plt <- f.move.sublist(ans.all.plt, ans.all$EXP)
      if (model.summ) 
        ans.all.plt$show <- f.show.con(ans.all.plt)
      f.plot.all(ans.all.plt, new.window = F)
      title(main = "Exponential model", font.main = 2.5)
      if (!ans.all$no.hill) {
        ans.all.plt <- ans.all
        if (ans.all$HILL$trend == FALSE) 
          ans.all.plt <- f.move.sublist(ans.all.plt, 
                                        ans.all$NULL.model)
        else ans.all.plt <- f.move.sublist(ans.all.plt, 
                                           ans.all$HILL)
        y.pos <- 95
        assign(".ypos", y.pos, immediate = T, pos = 1)
        ans.all.plt$yleg <- ""
        if (model.summ) 
          ans.all.plt$show <- f.show.con(ans.all.plt)
        f.plot.all(ans.all.plt, new.window = F)
        title(main = "Hill model", font.main = 2)
      }
    }
    if (ans.all$model.fam > 0) {
      if (HTML || WAPP) {
        new.window <- F
        f.graph.window(1, WAPP = WAPP, name.wapp = "exphill", 
                       plotprefix = ans.all$plotprefix, svg.plots = ans.all$svg.plots)
      }
      else {
        new.window <- T
        y.pos <- 0
        assign(".ypos", y.pos, immediate = T, pos = 1)
      }
      ans.all.plt <- ans.all
      if (ans.all$model.fam == 1) {
        ans.all.plt <- f.move.sublist(ans.all.plt, ans.all$EXP)
        ans.all.plt$heading <- ""
        if (model.summ) 
          ans.all.plt$show <- f.show.con(ans.all.plt)
        f.plot.all(ans.all.plt, new.window = new.window)
        title(main = "Exponential model", font.main = 2)
      }
      if (ans.all$model.fam == 2) {
        ans.all.plt <- f.move.sublist(ans.all.plt, ans.all$HILL)
        ans.all.plt$heading <- ""
        if (model.summ) 
          ans.all.plt$show <- f.show.con(ans.all.plt)
        f.plot.all(ans.all.plt, new.window = new.window)
        title(main = "Hill model", font.main = 2)
      }
    }
  }
  if (!ans.all$cont && ans.all$dtype != 3) {
    if (ans.all$gui && !ans.all$WAPP) {
      if (ans.all$model.fam > 0) {
        ans.all$model.ans <- ans.all$model.fam
        if (0) 
          if (ans.all$covar.no == 0) {
            if (ans.all$model.fam == 1) 
              ans.all$model.ans <- 8
            if (ans.all$model.fam == 2) 
              ans.all$model.ans <- 9
          }
        if (ans.all$covar.no > 0) {
          if (ans.all$model.fam == 1) 
            ans.all$model.ans <- 7
          if (ans.all$model.fam == 2) 
            ans.all$model.ans <- 8
        }
      }
    }
    if (ans.all$model.fam == 0) {
      if (!HTML) {
        f.graph.window(9, WAPP = WAPP, title = "quantal models", 
                       name.wapp = "setofmod", plotprefix = ans.all$plotprefix, 
                       svg.plots = ans.all$svg.plots)
      }
      count <- 0
      nr.models <- length(ans.all$model.list) - 2
      for (ii in 1:nr.models) {
        ans.all.plt <- ans.all
        name <- ans.all$model.list[ii]
        if (ans.all[[name]]$trend) 
          ans.all.plt <- f.move.sublist(ans.all.plt, 
                                        ans.all[[name]])
        else {
          ans.all.plt <- f.move.sublist(ans.all.plt, 
                                        ans.all$NULL.model)
          ans.all.plt$modelname <- ans.all[[name]]$modelname
        }
        ans.all.plt$heading <- ""
        if (ans.all$dtype == 2) 
          ans.all.plt$CI.plt <- FALSE
        count <- count + 1
        xleg <- ans.all.plt$xleg
        yleg <- ans.all.plt$yleg
        f.plot.all(ans.all.plt, new.window = FALSE)
        title(main = paste(ans.all.plt$modelname, "--"), 
              cex = 0.5)
      }
      if (ans.all$LVM) {
        ans.all.plt <- ans.all
        ans.all.plt$model.type <- 2
        if (1) 
          ans.all.plt <- f.move.sublist(ans.all.plt, 
                                        ans.all$EXP)
        else {
          ans.all.plt <- f.move.sublist(ans.all.plt, 
                                        ans.all$NULL.model)
          ans.all.plt$modelname <- ans.all[[name]]$modelname
        }
        ans.all.plt$heading <- ""
        if (ans.all$dtype == 2) 
          ans.all.plt$CI.plt <- FALSE
        count <- count + 1
        f.plot.all(ans.all.plt, new.window = FALSE)
        title(main = ans.all$EXP$modelname, cex = 0.5)
        ans.all.plt <- ans.all
        ans.all.plt$model.type <- 2
        if (1) 
          ans.all.plt <- f.move.sublist(ans.all.plt, 
                                        ans.all$HILL)
        else ans.all.plt <- f.move.sublist(ans.all.plt, 
                                           ans.all$NULL.model)
        ans.all.plt$heading <- ""
        if (ans.all$dtype == 2) 
          ans.all.plt$CI.plt <- FALSE
        count <- count + 1
        f.plot.all(ans.all.plt, new.window = FALSE)
        title(main = ans.all$HILL$modelname, cex = 0.5)
      }
      if (count < 9) {
        plot(1:3, 1:3, col = 0, yaxt = "n", xlab = "", 
             ylab = "", type = "n")
        text(1.2, 2.9, ans.all$odt.name, adj = 0)
        text(1.2, 2.5, "x-axis:", adj = 0)
        text(1.2, 2.3, xleg, adj = 0)
        text(1.2, 2.1, "y-axis:", adj = 0)
        text(1.2, 1.9, yleg, adj = 0)
        text(1.2, 1.2, paste("PROAST version", ans.all$PRversion), 
             adj = 0)
      }
      else mtext(ans.all$xleg, side = 1)
    }
    if (ans.all$model.fam > 0) {
      if (!HTML) 
        f.graph.window(1, WAPP = WAPP, name.wapp = "exphill", 
                       plotprefix = ans.all$plotprefix, svg.plots = ans.all$svg.plots)
      if (ans.all$model.ans %in% 1:9) {
        ii <- ans.all$model.ans
        ans.all.plt <- ans.all
        name <- ans.all$model.list[ii]
        ans.all.plt <- f.move.sublist(ans.all.plt, ans.all[[name]])
        if (ans.all$dtype == 2) 
          ans.all.plt$CI.plt <- FALSE
        if (model.summ) 
          ans.all.plt$show <- f.show.cat(ans.all.plt)
        f.plot.all(ans.all.plt, new.window = FALSE)
      }
    }
  }
  if (WAPP) 
    dev.off()
  if (exists("track")) 
    print("f.plot.gui:  END")
  return(invisible())
}

################
## f.show.con ##
################
f.show.con <- function (ans.all) 
{
  if (exists("track")) 
    print("f.show.con")
  with(ans.all, {
    nl <- "\n"
    tb <- "  "
    sp <- ": "
    if (length(ans.all$fct4) > 0) 
      nr.cc <- max(fct4)
    if (length(ans.all$fct5) > 0) 
      nr.dd <- max(fct5)
    if (fit.ans == 1) {
      loglik.txt <- "loglik"
    }
    if (fit.ans == 2) {
      loglik.txt <- "SP"
      loglik <- SP
      nr.var <- 0
    }
    if (length(MLE) < 22) {
      results <- c(round(loglik, 2), round(aic, 2), signif(MLE, 
                                                           4))
      txt <- c(loglik.txt, "AIC", text.par)
    }
    else {
      cat("\nATTENTION: Number of parameters too large for legend in plot\n")
      results <- c(round(loglik, 2), nr.var, nr.aa, nr.bb, 
                   nr.cc, nr.dd, signif(MLE[c(1, nr.var + 1, nr.var + 
                                                nr.aa + 1, nr.var + nr.aa + nr.bb + 1, nr.var + 
                                                nr.aa + nr.bb + nr.cc + 1)], 3))
      txt <- c(loglik.txt, "nr.var", "nr.aa", "nr.bb", 
               "nr.cc", "nr.dd", text.par[c(1, nr.var + 1, nr.var + 
                                              nr.aa + 1, nr.var + nr.aa + nr.bb + 1, nr.var + 
                                              nr.aa + nr.bb + nr.cc + 1)])
      results <- results[!is.na(results)]
      txt <- txt[!is.na(txt)]
    }
    show <- paste(" version:", PRversion[1])
    for (i in 1:length(txt)) show <- paste(show, nl, txt[i], 
                                           tb, results[i])
    if (!is.na(CED[1])) {
      if (nr.bb == 1) 
        fct2.txt <- NULL
      CES.txt <- "CES"
      if (NES.ans == 2 || model.ans == 47) 
        CES.txt <- "adjusted CES"
      CED.txt <- "CED"
      if (length(ans.all$ED50) != 0) 
        if (ED50) 
          CED.txt <- "ED50"
      else ED50 <- FALSE
      if ((CES != 0 || ED50) && !is.na(CED[1])) 
        if (!model.ans %in% c(12:16, 22:25, 46, 47)) {
          show <- paste(show, nl, CES.txt, tb, round(CES, 
                                                     3))
          CED.uniq <- unique(CED)
          for (ii in (1:nr.bb)) show <- paste(show, nl, 
                                              paste(CED.txt, fct2.txt[ii], sep = "-"), 
                                              signif(CED.uniq[ii], 4))
        }
      else {
        if (model.ans == 47) 
          show <- paste(show, nl, "adjusted CES", tb, 
                        round(CES, 4))
        else show <- paste(show, nl, CES.txt, tb, round(CES, 
                                                        3))
      }
      if (length(conf.int) > 0) 
        if (any(!is.na(conf.int))) {
          if (boot) 
            bt <- "bt"
          else bt <- ""
          BMDL.txt <- paste("CEDL", bt, sep = "")
          BMDU.txt <- paste("CEDU", bt, sep = "")
          if ((nr.aa != nr.bb) && (nr.bb > 1) && (nr.aa > 
                                                  1)) 
            for (ii in (1:nr.bb)) {
              show <- paste(show, nl, paste(BMDL.txt, 
                                            fct2.txt[ii], sep = "-"), tb, signif(conf.int[ii, 
                                                                                          1], 4))
              show <- paste(show, nl, paste(BMDU.txt, 
                                            fct2.txt[ii], sep = "-"), tb, signif(conf.int[ii, 
                                                                                          2], 4))
            }
          else if ((model.ans == 16) | (nr.bb == 1)) {
            show <- paste(show, nl, BMDL.txt, tb, signif(conf.int[1, 
                                                                  1], 4))
            show <- paste(show, nl, BMDU.txt, tb, signif(conf.int[1, 
                                                                  2], 4))
          }
          else if ((model.ans != 16) & (nr.aa > 1)) {
            for (ii in (1:nr.aa)) {
              show <- paste(show, nl, paste(BMDL.txt, 
                                            fct2.txt[ii], sep = "-"), tb, signif(conf.int[ii, 
                                                                                          1], 4))
              show <- paste(show, nl, paste(BMDU.txt, 
                                            fct2.txt[ii], sep = "-"), tb, signif(conf.int[ii, 
                                                                                          2], 4))
            }
          }
          else for (ii in (1:nr.bb)) {
            show <- paste(show, nl, paste(BMDL.txt, fct2.txt[ii], 
                                          sep = "-"), tb, signif(conf.int[ii, 1], 
                                                                 4))
            show <- paste(show, nl, paste(BMDU.txt, fct2.txt[ii], 
                                          sep = "-"), tb, signif(conf.int[ii, 2], 
                                                                 4))
          }
          if (nruns > 0) 
            show <- paste(show, nl, paste("nruns:", nruns))
        }
    }
    if (model.ans %in% c(12:15, 22:25, 46)) {
      CED.tmp <- MLE[(nr.var + nr.aa + 1):(nr.var + nr.aa + 
                                             nr.bb)]
      cc0 <- NA
      dd0 <- NA
      if (model.ans %in% c(13, 23)) 
        dd0 <- MLE[(nr.var + nr.aa + nr.bb + 1):(nr.var + 
                                                   nr.aa + nr.bb + nr.dd)]
      if (model.ans %in% c(14:15, 24:25, 46)) 
        cc0 <- MLE[(nr.var + nr.aa + nr.bb + 1):(nr.var + 
                                                   nr.aa + nr.bb + nr.cc)]
      if (model.ans %in% c(15, 25, 46)) 
        dd0 <- MLE[(nr.var + nr.aa + nr.bb + 2):(nr.var + 
                                                   nr.aa + nr.bb + 1 + nr.dd)]
      bb.tmp <- f.bb.con(model.ans, cc0, dd0, CED.tmp, 
                         CES, CED.ref = CED.tmp[ref.lev], ref.lev = ref.lev)
      if (length(MLE) < 10) 
        for (i in 1:length(bb.tmp)) show <- paste(show, 
                                                  "\n b: ", signif(bb.tmp[i], 4))
    }
    if (dtype %in% c(5, 15)) {
      show <- paste(show, nl, "inter.var:", signif(inter.var, 
                                                   2))
      show <- paste(show, nl, "intra.var:", signif(intra.var, 
                                                   2))
    }
    show <- paste(show, nl, "conv", sp, converged)
    show <- paste(show, nl, "scaling factor on x", sp, sf.x)
    if (any(Vdetlim > 0)) {
      if (length(unique(Vdetlim)) == 1) 
        show <- paste(show, nl, "detlim", sp, Vdetlim[1])
      else show <- paste(show, nl, "detlim", sp, min(Vdetlim, 
                                                     na.rm = T), "-", max(Vdetlim, na.rm = T))
    }
    if (!is.na(cens.up)) 
      show <- paste(show, nl, "cens.up", sp, cens.up)
    show <- paste(show, nl, "dtype", sp, dtype)
    if (!gui) {
      sel.tmp <- character(0)
      for (ii in 1:length(select.name)) sel.tmp <- paste(sel.tmp, 
                                                         "selected", sp, select.name[ii])
      if (length(sel.tmp) > 0) 
        show <- paste(show, nl, sel.tmp)
      else show <- paste(show, nl, "selected: all")
      rem.tmp <- character()
      if (is.na(remove.name)) 
        show <- paste(show, nl, "removed: NA")
      else if (remove.name == "none") 
        show <- paste(show, nl, "removed: none")
      else {
        for (ii in 1:length(remove.name)) rem.tmp <- paste(rem.tmp, 
                                                           " ", remove.name[ii])
        if (length(rem.tmp) >= 1) {
          show <- paste(show, nl, "removed:")
          show <- paste(show, rem.tmp)
        }
      }
    }
    if (quick.ans == 1) 
      show <- paste(show, nl, factor.name)
    if (length(covar.txt) > 1) 
      show <- paste(show, nl, paste("covariate:", covar.name, 
                                    "\n"))
    if (cond.ans != 3) 
      show <- paste(show, "\n fit cond:", cond.ans)
    if (exists("track")) 
      print("f.show.con:  END")
    return(show)
  })
}

##############
## f.bb.con ##
##############
f.bb.con <- function (model.ans, cc, dd, CED, CES, ref.lev = NA, CED.ref = NA) 
{
  if (exists("track2")) 
    print("f.bb.con")
  if (0) {
    print("f.bb.con")
    print(model.ans)
    print(CES)
  }
  CES.tmp <- CES
  if (model.ans == 12) {
    if (CES < 0) 
      CES.tmp <- CES
    bb <- logb(CES.tmp + 1)/CED
  }
  if (model.ans == 13) {
    if (CES < 0) 
      CES.tmp <- CES
    bb <- logb(CES.tmp + 1)/CED^dd
  }
  if (model.ans == 14) {
    if ((cc < 1)) 
      CES.tmp <- CES
    bb <- -logb((CES.tmp + 1 - cc)/(1 - cc))/CED
  }
  if (model.ans == 15) {
    if ((cc < 1)) 
      CES.tmp <- -abs(CES)
    dum <- (CES.tmp + 1 - cc)/(1 - cc)
    bb <- -log(dum)/(CED^dd)
  }
  if (model.ans == 22) {
    bb <- -CED * (1 + CES)/CES
  }
  if (model.ans == 23) {
    if (CES >= 0) 
      bb <- -CED/(CES/(1 + CES))^(1/dd)
    if (CES < 0) 
      bb <- CED/(-(CES/(1 + CES)))^(1/dd)
  }
  if (model.ans == 24) {
    bb <- CED/(CES/(cc - 1 - CES))
  }
  if (model.ans == 25) {
    if (cc == 0) 
      bb <- -CED * ((1 + CES)/CES)^(1/dd)
    else {
      CES <- abs(CES)
      if (cc >= 1) 
        bb <- CED/((CES/(cc - 1 - CES))^(1/dd))
      if (cc < 1) 
        bb <- CED/((-CES/(cc - 1 + CES))^(1/dd))
    }
  }
  if (model.ans == 46) {
    if (0) {
      print("f.bb.con")
      print(CED)
      print(CED.ref)
      print(ref.lev)
    }
    CED <- CED.ref/CED
    CED[ref.lev] <- CED.ref
    if ((cc < 1)) 
      CES.tmp <- -abs(CES)
    dum <- (CES.tmp + 1 - cc)/(1 - cc)
    bb <- -log(dum)/(CED^dd)
  }
  if (exists("track2")) 
    print("f.bb.con:END")
  return(bb)
}

################
## f.plot.all ##
################
f.plot.all <- function (ans.all, sep = F, bootstrap = F, new.window = TRUE) 
{
  if (exists("track")) 
    print("f.plot.all")
  if (ans.all$plot.type == 0) 
    return(invisible())
  if (ans.all$tans != 0) {
    ans.all$y.lim <- f.plot.cxt(ans.all, sep = sep)
    return(ans.all)
  }
  if (is.na(ans.all$xy.lim[1]) || ans.all$xy.lim[1] == 0) {
    x <- ans.all$x
    dum.contr <- min(x[x > 0])/5
    ans.all$xy.lim[1] <- dum.contr
  }
  if (0) 
    if (ans.all$model.ans == 46) {
      ans.all.saved <- ans.all
      ans.all$model.ans <- 15
      nr.aa <- max(ans.all$fct1)
      nr.bb <- max(ans.all$fct2)
      RPF.vec <- ans.all$regr.par[(nr.aa + 1):(nr.aa + 
                                                 nr.bb)]
      CED.ref <- RPF.vec[ans.all$ref.lev]
      CED.vec <- CED.ref/RPF.vec
      CED.vec[ans.all$ref.lev] <- CED.ref
      ans.all$regr.par[(nr.aa + 1):(nr.aa + nr.bb)] <- CED.vec
      ans.all <- f.pars(ans.all)
      if (ans.all$cont) 
        ans.all$CED <- CED.vec
      if (!ans.all$cont) 
        ans.all$CED.matr <- matrix(CED.vec, nrow = 1)
    }
  if (ans.all$dtype == 3 && !ans.all$combi.ans && (ans.all$plot.type %in% 
                                                   5:6)) {
    if (ans.all$quick.ans == 1) 
      ans.all$categ.ans <- menu(c("all categories in one plot", 
                                  "separate plots for each category"), title = paste("\n\nHow do you want to plot the categories?\n"))
    else ans.all$categ <- 1
    if (ans.all$nr.gr > 1) 
      sep <- TRUE
    if (ans.all$displ.ans > 0) {
      sep <- TRUE
      ans.all$fct2 <- ans.all$displ.fact
    }
    if (ans.all$categ.ans == 2) 
      ans.all$cex.1 <- 2
    ans.all$sep <- sep
  }
  if (sep) {
    f.plot.sep(ans.all)
    return()
  }
  with(ans.all, {
    if (cont) {
      if (new.window) 
        f.graph.window(1)
      if (plot.type %in% 10:11) {
        if (is.na(regr.resid[1])) 
          cat("\nno residuals available \n")
        else ans.all$yy <- regr.resid
        ans.all$CI.plt <- FALSE
      }
      ans.all <- f.plot.con(ans.all)
      if (model.ans != 11 && plot.type != 9) 
        f.lines.con(ans.all)
      if (plot.type < 9 && !is.na(CED[1]) && model.ans != 
          16 && fitted) {
        f.cedlines.con(ans.all)
      }
      f.mtext(ans.all)
    }
    if (!cont) {
      if (!combi.ans) {
        if (!bootstrap && new.window) 
          f.graph.window(1)
        if (dtype %in% c(4, 6, 84) && model.type == 2 && 
            plot.type %in% c(7:10)) {
          ans.all.plt <- f.model.bb(ans.all)
          out.list <- f.lines.cat(ans.all.plt)
          ans.all.plt$y.lim <- out.list$y.lim
          return(invisible())
        }
        else if (dtype %in% c(4, 6, 84) || (dtype == 
                                            2 && plot.type < 3)) {
          ans.all <- f.plot.frq(ans.all)
          if (!(model.ans == 14 && (model.type == 1))) {
            f.lines.frq(ans.all)
          }
          if (!is.na(ans.all$CED.matr[1])) 
            f.cedlines.bin(ans.all)
          if (dtype == 6 & plot.type == 8 & model.type == 
              1) {
            pi.expect <- f.expect.bin(model.ans, x, MLE[-1], 
                                      CES = CES, CES1 = CES1, CES2 = CES2, output = output)
            beta <- (alfa * (1 - pi.expect))/pi.expect
            fract.random <- rbeta(length(x), alfa, beta)
            plot(sort(fract.random), sort(kk/nn), cex = 0.65)
            abline(0, 1)
          }
        }
        else {
          if (model.type == 2) 
            ans.all.plt <- f.model.bb(ans.all)
          else ans.all.plt <- ans.all
          out.list <- f.lines.cat(ans.all.plt)
          ans.all$y.lim.plt <- out.list$y.lim
        }
        if (ans.all$categ.ans < 2) 
          f.mtext(ans.all)
      }
      if (combi.ans) {
        f.graph.window(4)
        ans.all.plt <- ans.all
        ans.all.plt$modelname <- ""
        ans.all.plt$plot.type <- 3 + ans.scale
        if (dtype %in% c(2, 3, 4, 6)) 
          ans.all.plt$ced.lines <- TRUE
        ans.all.plt$main <- ""
        out.list <- f.lines.cat(ans.all.plt)
        ans.all.plt$y.lim <- out.list$y.lim
        ans.all.plt$show <- f.show.cat(ans.all.plt)
        title(modelname, col = color[1], cex.main = 0.9, 
              font.main = 1)
        f.mtext(ans.all.plt)
        ans.all.plt$plot.type <- 9 + ans.scale
        ans.all.plt$modelname <- ""
        ans.all.plt$ced.lines <- F
        out.list <- f.lines.cat(ans.all.plt)
        ans.all.plt$y.lim <- out.list$y.lim
        if (!is.na(CED.matr[1])) {
          th.par.cum <- cumsum(th.par)
          CED.vect <- CED.matr[1, ]
          if (ans.scale == 1) {
            CED.vect <- log10(CED.vect)
            x <- log10(x[x > 0]/sf.x)
          }
          for (ii in 1:length(CED.vect)) lines(c(min(x), 
                                                 CED.vect[ii]), rep(th.par.cum[ii], 2), lty = 2)
          y.leg.2 <- paste(y.leg, " (latent var.uu)")
          if (dtype == 3 & nr.gr == 1) {
            CED.tmp <- round(CED.vect, rep(4, length(CED.vect)))
            CED.tmp <- CED.tmp[CED.tmp < max(x)]
            CED.left <- c(CED.tmp[1]/100, CED.tmp)
            CED.right <- c(CED.tmp, max(x))
            CED.mid <- (CED.left + CED.right)/2
            if (ans.scale == 1) 
              CED.mid <- 10^(CED.mid)
            uu <- f.expect.con(model.ans, CED.mid, regr.par, 
                               fct1 = 1, fct2 = 1, CES = CES, increase = increase)
            zz <- logb(uu)
            if (ans.scale == 1) 
              CED.mid <- log10(CED.mid)
            for (ii in 1:length(scores.orig)) text(CED.mid[ii], 
                                                   zz[ii], scores.orig[ii])
          }
        }
        th.par.cum <- cumsum(th.par)
        if (mode(CED.matr) == "NULL") 
          CED.vect <- 0
        else if (is.na(CED.matr[1])) 
          CED.vect <- 0
        else CED.vect <- CED.matr[1, ]
        if (dtype %in% 2:3) {
          ans.all.plt$plot.type <- 5 + ans.scale
        }
      }
    }
    if (model.ans == 46 && dtype == 3) 
      ans.all <- ans.all.saved
    if (exists("track")) 
      print("f.plot.all : END")
    return(ans.all)
  })
}

####################
## f.cedlines.con ##
####################
f.cedlines.con <- function (ans.all) 
{
  if (exists("track")) 
    print("f.cedlines.con")
  ans.all$cedes <- f.cedes.plt.con(ans.all)
  with(ans.all, {
    ES.x <- cedes$ES.x
    ES.y <- cedes$ES.y
    CED.x <- cedes$CED.x
    CED.y <- cedes$CED.y
    if (length(xans) > 1) 
      nr.gr <- 1
    for (ii in (1:nr.gr)) {
      lines(ES.x[ii, ], ES.y[ii, ], lty = 2)
      lines(CED.x[ii, ], CED.y[ii, ], lty = 2)
    }
    if (exists("track")) 
      print("f.cedlines.con:  END")
    return(invisible())
  })
}

#####################
## f.cedes.plt.con ##
#####################
f.cedes.plt.con <- function (ans.all) 
{
  if (exists("track")) 
    print("f.cedes.plt.con")
  if (is.null(ans.all$regr.par.matr)) 
    ans.all$regr.par.matr <- f.pars(ans.all)$regr.par.matr
  with(ans.all, {
    if (model.ans == 46) {
      CED.ref <- CED[ref.lev]
      CED <- CED.ref/CED
      CED[ref.lev] <- CED.ref
      regr.par[(nr.aa + 1):(nr.aa + nr.bb)] <- CED
      ans.all$model.ans <- 15
      model.ans <- 15
      ans.all$regr.par <- regr.par
      regr.par.matr <- f.pars(ans.all)$regr.par.matr
    }
    if (!is.matrix(regr.par.matr)) 
      regr.par.matr <- t(as.matrix(regr.par.matr))
    nr.subgr <- length(regr.par.matr[, 1])
    if (ans.m6.sd == 2) {
      regr.par <- regr.par.matr[1, ]
      GCES <- CES
      qq <- regr.par[4]
      if (model.ans == 6) {
        if (is.na(CED)) 
          cat("\n\nATTENTION: CED not yet calculated ......  \n\n")
        CED.0 <- CED
      }
      if (model.ans == 47) 
        CED.0 <- regr.par[3]
    }
    if (nr.bb == 1) 
      CED <- rep(CED, nr.aa)
    dum.contr <- xy.lim[1]
    low.y <- xy.lim[4]
    if (is.na(low.y)) 
      low.y <- min(y)
    if (low.y == 9999) 
      low.y <- min(y)
    ES.x <- numeric()
    ES.y <- numeric()
    CED.x <- numeric()
    CED.y <- numeric()
    for (jj in 1:nr.subgr) {
      if (ans.m6.sd == 2) {
        regr.par <- regr.par.matr[jj, ]
        sign.q.0 <- sign.q[jj]
        CES.0 <- exp(qq * sqrt(regr.par[1]))^GCES - 1
        ES.0 <- exp(log(regr.par[2]) + sign.q.0 * log(CES.0 + 
                                                        1))
      }
      else {
        CED.0 <- CED[jj]
        CES.0 <- CES
        ES.0 <- f.expect.con(model.ans, CED.0, regr.par.matr[jj, 
                                                             ], fct1 = 1, fct2 = 1, fct5 = 1, CES = CES.0, 
                             increase = increase, ref.lev = ref.lev)
      }
      if (plot.type %in% c(1, 2, 5)) 
        ES.y <- rbind(ES.y, rep(ES.0, 2))
      if (plot.type %in% c(3, 4, 6)) 
        ES.y <- rbind(ES.y, rep(log10(ES.0), 2))
      if (plot.type %in% 7:8) 
        ES.y <- rbind(ES.y, rep(sqrt(ES.0), 2))
      if (plot.type %in% c(1, 3, 7)) 
        ES.x <- rbind(ES.x, c(min(x), CED.0))
      if (plot.type %in% c(2, 4, 8)) 
        ES.x <- rbind(ES.x, c(log10(dum.contr), log10(CED.0)))
      if (plot.type %in% c(5, 6)) 
        ES.x <- rbind(ES.x, c(sqrt(dum.contr), sqrt(CED.0)))
      if (plot.type %in% c(1, 2, 5)) 
        CED.y <- rbind(CED.y, c(low.y, ES.0))
      if (plot.type %in% c(3, 4, 6)) 
        CED.y <- rbind(CED.y, c(log10(low.y), log10(ES.0)))
      if (plot.type %in% c(7, 8)) 
        CED.y <- rbind(CED.y, c(log10(low.y), sqrt(ES.0)))
      if (plot.type %in% c(1, 3, 7)) 
        CED.x <- rbind(CED.x, rep(CED.0, 2))
      if (plot.type %in% c(2, 4, 8)) 
        CED.x <- rbind(CED.x, rep(log10(CED.0), 2))
      if (plot.type %in% c(5, 6)) 
        CED.x <- rbind(CED.x, rep(sqrt(CED.0), 2))
    }
    ES.x <- matrix(ES.x, ncol = 2)
    ES.y <- matrix(ES.y, ncol = 2)
    CED.x <- matrix(CED.x, ncol = 2)
    CED.y <- matrix(CED.y, ncol = 2)
    if (0) 
      if (plot.type %in% c(3, 4, 6) && !is.na(xy.lim[6])) {
        CED.y <- CED.y + xy.lim[6]
        ES.y <- ES.y + xy.lim[6]
      }
    cedes.lst <- list()
    cedes.lst$CED.x = CED.x
    cedes.lst$CED.y = CED.y
    cedes.lst$ES.x = ES.x
    cedes.lst$ES.y = ES.y
    if (0) {
      print("f.cedes.plt.con")
      print(CED.x)
      print(CED.y)
      print(ES.x)
      print(ES.y)
    }
    if (exists("track")) 
      print("f.cedes.plt.con: END")
    return(cedes.lst)
  })
}

#############
## f.mtext ##
#############
f.mtext <- function (ans.all) 
{
  if (exists("track")) 
    print("f.mtext")
  modelname <- ans.all$modelname
  show <- ans.all$show
  if (show != "") {
    if (length(ans.all$y.lim.plt) > 0) 
      y.lim <- ans.all$y.lim.plt
    else y.lim <- ans.all$xy.lim[4:5]
    loc.txt <- y.lim[2] + 0.05 * (y.lim[2] - y.lim[1])
    mtext(show, 4, 1, adj = 0, padj = 1, las = 1, at = loc.txt, 
          cex = 0.7)
    if (ans.all$output) 
      title(main = paste("\n\n", modelname), cex.main = 1.2, 
            font.main = 1)
  }
  if (exists("track")) 
    print("f.mtext:  END")
  return(invisible())
}

#####################
## f.store.results ##
#####################
f.store.results <- function (ans.all, store.name = 0, add.ext = 0) 
{
  if (exists("track")) 
    print("f.store.results")
  if (ans.all$quick.ans == 1) {
    ans.all$EXP <- NULL
    ans.all$HILL <- NULL
  }
  if (store.name == 0) 
    store.name <- "proastTemp"
    #store.name <- eval(as.character(parse(prompt = "\n\ngive name for file to store results (or type 0 if none)  >   ")))
  name.brief <- store.name
  if (store.name != 0 & add.ext != 0) 
    store.name <- paste(store.name, add.ext, sep = ".")
  ans.all$date <- date()
  if (store.name != 0) {
    f.assign(eval(store.name), ans.all)
  }
  return(name.brief)
}

###############
## f.CED.all ##
###############
f.CED.all <- function (exp.out, hill.out, CED.all, y.leg) 
{
  if (exists("track")) 
    print("f.CED.all")
  CED.matr <- CED.all$CED.matr
  endpoints <- CED.all$endpoints
  Mces <- CED.all$Mces
  Msd <- CED.all$Msd
  covar.txt <- CED.all$covar.txt
  if (length(exp.out) > 0) {
    nr.aa <- max(exp.out$fct1)
    nr.bb <- max(exp.out$fct2)
    CED <- exp.out$CED
    conf.int <- exp.out$conf.int
    if (length(covar.txt) > 2) {
      cat("\nATTENTION: plot of all CED CIs not implemented for more than two subgroups\n")
      return(CED.all)
    }
    else if (length(conf.int[, 1]) == 1 && length(covar.txt) > 
             1) {
      confint.tmp <- conf.int[1, ]
      conf.int <- matrix(ncol = 2, nrow = length(covar.txt))
      for (ii in 1:length(covar.txt)) conf.int[ii, ] <- confint.tmp
    }
    if (!is.na(CED[1])) {
      col <- 1
      CED.tmp <- rep(NA, length(CED.matr[1, ]))
      for (kk in 1:length(covar.txt)) {
        CED.tmp[(col):(col + 1)] <- signif(conf.int[kk, 
                                                    ], 3)
        col <- col + 2
      }
      CED.matr <- rbind(CED.matr, CED.tmp)
      endpoints <- c(endpoints, y.leg)
    }
    ces.tmp <- c(y.leg, round(exp.out$CES, 3))
    sd.tmp <- c(y.leg, round(exp.out$sd, 4))
  }
  if (length(hill.out) > 0) {
    nr.aa <- max(hill.out$fct1)
    nr.bb <- max(hill.out$fct2)
    CED <- hill.out$CED
    conf.int <- hill.out$conf.int
    if (length(conf.int[, 1]) == 1 && length(covar.txt) > 
        1) {
      confint.tmp <- conf.int[1, ]
      conf.int <- matrix(ncol = 2, nrow = length(covar.txt))
      for (ii in 1:length(covar.txt)) conf.int[ii, ] <- confint.tmp
    }
    if (!is.na(CED[1])) {
      col <- 1
      CED.tmp <- rep(NA, length(CED.matr[1, ]))
      for (kk in 1:length(covar.txt)) {
        CED.tmp[(col):(col + 1)] <- signif(conf.int[kk, 
                                                    ], 3)
        col <- col + 2
      }
      CED.matr <- rbind(CED.matr, CED.tmp)
      endpoints <- c(endpoints, y.leg)
    }
    ces.tmp <- c(ces.tmp, round(hill.out$CES, 3))
    sd.tmp <- c(sd.tmp, round(hill.out$sd, 4))
  }
  Mces <- rbind(Mces, ces.tmp)
  Msd <- rbind(Msd, sd.tmp)
  CED.all$endpoints <- endpoints
  CED.all$CED.matr <- CED.matr
  CED.all$CES.0 <- exp.out$CES.0
  CED.all$Mces <- Mces
  CED.all$Msd <- Msd
  if (exists("track")) 
    print("f.CED.all:   END")
  return(CED.all)
}

##################
## f.alert.full ##
##################
f.alert.full <- function () 
{
  #"\nATTENTION: the AIC of the best model (minimum AIC) is more than two units larger than that of the full model.\nCheck the following options:\n1. the data consist of subgroups not taken into account in the model (e.g. various studies, or two sexes)\n2. the data contain litter effects not taken into account\n3. the (mean) responses shows scatter greater than expected based on random error (the fitted curve does not hit all the CIs)\n4. the response in the top dose group deviates substantially (note the CI around the response) from the fitted model \nIn most cases one of these explanations applies. The following actions may be considered:\n1. use the factor defining the subgroups as a covariate and re-analyse the data\n2. re-analyse the data with litter effects taken into account\n3. the greater scatter in (mean) responses will result in a wider BMD CI; no further action needed\n4. consider to leave out the top dose\n\n"
}

