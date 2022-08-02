suppressMessages(require('parallel'))
suppressMessages(require('glmnet'))
cores = detectCores()/8


### directory address for input data ###
dirs = list(
  wbe_dir = './input/WBE_datasets',
  expression_dir = './input/GTEx_v7',
  genotype_dir = './input/Genotype_datasets',
  input_info = './input/input_info.Rdata'
)


genpath = function(d,n) {
  file.path(d,paste0(n,'.Rdata'))
}  


scaling = function(mat, means, stds) {
  for (p in colnames(mat)) {
    mat[,p] = (mat[,p]-means[p])/stds[p]; return(mat)
  }
}


data_normalization = function(x_train, x_test, y_train, y_test) {
    input = vector('list',4)
    names(input) = c('x_train','x_test','y_train','y_test')
    
    y_mean = mean(y_train) 
    y_sd = sqrt(var(y_train))

    x_means = apply(x_train, MARGIN = 2, mean)
    x_stds = apply(x_train, MARGIN = 2, function(x) {sqrt(var(x))})

    input$y_train = (y_train - y_mean)/y_sd
    input$y_test = (y_test - y_mean)/y_sd
    input$x_train = scaling(x_train, x_means, x_stds)
    input$x_test = scaling(x_test, x_means, x_stds)
    
    return(input)
}


# empirically calculate the variance of the prediction outcome
# for the weights in IVW model
variance_estimator = function(x_train, y_train, nfolds, alpha, lambda) {
  
    N = length(y_train)
    foldids = sample(rep(1:nfolds, length.out = N))
    var_vec = array(0L, nfolds)
    
    for (i in seq(nfolds)) {
      which = (foldids == i) 
      fit = glmnet(x_train[!which,], y_train[!which], alpha = alpha, lambda = lambda)
      var_vec[i] = var(y_train[which] - predict(fit, newx = x_train[which,]))
    }
    
    predvar = mean(var_vec)
    
    return (predvar)
}


# apply IVW to combine the outcomes of GEN model(test) and WBE model(test)
estim_ivw = function(x, w, y) {
  ivw_w = 1/w
  ivw_W = ivw_w/sum(ivw_w)
  
  return(as.numeric(cor(y, (x%*%ivw_W))**2))
}



# regularized Linear Regression for gene expression prediction
# we decided to use Ridge regression as our regularized regression method (alpha = 0)
RLR = function(x_train, x_test, y_train, y_test, IVW = F, MERGED = F, alpha = 0, lambda_set = (10^seq(-3,3,length.out=50)) ) {
  
    res = tryCatch({
        input = data_normalization(x_train, x_test, y_train, y_test)
        cvfit = cv.glmnet(input$x_train, input$y_train, alpha = alpha, lambda = lambda_set)
        
        # predict expression level for a given gene
        pred_test = predict(cvfit, s = cvfit$lambda.min, newx = input$x_test)
        if (var(pred_test) == 0) {
            return(list(completed = FALSE, m = 'the prediction has zero variance'))
        }
        
        test_r2 = cor(as.vector(input$y_test), as.vector(pred_test))**2
        
        
        # calculate weights (inverse variance) for IVW model
        if (IVW) {
          pred_var = variance_estimator(rbind(input$x_train, input$x_test), c(input$y_train, input$y_test), 10, alpha, cvfit$lambda.min)
        } else {pred_var = NA}
        
        # examine the presence of genotype predictor with non-zero effect for MERGED model
        if (MERGED) {
          
          # we use LASSO for examining the presence of genotype predictor with non-zero effect
          efit = cv.glmnet(input$x_train, input$y_train, alpha = 1, lambda = lambda_set)
          beta <- as.vector(coef(efit, s = efit$lambda.min))[-1]
          nonzero <- length(which(abs(beta) > 0))
          
          if (nonzero > 0) {
            eqtl_set = rownames(coef(efit, s = efit$lambda.min))[-1]
            eqtl_set = eqtl_set[which(abs(beta) > 0)]
          } else {
            eqtl_set = NULL
          }
        } else {
          nonzero = NA
          eqtl_set = NA
        }
        
        return(list(completed = TRUE, y_test = input$y_test, pred_test = pred_test, test_r2 = test_r2, pred_var = pred_var, nonzero = nonzero, eqtl_set = eqtl_set))
        
    }, error = function(e) {list(completed = FALSE, m = e)}
    )
    
    res
}


run_models = function(y_train, y_val, y_test, wbe_train, wbe_val, wbe_test, gen_train, gen_val, gen_test) {
  
    res = tryCatch(
    {
        fit = test_r2 = weight_var = c()
        test_y = pred_test = list()
        m = c()

        # WBE model (val)
        # test set: validation set
        imputed = RLR(x_train = wbe_train, x_test = wbe_val, y_train = y_train, y_test = y_val, IVW = T, MERGED = F)
        if (imputed$completed) {
            test_r2['WBE_val'] = imputed$test_r2
            weight_var['WBE_val'] = imputed$pred_var
            fit['WBE_val'] = T
        } else {fit['WBE_val'] = F; test_r2['WBE_val'] = NA; m['WBE_val'] = imputed$m}
        remove(imputed)

        # WBE model (test)
        # test set: test set
        imputed = RLR(x_train = rbind(wbe_train, wbe_val), x_test = wbe_test, y_train = c(y_train, y_val), y_test = y_test, IVW = F, MERGED = F)
        if (imputed$completed) {
            pred_test[['WBE_test']] = imputed$pred_test
            test_y[['WBE_test']] = imputed$y_test
            test_r2['WBE_test'] = imputed$test_r2
            fit['WBE_test'] = T
        } else {fit['WBE_test'] = F; test_r2['WBE_test'] = NA; m['WBE_test'] = imputed$m}
        remove(imputed)

        # GEN model (val)
        # test set: validation set 
        imputed = RLR(x_train = gen_train, x_test = gen_val, y_train = y_train, y_test = y_val, IVW = T, MERGED = T)
        if (imputed$completed) {
            test_r2['GEN_val'] = imputed$test_r2
            weight_var['GEN_val'] = imputed$pred_var
            eqtl_N = imputed$nonzero
            eqtl_set = imputed$eqtl_set
            fit['GEN_val'] = T
        } else {fit['GEN_val'] = F; test_r2['GEN_val'] = NA; m['GEN_val'] = imputed$m}
        remove(imputed)

        # GEN model (test)
        # test set: test set 
        imputed = RLR(x_train = rbind(gen_train, gen_val), x_test = gen_test, y_train = c(y_train, y_val), y_test = y_test, IVW = F, MERGED = F)
        if (imputed$completed) {
            pred_test[['GEN_test']] = imputed$pred_test
            test_y[['GEN_test']] = imputed$y_test
            test_r2['GEN_test'] = imputed$test_r2
            fit['GEN_test'] = T
        } else {fit['GEN_test'] = F; test_r2['GEN_test'] = NA; m['GEN_test'] = imputed$m}
        remove(imputed)
            
          
        #####  combined models  #####
        
        # MERGED model 
        # impute gene expression using the merged dataset only when at least 1 genotype predictor with nonzero effect exists
        if (!is.na(eqtl_N) & eqtl_N > 0) {
          
            tag = 'MERGED_test'
            gen_train = as.matrix(gen_train[,eqtl_set])
            gen_val = as.matrix(gen_val[,eqtl_set])
            gen_test = as.matrix(gen_test[,eqtl_set])
          
            if(eqtl_N == 1) {
              colnames(gen_train) = eqtl_set
              colnames(gen_val) = eqtl_set
              colnames(gen_test) = eqtl_set
            }
          
            imputed = RLR(x_train = cbind(rbind(wbe_train, wbe_val), rbind(gen_train, gen_val)), x_test = cbind(wbe_test, gen_test), y_train = c(y_train, y_val), y_test = y_test, IVW = F, MERGED = F)
            if(imputed$completed){
              pred_test[[tag]] = imputed$pred_test
              test_y[[tag]] = imputed$y_test
              test_r2[tag] = imputed$test_r2
              fit[tag] = T
            } else {fit[tag] = F; test_r2[tag] = NA; m[tag] = imputed$m}
            remove(imputed)
          
        } else {
            tag = 'MERGED_test'
            fit[tag] = F; test_r2[tag] = test_r2['WBE_test']; m[tag] = "The gene has no eQTL."
        }
        
        
        # IVW model
        # combine the outcomes of GEN model(test) and WBE model(test) only when GEN model(val) outperformed WBE model(val)
        if (test_r2['GEN_val'] > test_r2['WBE_val']) {
          tag = 'IVW_test'
          ivw = T
          
          y = y_test
          x = cbind(pred_test[['WBE_test']], pred_test[['GEN_test']])
          colnames(x) = c('WBE', 'GEN')
          w = c(weight_var['WBE_val'], weight_var['GEN_val']) 
          test_r2[tag] = estim_ivw(x, w, y)
          
        } else {
          tag = 'IVW_test'
          fit[tag] = F; test_r2[tag] = test_r2['WBE_test']; ivw = F; m[tag] = "WBE(val) > GEN(val)."
        }
        
        ############################# 
        
        
        list(completed = TRUE, m = m, fit = fit, y = test_y, pred = pred_test, r2 = test_r2, ivw = ivw, eqtl_N = eqtl_N)
        
    }, error = function(e) list(completed = FALSE, m = e) 
    )
    res 
}

run_analysis <- function(g) {
  
    res = tryCatch(
        {
          # load genotype data 
          f = genpath(dirs$genotype_dir, g)
          if (!(file.exists(f))) return(list(completed = FALSE, m = 'no genotype file'))
          load(f) ## This file contains genotypes of GTEx XXX samples
          
          # generate model inputs
          GEN = as.matrix(X[,apply(as.matrix(X),2,FUN = function(x) {s = x%%1 == 0; min(sum(x[s]), sum(2-x[s])) > (length(x) * 0.1)})]); remove(X);
          gen_selected = abs(cor(Y[trains,g], GEN[trains,])) > 0.1
          wbe_selected = abs(cor(Y[trains,g], WBE[trains,])) > 0.2
          if(sum(gen_selected) < 10) return(list(completed = FALSE, m = 'sum(# of genotype predictors selected) < 10'))
          if(sum(wbe_selected) < 10) return(list(completed = FALSE, m = 'sum(# of WBE predictors selected) < 10'))

          # run analysis
          res = run_models(y_train = Y[trains,g], y_val = Y[vals,g], y_test = Y[tests,g],
                           wbe_train = WBE[trains, wbe_selected], wbe_val = WBE[vals, wbe_selected], wbe_test = WBE[tests, wbe_selected],
                           gen_train = GEN[trains, gen_selected], gen_val = GEN[vals, gen_selected], gen_test = GEN[tests, gen_selected]
                          )
        }, error = function(e) list(completed = FALSE, m = paste0(g,': ',e))
    )
    
    res
}



# tissue-specific gene expression imputation
GEX_impute = function(t){
  
    cat('Target tissue: ', t, '\n')
  
    trains = sample_ids_info$train[[t]] 
    vals = sample_ids_info$validation[[t]] 
    tests = sample_ids_info$test[[t]]

    load(genpath(dirs$expression_dir, t)) ## Load target expression data (Y)
    Y = tissue_log_tpm_plusone; remove(tissue_log_tpm_plusone)
    cat('Total number of target genes: ', ncol(Y),'\n')
    cat('Total sample counts: ', nrow(Y), '\n')
   
    load(genpath(dirs$wbe_dir, t)) ## Load WBE expression data (WBE)
    WBE = wbe_log_tpm_plusone; remove(wbe_log_tpm_plusone)
    
    genes_to_be_analyzed = colnames(Y)

    cluster <- makeCluster(cores)
    clusterEvalQ(cluster, library(glmnet))
    clusterExport(cl = cluster, varlist = c('Y', 'WBE', 'trains', 'vals', 'tests', 'run_analysis', 'run_models', 'RLR', 'variance_estimator', 'estim_ivw', 'genpath', 'scaling', 'data_normalization', 'dirs'), envir = environment())

    
    results = parLapply(cluster, fun = run_analysis, X = genes_to_be_analyzed)
    names(results) = genes_to_be_analyzed
    
    
    # summarize the imputation results (Mean R2)
    available_genes = c()
    for (g in names(results)){
      if (results[[g]]$completed) available_genes = c(available_genes, g)
    }
    n_genes = length(available_genes)
    
    eqtl_p = ivw_p = 0
    wbe_r2 = gen_r2 = merged_r2 = ivw_r2 = rep(0, n_genes)
    names(wbe_r2) = names(gen_r2) = names(merged_r2) = names(ivw_r2) = available_genes
    
    for (g in available_genes) {
      wbe_r2[g] = results[[g]]$r2['WBE_test'] 
      gen_r2[g] = results[[g]]$r2['GEN_test'] 
      merged_r2[g] = results[[g]]$r2['MERGED_test']
      ivw_r2[g] = results[[g]]$r2['IVW_test']
      
      if (results[[g]]$eqtl_N > 0) eqtl_p = eqtl_p + 1
      if (results[[g]]$ivw) ivw_p = ivw_p + 1
    }
    
    eqtl_p = eqtl_p/n_genes
    ivw_p = ivw_p/n_genes
    
    mean_r2 = c(mean(wbe_r2, na.rm = T), mean(gen_r2, na.rm = T), mean(merged_r2, na.rm = T), mean(ivw_r2, na.rm = T), eqtl_p, ivw_p)
    names(mean_r2) = c('WBE', 'GEN', 'MERGED', 'IVW', 'MERGED_proportion', 'IVW_proportion')
    cat('---------------------------------------------\n')

    remove(Y, WBE, results, genes_to_be_analyzed)
    stopCluster(cluster)
    return(mean_r2)
    
}


load(dirs$input_info)
cat('Tissue-specific gene expression imputation using Genotype and WBE as predictors begins ... ...\n\n\n')
result.summary <- t(sapply(tissues[46:47], FUN = GEX_impute))
cat('Completed ... ...\n\n\n\n')
print(result.summary)

tissue_wise_avg = apply(result.summary, 2, mean)
cat(">> Tissue-wise Average:\n")
print(tissue_wise_avg)

save(result.summary, file = 'output.Rdata')
