sqtl.seeker.modified <- function (tre.df, genotype.f, gene.loc, genic.window = 5000, 
          min.nb.ext.scores = 1000, nb.perm.max = 1e+06, nb.perm.max.svQTL = 10000, 
          svQTL = FALSE, approx = TRUE, verbose = TRUE) 
{
  . = nb.groups = snpId = NULL
  check.genotype <- function(geno.df, tre.df) {
    apply(geno.df, 1, function(geno.snp) {
      if (sum(as.numeric(geno.snp) == -1) > 2) {
        return("Missing genotype")
      }
      geno.snp.t = table(geno.snp[geno.snp > -1])
      if (sum(geno.snp.t >= 2) < 5) {
        return("One group of >5 samples")
      }
      nb.diff.pts = sapply(names(geno.snp.t)[geno.snp.t > 
                                               1], function(geno.i) {
                                                 sQTLseekeR:::nbDiffPt(tre.df[, which(geno.snp == geno.i)])
                                               })
      if (sum(nb.diff.pts >= 2) < 5) {
        return("One group of >5 different splicing")
      }
      return("PASS")
    })
  }
  analyze.gene.f <- function(tre.gene) {
    if (verbose) 
      message(tre.gene$geneId[1])
      gr.gene = with(gene.loc[which(gene.loc$geneId == tre.gene$geneId[1]), 
                            ], GenomicRanges::GRanges(chr, IRanges::IRanges(start, 
                                                                            end)))
    if (genic.window > 0) {
      gr.gene = GenomicRanges::resize(gr.gene, GenomicRanges::width(gr.gene) + 
                                        2 * genic.window, fix = "center")
    }
    if (length(gr.gene) > 0) {
      #print("Check if tre.gene is na:\n")
      tre.gene = tre.gene[, !is.na(tre.gene[1, ])]
      if ( dim(tre.gene)[2] < 3 ) {
        message("Not enough samples")
      } else {
      #print(tre.gene)
      #print("Read table genotype f as genotype headers:\n")
      genotype.headers = as.character(utils::read.table(genotype.f, 
                                                        as.is = TRUE, nrows = 1))
      #print(genotype.headers)
      #print("Colnames tre.gene:\n")
      #print(colnames(tre.gene))
      #print("Checking com samples :\n")
      com.samples = intersect(colnames(tre.gene), genotype.headers)
      #print(com.samples)
      #print("Calculate distance:\n")
      tre.dist = sQTLseekeR:::hellingerDist(tre.gene[, com.samples])
      #print(tre.dist)
      res.df = data.frame()
      #print("Gr.gene.spl is gr.gene is:\n")
      gr.gene.spl = gr.gene
      #print(gr.gene.spl)
      if (any(GenomicRanges::width(gr.gene) > 20000)) {
        #print("Inside if any genomic ranges width gr gene > 20000\n")
        gr.gene.spl = gr.gene[which(GenomicRanges::width(gr.gene) <= 
                                      20000)]
        #print(gr.gene.spl)
        #print("For each unique genomicragne width gr.gene bigger than 20000:\n")
        for (ii in unique(which(GenomicRanges::width(gr.gene) > 
                                20000))) {
          #print("Pos.breaks:\n")
          pos.breaks = unique(round(seq(GenomicRanges::start(gr.gene[ii]), 
                                        GenomicRanges::end(gr.gene[ii]), length.out = floor(GenomicRanges::width(gr.gene[ii])/10000) + 
                                          1)))
          #print(pos.breaks)
          #print("lenght pos.breaks")
          #print(length(pos.breaks))
          
          #print("gr.gene.spl.ii:\n")
          gr.gene.spl.ii = rep(gr.gene[ii], length(pos.breaks) - 
                                 1)
          #print(gr.gene.spl.ii)
          #print("Start:\n")
          GenomicRanges::start(gr.gene.spl.ii) = pos.breaks[-length(pos.breaks)]
          #print(GenomicRanges::start(gr.gene.spl.ii))
          #print("posbreaks length = pos.breaks length + 1\n")
          pos.breaks[length(pos.breaks)] = pos.breaks[length(pos.breaks)] + 
            1
          #print(pos.breaks[length(pos.breaks)])
          #print("End:\n")
          GenomicRanges::end(gr.gene.spl.ii) = pos.breaks[-1] - 
            1
          #print(GenomicRanges::end(gr.gene.spl.ii))
          #print("Gr.gene.spl:\n")
          gr.gene.spl = c(gr.gene.spl, gr.gene.spl.ii)
          #print(gr.gene.spl)
        }
      }
      res.df = lapply(1:length(gr.gene.spl), function(ii) {
        res.range = data.frame()
        if (verbose) {
          message("  Sub-range ", ii)
        }
        #print("Gr.gene.spl[ii]:\n")
        #print(gr.gene.spl[ii])
        
        #print(GenomicRanges::start(gr.gene.spl[ii]))
        if ( GenomicRanges::start(gr.gene.spl[ii]) >= 0 ) {
          genotype.gene = read.bedix(genotype.f, gr.gene.spl[ii])
          #print("Genotype.gene:\n")
          #print(genotype.gene)
          if (verbose & is.null(genotype.gene)) {
            message("    No SNPs in the genomic range.")
          }
          if (!is.null(genotype.gene)) {
            #print("Colnames genotype gene:\n")
            #print(colnames(genotype.gene))
            
            #print("Com.samples:\n")
            #print(com.samples)
            snps.to.keep = check.genotype(genotype.gene[, 
                                                        com.samples], tre.gene[, com.samples])
            #print("Snps.to.keep \n")
            #print(snps.to.keep)
            if (verbose) {
              snps.to.keep.t = table(snps.to.keep)
              message("    ", paste(names(snps.to.keep.t), 
                                    snps.to.keep.t, sep = ":", collapse = ", "))
            }
            if (any(snps.to.keep == "PASS")) {
              #print("check if gorup of 5 samples enters here")
              genotype.gene = genotype.gene[snps.to.keep == 
                                              "PASS", ]
              res.range = dplyr::do(dplyr::group_by(genotype.gene, 
                                                    snpId), sQTLseekeR:::compFscore(., tre.dist, tre.gene, 
                                                                                    svQTL = svQTL))
            }
          }
          #print("res.range is: \n")
          #print(res.range)
          return(res.range)
        }
      })
      range.done = which(unlist(lapply(res.df, nrow)) > 
                           0)
      if (length(range.done) > 0) {
        res.df = res.df[range.done]
        res.df = do.call(rbind, res.df)
        res.df = dplyr::do(dplyr::group_by(res.df, nb.groups), 
                           sQTLseekeR:::compPvalue(., tre.dist, approx = approx, min.nb.ext.scores = min.nb.ext.scores, 
                                      nb.perm.max = nb.perm.max))
        if (svQTL) {
          res.df = dplyr::do(dplyr::group_by(res.df, 
                                             nb.groups), sQTLseekeR:::compPvalue(., tre.dist, svQTL = TRUE, 
                                                                    min.nb.ext.scores = min.nb.ext.scores, nb.perm.max = nb.perm.max.svQTL))
        }
        return(data.frame(done = TRUE, res.df))
      }
    }
    }
    else {
      if (verbose) {
        warning("Issue with the gene location.")
      }
    }
    return(data.frame(done = FALSE))
  }
  ret.df = lapply(unique(tre.df$geneId), function(gene.i) {
    #print("Gene.i:\n")
    #print(gene.i)
    #print("Which gene is in tre.df:\n")
    #print(which(tre.df$geneId == gene.i))
    df = tre.df[which(tre.df$geneId == gene.i), ]
    #print(df)
    data.frame(geneId = gene.i, analyze.gene.f(df))
  })
  done = which(unlist(lapply(ret.df, ncol)) > 2)
  if (length(done) > 0) {
    ret.df = ret.df[done]
    ret.df = do.call(rbind, ret.df)
    ret.df$done = NULL
    return(ret.df)
  }
  else {
    return(NULL)
  }
}