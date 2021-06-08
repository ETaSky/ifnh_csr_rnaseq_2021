# Recycle bin for codes
# tpm
ft_1 <- ft_0 %>% select(-c(Chr:Strand)) %>% pivot_longer(cols = -c(1:2), names_to = "sample_id", values_to = "raw_ct") %>% filter(raw_ct>0) %>% mutate(sample_id = gsub("data\\/mapped_STAR\\/Aligned\\.out\\.bam\\:", "", sample_id), sample_id = gsub("-", "\\.", sample_id), l_ef = pmax(0, Length - avg_fragment_len) + 1, ct_rate = raw_ct/l_ef) %>% group_by(sample_id) %>% mutate(tpm = ct_rate/sum(ct_rate)*1E6) %>% ungroup
# vocano plot
plot_volcano=function(dat){
    g = ggplot(data = dat, aes(x=logFC, y = -log10(adj.P.Val), col=diffexpressed, label=difflabel)) + 
        geom_point() + 
        #geom_vline(xintercept=c(-1, 1), col="red") + 
        geom_hline(yintercept=-log10(0.1), col="red") + 
        geom_text_repel() + 
        scale_color_manual(values = c(No="gray70", Yes="red"), name = "differentiated\nexpressed") + 
        labs(x = "Log2 Fold Change", y = "-log10(pValue)") + 
        theme_minimal()
    g
}

# limma codes
gene_to_keep = ft_wk %>% mutate(Flag = tpm>=10) %>% group_by(Geneid) %>% summarise(N_Flag = sum(Flag)) %>% filter(N_Flag>0.2*nsamples) %>% pull(Geneid)

    dat_mat = ft_wk %>% filter(Geneid %in% gene_to_keep) %>% pivot_wider(id_cols = Geneid, names_from = sample_id, values_from = tpm_log2, values_fill = log2(3)) %>% column_to_rownames(var = "Geneid")
    
    grps = mt_wk$BirthMode[match(colnames(dat_mat), mt_wk$TubeID)]

    # dat_mat = ft_wk %>% pivot_wider(id_cols = Geneid, names_from = sample_id, values_from = raw_ct, values_fill = 0) %>% column_to_rownames(var = "Geneid")
    # x = DGEList(counts = dat_mat, genes = ft_0[ft_0$Geneid %in% rownames(dat_mat), c("Geneid", "Length")])
    # 
    # grps = mt_wk$BirthMode[match(colnames(dat_mat), mt_wk$TubeID)]
    # x$samples$group = grps
    # 
    # cpm <- cpm(x)
    # lcpm <- cpm(x, log = T)
    # 
    # L <- mean(x$samples$lib.size) * 1e-6
    # M <- median(x$samples$lib.size) * 1e-6
    # c(L, M)
    # 
    # summary(lcpm)
    # 
    # keep.exprs <- filterByExpr(x, group = grps)
    # x <- x[keep.exprs,, keep.lib.sizes=FALSE]
    # 
    # lcpm.cutoff <- log2(10/M + 2/L)
    # col <- brewer.pal(nsamples, "Paired")
    # par(mfrow=c(1,2))
    # plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
    # title(main="A. Raw data", xlab="Log-cpm")
    # abline(v=lcpm.cutoff, lty=3)
    # for (i in 2:nsamples){
    #     den <- density(lcpm[,i])
    #     lines(den$x, den$y, col=col[i], lwd=2)
    # }
    # 
    # lcpm <- cpm(x, log=TRUE)
    # plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
    # title(main="B. Filtered data", xlab="Log-cpm")
    # abline(v=lcpm.cutoff, lty=3)
    # for (i in 2:nsamples){
    #     den <- density(lcpm[,i])
    #     lines(den$x, den$y, col=col[i], lwd=2)
    # }
    # 
    # x <- calcNormFactors(x, method = "TMM")
    # 
    # lcpm <- cpm(x, log=TRUE)
    # par(mfrow=c(1,2))
    # 
    # cols = BM_colors[grps]
    # plotMDS(lcpm, labels = grps, col = cols)
    # 
    # design <- model.matrix(~ 0 + as.factor(grps))
    # colnames(design) = c("CS", "VF")
    # contrast.matrix = makeContrasts(CS-VF, levels = design)
    # 
    # par(mfrow=c(1,2))
    # v <- voom(x, design, plot=TRUE)
    # v
    # 
    # vfit <- lmFit(v, design)
    # vfit <- contrasts.fit(vfit, contrasts=contrast.matrix)
    # efit <- eBayes(vfit)
    # plotSA(efit, main="Final model: Mean-variance trend")
    # 
    # summary(decideTests(efit))
    
    cols = BM_colors[grps]
    plotMDS(dat_mat, labels = grps, col = cols)

    design <- model.matrix(~ 0 + as.factor(grps))
    colnames(design) = c("CS", "VF")
    rownames(design) = colnames(dat_mat)
    contrast.matrix = makeContrasts(CS-VF, levels = design)
    
    fit <- lmFit(dat_mat, design)
    fit2 <- contrasts.fit(fit, contrasts = contrast.matrix)
    efit <- eBayes(fit2)
    summary(decideTests(efit, adjust.method = "fdr", p.value = 0.10))

    de_gene = topTable(efit, coef = "CS - VF", adjust.method = "fdr", number = Inf) %>% rownames_to_column(var = "gene") %>% mutate(diffexpressed = ifelse(adj.P.Val<=0.1, "Yes", "No"), difflabel = ifelse(diffexpressed=="Yes", gene, NA))
    plot_volcano(de_gene)