# ===========================================================================
#
# gds2pgen.r: Format Conversion from PLINK2 PGEN to GDS
#
# Copyright (C) 2025    Xiuwen Zheng (zhengx@u.washington.edu)
#
# This is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License Version 3 as
# published by the Free Software Foundation.
#
# gds2pgen is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with gds2pgen.
# If not, see <http://www.gnu.org/licenses/>.


#############################################################
# Internal functions
#

.cat <- function(...) cat(..., "\n", sep="")

tm <- function() strftime(Sys.time(), "%Y-%m-%d %H:%M:%S")


#############################################################
# Format conversion from BGEN to GDS
#
seqReadPVAR <- function(pvar, sel=NULL)
{
    # check
    error <- TRUE
    if (is.character(pvar))
    {
        pvar <- NewPvar(pvar)
        on.exit(ClosePvar(pvar))
        error <- FALSE
    }
    if (error && is.list(pvar))
    {
        stopifnot(is.character(pvar$class), pvar$class=="pvar")
        stopifnot(is(pvar$pvar, "externalptr"))
        error <- FALSE
    }
    if (error)
        stop("'pvar' should be a file name or an object returned from pgenlibr::NewPvar().")
    # selection
    if (is.null(sel))
        sel <- seq_len(GetVariantCt(pvar))
    else if (is.logical(sel))
        sel <- which(sel)
    else if (!is.numeric(sel))
        stop("'sel' should be NULL, a logical vector or a numeric vector.")

    # read
    data.frame(
        chrom  = vapply(sel, function(i) GetVariantChrom(pvar, i), ""),
        pos    = vapply(sel, function(i) GetVariantPos(pvar, i), 0L),
        allele = vapply(sel, function(i)
            paste(vapply(seq_len(GetAlleleCt(pvar, i)), function(j)
                GetAlleleCode(pvar, i, j), ""), collapse=","), ""),
        rsid   = vapply(sel, function(i) GetVariantId(pvar, i), "")
    )
}


#############################################################
# Format conversion from PGEN to GDS
#
seqPGEN2GDS <- function(pgen.fn, pvar.fn, psam.fn, out.gdsfn,
    compress.geno="LZMA_RA", compress.annotation="LZMA_RA",
    start=1L, count=NA_integer_, ignore.chr.prefix=c("chr", "0"),
    include.pheno=TRUE, optimize=TRUE, digest=TRUE, parallel=FALSE,
    verbose=TRUE)
{
    # check
    stopifnot(is.character(pgen.fn), length(pgen.fn)==1L)
    if (missing(pvar.fn) && missing(psam.fn))
    {
        fn <- gsub("\\.pgen$", "", pgen.fn, ignore.case=TRUE)
        pvar.fn <- paste0(fn, ".pvar")
        psam.fn <- paste0(fn, ".psam")
        if (!grepl("\\.pgen$", pgen.fn, ignore.case=TRUE))
            pgen.fn <- paste0(fn, ".pgen")
    }
    stopifnot(is.character(pvar.fn), length(pvar.fn)==1L)
    stopifnot(is.character(psam.fn), length(psam.fn)==1L)
    stopifnot(is.character(out.gdsfn), length(out.gdsfn)==1L)
    stopifnot(is.character(compress.geno), length(compress.geno)==1L)
    stopifnot(is.character(compress.annotation), length(compress.annotation)==1L)
    stopifnot(is.numeric(start), length(start)==1L)
    stopifnot(is.numeric(count), length(count)==1L)
    stopifnot(is.character(ignore.chr.prefix))
    stopifnot(is.logical(optimize), length(optimize)==1L)
    stopifnot(is.logical(include.pheno) | is.character(include.pheno))

    # open pgen file
    if (verbose)
    {
        .cat(date())
        .cat("PLINK2 PGEN to SeqArray GDS:")
        .cat("    pgen file: ", SeqArray:::.pretty_size(file.size(pgen.fn)))
        .cat("        ", pgen.fn)
        .cat("    pvar file: ", SeqArray:::.pretty_size(file.size(pvar.fn)))
        .cat("        ", pvar.fn)
    }
    pgen.fn <- normalizePath(pgen.fn, mustWork=FALSE)
    pvar.fn <- normalizePath(pvar.fn, mustWork=FALSE)
    psam.fn <- normalizePath(psam.fn, mustWork=FALSE)
    pvar <- pgen <- NULL
    pvar <- pgenlibr::NewPvar(pvar.fn)
    pgen <- pgenlibr::NewPgen(pgen.fn, pvar=pvar)
    on.exit({
    	if (!is.null(pvar)) ClosePgen(pgen)
    	if (!is.null(pgen)) ClosePvar(pvar)
    })
    nvar <- GetVariantCt(pgen)
    nsamp <- GetRawSampleCt(pgen)

    if (is.na(start) || start<1L) start <- 1L
    if (!is.integer(start)) start <- as.integer(start)
	count <- as.integer(count)
    if (is.na(count) || count<1L)
    {
    	count <- as.integer(nvar - start + 1L)
    } else {
    	count <- min(start + count, nvar + 1L) - start
    }
    last <- start + count

    # read psam file
    if (verbose)
    {
        .cat("    psam file: ", SeqArray:::.pretty_size(file.size(psam.fn)))
        .cat("        ", psam.fn)
        .cat("    # of samples: ", nsamp)
        .cat("    # of variants: ", nvar)
    }
    fam <- read.table(psam.fn, header=TRUE, comment.char="",
        stringsAsFactors=FALSE)
    s <- names(fam)
    s[s=="X.IID"] <- "IID"; s[s=="X.FID"] <- "FID"
    names(fam) <- s
    if (anyDuplicated(fam$IID) == 0L)
    {
        sample.id <- fam$IID
    } else {
        stop("Sample IDs in PLINK2 PGEN are not unique (IID column)!")
    }

    if (verbose)
    {
        cat("    Output:\n        ", out.gdsfn, "\n", sep="")
        if (start!=1L || count!=nvar)
            cat("        (starting from ", start, ", count: ", count, "\n", sep="")
    }

    # the number of parallel tasks
    pnum <- SeqArray:::.NumParallel(parallel)
    if (pnum > 1L)
    {
        if (count >= pnum)
        {
            fn <- sub("^([^.]*).*", "\\1", basename(out.gdsfn))
            psplit <- SeqArray:::.file_split(count, pnum, start)
            # need unique temporary file names
            ptmpfn <- character()
            while (length(ptmpfn) < pnum)
            {
                s <- tempfile(pattern=sprintf("%s_tmp%02d_",
                    fn, length(ptmpfn)+1L), tmpdir=dirname(out.gdsfn))
                file.create(s)
                if (!(s %in% ptmpfn)) ptmpfn <- c(ptmpfn, s)
            }
            if (verbose)
            {
                .cat("    output to path: ", file.path(dirname(out.gdsfn), ""))
                cat(sprintf("    writing to %d files:\n", pnum))
                cat(sprintf("        %s [%s..%s]\n", basename(ptmpfn),
                    SeqArray:::.pretty(psplit[[1L]]),
                    SeqArray:::.pretty(psplit[[1L]] + psplit[[2L]] - 1L)),
                    sep="")
                flush.console()
            }

            # conversion in parallel
            seqParallel(parallel, NULL,
                FUN = function(pgen.fn, pvar.fn, psam.fn,
                    compress.geno, compress.annotation, ptmpfn, psplit)
                {
                    # the process id, starting from one
                    i <- SeqArray:::process_index
                    gds2pgen::seqPGEN2GDS(pgen.fn, pvar.fn, psam.fn, ptmpfn[i],
                        compress.geno = compress.geno,
                        compress.annotation = compress.annotation,
                        start = psplit[[1L]][i], count = psplit[[2L]][i],
                        optimize = FALSE, digest = FALSE, parallel = FALSE,
                        verbose = FALSE
                    )
                }, split = "none",
                pgen.fn=pgen.fn, pvar.fn=pvar.fn, psam.fn=psam.fn,
                compress.geno=compress.geno,
                compress.annotation=compress.annotation,
                ptmpfn=ptmpfn, psplit=psplit
            )
            if (verbose)
            {
                cat("    done splitting (", date(), ").\n", sep="")
                cat("    --------\n")
            }

        } else {
            pnum <- 1L
            message("No use of parallel environment!")
        }
    }

    # create GDS file
    dstfile <- createfn.gds(out.gdsfn)
    on.exit({
        if (!is.null(dstfile)) closefn.gds(dstfile)
    }, add=TRUE)
    put.attr.gdsn(dstfile$root, "FileFormat", "SEQ_ARRAY")
    put.attr.gdsn(dstfile$root, "FileVersion", "v1.0")

    n <- addfolder.gdsn(dstfile, "description")
    put.attr.gdsn(n, "source.format", "PLINK2 PGEN Format")

    # add sample.id
    if (verbose) cat("    sample.id  ")
    n <- add.gdsn(dstfile, "sample.id", sample.id, compress=compress.annotation,
        closezip=TRUE)
    SeqArray:::.DigestCode(n, digest, verbose, FALSE)

    # add variant.id
    if (verbose) cat("    variant.id  ")
    n <- add.gdsn(dstfile, "variant.id", seq.int(start, length.out=count),
        compress=compress.annotation, closezip=TRUE)
    SeqArray:::.DigestCode(n, digest, verbose, FALSE)

    # add chromosome, position, allele
    if (verbose) cat("    chromosome, position, allele, rsid ...\n")
    v <- seqReadPVAR(pvar, seq.int(start, length.out=count))
    allele_max_cnt <- GetMaxAlleleCt(pvar)
    if (verbose)
        .cat("        maximum allele count per variant: ", allele_max_cnt)
    ## chrom
    if (length(ignore.chr.prefix))
    {
        s <- paste0("^(", paste(ignore.chr.prefix, collapse="|"), ")")
        v$chrom <- gsub(s, "", v$chrom)
    }
    n <- add.gdsn(dstfile, "chromosome", v$chrom, compress=compress.annotation,
        closezip=TRUE)
    if (verbose) cat("        ")
    SeqArray:::.DigestCode(n, digest, verbose, FALSE)
    ## pos
    n <- add.gdsn(dstfile, "position", v$pos, compress=compress.annotation,
        closezip=TRUE)
    if (verbose) cat("        ")
    SeqArray:::.DigestCode(n, digest, verbose, FALSE)
    ## allele
    n <- add.gdsn(dstfile, "allele", v$allele, compress=compress.annotation,
        closezip=TRUE)
    if (verbose) cat("        ")
    SeqArray:::.DigestCode(n, digest, verbose, FALSE)

    # add a folder for genotypes & phase
    ngen <- addfolder.gdsn(dstfile, "genotype")
    put.attr.gdsn(ngen, "VariableName", "GT")
    put.attr.gdsn(ngen, "Description", "Genotype")
    npha <- addfolder.gdsn(dstfile, "phase")
    # add annotation folder
    nann <- addfolder.gdsn(dstfile, "annotation")
    # add annotation/id
    n <- add.gdsn(nann, "id", v$rsid, compress=compress.annotation,
        closezip=TRUE)
    if (verbose) cat("        ")
    SeqArray:::.DigestCode(n, digest, verbose, FALSE)
    # remove the variant v to reduce memory usage
    remove(v)
    # RLE-coded chromosome
    SeqArray:::.optim_chrom(dstfile)

    # add nodes for genotypes
    if (verbose)
        .cat("    genotype [", tm(), "] ...")
    n_g <- add.gdsn(ngen, "data", storage="bit2", valdim=c(2L, nsamp, 0L),
        compress=compress.geno)
    n_i <- add.gdsn(ngen, "@data", storage="uint8", compress=compress.annotation,
        visible=FALSE)
    n_p <- add.gdsn(npha, "data", storage="bit1", valdim=c(nsamp, 0L),
        compress=compress.geno)

    # progress bar
    progfile <- file(paste0(out.gdsfn, ".progress"), "wt")
    flush(progfile)
    on.exit({
        close(progfile)
        unlink(paste0(out.gdsfn, ".progress"), force=TRUE)
    }, add=TRUE)

    if (pnum <= 1L)
    {
        # process genotypes
        buf <- IntAlleleCodeBuf(pgen)  # an integer matrix
        buf2 <- BoolBuf(pgen)  # a logical vector
        ii <- integer(1L)
        read_fc <- quote(ReadAlleles(pgen, buf, ii, buf2))
        # call C
        .Call(SEQ_PGEN_Allele_Import, read_fc, buf, ii, buf2, new.env(),
            dstfile$root, start, count, progfile, verbose)
        # close the nodes
        remove(read_fc, buf, ii)
        readmode.gdsn(n_g)
        if (verbose) cat("      ")
        SeqArray:::.DigestCode(n_g, digest, verbose)
        readmode.gdsn(n_i)
        SeqArray:::.DigestCode(n_i, digest, FALSE)

    } else {
        varnm <- c("genotype/data", "genotype/@data", "phase/data")
        # "phase/data"
        # open all temporary files
        for (fn in ptmpfn)
        {
            if (verbose)
                cat("        adding", sQuote(basename(fn)))
            # open the gds file
            tmpgds <- openfn.gds(fn)
            # merge variables
            for (nm in varnm)
                append.gdsn(index.gdsn(dstfile, nm), index.gdsn(tmpgds, nm))
            # close the file
            closefn.gds(tmpgds)
            if (verbose) .cat(" [", tm(), " done]")
        }

        # remove temporary files
        unlink(ptmpfn, force=TRUE)
    }

    # additional nodes for genotype
    n <- add.gdsn(ngen, "extra.index", storage="int32", valdim=c(3L,0L),
        compress=compress.geno, closezip=TRUE)
    put.attr.gdsn(n, "R.colnames",
        c("sample.index", "variant.index", "length"))
    add.gdsn(ngen, "extra", storage="int16", compress=compress.geno, closezip=TRUE)
    # additional nodes for phase data
    n <- add.gdsn(npha, "extra.index", storage="int32", valdim=c(3L,0L),
        compress=compress.geno, closezip=TRUE)
    put.attr.gdsn(n, "R.colnames",
        c("sample.index", "variant.index", "length"))
    add.gdsn(npha, "extra", storage="bit1", compress=compress.geno, closezip=TRUE)

    # add annotation/qual
    n <- add.gdsn(nann, "qual", storage="float", compress=compress.annotation)
    SeqArray:::.append_rep_gds(n, NaN, nvar)
    readmode.gdsn(n)
    if (verbose) cat("    annotation/qual")
    SeqArray:::.DigestCode(n, digest, verbose)

    # add filter
    n <- add.gdsn(nann, "filter", storage="int32", compress=compress.annotation)
    SeqArray:::.append_rep_gds(n, as.raw(1L), nvar)
    readmode.gdsn(n)
    put.attr.gdsn(n, "R.class", "factor")
    put.attr.gdsn(n, "R.levels", c("PASS"))
    put.attr.gdsn(n, "Description", c("All filters passed"))
    if (verbose) cat("    annotation/filter")
    SeqArray:::.DigestCode(n, digest, verbose)

    # add the INFO field
    addfolder.gdsn(nann, "info")
    # add the FORMAT field
    addfolder.gdsn(nann, "format")

    # optimize access efficiency
    if (verbose)
        cat("Done.\n", date(), "\n", sep="")
    closefn.gds(dstfile)
    dstfile <- NULL
    if (optimize)
    {
        if (verbose)
            cat("Optimize the access efficiency ...\n")
        cleanup.gds(out.gdsfn, verbose=verbose)
        if (verbose) cat(date(), "\n", sep="")
    }

    # output
    invisible(normalizePath(out.gdsfn))
}

