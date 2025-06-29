# ===========================================================================
#
# pgen2gds.r: Format Conversion from PLINK2 PGEN to GDS
#
# Copyright (C) 2025    Xiuwen Zheng (zhengx@u.washington.edu)
#
# This is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License Version 3 as
# published by the Free Software Foundation.
#
# pgen2gds is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with pgen2gds.
# If not, see <http://www.gnu.org/licenses/>.


#############################################################
# Internal functions
#

.cat <- function(...) cat(..., "\n", sep="")

tm <- function() strftime(Sys.time(), "%Y-%m-%d %H:%M:%S")

pretty_size <- function(x)
{
    stopifnot(is.numeric(x), length(x)==1L)
    if (is.na(x) || !is.finite(x))
        "NA"
    else if (x >= 1024^4)
        sprintf("%.1fT", x/1024^4)
    else if (x >= 1024^3)
        sprintf("%.1fG", x/1024^3)
    else if (x >= 1024^2)
        sprintf("%.1fM", x/1024^2)
    else if (x >= 1024)
        sprintf("%.1fK", x/1024)
    else
        sprintf("%gB", x)
}

# append to a gds node using the data return from a user-defined function
append_fc_gds <- function(nd, start, count, type, fc, fc2=identity)
{
    last <- start + count
    i <- start
    while (i < last)
    {
        cnt <- min(i+10000L, last) - i
        s <- vapply(seq.int(i, length.out=cnt), fc, type)
        append.gdsn(nd, fc2(s))
        i <- i + cnt
    }
    readmode.gdsn(nd)
    invisible()
}


#############################################################
# Read PLINK2 pvar file and return a data.frame
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
    optimize=TRUE, digest=TRUE, parallel=FALSE, verbose=TRUE)
{
    # check
    stopifnot(is.character(pgen.fn), length(pgen.fn)==1L)
    if (missing(pvar.fn) && missing(psam.fn))
    {
        fn <- gsub("\\.pgen$", "", pgen.fn, ignore.case=TRUE)
        # set pvar file
        pvar.fn <- paste0(fn, ".pvar.zst")
        if (!file.exists(pvar.fn))
        {
            pvar.fn <- paste0(fn, ".pvar.gz")
            if (!file.exists(pvar.fn))
                pvar.fn <- paste0(fn, ".pvar")
        }
        # set psam file
        psam.fn <- paste0(fn, ".psam.zst")
        if (!file.exists(psam.fn))
        {
            psam.fn <- paste0(fn, ".psam.gz")
            if (!file.exists(psam.fn))
                psam.fn <- paste0(fn, ".psam")
        }
        # check pgen
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

    # check files
    if (!file.exists(pgen.fn)) stop("No file ", sQuote(pgen.fn), ".")
    if (!file.exists(pvar.fn)) stop("No file ", sQuote(pvar.fn), ".")
    if (!file.exists(psam.fn)) stop("No file ", sQuote(psam.fn), ".")

    # open pgen file
    if (verbose)
    {
        .cat("##< ", tm())
        .cat("PLINK2 PGEN to SeqArray GDS:")
        .cat("    pgen file (", pretty_size(file.size(pgen.fn)), "):")
        .cat("        ", pgen.fn)
        .cat("    pvar file (", pretty_size(file.size(pvar.fn)), "):")
        .cat("        ", pvar.fn)
    }
    pgen.fn <- normalizePath(pgen.fn, mustWork=FALSE)
    pvar.fn <- normalizePath(pvar.fn, mustWork=FALSE)
    pvar <- pgen <- NULL
    pvar <- pgenlibr::NewPvar(pvar.fn)
    pgen <- pgenlibr::NewPgen(pgen.fn, pvar=pvar)
    on.exit({
    	if (!is.null(pgen)) ClosePgen(pgen)
    	if (!is.null(pvar)) ClosePvar(pvar)
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

    # read psam file
    if (verbose)
    {
        .cat("    psam file (", pretty_size(file.size(psam.fn)), "):")
        .cat("        ", psam.fn)
        .cat("        reading ...")
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
    if (nsamp != nrow(fam))
        stop("Inconsistent number of rows in psam file.")

    if (verbose)
    {
        .cat("    # of samples: ", nsamp)
        .cat("    # of variants: ", nvar)
        cat("    Output:\n        ", out.gdsfn, "\n", sep="")
        if (start!=1L || count!=nvar)
            cat("        (starting from ", start, ", count: ", count, ")\n", sep="")
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

            # reset memory
            ClosePgen(pgen); pgen <- NULL
            ClosePvar(pvar); pvar <- NULL
            gc(FALSE, reset=TRUE, full=TRUE)

            # conversion in parallel
            seqParallel(parallel, NULL,
                FUN = function(pgen.fn, pvar.fn, psam.fn,
                    compress.geno, compress.annotation, ptmpfn, psplit)
                {
                    # the process id, starting from one
                    i <- SeqArray:::process_index
                    tryCatch(
                    {
                        pgen2gds::seqPGEN2GDS(pgen.fn, pvar.fn, psam.fn,
                            ptmpfn[i], compress.geno=compress.geno,
                            compress.annotation=compress.annotation,
                            start=psplit[[1L]][i], count=psplit[[2L]][i],
                            optimize=FALSE, digest=FALSE, parallel=FALSE,
                            verbose=FALSE)
                    }, error = function(e) {
                        # capture full traceback
                        trace <- capture.output({
                            cat("Error: ", e$message, "\n", sep="")
                            traceback()
                        })
                        con <- file(paste0(ptmpfn[i], ".progress.txt"), open="at")
                        writeLines(trace, con)
                        close(con)
                        stop(e$message)
                    })
                    invisible()  # return
                }, split = "none",
                pgen.fn=pgen.fn, pvar.fn=pvar.fn, psam.fn=psam.fn,
                compress.geno=compress.geno,
                compress.annotation=compress.annotation,
                ptmpfn=ptmpfn, psplit=psplit
            )
            if (verbose)
            {
                cat("    done splitting (", tm(), ")\n", sep="")
                cat("    --------\n")
            }

            # reopen the file
            pvar <- pgenlibr::NewPvar(pvar.fn)

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

    # add chromosome
    if (verbose) cat("    chromosome  ")
    if (length(ignore.chr.prefix))
    {
        ignore.chr.prefix <- paste0("^(",
            paste(ignore.chr.prefix, collapse="|"), ")")
    }
    n <- add.gdsn(dstfile, "chromosome", storage="string",
        compress=compress.annotation)
    append_fc_gds(n, start, count, "", function(i) GetVariantChrom(pvar, i),
        function(s) {
            if (length(ignore.chr.prefix)) s <- gsub(ignore.chr.prefix, "", s)
            s
        })
    SeqArray:::.DigestCode(n, digest, verbose, FALSE)
    # RLE-coded chromosome
    SeqArray:::.optim_chrom(dstfile)

    # add position
    if (verbose) cat("    position  ")
    n <- add.gdsn(dstfile, "position", storage="int32",
        compress=compress.annotation)
    append_fc_gds(n, start, count, 0L, function(i) GetVariantPos(pvar, i))
    SeqArray:::.DigestCode(n, digest, verbose, FALSE)

    # add allele
    if (verbose) cat("    allele  ")
    n <- add.gdsn(dstfile, "allele", storage="string",
        compress=compress.annotation)
    append_fc_gds(n, start, count, "", function(i) {
            paste(vapply(seq_len(GetAlleleCt(pvar, i)), function(j)
                GetAlleleCode(pvar, i, j), ""), collapse=",")
        })
    SeqArray:::.DigestCode(n, digest, verbose, FALSE)

    # add a folder for genotypes & phase
    ngen <- addfolder.gdsn(dstfile, "genotype")
    put.attr.gdsn(ngen, "VariableName", "GT")
    put.attr.gdsn(ngen, "Description", "Genotype")
    npha <- addfolder.gdsn(dstfile, "phase")
    # add annotation folder
    nann <- addfolder.gdsn(dstfile, "annotation")
    # add annotation/id (rsid)
    if (verbose) cat("    annotation/id  ")
    n <- add.gdsn(nann, "id", storage="string", compress=compress.annotation)
    append_fc_gds(n, start, count, "", function(i) GetVariantId(pvar, i))
    SeqArray:::.DigestCode(n, digest, verbose, FALSE)

    # add nodes for genotypes
    if (verbose)
        .cat("    genotype [", tm(), "] ...")
    n_g <- add.gdsn(ngen, "data", storage="bit2", valdim=c(2L, nsamp, 0L),
        compress=compress.geno)
    n_i <- add.gdsn(ngen, "@data", storage="uint8", compress=compress.annotation,
        visible=FALSE)
    n_p <- add.gdsn(npha, "data", storage="bit1", valdim=c(nsamp, 0L),
        compress=compress.geno)

    # progress information
    progfilename <- paste0(out.gdsfn, ".progress.txt")
    progfile <- file(progfilename, "wt")
    progfile_to_rm <- FALSE
    writeLines(paste("#", tm()), progfile)
    writeLines(paste("Input:", basename(pgen.fn)), progfile)
    if (start!=1L || count!=nvar)
        writeLines(paste0("    start: ", start, ", count: ", count), progfile)
    writeLines(paste("Output:", basename(out.gdsfn)), progfile)
    flush(progfile)
    on.exit({
        close(progfile)
        if (progfile_to_rm) unlink(progfilename, force=TRUE)
    }, add=TRUE)

    if (pnum <= 1L)
    {
        # process genotypes
        buf <- IntBuf(pgen)  # an integer vector
        ii <- integer(1L)    # variant index
        ia <- integer(1L)    # allele index
        read_gt_fc <- quote(ReadHardcalls(pgen, buf, ii, ia))
        allele_num_fc <- quote(GetAlleleCt(pvar, ii))
        # call C
        .Call(SEQ_PGEN_Geno_Import, read_gt_fc, allele_num_fc, buf, ii, ia,
            new.env(), dstfile$root, start, count, progfile, verbose)
        # remove temporary variables
        remove(read_gt_fc, allele_num_fc, buf, ii, ia)
        # close the nodes
        readmode.gdsn(n_g)
        if (verbose) cat("      ")
        SeqArray:::.DigestCode(n_g, digest, verbose)
        readmode.gdsn(n_i)
        SeqArray:::.DigestCode(n_i, digest, FALSE)
        SeqArray:::.append_rep_gds(n_p, raw(1L), as.double(count)*nsamp)
        readmode.gdsn(n_p)
        SeqArray:::.DigestCode(n_p, digest, FALSE)

    } else {
        varnm <- c("genotype/data", "genotype/@data", "phase/data")
        # "phase/data"
        # open all temporary files
        for (fn in ptmpfn)
        {
            if (verbose)
                cat("        adding", sQuote(basename(fn)))
            writeLines(paste("Adding", basename(fn)), progfile)
            flush(progfile)
            # open the gds file
            tmpgds <- openfn.gds(fn)
            n <- objdesp.gdsn(index.gdsn(tmpgds, "variant.id"))$dim
            # merge variables
            for (nm in varnm)
                append.gdsn(index.gdsn(dstfile, nm), index.gdsn(tmpgds, nm))
            # close the file
            closefn.gds(tmpgds)
            if (verbose) .cat(" [", tm(), " done]")
            writeLines(paste0("    [", tm(), ", ",
                prettyNum(n, big.mark=",", scientific=FALSE),
                " variants added]"), progfile)
            flush(progfile)
        }

        # remove temporary files
        unlink(ptmpfn, force=TRUE)
        writeLines(paste("Done. #", tm()), progfile)
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
        if (optimize) .cat("Done.  # ", tm()) else cat("Done.\n")
    closefn.gds(dstfile)
    dstfile <- NULL
    progfile_to_rm <- TRUE
    if (optimize)
    {
        if (verbose)
            cat("Optimize the access efficiency ...\n")
        cleanup.gds(out.gdsfn, verbose=verbose)
    }
    if (verbose) .cat("##> ", tm())

    # output
    invisible(normalizePath(out.gdsfn))
}

