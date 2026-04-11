# test-pgen2gds.R — tests for pgen2gds package

pgen_fn <- system.file("extdata", "plink2_gen.pgen", package = "pgen2gds")
pvar_fn <- system.file("extdata", "plink2_gen.pvar", package = "pgen2gds")
psam_fn <- system.file("extdata", "plink2_gen.psam", package = "pgen2gds")


# ---- seqReadPVAR -----------------------------------------------------------

test_that("seqReadPVAR reads pvar file correctly",
{
    df <- seqReadPVAR(pvar_fn)
    expect_s3_class(df, "data.frame")
    expect_named(df, c("chrom", "pos", "allele", "rsid"))
    expect_true(nrow(df) > 0L)
    expect_type(df$pos, "integer")
    expect_type(df$chrom, "character")
    expect_type(df$allele, "character")
    expect_type(df$rsid, "character")
    expect_true(all(nchar(df$allele) > 0L))
})

test_that("seqReadPVAR with logical sel returns subset",
{
    df_all <- seqReadPVAR(pvar_fn)
    sel <- rep(FALSE, nrow(df_all))
    sel[1:5] <- TRUE
    df_sub <- seqReadPVAR(pvar_fn, sel = sel)
    expect_equal(nrow(df_sub), 5L)
    expect_equal(df_sub$pos, df_all$pos[1:5])
    expect_equal(df_sub$chrom, df_all$chrom[1:5])
})

test_that("seqReadPVAR with numeric sel returns subset",
{
    df_all <- seqReadPVAR(pvar_fn)
    df_sub <- seqReadPVAR(pvar_fn, sel = 1:3)
    expect_equal(nrow(df_sub), 3L)
    expect_equal(df_sub$pos, df_all$pos[1:3])
    expect_equal(df_sub$rsid, df_all$rsid[1:3])
})

test_that("seqReadPVAR rejects invalid input types",
{
    expect_error(seqReadPVAR(42), "'pvar' should be a file name")
    expect_error(seqReadPVAR(NULL), "'pvar' should be a file name")
})

test_that("seqReadPVAR rejects invalid sel types",
{
    expect_error(seqReadPVAR(pvar_fn, sel = "bad"), "'sel' should be")
})


# ---- seqPGEN2GDS basic conversion ------------------------------------------

test_that("seqPGEN2GDS creates a valid GDS file",
{
    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE),
        add=TRUE)

    result <- seqPGEN2GDS(pgen_fn, out.gdsfn=out, optimize=FALSE,
        verbose=FALSE)
    expect_true(file.exists(out))
    expect_equal(normalizePath(result), normalizePath(out))

    gds <- SeqArray::seqOpen(out)
    on.exit(SeqArray::seqClose(gds), add=TRUE)

    samp <- SeqArray::seqGetData(gds, "sample.id")
    vars <- SeqArray::seqGetData(gds, "variant.id")
    chrom <- SeqArray::seqGetData(gds, "chromosome")
    pos <- SeqArray::seqGetData(gds, "position")
    allele <- SeqArray::seqGetData(gds, "allele")

    expect_true(length(samp) > 0L)
    expect_true(length(vars) > 0L)
    expect_equal(length(chrom), length(vars))
    expect_equal(length(pos), length(vars))
    expect_equal(length(allele), length(vars))
})

test_that("seqPGEN2GDS auto-derives pvar/psam from pgen path",
{
    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE),
        add=TRUE)

    # only pass pgen.fn, let pvar and psam be auto-derived
    result <- seqPGEN2GDS(pgen_fn, out.gdsfn=out, optimize=FALSE,
        verbose=FALSE)
    expect_true(file.exists(out))

    gds <- SeqArray::seqOpen(out)
    on.exit(SeqArray::seqClose(gds), add=TRUE)
    expect_true(length(SeqArray::seqGetData(gds, "variant.id")) > 0L)
})


# ---- variant.sel ------------------------------------------------------------

test_that("seqPGEN2GDS with numeric variant.sel",
{
    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE),
        add=TRUE)

    seqPGEN2GDS(pgen_fn, out.gdsfn=out, variant.sel = 1:10,
        optimize=FALSE, verbose=FALSE)

    gds <- SeqArray::seqOpen(out)
    on.exit(SeqArray::seqClose(gds), add=TRUE)
    expect_equal(length(SeqArray::seqGetData(gds, "variant.id")), 10L)
})

test_that("seqPGEN2GDS with logical variant.sel",
{
    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE),
        add=TRUE)

    # get total variant count from pvar
    df <- seqReadPVAR(pvar_fn)
    sel <- rep(FALSE, nrow(df))
    sel[c(1, 5, 10)] <- TRUE

    seqPGEN2GDS(pgen_fn, out.gdsfn=out, variant.sel = sel,
        optimize=FALSE, verbose=FALSE)

    gds <- SeqArray::seqOpen(out)
    on.exit(SeqArray::seqClose(gds), add=TRUE)
    expect_equal(length(SeqArray::seqGetData(gds, "variant.id")), 3L)
})


# ---- start / count ---------------------------------------------------------

test_that("seqPGEN2GDS with start and count",
{
    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE),
        add=TRUE)

    seqPGEN2GDS(pgen_fn, out.gdsfn=out, start = 5L, count = 10L,
        optimize=FALSE, verbose=FALSE)

    gds <- SeqArray::seqOpen(out)
    on.exit(SeqArray::seqClose(gds), add=TRUE)

    vid <- SeqArray::seqGetData(gds, "variant.id")
    expect_equal(length(vid), 10L)
    expect_equal(vid[1], 5L)
})

test_that("seqPGEN2GDS with start only (count defaults to remainder)",
{
    out_all <- tempfile(fileext=".gds")
    out_sub <- tempfile(fileext=".gds")
    on.exit(unlink(c(out_all, out_sub,
        paste0(out_all, ".progress.txt"),
        paste0(out_sub, ".progress.txt")), force=TRUE), add=TRUE)

    seqPGEN2GDS(pgen_fn, out.gdsfn = out_all, optimize=FALSE, verbose=FALSE)
    gds_all <- SeqArray::seqOpen(out_all)
    on.exit(SeqArray::seqClose(gds_all), add=TRUE)
    total <- length(SeqArray::seqGetData(gds_all, "variant.id"))

    seqPGEN2GDS(pgen_fn, out.gdsfn = out_sub, start = 11L,
        optimize=FALSE, verbose=FALSE)
    gds_sub <- SeqArray::seqOpen(out_sub)
    on.exit(SeqArray::seqClose(gds_sub), add=TRUE)
    expect_equal(length(SeqArray::seqGetData(gds_sub, "variant.id")),
        total - 10L)
})


# ---- chromosome prefix stripping -------------------------------------------

test_that("seqPGEN2GDS strips chr prefix by default",
{
    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE),
        add=TRUE)

    seqPGEN2GDS(pgen_fn, out.gdsfn=out, optimize=FALSE, verbose=FALSE)

    gds <- SeqArray::seqOpen(out)
    on.exit(SeqArray::seqClose(gds), add=TRUE)
    chrom <- SeqArray::seqGetData(gds, "chromosome")
    # none should start with "chr" after stripping
    expect_false(any(grepl("^chr", chrom, ignore.case=TRUE)))
})

test_that("seqPGEN2GDS preserves chr prefix when ignore.chr.prefix is empty",
{
    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE),
        add=TRUE)

    seqPGEN2GDS(pgen_fn, out.gdsfn=out, ignore.chr.prefix = character(0),
        optimize=FALSE, verbose=FALSE)

    gds <- SeqArray::seqOpen(out)
    on.exit(SeqArray::seqClose(gds), add=TRUE)
    chrom <- SeqArray::seqGetData(gds, "chromosome")
    # read raw chroms from pvar for comparison
    df <- seqReadPVAR(pvar_fn)
    expect_equal(chrom, df$chrom)
})


# ---- genotype data consistency ---------------------------------------------

test_that("seqPGEN2GDS genotype data has correct dimensions",
{
    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE),
        add=TRUE)

    seqPGEN2GDS(pgen_fn, out.gdsfn=out, optimize=FALSE, verbose=FALSE)

    gds <- SeqArray::seqOpen(out)
    on.exit(SeqArray::seqClose(gds), add=TRUE)

    nsamp <- length(SeqArray::seqGetData(gds, "sample.id"))
    nvar <- length(SeqArray::seqGetData(gds, "variant.id"))

    # genotypes should be retrievable
    SeqArray::seqSetFilter(gds, variant.sel = 1:min(5, nvar), verbose=FALSE)
    geno <- SeqArray::seqGetData(gds, "genotype")
    expect_true(is.array(geno))
    # SeqArray genotype: ploidy x nsamp x nvar_selected
    expect_equal(dim(geno)[2], nsamp)
})


# ---- annotation nodes ------------------------------------------------------

test_that("seqPGEN2GDS creates annotation/id with correct length",
{
    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE),
        add=TRUE)

    seqPGEN2GDS(pgen_fn, out.gdsfn=out, variant.sel = 1:20,
        optimize=FALSE, verbose=FALSE)

    gds <- SeqArray::seqOpen(out)
    on.exit(SeqArray::seqClose(gds), add=TRUE)
    rsid <- SeqArray::seqGetData(gds, "annotation/id")
    expect_equal(length(rsid), 20L)
    expect_type(rsid, "character")
})


# ---- error handling ---------------------------------------------------------

test_that("seqPGEN2GDS rejects missing files",
{
    expect_error(
        seqPGEN2GDS("nonexistent.pgen", "nonexistent.pvar",
            "nonexistent.psam", out.gdsfn = "out.gds", verbose=FALSE),
        "No file")
})

test_that("seqPGEN2GDS rejects empty variant.sel",
{
    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE),
        add=TRUE)
    expect_error(
        seqPGEN2GDS(pgen_fn, out.gdsfn=out, variant.sel = integer(0),
            verbose=FALSE),
        "no selected variant")
})

test_that("seqPGEN2GDS rejects variant.sel with NA",
{
    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE),
        add=TRUE)
    expect_error(
        seqPGEN2GDS(pgen_fn, out.gdsfn=out,
            variant.sel = c(1L, NA_integer_, 3L), verbose=FALSE),
        "should not have NA")
})

test_that("seqPGEN2GDS rejects variant.sel with duplicates",
{
    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE),
        add=TRUE)
    expect_error(
        seqPGEN2GDS(pgen_fn, out.gdsfn=out,
            variant.sel = c(1L, 2L, 2L, 3L), verbose=FALSE),
        "should not have any duplicate")
})

test_that("seqPGEN2GDS rejects out-of-range variant.sel",
{
    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE),
        add=TRUE)
    expect_error(
        seqPGEN2GDS(pgen_fn, out.gdsfn=out,
            variant.sel = c(0L, 1L), verbose=FALSE),
        "should be between")
})


# ---- sample.sel -------------------------------------------------------------

test_that("seqPGEN2GDS with numeric sample.sel",
{
    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE),
        add=TRUE)

    seqPGEN2GDS(pgen_fn, out.gdsfn=out, sample.sel = 1:10,
        optimize=FALSE, verbose=FALSE)

    gds <- SeqArray::seqOpen(out)
    on.exit(SeqArray::seqClose(gds), add=TRUE)
    samp <- SeqArray::seqGetData(gds, "sample.id")
    expect_equal(length(samp), 10L)
})

test_that("seqPGEN2GDS with logical sample.sel",
{
    # get total sample count from a full conversion
    out_all <- tempfile(fileext=".gds")
    on.exit(unlink(c(out_all, paste0(out_all, ".progress.txt")), force=TRUE),
        add=TRUE)
    seqPGEN2GDS(pgen_fn, out.gdsfn = out_all, optimize=FALSE, verbose=FALSE)
    gds_all <- SeqArray::seqOpen(out_all)
    on.exit(SeqArray::seqClose(gds_all), add=TRUE)
    all_samp <- SeqArray::seqGetData(gds_all, "sample.id")
    nsamp_total <- length(all_samp)

    sel <- rep(FALSE, nsamp_total)
    sel[c(1, 5, 10)] <- TRUE

    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE),
        add=TRUE)
    seqPGEN2GDS(pgen_fn, out.gdsfn=out, sample.sel = sel,
        optimize=FALSE, verbose=FALSE)

    gds <- SeqArray::seqOpen(out)
    on.exit(SeqArray::seqClose(gds), add=TRUE)
    samp <- SeqArray::seqGetData(gds, "sample.id")
    expect_equal(length(samp), 3L)
    expect_equal(samp, all_samp[c(1, 5, 10)])
})

test_that("seqPGEN2GDS sample.sel produces correct sample IDs",
{
    out_all <- tempfile(fileext=".gds")
    out_sub <- tempfile(fileext=".gds")
    on.exit(unlink(c(out_all, out_sub,
        paste0(out_all, ".progress.txt"),
        paste0(out_sub, ".progress.txt")), force=TRUE), add=TRUE)

    seqPGEN2GDS(pgen_fn, out.gdsfn = out_all, optimize=FALSE,
        verbose=FALSE)
    gds_all <- SeqArray::seqOpen(out_all)
    on.exit(SeqArray::seqClose(gds_all), add=TRUE)
    all_samp <- SeqArray::seqGetData(gds_all, "sample.id")

    idx <- c(3L, 7L, 15L, 20L)
    seqPGEN2GDS(pgen_fn, out.gdsfn = out_sub, sample.sel = idx,
        optimize=FALSE, verbose=FALSE)
    gds_sub <- SeqArray::seqOpen(out_sub)
    on.exit(SeqArray::seqClose(gds_sub), add=TRUE)
    sub_samp <- SeqArray::seqGetData(gds_sub, "sample.id")
    expect_equal(sub_samp, all_samp[idx])
})

test_that("seqPGEN2GDS sample.sel genotype dimensions are correct",
{
    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE),
        add=TRUE)

    seqPGEN2GDS(pgen_fn, out.gdsfn=out, sample.sel = 1:5,
        variant.sel = 1:10, optimize=FALSE, verbose=FALSE)

    gds <- SeqArray::seqOpen(out)
    on.exit(SeqArray::seqClose(gds), add=TRUE)
    expect_equal(length(SeqArray::seqGetData(gds, "sample.id")), 5L)
    expect_equal(length(SeqArray::seqGetData(gds, "variant.id")), 10L)

    SeqArray::seqSetFilter(gds, variant.sel = 1:5, verbose=FALSE)
    geno <- SeqArray::seqGetData(gds, "genotype")
    # ploidy x nsamp x nvar
    expect_equal(dim(geno)[2], 5L)
})

test_that("seqPGEN2GDS rejects empty sample.sel",
{
    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE),
        add=TRUE)
    expect_error(
        seqPGEN2GDS(pgen_fn, out.gdsfn=out, sample.sel = integer(0),
            verbose=FALSE),
        "no selected sample")
})

test_that("seqPGEN2GDS rejects sample.sel with NA",
{
    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE),
        add=TRUE)
    expect_error(
        seqPGEN2GDS(pgen_fn, out.gdsfn=out,
            sample.sel = c(1L, NA_integer_, 3L), verbose=FALSE),
        "should not have NA")
})

test_that("seqPGEN2GDS rejects sample.sel with duplicates",
{
    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE),
        add=TRUE)
    expect_error(
        seqPGEN2GDS(pgen_fn, out.gdsfn=out,
            sample.sel = c(1L, 2L, 2L, 3L), verbose=FALSE),
        "should not have any duplicate")
})

test_that("seqPGEN2GDS rejects out-of-range sample.sel",
{
    out <- tempfile(fileext=".gds")
    on.exit(unlink(c(out, paste0(out, ".progress.txt")), force=TRUE),
        add=TRUE)
    expect_error(
        seqPGEN2GDS(pgen_fn, out.gdsfn=out,
            sample.sel = c(0L, 1L), verbose=FALSE),
        "should be between")
})
