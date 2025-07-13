\name{seqPGEN2GDS}
\alias{seqPGEN2GDS}
\title{Reformat PLINK2 PGEN files}
\description{
    Reformats PLINK2 pgen files to GDS format.
}

\usage{
seqPGEN2GDS(pgen.fn, pvar.fn, psam.fn, out.gdsfn,
    compress.geno="LZMA_RA", compress.annotation="LZMA_RA",
    start=1L, count=NA_integer_, ignore.chr.prefix=c("chr", "0"),
    reference=NULL, optimize=TRUE, digest=TRUE, parallel=FALSE, verbose=TRUE)
}

\arguments{
    \item{pgen.fn}{a file name for the pgen file}
    \item{pvar.fn}{a file name for the pvar file, or missing}
    \item{psam.fn}{a file name for the psam file, or missing}
    \item{out.gdsfn}{the file name of output GDS file}
    \item{compress.geno}{the compression method for "genotype"; optional
        values are defined in the function \code{add.gdsn}}
    \item{compress.annotation}{the compression method for the GDS variables,
        except "genotype"; optional values are defined in the function
        \code{add.gdsn}}
    \item{ignore.chr.prefix}{a vector of character, indicating the prefix of
        chromosome which should be ignored, e.g., \code{"chr"}; it is not
        case-sensitive}
    \item{start}{the starting variant if importing part of the pgen file}
    \item{count}{the maximum count of variant if importing part of the pgen
        file, -1 indicates importing to the end}
    \item{reference}{genome reference, like "GRCh37", "GRCh38"; it is not
        specified if \code{reference=NULL}}
    \item{optimize}{if \code{TRUE}, optimize the access efficiency by calling
        \code{\link{cleanup.gds}}}
    \item{digest}{a logical value (TRUE/FALSE) or a character ("md5", "sha1",
        "sha256", "sha384" or "sha512"); add md5 hash codes to the GDS file
        if TRUE or a digest algorithm is specified}
    \item{parallel}{\code{FALSE} (serial processing), \code{TRUE} (parallel
        processing), a numeric value indicating the number of cores, or a
        cluster object for parallel processing; \code{parallel} is passed to
        the argument \code{cl} in \code{\link{seqParallel}}, see
        \code{\link{seqParallel}} for more details}
    \item{verbose}{if \code{TRUE}, show information}
}

\value{
    Return the file name of SeqArray GDS file with an absolute path.
}

\references{
    \url{https://bitbucket.org/gavinband/bgen},
    \url{http://bgenformat.org}
}
\author{Xiuwen Zheng}
\seealso{
    \code{\link{seqReadPVAR}}
}

\examples{
pgen_fn <- system.file("extdata", "plink2_gen.pgen", package="pgen2gds")

seqPGEN2GDS(pgen_fn, out.gdsfn="test.gds")

# delete the temporary file
unlink("test.gds", force=TRUE)
}

\keyword{GDS}
\keyword{PLINK2}
\keyword{utilities}
