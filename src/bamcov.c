/*  modified version of coverage.c -- samtools coverage subcommand

    Copyright (C) 2018,2019 Florian Breitwieser
    Portions copyright (C) 2019-2021 Genome Research Ltd.

    Author: Florian P Breitwieser <florian.bw@gmail.com>

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE.  */

/* This program calculates coverage from multiple BAMs
 * simultaneously, to achieve random access and to use the BED interface.
 * To compile this program separately, you may:
 *
 * gcc -g -Wall -O2 -I. -I../htslib -o bamcov  bamcov.c sam_opts.c sam_utils.c  ../htslib/libhts.a  -lz  -lbz2  -lcurl -llzma -lm  -lpthread -lcrypto
 * 
 */

// C headers
//#include <config.h>

#include <ctype.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>  // variadic functions
#include <limits.h>  // INT_MAX
#include <math.h>    // round
#include <stdbool.h>
#include <string.h>
#include <unistd.h>
#include <assert.h>

#include <sys/ioctl.h>

#include "htslib/sam.h"
#include "htslib/hts.h"
#include "samtools.h"
#include "sam_opts.h"
#include <sys/stat.h>

const char *VERSION = "0.1";

typedef struct {  // auxiliary data structure to hold stats on coverage
    unsigned long long n_covered_bases;
    unsigned long long summed_coverage;
    unsigned long long summed_baseQ;
    unsigned long long summed_mapQ;
    unsigned int n_reads;
    unsigned int n_selected_reads;
    bool covered;
    hts_pos_t beg;
    hts_pos_t end;
    int64_t bin_width;
} stats_aux_t;

typedef struct {  // auxiliary data structure to hold a BAM file
    samFile *fp;     // file handle
    sam_hdr_t *hdr;  // file header
    hts_itr_t *iter; // iterator to a region - NULL for us by default
    int min_mapQ;    // mapQ filter
    int min_len;     // length filter
    int fail_flags;
    int required_flags;
    stats_aux_t *stats;
} bam_aux_t;

static int is_url(const char *s)
{
    static const char uri_scheme_chars[] =
        "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+.-";
    return s[strspn(s, uri_scheme_chars)] == ':';
}

#define MAX_PATH_LEN 1024
int read_file_list(const char *file_list,int *n,char **argv[])
{
    char buf[MAX_PATH_LEN];
    int len, nfiles = 0;
    char **files = NULL;
    struct stat sb;

    *n = 0;
    *argv = NULL;

    FILE *fh = fopen(file_list,"r");
    if ( !fh )
    {
        fprintf(stderr,"%s: %s\n", file_list,strerror(errno));
        return 1;
    }

    files = calloc(nfiles,sizeof(char*));
    nfiles = 0;
    while ( fgets(buf,MAX_PATH_LEN,fh) )
    {
        // allow empty lines and trailing spaces
        len = strlen(buf);
        while ( len>0 && isspace(buf[len-1]) ) len--;
        if ( !len ) continue;

        // check sanity of the file list
        buf[len] = 0;
        if (! (is_url(buf) || stat(buf, &sb) == 0))
        {
            // no such file, check if it is safe to print its name
            int i, safe_to_print = 1;
            for (i=0; i<len; i++)
                if (!isprint(buf[i])) { safe_to_print = 0; break; }
            if ( safe_to_print )
                fprintf(stderr,"The file list \"%s\" appears broken, could not locate: %s\n", file_list,buf);
            else
                fprintf(stderr,"Does the file \"%s\" really contain a list of files and do all exist?\n", file_list);
            return 1;
        }

        nfiles++;
        files = realloc(files,nfiles*sizeof(char*));
        files[nfiles-1] = strdup(buf);
    }
    fclose(fh);
    if ( !nfiles )
    {
        fprintf(stderr,"No files read from %s\n", file_list);
        return 1;
    }
    *argv = files;
    *n    = nfiles;
    return 0;
}


static int usage() {
    fprintf(stdout, "Usage: bamcov [options] in1.bam [in2.bam [...]]\n\n"
            "Input options:\n"
            "  -b, --bam-list FILE     list of input BAM filenames, one per line\n"
            "  -l, --min-read-len INT  ignore reads shorter than INT bp [0]\n"
            "  -q, --min-MQ INT        mapping quality threshold [0]\n"
            "  -Q, --min-BQ INT        base quality threshold [0]\n"
            "  --rf <int|str>          required flags: skip reads with mask bits unset []\n"
            "  --ff <int|str>          filter flags: skip reads with mask bits set \n"
            "                                      [UNMAP,SECONDARY,QCFAIL,DUP]\n"
            "  -d, --depth INT         maximum allowed coverage depth [1000000].\n"
            "                          If 0, depth is set to the maximum integer value,\n"
            "                          effectively removing any depth limit.\n"
            "Output options:\n"
            "  -o, --output FILE       write output to FILE [stdout]\n"
            "  -w, --win INT           size in bp of each bin [100]\n"
            "  -r, --region REG        show specified region. Format: chr:start-end. \n"
            "  -h, --help              help (this page)\n");

    fprintf(stdout, "\nGeneric options:\n");
    sam_global_opt_help(stdout, "-.--.--.");

    return EXIT_SUCCESS;
}


// read one alignment from one BAM file
static int read_bam(void *data, bam1_t *b) {
    bam_aux_t *aux = (bam_aux_t*)data; // data in fact is a pointer to an auxiliary structure
    int nref = sam_hdr_nref(aux->hdr);
    int ret;
    while (1) {
        if((ret = aux->iter? sam_itr_next(aux->fp, aux->iter, b) : sam_read1(aux->fp, aux->hdr, b)) < 0) break;
        if (b->core.tid >= 0 && b->core.tid < nref)
            aux->stats[b->core.tid].n_reads++;

        if ( aux->fail_flags && (b->core.flag & aux->fail_flags) ) continue;
        if ( aux->required_flags && !(b->core.flag & aux->required_flags) ) continue;
        if ( b->core.qual < aux->min_mapQ ) continue;
        if ( aux->min_len && bam_cigar2qlen(b->core.n_cigar, bam_get_cigar(b)) < aux->min_len ) continue;
        if (b->core.tid >= 0 && b->core.tid < nref) {
            aux->stats[b->core.tid].n_selected_reads++;
            aux->stats[b->core.tid].summed_mapQ += b->core.qual;
        }
        break;
    }
    return ret;
}

void print_hist(FILE *file_out, const sam_hdr_t *h, const stats_aux_t *stats, int tid, const uint32_t *hist,
        const int hist_size, const bool full_utf) {
    int i;

    for (i = 0; i < hist_size ; ++i) {
      fprintf(file_out,"%s\t%li\t%li\t%i\n",  
          sam_hdr_tid2name(h, tid), i*stats[tid].bin_width + stats[tid].beg,(i+1)*stats[tid].bin_width  + stats[tid].beg, hist[i]);
    }

}

int main_coverage(int argc, char *argv[]) {
    int status = EXIT_SUCCESS;

    int ret, tid = -1, old_tid = -1, pos, i, j;

    int max_depth = 1000000;
    int opt_min_baseQ = 0;
    int opt_min_mapQ = 0;
    int opt_min_len = 0;
    int opt_win = 100;
    char *opt_output_file = NULL;
    bam_aux_t **data = NULL;
    bam_mplp_t mplp = NULL;
    const bam_pileup1_t **plp = NULL;
    uint32_t *hist = NULL;
    stats_aux_t *stats = NULL;
    char *opt_reg = 0; // specified region
    char *opt_file_list = NULL;
    int n_bam_files = 0;
    char **fn = NULL;
    int fail_flags = (BAM_FUNMAP | BAM_FSECONDARY | BAM_FQCFAIL | BAM_FDUP); // Default fail flags
    int required_flags = 0;

    int *n_plp = NULL;
    sam_hdr_t *h = NULL; // BAM header of the 1st input

    bool opt_full_utf = true;

    FILE *file_out = stdout;

    sam_global_args ga = SAM_GLOBAL_ARGS_INIT;
    static const struct option lopts[] = {
        SAM_OPT_GLOBAL_OPTIONS('-', 0, '-', '-', 0, '-'),
        {"rf", required_argument, NULL, 1}, // require flag
        {"ff", required_argument, NULL, 2}, // filter flag
        {"incl-flags", required_argument, NULL, 1}, // require flag
        {"excl-flags", required_argument, NULL, 2}, // filter flag
        {"bam-list", required_argument, NULL, 'b'},
        {"min-read-len", required_argument, NULL, 'l'},
        {"min-MQ", required_argument, NULL, 'q'},
        {"min-mq", required_argument, NULL, 'q'},
        {"min-BQ", required_argument, NULL, 'Q'},
        {"min-bq", required_argument, NULL, 'Q'},
        {"output", required_argument, NULL, 'o'},
        {"win", required_argument, NULL, 'w'},
        {"region", required_argument, NULL, 'r'},
        {"help", no_argument, NULL, 'h'},
        {"depth", required_argument, NULL, 'd'},
        { NULL, 0, NULL, 0 }
    };

    // parse the command line
    int c;
    opterr = 0;
    while ((c = getopt_long(argc, argv, "o:l:q:Q:hw:r:b:d:", lopts, NULL)) != -1) {
        switch (c) {
            case 1:
                if ((required_flags = bam_str2flag(optarg)) < 0) {
                    fprintf(stderr,"Could not parse --rf %s\n", optarg); return EXIT_FAILURE;
                }; break;
            case 2:
                if ((fail_flags = bam_str2flag(optarg)) < 0) {
                    fprintf(stderr,"Could not parse --ff %s\n", optarg); return EXIT_FAILURE;
                }; break;
            case 'o': opt_output_file = optarg; break;
            case 'l': opt_min_len = atoi(optarg); break;
            case 'q': opt_min_mapQ = atoi(optarg); break;
            case 'Q': opt_min_baseQ = atoi(optarg); break;
            case 'd': max_depth = atoi(optarg); break; // maximum coverage depth
            case 'w': opt_win = atoi(optarg); break;
            case 'r': opt_reg = optarg; break;   // parsing a region requires a BAM header (strdup unnecessary)
            case 'b': opt_file_list = optarg; break;
            case 'h': return usage();
            default:  if (parse_sam_global_opt(c, optarg, lopts, &ga) == 0) break;
                          /* else fall-through */
            case '?':
                if (optopt != '?') {  // '-?' appeared on command line
                    if (optopt) { // Bad short option
                        print_error("coverage", "invalid option -- '%c'", optopt);
                    } else { // Bad long option
                        // Do our best.  There is no good solution to finding
                        // out what the bad option was.
                        // See, e.g. https://stackoverflow.com/questions/2723888/where-does-getopt-long-store-an-unrecognized-option
                        if (optind > 0 && strncmp(argv[optind - 1], "--", 2) == 0) {
                            print_error("coverage", "unrecognised option '%s'",
                                        argv[optind - 1]);
                        }
                    }
                }
                return usage();
        }
    }
    if (optind == argc && !opt_file_list)
        return usage();

    // output file provided by user
    if (opt_output_file != NULL && strcmp(opt_output_file,"-")!=0) {
        file_out = fopen( opt_output_file, "w" );
        if (file_out == NULL) {
            print_error_errno("coverage", "Cannot open \"%s\" for writing.", opt_output_file);
            return EXIT_FAILURE;
        }
    }

    // Open all BAM files
    if (opt_file_list) {
        // Read file names from opt_file_list into argv, and record the number of files in n_bam_files
        if (read_file_list(opt_file_list, &n_bam_files, &fn)) {
            print_error_errno("coverage", "Cannot open file list \"%s\".", opt_file_list);
            return EXIT_FAILURE;
        }
        argv = fn;
        optind = 0;
    } else {
        n_bam_files = argc - optind; // the number of BAMs on the command line
    }

    data = (bam_aux_t **)calloc(n_bam_files, sizeof(bam_aux_t*)); // data[i] for the i-th BAM file
    if (!data) {
        print_error_errno("coverage", "Failed to allocate memory");
        status = EXIT_FAILURE;
        goto coverage_end;
    }

    for (i = 0; i < n_bam_files; ++i) {
        int rf;
        data[i] = (bam_aux_t *) calloc(1, sizeof(bam_aux_t));
        if (!data[i]) {
            print_error_errno("coverage", "Failed to allocate memory");
            status = EXIT_FAILURE;
            goto coverage_end;
        }
        data[i]->fp = sam_open_format(argv[optind+i], "r", &ga.in); // open BAM

        if (data[i]->fp == NULL) {
            print_error_errno("coverage", "Could not open \"%s\"", argv[optind+i]);
            status = EXIT_FAILURE;
            goto coverage_end;
        }
        rf = SAM_FLAG | SAM_RNAME | SAM_POS | SAM_MAPQ | SAM_CIGAR | SAM_SEQ;
        if (opt_min_baseQ) rf |= SAM_QUAL;

        // Set CRAM options on file handle - returns 0 on success
        if (hts_set_opt(data[i]->fp, CRAM_OPT_REQUIRED_FIELDS, rf)) {
            print_error("coverage", "Failed to set CRAM_OPT_REQUIRED_FIELDS value");
            status = EXIT_FAILURE;
            goto coverage_end;
        }
        if (hts_set_opt(data[i]->fp, CRAM_OPT_DECODE_MD, 0)) {
            print_error("coverage", "Failed to set CRAM_OPT_DECODE_MD value");
            status = EXIT_FAILURE;
            goto coverage_end;
        }
        data[i]->min_mapQ = opt_min_mapQ;            // set the mapQ filter
        data[i]->min_len  = opt_min_len;             // set the qlen filter
        data[i]->hdr = sam_hdr_read(data[i]->fp);    // read the BAM header
        data[i]->fail_flags = fail_flags;
        data[i]->required_flags = required_flags;
        if (data[i]->hdr == NULL) {
            print_error_errno("coverage", "Could not read header for \"%s\"", argv[optind+i]);
            status = EXIT_FAILURE;
            goto coverage_end;
        }

        // Lookup region if specified
        if (opt_reg) { // if a region is specified
            hts_idx_t *idx = sam_index_load(data[i]->fp, argv[optind+i]);  // load the index
            if (idx == NULL) {
                print_error_errno("coverage", "Failed to load index for \"%s\"", argv[optind+i]);
                status = EXIT_FAILURE;
                goto coverage_end;
            }
            data[i]->iter = sam_itr_querys(idx, data[i]->hdr, opt_reg); // set the iterator
            hts_idx_destroy(idx); // the index is not needed any more; free the memory
            if (data[i]->iter == NULL) {
                print_error("coverage", "Failed to parse region \"%s\". Check the region format or region name presence in the file \"%s\"", opt_reg, argv[optind+i]);
                status = EXIT_FAILURE;
                goto coverage_end;
            }
        }
    }

    h = data[0]->hdr; // easy access to the header of the 1st BAM
    int n_targets = sam_hdr_nref(h);
    stats = calloc(n_targets, sizeof(stats_aux_t));
    if (!stats) {
        print_error_errno("coverage", "Failed to allocate memory");
        status = EXIT_FAILURE;
        goto coverage_end;
    }

    int64_t n_bins =0;
    if (opt_reg) {
        stats_aux_t *s = stats + data[0]->iter->tid;
        s->beg = data[0]->iter->beg; // and to the parsed region coordinates
        s->end = data[0]->iter->end;
        if (s->end == HTS_POS_MAX) {
            s->end = sam_hdr_tid2len(h, data[0]->iter->tid);
        }
        n_bins = ceil( 1.0*(s->end - s->beg)/opt_win);
        s->bin_width = opt_win;
    }

    for (i=0; i<n_bam_files; i++)
        data[i]->stats = stats;

    int64_t current_bin = 0;

    // the core multi-pileup loop
    mplp = bam_mplp_init(n_bam_files, read_bam, (void**)data); // initialization
    if (max_depth > 0)
        bam_mplp_set_maxcnt(mplp, max_depth);  // set maximum coverage depth
    else if (!max_depth)
        bam_mplp_set_maxcnt(mplp, INT_MAX);


    // Extra info for histogram and coverage counting
    n_plp = (int*) calloc(n_bam_files, sizeof(int*)); // n_plp[i] is the number of covering reads from the i-th BAM
    plp = (const bam_pileup1_t**) calloc(n_bam_files, sizeof(bam_pileup1_t*)); // plp[i] points to the array of covering reads (internal in mplp)
    if ( !n_plp || !plp) {
        print_error_errno("coverage", "Failed to allocate memory");
        status = EXIT_FAILURE;
        goto coverage_end;
    }
    while ((ret=bam_mplp_auto(mplp, &tid, &pos, n_plp, plp)) > 0) { // come to the next covered position

        if (tid != old_tid) { // Next target sequence

            stats[tid].covered = true;
            if (!opt_reg)
                stats[tid].end = sam_hdr_tid2len(h, tid);

            if (old_tid >= 0) {
                print_hist(file_out, h, stats, old_tid, hist, n_bins, opt_full_utf);
                if (hist) free(hist);
            }

            if (!opt_reg)
                stats[tid].end = sam_hdr_tid2len(h, tid);

            n_bins = ceil( 1.0*(stats[tid].end-stats[tid].beg )/opt_win);
            stats[tid].bin_width = opt_win;
            hist = (uint32_t*) calloc(n_bins, sizeof(uint32_t));
            if ( !hist ) {
                print_error_errno("coverage", "Failed to allocate memory");
                status = EXIT_FAILURE;
                goto coverage_end;
            }
            old_tid = tid;
        }
        if (pos < stats[tid].beg || pos >= stats[tid].end) continue; // out of range; skip
        if (tid >= n_targets) continue;     // diff number of @SQ lines per file?

        current_bin = (pos - stats[tid].beg) / stats[tid].bin_width;

        bool count_base = false;
        for (i = 0; i < n_bam_files; ++i) { // base level filters have to go here
            int depth_at_pos = n_plp[i];
            for (j = 0; j < n_plp[i]; ++j) {
                const bam_pileup1_t *p = plp[i] + j; // DON'T modify plp[][] unless you really know

                if (p->is_del || p->is_refskip) --depth_at_pos; // having dels or refskips at tid:pos
                else if (p->qpos < p->b->core.l_qseq &&
                        bam_get_qual(p->b)[p->qpos] < opt_min_baseQ) --depth_at_pos; // low base quality
                else
                    stats[tid].summed_baseQ += bam_get_qual(p->b)[p->qpos];
            }
            if (current_bin < n_bins)
              hist[current_bin] += depth_at_pos;
          
            if (depth_at_pos > 0) {
                count_base = true;
                stats[tid].summed_coverage += depth_at_pos;
            }
        }
        if (count_base) {
            stats[tid].n_covered_bases++;
        }
    }

    if (tid == -1 && opt_reg && *opt_reg != '*') {
        // Region specified but no data covering it.
        tid = data[0]->iter->tid;
        n_bins = ceil( 1.0*(stats[tid].end-stats[tid].beg )/opt_win);
        stats[tid].bin_width = opt_win;
        hist = (uint32_t*) calloc(n_bins, sizeof(uint32_t));
    }

    if (tid < n_targets && tid >=0) {
        print_hist(file_out, h, stats, tid, hist, n_bins, opt_full_utf);
    }

    if (ret < 0) status = EXIT_FAILURE;

coverage_end:
    if (n_plp) free(n_plp);
    if (plp) free(plp);
    if (mplp) bam_mplp_destroy(mplp);

    if (hist) free(hist);
    if (stats) free(stats);

    // Close files and free data structures
    if (!(file_out == stdout || fclose(file_out) == 0)) {
        if (status == EXIT_SUCCESS) {
            print_error_errno("coverage", "error on closing \"%s\"",
                    (opt_output_file && strcmp(opt_output_file, "-") != 0?
                     opt_output_file : "stdout"));
            status = EXIT_FAILURE;
        }
    }

    if (data) {
        for (i = 0; i < n_bam_files && data[i]; ++i) {
            sam_hdr_destroy(data[i]->hdr);
            if (data[i]->fp) sam_close(data[i]->fp);
            hts_itr_destroy(data[i]->iter);
            free(data[i]);
        }
        free(data);
    }

    if (opt_file_list && fn) {
        for (i = 0; i < n_bam_files; ++i)
            free(fn[i]);
        free(fn);
    }
    sam_global_args_free(&ga);

    return status;
}

int main(int argc, char *argv[]) {
    return main_coverage(argc, argv);
}
