#ifndef BEDTK_LITE_H_
#define BEDTK_LITE_H_

#include "cgranges.h"
#include "kseq.h"
#include "ketopt.h"
KSTREAM_INIT(gzFile, gzread, 0x10000)

cgranges_t *read_bed3(const char *fn);
char *parse_bed3b(char *s, int32_t *st_, int32_t *en_, char **r);

/*
static int x(const char *fp) {
  gzFile f = gzopen(fp, "r");
  if (f == NULL) {
      fprintf(stderr, "\nFailed to open file '%s'!\n", fp);
      return 1;
  }

  kstring_t str;
  uint64_t total = 0;
  kstream_t *ks = ks_init(f);
  char *contig, *rest;
  int32_t start, end;
  while (ks_getuntil(ks, KS_SEP_LINE, &str, 0) >= 0) {
    total++;

    contig = parse_bed3b(str.s, &start, &end, &rest);
    if (contig == NULL) {
        fprintf(stderr, "\nContig not found!\n");
        return 1;
    }
  }
  return 0;
}
*/

#endif
