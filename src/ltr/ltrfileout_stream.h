/*
  Copyright (c) 2008 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2008 Center for Bioinformatics, University of Hamburg

  Permission to use, copy, modify, and distribute this software for any
  purpose with or without fee is hereby granted, provided that the above
  copyright notice and this permission notice appear in all copies.

  THE SOFTWARE IS PROVIDED "AS IS" AND THE AUTHOR DISCLAIMS ALL WARRANTIES
  WITH REGARD TO THIS SOFTWARE INCLUDING ALL IMPLIED WARRANTIES OF
  MERCHANTABILITY AND FITNESS. IN NO EVENT SHALL THE AUTHOR BE LIABLE FOR
  ANY SPECIAL, DIRECT, INDIRECT, OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES
  WHATSOEVER RESULTING FROM LOSS OF USE, DATA OR PROFITS, WHETHER IN AN
  ACTION OF CONTRACT, NEGLIGENCE OR OTHER TORTIOUS ACTION, ARISING OUT OF
  OR IN CONNECTION WITH THE USE OR PERFORMANCE OF THIS SOFTWARE.
*/

#ifndef LTRFILEOUT_STREAM_H
#define LTRFILEOUT_STREAM_H

#include "core/bioseq.h"
#include "extended/node_stream_api.h"
#include "extended/region_mapping.h"
#include "ltr/ltrdigest_def.h"

/* implements the ``genome_stream'' interface */
typedef struct GtLTRFileOutStream GtLTRFileOutStream;

const GtNodeStreamClass* gt_ltr_fileout_stream_class(void);

GtNodeStream* gt_ltr_fileout_stream_new(GtNodeStream *in_stream,
                                        int tests_to_run,
                                        GtRegionMapping *regionmapping,
                                        char *file_prefix,
                                        GtPPTOptions *ppt_opts,
                                        GtPBSOptions *pbs_opts,
#ifdef HAVE_HMMER
                                        GtPdomOptions *pdom_opts,
#endif
                                        const char *trnafilename,
                                        const char *seqfilename,
                                        const char *gfffilename,
                                        unsigned int seqnamelen,
                                        GtError *e);

#endif
