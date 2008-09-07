/*
  Copyright (c) 2007-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2007-2008 Center for Bioinformatics, University of Hamburg

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

#ifndef REGION_MAPPING_H
#define REGION_MAPPING_H

#include "core/str.h"

/* maps a sequence-region to a sequence file */
typedef struct RegionMapping RegionMapping;

RegionMapping* region_mapping_new_mapping(GT_Str *mapping_filename, GT_Error*);
RegionMapping* region_mapping_new_seqfile(GT_Str *sequence_filename);
RegionMapping* region_mapping_ref(RegionMapping*);
int            region_mapping_get_raw_sequence(RegionMapping*, const char**,
                                               GT_Str *seqid, GT_Error*);
int            region_mapping_get_raw_sequence_length(RegionMapping*,
                                                      unsigned long*,
                                                      GT_Str *seqid, GT_Error*);
void           region_mapping_delete(RegionMapping*);

#endif
