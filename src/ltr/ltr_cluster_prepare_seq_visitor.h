/*
  Copyright (c) 2012 Sascha Steinbiss <steinbiss@zbh.uni-hamburg.de>
  Copyright (c) 2012 Center for Bioinformatics, University of Hamburg

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

#ifndef LTR_CLUSTER_PREPARE_SEQ_VISITOR_H
#define LTR_CLUSTER_PREPARE_SEQ_VISITOR_H

/* implements the ''node visitor'' interface */
typedef struct GtLTRClusterPrepareSeqVisitor GtLTRClusterPrepareSeqVisitor;

#include "core/array.h"
#include "core/str_array.h"
#include "extended/node_visitor.h"
#include "extended/region_mapping_api.h"

GtNodeVisitor* gt_ltr_cluster_prepare_seq_visitor_new(GtRegionMapping *rm,
                                                      GtError *err);
GtHashmap*     gt_ltr_cluster_prepare_seq_visitor_get_encseqs(
                                              GtLTRClusterPrepareSeqVisitor *v);
GtStrArray*    gt_ltr_cluster_prepare_seq_visitor_get_features(
                                              GtLTRClusterPrepareSeqVisitor *v);
GtHashmap*     gt_ltr_cluster_prepare_seq_visitor_get_id_to_feat_hash(
                                              GtLTRClusterPrepareSeqVisitor *v);
GtHashmap*     gt_ltr_cluster_prepare_seq_visitor_get_feat_to_id_hash(
                                              GtLTRClusterPrepareSeqVisitor *v);

#endif
