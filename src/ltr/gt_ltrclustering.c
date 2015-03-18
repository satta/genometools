/*
  Copyright (c) 2011 Sascha Kastens <mail@skastens.de>
  Copyright (c) 2014 Sascha Steinbiss <sascha@steinbiss.name>
  Copyright (c) 2011 Center for Bioinformatics, University of Hamburg

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

#include "core/encseq.h"
#include "core/ma.h"
#include "core/output_file_api.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "extended/gff3_in_stream.h"
#include "extended/gff3_out_stream_api.h"
#include "extended/seqid2file.h"
#include "ltr/ltr_cluster_stream.h"
#include "ltr/ltr_classify_stream.h"
#include "ltr/gt_ltrclustering.h"

typedef struct {
  GtFile *outfp;
  GtOutputFileInfo *ofi;
  GtUword psmall, plarge;
  double xdrop, identity;
  int wordsize;
  GtStrArray *feats;
  GtSeqid2FileInfo *s2fi;
} GtLTRClusteringArguments;

static void* gt_ltrclustering_arguments_new(void)
{
  GtLTRClusteringArguments *arguments = gt_calloc((size_t) 1,
                                                  sizeof (*arguments));
  arguments->ofi = gt_output_file_info_new();
  arguments->s2fi = gt_seqid2file_info_new();
  arguments->feats = gt_str_array_new();
  return arguments;
}

static void gt_ltrclustering_arguments_delete(void *tool_arguments)
{
  GtLTRClusteringArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_file_delete(arguments->outfp);
  gt_str_array_delete(arguments->feats);
  gt_seqid2file_info_delete(arguments->s2fi);
  gt_output_file_info_delete(arguments->ofi);
  gt_free(arguments);
}

static GtOptionParser* gt_ltrclustering_option_parser_new(void *tool_arguments)
{
  GtLTRClusteringArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option...] [GFF3_file ...]",
                            "Cluster features of LTRs.");

  /* -psmall */
  option = gt_option_new_uword_min_max("psmall", "percentage of"
                                       " the smaller sequence a match needs to"
                                       " cover in order to cluster the two"
                                       " sequences of the match",
                                       &arguments->psmall, 80UL, 0, 100UL);

  gt_option_parser_add_option(op, option);

  /* -plarge */
  option = gt_option_new_uword_min_max("plarge", "percentage of"
                                       " the larger sequence a match needs to"
                                       " cover in order to cluster the two"
                                       " sequences of the match",
                                       &arguments->plarge, 20UL, 0, 100UL);
  gt_option_parser_add_option(op, option);

  /* -usefeatures */
  option = gt_option_new_string_array("usefeatures", "list of features used "
                                      "for classification (e.g. lLTR, rLTR, "
                                      "etc.), terminate with '--'",
                                      arguments->feats);
  gt_option_parser_add_option(op, option);

  gt_seqid2file_register_options(op, arguments->s2fi);
  gt_output_file_info_register_options(arguments->ofi, op, &arguments->outfp);

  return op;
}

static int gt_ltrclustering_runner(int argc, const char **argv,
                                       int parsed_args, void *tool_arguments,
                                       GtError *err)
{
  GtLTRClusteringArguments *arguments = tool_arguments;
  GtNodeStream *gff3_in_stream      = NULL,
               *last_stream         = NULL,
               *ltr_cluster_stream  = NULL,
               *ltr_classify_stream = NULL,
               *gff3_out_stream     = NULL;
  GtRegionMapping *rm = NULL;
  int had_err = 0, arg = parsed_args;
  gt_error_check(err);
  gt_assert(arguments);
  arg++;

  /* add region mapping if given */
  if (gt_seqid2file_option_used(arguments->s2fi)) {
    rm = gt_seqid2file_region_mapping_new(arguments->s2fi, err);
    if (!rm)
      had_err = -1;
  }

  if (!had_err) {
    last_stream = gff3_in_stream = gt_gff3_in_stream_new_unsorted(argc - arg,
                                                                  argv + arg);
    gt_assert(gff3_in_stream);
    last_stream = ltr_cluster_stream = gt_ltr_cluster_stream_new(last_stream,
                                                             rm,
                                                             GT_UNDEF_INT,
                                                             GT_UNDEF_INT,
                                                             GT_UNDEF_INT,
                                                             GT_UNDEF_INT,
                                                             GT_UNDEF_INT,
                                                             GT_UNDEF_INT,
                                                             GT_UNDEF_INT,
                                                             10,
                                                             GT_UNDEF_INT,
                                                             GT_UNDEF_INT,
                                                             arguments->plarge,
                                                             arguments->psmall,
                                                             NULL,
                                                             err);
    if (!ltr_cluster_stream)
      had_err = -1;
  }
  if (!had_err) {
    last_stream = ltr_classify_stream = gt_ltr_classify_stream_new(last_stream,
                                                                   NULL,
                                                                   NULL,
                                                                   NULL,
                                                                   NULL,
                                                                   err);
    if (!ltr_classify_stream)
      had_err = -1;
  }
  if (!had_err) {
      last_stream = gff3_out_stream = gt_gff3_out_stream_new(last_stream,
                                                             arguments->outfp);
  }
  if (!had_err)
    had_err = gt_node_stream_pull(last_stream, err);

  gt_node_stream_delete(ltr_classify_stream);
  gt_node_stream_delete(ltr_cluster_stream);
  gt_node_stream_delete(gff3_in_stream);
  gt_node_stream_delete(gff3_out_stream);
  gt_region_mapping_delete(rm);

  return had_err;
}

GtTool* gt_ltrclustering(void)
{
  return gt_tool_new(gt_ltrclustering_arguments_new,
                  gt_ltrclustering_arguments_delete,
                  gt_ltrclustering_option_parser_new,
                  NULL,
                  gt_ltrclustering_runner);
}
