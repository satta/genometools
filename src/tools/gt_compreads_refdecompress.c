/*
  Copyright (c) 2011 Joachim Bonnet <joachim.bonnet@studium.uni-hamburg.de>
  Copyright (c) 2012 Dirk Willrodt <willrodt@zbh.uni-hamburg.de>
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

#include <string.h>
#include "core/alphabet_api.h"
#include "core/basename_api.h"
#include "core/fa.h"
#include "core/fileutils_api.h"
#include "core/ma.h"
#include "core/mathsupport.h"
#include "core/range_api.h"
#include "core/showtime.h"
#include "core/undef_api.h"
#include "core/unused_api.h"
#include "core/xansi_api.h"
#include "extended/rcr.h"
#include "tools/gt_compreads_refdecompress.h"

typedef struct {
  GtStr  *file,
         *ref,
         *name;
  bool verbose,
       qnames;
} GtCsrRcrDecodeArguments;

static void* gt_compreads_refdecompress_arguments_new(void)
{
  GtCsrRcrDecodeArguments *arguments = gt_calloc((size_t) 1, sizeof *arguments);
  arguments->file = gt_str_new();
  arguments->ref = gt_str_new();
  arguments->name = gt_str_new();

  return arguments;
}

static void gt_compreads_refdecompress_arguments_delete(void *tool_arguments)
{
  GtCsrRcrDecodeArguments *arguments = tool_arguments;
  if (!arguments) return;
  gt_str_delete(arguments->file);
  gt_str_delete(arguments->name);
  gt_str_delete(arguments->ref);
  gt_free(arguments);
}

static GtOptionParser* gt_compreads_refdecompress_option_parser_new(void *tool_arguments)
{
  GtCsrRcrDecodeArguments *arguments = tool_arguments;
  GtOptionParser *op;
  GtOption *option;
  gt_assert(arguments);

  /* init */
  op = gt_option_parser_new("[option ...] (-file file)",
                         "Decodes a given RCR (Reference Compressed Reads).");

  option = gt_option_new_bool("v", "be verbose",
                              &arguments->verbose, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_bool("qnames", "decode read names, default is "
                              "to just number them",
                              &arguments->qnames, false);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_string("ref", "Index file (generated by the gt encseq"
                                " tool) for reference genome.",
                                arguments->ref, NULL);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_string("rcr", "specify base name of file containing"
                                " RCR.",
                                arguments->file, NULL);
  gt_option_is_mandatory(option);
  gt_option_parser_add_option(op, option);

  option = gt_option_new_string("name", "specify base name for decoded RCR"
                                " (suffix will be \".rcr.decoded\")",
                                arguments->name, NULL);
  gt_option_parser_add_option(op, option);

  gt_option_parser_set_min_max_args(op, 0U, 0U);
  return op;
}

static int gt_compreads_refdecompress_runner(GT_UNUSED int argc,
                                    GT_UNUSED const char **argv,
                                    GT_UNUSED int parsed_args,
                                    void *tool_arguments, GtError *err)
{
  GtCsrRcrDecodeArguments *arguments = tool_arguments;
  int had_err = 0;
  GtRcrDecoder *rcrd = NULL;
  GtTimer *timer = NULL;
  GtEncseq *encseq = NULL;
  GtEncseqLoader *el = NULL;

  gt_error_check(err);
  gt_assert(arguments);

  if (gt_showtime_enabled()) {
    timer = gt_timer_new_with_progress_description("start");
    gt_timer_start(timer);
    gt_assert(timer);
  }

  if (gt_str_length(arguments->name) == 0) {
    char *basename = gt_basename(gt_str_get(arguments->file));
    gt_str_set(arguments->name, basename);
    gt_free(basename);
  }
  if (timer != NULL)
    gt_timer_show_progress(timer, "load encseq", stdout);
  el = gt_encseq_loader_new();
  gt_encseq_loader_enable_autosupport(el);
  gt_encseq_loader_require_description_support(el);
  encseq = gt_encseq_loader_load(el, gt_str_get(arguments->ref), err);
  gt_encseq_loader_delete(el);
  if (!encseq) {
    gt_error_set(err, "Could not load GtEncseq %s", gt_str_get(arguments->ref));
    had_err = -1;
  }

  if (!had_err) {
    if (!gt_alphabet_is_dna(gt_encseq_alphabet(encseq))) {
      gt_error_set(err, "Alphabet in %s has to be DNA",
                   gt_str_get(arguments->ref));
      had_err = -1;
    }
  }

  if (!had_err) {
    rcrd = gt_rcr_decoder_new(gt_str_get(arguments->file), encseq, timer,
                              err);
    if (rcrd == NULL)
      had_err = -1;
  }
  if (!had_err &&
      arguments->qnames)
    had_err = gt_rcr_decoder_enable_description_support(rcrd, err);
  if (!had_err)
    had_err = gt_rcr_decoder_decode(rcrd,
                                    gt_str_get(arguments->name), timer, err);
  if (timer != NULL) {
    gt_timer_show_progress_final(timer, stdout);
    gt_timer_delete(timer);
  }
  gt_rcr_decoder_delete(rcrd);
  gt_encseq_delete(encseq);
  return had_err;
}

GtTool* gt_compreads_refdecompress(void)
{
  return gt_tool_new(gt_compreads_refdecompress_arguments_new,
                     gt_compreads_refdecompress_arguments_delete,
                     gt_compreads_refdecompress_option_parser_new,
                     NULL,
                     gt_compreads_refdecompress_runner);
}
