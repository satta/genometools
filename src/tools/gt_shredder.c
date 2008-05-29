/*
  Copyright (c) 2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
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

#include "libgtcore/bioseq.h"
#include "libgtcore/fasta.h"
#include "libgtcore/option.h"
#include "libgtcore/ma.h"
#include "libgtcore/unused.h"
#include "libgtext/shredder.h"
#include "tools/gt_shredder.h"

typedef struct {
  unsigned long coverage,
                minlength,
                maxlength,
                overlap;
} ShredderArguments;

static void* gt_shredder_arguments_new(void)
{
  return ma_calloc(1, sizeof (ShredderArguments));
}

static void gt_shredder_arguments_delete(void *tool_arguments)
{
  ShredderArguments *arguments = tool_arguments;
  if (!arguments) return;
  ma_free(arguments);
}

static OptionParser* gt_shredder_option_parser_new(void *tool_arguments)
{
  ShredderArguments *arguments = tool_arguments;
  OptionParser *op;
  Option *o;
  assert(arguments);
  op = option_parser_new("[option ...] sequence_file",
                         "Shredder sequence_file into consecutive pieces of "
                         "random length.");
  o = option_new_ulong_min("coverage", "Set the number of times the "
                           "sequence_file is shreddered", &arguments->coverage,
                           1, 1);
  option_parser_add_option(op, o);
  o = option_new_ulong("minlength", "Set the minimum length of the shreddered "
                       "fragments", &arguments->minlength, 300);
  option_parser_add_option(op, o);
  o = option_new_ulong("maxlength", "Set the maximum length of the shreddered "
                       "fragments", &arguments->maxlength, 700);
  option_parser_add_option(op, o);
  o = option_new_ulong("overlap", "Set the overlap between of consecutive "
                       "pieces", &arguments->overlap, 0);
  option_parser_add_option(op, o);
  option_parser_set_min_max_args(op, 1, 1);
  return op;
}

static int gt_shredder_runner(UNUSED int argc, const char **argv,
                              int parsed_args, void *tool_arguments, Error *err)
{
  ShredderArguments *arguments = tool_arguments;
  unsigned long i;
  Bioseq *bioseq;
  int had_err = 0;
  Str *desc;

  error_check(err);
  assert(arguments);

  /* init */
  desc = str_new();
  bioseq = bioseq_new(argv[parsed_args], err);
  if (!bioseq)
    had_err = -1;

  /* shredder */
  for (i = 0; !had_err && i < arguments->coverage; i++) {
    Shredder *shredder;
    unsigned long fragment_length;
    const char *fragment;
    shredder = shredder_new(bioseq, arguments->minlength, arguments->maxlength);
    shredder_set_overlap(shredder, arguments->overlap);
    while ((fragment = shredder_shred(shredder, &fragment_length, desc))) {
      str_append_cstr(desc, " [shreddered fragment]");
      fasta_show_entry(str_get(desc), fragment, fragment_length, 0);
      str_reset(desc);
    }
    shredder_delete(shredder);
  }

  /* free */
  bioseq_delete(bioseq);
  str_delete(desc);

  return had_err;
}

Tool* gt_shredder(void)
{
  return tool_new(gt_shredder_arguments_new,
                  gt_shredder_arguments_delete,
                  gt_shredder_option_parser_new,
                  NULL,
                  gt_shredder_runner);
}
