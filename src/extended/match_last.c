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

#include <stdlib.h>
#include <string.h>
#include "core/assert_api.h"
#include "core/ma.h"
#include "core/undef_api.h"
#include "extended/match_last.h"
#include "extended/match_rep.h"

struct GtMatchLAST {
  GtMatch parent_instance;
  unsigned long seqno1,
                seqno2,
                score;
};

const GtMatchClass* gt_match_last_class()
{
  static const GtMatchClass *matchc = NULL;
  if (!matchc)
    matchc = gt_match_class_new(sizeof (GtMatchLAST),
                                NULL,
                                NULL);
  return matchc;
}

GtMatch* gt_match_last_new(const char *seqid1,
                           const char *seqid2,
                           unsigned long score,
                           unsigned long seqno1,
                           unsigned long seqno2,
                           unsigned long start_seq1,
                           unsigned long start_seq2,
                           unsigned long end_seq1,
                           unsigned long end_seq2)
{
  GtMatch *match;
  GtMatchLAST *matchs;
  match = gt_match_create(gt_match_last_class(), start_seq1, end_seq1,
                          start_seq2, end_seq2, seqid1, seqid2);
  matchs = gt_match_last_cast(match);
  matchs->seqno1 = seqno1;
  matchs->seqno2 = seqno2;
  matchs->score = score;
  return match;
}

unsigned long gt_match_last_get_seqno1(const GtMatchLAST *ms)
{
  gt_assert(ms);
  return ms->seqno1;
}

unsigned long gt_match_last_get_seqno2(const GtMatchLAST *ms)
{
  gt_assert(ms);
  return ms->seqno2;
}

unsigned long gt_match_last_get_score(const GtMatchLAST *ms)
{
  gt_assert(ms);
  return ms->score;
}
