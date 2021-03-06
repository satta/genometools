/*
  Copyright (c) 2011-2012 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2011-2012 Center for Bioinformatics, University of Hamburg

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

/*
  THIS FILE IS GENERATED by
  scripts/gen-esa-bottomup.rb
  --key spmvar.
  DO NOT EDIT.
*/

#include <limits.h>
#include "core/ma_api.h"
#include "esa-seqread.h"
#include "seqnumrelpos.h"

static void initBUinfo_spmvar(GtBUinfo_spmvar *,
                              GtBUstate_spmvar *);

static void freeBUinfo_spmvar(GtBUinfo_spmvar *,
                              GtBUstate_spmvar *);

static int processleafedge_spmvar(bool,
    GtUword,
    GtBUinfo_spmvar *,
    GtUword, /* seqnum */
    GtUword, /* relpos */
    GtBUstate_spmvar *,
    GtError *err);

static int processbranchingedge_spmvar(bool firstsucc,
    GtUword,
    GtBUinfo_spmvar *,
    GtUword,
    GtUword,
    GtBUinfo_spmvar *,
    GtBUstate_spmvar *,
    GtError *);

static int processlcpinterval_spmvar(GtUword,
    GtBUinfo_spmvar *,
    GtBUstate_spmvar *,
    GtError *err);

#define TOP_ESA_BOTTOMUP_spmvar\
        stack->spaceGtBUItvinfo[stack->nextfreeGtBUItvinfo-1]

#define POP_ESA_BOTTOMUP_spmvar\
        (stack->spaceGtBUItvinfo + (--stack->nextfreeGtBUItvinfo))

#define PUSH_ESA_BOTTOMUP_spmvar(LCP,LB)\
        if (stack->nextfreeGtBUItvinfo >= stack->allocatedGtBUItvinfo)\
        {\
          gt_assert(stack->nextfreeGtBUItvinfo ==\
                    stack->allocatedGtBUItvinfo);\
          stack->spaceGtBUItvinfo\
            = allocateBUstack_spmvar(stack->spaceGtBUItvinfo,\
                              stack->allocatedGtBUItvinfo,\
                              stack->allocatedGtBUItvinfo+incrementstacksize,\
                              bustate);\
          stack->allocatedGtBUItvinfo += incrementstacksize;\
        }\
        gt_assert(stack->spaceGtBUItvinfo != NULL);\
        stack->spaceGtBUItvinfo[stack->nextfreeGtBUItvinfo].lcp = LCP;\
        stack->spaceGtBUItvinfo[stack->nextfreeGtBUItvinfo].lb = LB;\
        stack->spaceGtBUItvinfo[stack->nextfreeGtBUItvinfo++].rb = ULONG_MAX

typedef struct
{
  GtUword lcp, lb, rb;
  GtBUinfo_spmvar info;
} GtBUItvinfo_spmvar;

typedef struct
{
  GtBUItvinfo_spmvar *spaceGtBUItvinfo;
  GtUword allocatedGtBUItvinfo,
                nextfreeGtBUItvinfo;
} GtArrayGtBUItvinfo_spmvar;

GtArrayGtBUItvinfo_spmvar *gt_GtArrayGtBUItvinfo_new_spmvar(void)
{
  GtArrayGtBUItvinfo_spmvar *stack = gt_malloc(sizeof (*stack));

  stack->spaceGtBUItvinfo = NULL;
  stack->allocatedGtBUItvinfo = stack->nextfreeGtBUItvinfo = 0;
  return stack;
}

void gt_GtArrayGtBUItvinfo_delete_spmvar(
                           GtArrayGtBUItvinfo_spmvar *stack,
                           GtBUstate_spmvar *state)
{
  GtUword idx;

  for (idx=0; idx<stack->allocatedGtBUItvinfo; idx++)
  {
    freeBUinfo_spmvar(&stack->spaceGtBUItvinfo[idx].info,state);
  }
  gt_free(stack->spaceGtBUItvinfo);
  gt_free(stack);
}

static GtBUItvinfo_spmvar *allocateBUstack_spmvar(
                                   GtBUItvinfo_spmvar *ptr,
                                   GtUword currentallocated,
                                   GtUword allocated,
                                   GtBUstate_spmvar *state)
{
  GtUword idx;
  GtBUItvinfo_spmvar *itvinfo;

  itvinfo = gt_realloc(ptr,sizeof (*itvinfo) * allocated);
  gt_assert(allocated > currentallocated);
  for (idx=currentallocated; idx<allocated; idx++)
  {
    initBUinfo_spmvar(&itvinfo[idx].info,state);
  }
  gt_assert(itvinfo != NULL);
  return itvinfo;
}

static int gt_esa_bottomup_RAM_spmvar(const GtUword *bucketofsuffixes,
                        const uint16_t *lcptab_bucket,
                        GtUword numberofsuffixes,
                        GtArrayGtBUItvinfo_spmvar *stack,
                        GtBUstate_spmvar *bustate,
                        const GtSeqnumrelpos *snrp,
                        GtError *err)
{
  const GtUword incrementstacksize = 32UL;
  GtUword lcpvalue,
                previoussuffix,
                idx;
  GtBUItvinfo_spmvar *lastinterval = NULL;
  bool haserr = false, firstedge, firstedgefromroot;

  PUSH_ESA_BOTTOMUP_spmvar(0,0);
  firstedgefromroot = true;
  gt_assert (numberofsuffixes > 0);
  for (idx = 0; !haserr && idx < numberofsuffixes-1; idx++)
  {
    lcpvalue = (GtUword) lcptab_bucket[idx+1];
    previoussuffix = bucketofsuffixes[idx];
    gt_assert(stack->nextfreeGtBUItvinfo > 0);
    if (lcpvalue <= TOP_ESA_BOTTOMUP_spmvar.lcp)
    {
      if (TOP_ESA_BOTTOMUP_spmvar.lcp > 0 || !firstedgefromroot)
      {
        firstedge = false;
      } else
      {
        firstedge = true;
        firstedgefromroot = false;
      }
      if (processleafedge_spmvar(firstedge,
                          TOP_ESA_BOTTOMUP_spmvar.lcp,
                          &TOP_ESA_BOTTOMUP_spmvar.info,
                          gt_seqnumrelpos_decode_seqnum(snrp,previoussuffix),
                          gt_seqnumrelpos_decode_relpos(snrp,previoussuffix),
                          bustate,
                          err) != 0)
      {
        haserr = true;
      }
    }
    gt_assert(lastinterval == NULL);
    while (!haserr && lcpvalue < TOP_ESA_BOTTOMUP_spmvar.lcp)
    {
      lastinterval = POP_ESA_BOTTOMUP_spmvar;
      lastinterval->rb = idx;
      if (processlcpinterval_spmvar(lastinterval->lcp,
                             &lastinterval->info,
                             bustate,
                             err) != 0)
      {
        haserr = true;
      }
      if (lcpvalue <= TOP_ESA_BOTTOMUP_spmvar.lcp)
      {
        if (TOP_ESA_BOTTOMUP_spmvar.lcp > 0 || !firstedgefromroot)
        {
          firstedge = false;
        } else
        {
          firstedge = true;
          firstedgefromroot = false;
        }
        if (processbranchingedge_spmvar(firstedge,
               TOP_ESA_BOTTOMUP_spmvar.lcp,
               &TOP_ESA_BOTTOMUP_spmvar.info,
               lastinterval->lcp,
               lastinterval->rb - lastinterval->lb + 1,
               &lastinterval->info,
               bustate,
               err) != 0)
        {
          haserr = true;
        }
        lastinterval = NULL;
      }
    }
    if (!haserr && lcpvalue > TOP_ESA_BOTTOMUP_spmvar.lcp)
    {
      if (lastinterval != NULL)
      {
        GtUword lastintervallb = lastinterval->lb;
        GtUword lastintervallcp = lastinterval->lcp,
              lastintervalrb = lastinterval->rb;
        PUSH_ESA_BOTTOMUP_spmvar(lcpvalue,lastintervallb);
        if (processbranchingedge_spmvar(true,
                       TOP_ESA_BOTTOMUP_spmvar.lcp,
                       &TOP_ESA_BOTTOMUP_spmvar.info,
                       lastintervallcp,
                       lastintervalrb - lastintervallb + 1,
                       NULL,
                       bustate,
                       err) != 0)
        {
          haserr = true;
        }
        lastinterval = NULL;
      } else
      {
        PUSH_ESA_BOTTOMUP_spmvar(lcpvalue,idx);
        if (processleafedge_spmvar(true,
                            TOP_ESA_BOTTOMUP_spmvar.lcp,
                            &TOP_ESA_BOTTOMUP_spmvar.info,
                            gt_seqnumrelpos_decode_seqnum(snrp,previoussuffix),
                            gt_seqnumrelpos_decode_relpos(snrp,previoussuffix),
                            bustate,
                            err) != 0)
        {
          haserr = true;
        }
      }
    }
  }
  gt_assert(stack->nextfreeGtBUItvinfo > 0);
  if (!haserr && TOP_ESA_BOTTOMUP_spmvar.lcp > 0)
  {
    GtUword lastsuftabvalue = bucketofsuffixes[numberofsuffixes-1];
    if (processleafedge_spmvar(false,
                        TOP_ESA_BOTTOMUP_spmvar.lcp,
                        &TOP_ESA_BOTTOMUP_spmvar.info,
                        gt_seqnumrelpos_decode_seqnum(snrp,lastsuftabvalue),
                        gt_seqnumrelpos_decode_relpos(snrp,lastsuftabvalue),
                        bustate,
                        err) != 0)
    {
      haserr = true;
    } else
    {
      TOP_ESA_BOTTOMUP_spmvar.rb = idx;
      if (processlcpinterval_spmvar(TOP_ESA_BOTTOMUP_spmvar.lcp,
                             &TOP_ESA_BOTTOMUP_spmvar.info,
                             bustate,
                             err) != 0)
      {
        haserr = true;
      }
    }
  }
  stack->nextfreeGtBUItvinfo = 0; /* empty the stack */
  return haserr ? -1 : 0;
}
