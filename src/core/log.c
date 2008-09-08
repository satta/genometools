/*
  Copyright (c) 2006-2008 Gordon Gremme <gremme@zbh.uni-hamburg.de>
  Copyright (c) 2006-2008 Center for Bioinformatics, University of Hamburg

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

#include <stdarg.h>
#include <stdio.h>
#include "core/log.h"
#include "core/xansi.h"

static bool logging = false;
static FILE *logfp = NULL;

void gt_log_enable(void)
{
  logfp = stderr;
  logging = true;
}

bool gt_log_enabled(void)
{
  return logging;
}

void gt_log_log(const char *format, ...)
{
  va_list ap;
  if (!logging) return;
  va_start(ap, format);
  gt_log_vlog(format, ap);
  va_end(ap);
}

void gt_log_vlog(const char *format, va_list ap)
{
  if (!logging) return;
  assert(logfp);
  fprintf(logfp, "debug: ");
  (void) vfprintf(logfp, format, ap);
  (void) putc('\n', logfp);
}

FILE* gt_log_fp(void)
{
  assert(logging);
  return logfp;
}

void gt_log_set_fp(FILE *fp)
{
  assert(logging);
  logfp = fp;
}
