/*
  Copyright (c) 2007 Stefan Kurtz <kurtz@zbh.uni-hamburg.de>
  Copyright (c) 2007 Center for Bioinformatics, University of Hamburg
  See LICENSE file or http://genometools.org/license.html for license details.
*/

#include <stdlib.h>
#include <string.h>
#include <stdbool.h>
#include "libgtcore/env.h"
#include "types.h"
#include "arraydef.h"
#include "spacedef.h"
#include "alphadef.h"
#include "chardef.h"
#include "addnextchar.h"
#include "encseq-def.h"
#include "intcode-def.h"
#include "sarr-def.h"

#include "kmer2string.pr"
#include "compfilenm.pr"
#include "mappedstr.pr"
#include "alphabet.pr"
#include "sfxmap.pr"

static Codetype qgram2codefillspecial(unsigned int numofchars,
                                      unsigned int kmersize,
                                      const Encodedsequence *encseq,
                                      Uint64 startpos,
                                      Uint64 totallength)
{
  Codetype integercode;
  Uint64 pos;
  bool foundspecial;
  Uchar cc;

  if (startpos >= totallength)
  {
    integercode = numofchars - 1;
    foundspecial = true;
  } else
  {
    cc = getencodedchar64(encseq,startpos);
    if (ISSPECIAL(cc))
    {
      integercode = numofchars - 1;
      foundspecial = true;
    } else
    {
      integercode = (Codetype) cc;
      foundspecial = false;
    }
  }
  for (pos = startpos + (Uint64) 1; pos< startpos + (Uint64) kmersize; pos++)
  {
    if (foundspecial)
    {
      ADDNEXTCHAR(integercode,numofchars-1,numofchars);
    } else
    {
      if (pos >= totallength)
      {
        ADDNEXTCHAR(integercode,numofchars-1,numofchars);
        foundspecial = true;
      } else
      {
        cc = getencodedchar64(encseq,pos);
        if (ISSPECIAL(cc))
        {
          ADDNEXTCHAR(integercode,numofchars-1,numofchars);
          foundspecial = true;
        } else
        {
          ADDNEXTCHAR(integercode,cc,numofchars);
        }
      }
    }
  }
  return integercode;
}

DECLAREARRAYSTRUCT(Codetype);

static void outkmeroccurrence(void *processinfo,
                              Codetype code,
                              /*@unused@*/ Uint64 position,
                              /*@unused@*/ const Firstspecialpos
                                                 *firstspecialposition,
                              Env *env)
{
  ArrayCodetype *codelist = (ArrayCodetype *) processinfo;

  STOREINARRAY(codelist,Codetype,1024,code);
}

static void collectkmercode(ArrayCodetype *codelist,
                            const Encodedsequence *encseq,
                            unsigned int kmersize,
                            unsigned int numofchars,
                            Uint64 stringtotallength,
                            Env *env)
{
  Uint64 offset;
  Codetype code;

  for (offset=0; offset<=stringtotallength; offset++)
  {
    code = qgram2codefillspecial(numofchars,
                                 kmersize,
                                 encseq,
                                 offset,
                                 stringtotallength);
    STOREINARRAY(codelist,Codetype,1024,code);
  }
}

static int comparecodelists(const ArrayCodetype *codeliststream,
                            const ArrayCodetype *codeliststring,
                            unsigned int kmersize,
                            unsigned int numofchars,
                            const char *characters,
                            Env *env)
{
  Uint i;
  char buffer1[64+1], buffer2[64+1];

  if (codeliststream->nextfreeCodetype != codeliststring->nextfreeCodetype)
  {
    env_error_set(env,
                  "length codeliststream= %lu != %lu =length codeliststring",
                  (Showuint) codeliststream->nextfreeCodetype,
                  (Showuint) codeliststring->nextfreeCodetype);
    return -1;
  }
  for (i=0; i<codeliststream->nextfreeCodetype; i++)
  {
    if (codeliststream->spaceCodetype[i] != codeliststring->spaceCodetype[i])
    {
      kmercode2string(buffer1,
                      codeliststream->spaceCodetype[i],
                      numofchars,
                      kmersize,
                      characters);
      kmercode2string(buffer2,
                      codeliststring->spaceCodetype[i],
                      numofchars,
                      kmersize,
                      characters);
      env_error_set(env,
                    "codeliststream[%lu] = %lu != %lu = "
                    "codeliststring[%lu]\n%s != %s",
                    (Showuint) i,
                    (Showuint) codeliststream->spaceCodetype[i],
                    (Showuint) codeliststring->spaceCodetype[i],
                    (Showuint) i,
                    buffer1,
                    buffer2);
      return -1;
    }
  }
  return 0;
}

static int verifycodelists(const Encodedsequence *encseq,
                           const Uchar *characters,
                           unsigned int kmersize,
                           unsigned int numofchars,
                           Uint64 stringtotallength,
                           const ArrayCodetype *codeliststream,
                           Env *env)
{
  bool haserr = false;
  ArrayCodetype codeliststring;

  INITARRAY(&codeliststring,Codetype);
  collectkmercode(&codeliststring,
                  encseq,
                  kmersize,
                  numofchars,
                  stringtotallength,
                  env);
  if (comparecodelists(codeliststream,
                       &codeliststring,
                       kmersize,
                       numofchars,
                       (const char *) characters,
                       env) != 0)
  {
    haserr = true;
  }
  FREEARRAY(&codeliststring,Codetype);
  return haserr ? -1 : 0;
}

int verifymappedstr(const Suffixarray *suffixarray,Env *env)
{
  unsigned int numofchars;
  ArrayCodetype codeliststream;
  bool haserr = false;

  numofchars = getnumofcharsAlphabet(suffixarray->alpha);
  INITARRAY(&codeliststream,Codetype);
  if (getfastastreamkmers((const char **) suffixarray->filenametab,
                          suffixarray->numoffiles,
                          outkmeroccurrence,
                          &codeliststream,
                          numofchars,
                          suffixarray->prefixlength,
                          getsymbolmapAlphabet(suffixarray->alpha),
                          env) != 0)
  {
    haserr = true;
  }
  if (!haserr)
  {
    if (verifycodelists(suffixarray->encseq,
                        getcharactersAlphabet(suffixarray->alpha),
                        suffixarray->prefixlength,
                        numofchars,
                        getencseqtotallength(suffixarray->encseq),
                        &codeliststream,
                        env) != 0)
    {
      haserr = true;
    }
  }
  FREEARRAY(&codeliststream,Codetype);
  return haserr ? -1 : 0;
}
