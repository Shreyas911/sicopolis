/*
 * TAPENADE Automatic Differentiation Engine
 * Copyright (C) 1999-2024 Inria
 * See the LICENSE.md file in the project root for more information.
 *
 * 
 * This package gathers functions used to profile a piece of (adjoint differentiated) code in terms of memory and time. 
 */

// Compile with with -D_ADPROFILETRACE=<num> to trace profile computations about the checkpoint
//  whose ckpRank (from 1 up) is <num>, or about all checkpoints if passing -D_ADPROFILETRACE=-1

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <signal.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include <unistd.h>
#include <stdint.h>
#include <inttypes.h>

#include <assert.h>

#include "adProfile.h"
#include "adStack.h"
/*
 *Profiling will work as follows
 * SNPWRITE: 
 * -> create new record with XXX Some information
 * -> Caller record id is set to a static variable
 * -> Save current stack size in current record
 * BEGIN ADV: 
 * -> Save current stack size in current record
 * -> Compute snapshot size (diff from both current stack sizes)
 * -> Update current caller id to this record
 * END ADV: 
 * -> Add current call to the "wait_for_bwd" stack 
 * 
 */

typedef uint64_t stack_size_t;
typedef int64_t  delta_size_t;

typedef clock_t cycle_time_t;
cycle_time_t mytime_() {
        return clock();
}

/** The maximum number of code checkpoint locations that we can handle */
#define MAX_NB_CHECKPOINTS 2000

/** Elementary cost/benefit data of the decision to NOT checkpoint a given static checkpoint location. */
typedef struct _CostBenefit{
  unsigned int ckpRank; // Rank (>=1) of the static checkpoint location this cost/benefit applies to. 
  cycle_time_t deltaT; // The (always >=0) spared runtime benefit of NOT checkpointing this location.
  delta_size_t deltaPM; // The extra peak memory cost of NOT checkpointing this location.
  delta_size_t deltaPMTurn; // Same as deltaPM, but restricted to the peak memory at the next turn.
  struct _CostBenefit *next; // Following CostBenefit infos, about other checkpoint locations.
} CostBenefit;

/** One element in a list of runtime checkpoints, each with some time/peak memory cost data.
 * The checkpoint locations in the list are all directly contained in the same parent runtime checkpoint.
 * The list is ordered from the most recent (downstream) checkpoint till the most upstream one. */
typedef struct _LevelCkpStack {
  unsigned int ckpRank; // The rank (>=1) of the static checkpoint that is now encountered
  cycle_time_t ckpT;
  cycle_time_t basisCkpT;
#ifdef _ADPROFILETRACE
  stack_size_t sizeBeforeSnp;
#endif
  stack_size_t sizeAfterSnp;
  struct _LevelCkpStack *previous;
} LevelCkpStack;

/** One element in the list of nested runtime checkpoints at some present time during the profiling run.
 * The list is ordered from the deepest runtime checkpoint to the topmost enclosing runtime checkpoint. */
//TODO meilleur nom que CkpStack (stack deja pris!)
typedef struct _CkpStack {
  LevelCkpStack *levelCheckpoints;
  CostBenefit *costBenefits;
  stack_size_t stackPeak; // Peak stack size reached during reversal of this level.
  stack_size_t stackPeakTurn; // Peak stack size reached at the next turn.
  LevelCkpStack *repeatLevel; // for management of start/reset/end Repeat
  struct _CkpStack *above;
  struct _CkpStack *below;
} CkpStack;

/** The main global stack of active checkpoints at some present time during the profiling run. */
CkpStack initialCkpStack = {.levelCheckpoints=NULL, .costBenefits=NULL,
                            .stackPeak=0, .stackPeakTurn=0, .repeatLevel=NULL,
                            .above=NULL, .below=NULL};
CkpStack *curCheckpoints = &initialCkpStack;

/************************* UTILITIES ************************/ 

/** Describes one (static) checkpoint location in the code source, with a unique ckpRank attached to it.
 * A (yes/no) checkpointing decision applies to one static checkpoint location globally. */
typedef struct {
  char* fname; // Name of the file in which the checkpoint appears
  unsigned int line_nb; // Line number in the original code where the checkpoint appears
  char* callee_name; // Name of the procedure being called (when checkpoint is on a call)
  unsigned int occurrences;
  unsigned int ckpRank; // Unique integer rank (from 1 up) attached to this checkpoint location.
} CkpLocation;

/** The array of all static checkpoint locations encountered so far */
CkpLocation* allCkpLocations[MAX_NB_CHECKPOINTS] = {};

/** Get or create the unique ckpRank (>=1) associated to the given static checkpoint location,
 * which is the combination of the file name and the line number in this file at which this
 * checkpoint is located. When this checkpoint is at a call site, callname should be passed
 * the procedure name, otherwise NULL. */
unsigned int getCheckpointLocationRank(char* callname, char* filename, unsigned int lineno) {
  unsigned int i = 1;
  int found = 0;
  while (i<MAX_NB_CHECKPOINTS && !found) {
    if (allCkpLocations[i]==NULL) {
      allCkpLocations[i] = (CkpLocation *)malloc(sizeof(CkpLocation));
      allCkpLocations[i]->ckpRank = i;
      allCkpLocations[i]->fname = (char*)malloc((1+strlen(filename))*sizeof(char));
      strcpy(allCkpLocations[i]->fname, filename);
      allCkpLocations[i]->line_nb = lineno;
      if (callname) {
        allCkpLocations[i]->callee_name = (char*)malloc((1+strlen(callname))*sizeof(char));
        strcpy(allCkpLocations[i]->callee_name, callname);
      } else {
        allCkpLocations[i]->callee_name = NULL;
      }
      allCkpLocations[i]->occurrences = 0;
      found = 1;
    } else if (allCkpLocations[i]->line_nb==lineno
               && 0==strcmp(allCkpLocations[i]->fname, filename)
               && 0==strcmp(allCkpLocations[i]->callee_name, callname)) {
      found = 1;
    }
    if (!found) ++i;
  }
  return i;
}

LevelCkpStack* stackNewLevelCkp(unsigned int ckpRank, LevelCkpStack* levelCheckpoints) {
  LevelCkpStack* additionalLevelCkp = (LevelCkpStack*)malloc(sizeof(LevelCkpStack));
  additionalLevelCkp->ckpRank = ckpRank;
  additionalLevelCkp->ckpT = 0;
  additionalLevelCkp->basisCkpT = 0;
#ifdef _ADPROFILETRACE
  additionalLevelCkp->sizeBeforeSnp = 0;
#endif
  additionalLevelCkp->sizeAfterSnp = 0;
  additionalLevelCkp->previous = levelCheckpoints;
  return additionalLevelCkp;
}

void releaseLevelCkp(CkpStack* checkpoints) {
  LevelCkpStack* levelCkps = checkpoints->levelCheckpoints;
  checkpoints->levelCheckpoints = levelCkps->previous;
  if (checkpoints->repeatLevel==NULL) { // Don't free if we are in Repeat mode!
    free(levelCkps);
  }
}

CkpStack* openNewCkp(CkpStack* checkpoints) {
  CkpStack* newCkp = checkpoints->above;
  if (!newCkp) {
    newCkp = (CkpStack*)malloc(sizeof(CkpStack));
    checkpoints->above = newCkp;
    newCkp->repeatLevel = NULL;
    newCkp->below = checkpoints;
    newCkp->above = NULL;
  }
  newCkp->levelCheckpoints = NULL;
  newCkp->costBenefits = NULL;
  newCkp->stackPeak = 0;
  newCkp->stackPeakTurn = 0;
  return newCkp;
}

CkpStack* closeCkp(CkpStack* checkpoints) {
  return checkpoints->below;
}

/** Advances in the given (address of a) chain of CostBenefit's, till finding a
 * CostBenefit about the checkpoint indexed as "ckpRank". This chain is ordered by growing ckpRank's.
 * If no such CostBenefit is found, creates and inserts a new, empty one, at the correct location.
 * This insertion modifies the given chain of CostBenefit's by side-effect.
 * Returns the address of the sub-chain of CostBenefit that starts with the one found or created. */
CostBenefit** getSetCostBenefit(CostBenefit **toCostBenefits, unsigned int ckpRank) {
  while (*toCostBenefits!=NULL && (*toCostBenefits)->ckpRank < ckpRank) {
    toCostBenefits = &((*toCostBenefits)->next);
  }
  if (*toCostBenefits==NULL || (*toCostBenefits)->ckpRank > ckpRank) {
    CostBenefit* additionalCostBenefits = (CostBenefit*)malloc(sizeof(CostBenefit));
    additionalCostBenefits->ckpRank = ckpRank;
    additionalCostBenefits->deltaT = 0;
    additionalCostBenefits->deltaPM = 0;
    additionalCostBenefits->deltaPMTurn = 0;
    additionalCostBenefits->next = *toCostBenefits;
    *toCostBenefits = additionalCostBenefits;
  }
  return toCostBenefits;
}

#ifdef _ADPROFILETRACE

int depthCkp(CkpStack* checkpoints) {
  int i = 0;
  while (checkpoints) {
    ++i;
    checkpoints = checkpoints->below;
  }
  return i;
}

int costBenefitsLength(CostBenefit* costBenefits) {
  int i = 0;
  while (costBenefits) {
    ++i;
    costBenefits = costBenefits->next;
  }
  return i;
}

void dumpCostBenefits(CostBenefit *costBenefits) {
  while (costBenefits) {
    printf(" %u:(T:%"PRIu64";M:%"PRId64"; turn:%"PRId64")",
           costBenefits->ckpRank, costBenefits->deltaT, costBenefits->deltaPM, costBenefits->deltaPMTurn);
    costBenefits = costBenefits->next;
  }
}

#endif

/****************** MAIN PUBLISHED PRIMITIVES ****************/

void adProfileAdj_SNPWrite(char* callname, char* filename, unsigned int lineno) {
  unsigned int ckpRank = getCheckpointLocationRank(callname, filename, lineno);
  curCheckpoints->levelCheckpoints = stackNewLevelCkp(ckpRank, curCheckpoints->levelCheckpoints);
#ifdef _ADPROFILETRACE
  curCheckpoints->levelCheckpoints->sizeBeforeSnp = adStack_getCurrentStackSize();
#endif
  curCheckpoints->levelCheckpoints->basisCkpT = mytime_();
}

void adProfileAdj_beginAdvance(char* callname, char* filename, unsigned int lineno) {
  assert(getCheckpointLocationRank(callname, filename, lineno)==curCheckpoints->levelCheckpoints->ckpRank);
  curCheckpoints->levelCheckpoints->sizeAfterSnp = adStack_getCurrentStackSize();
}

void adProfileAdj_endAdvance(char* callname, char* filename, unsigned int lineno) {
  curCheckpoints->levelCheckpoints->ckpT = mytime_() - curCheckpoints->levelCheckpoints->basisCkpT;
  assert(getCheckpointLocationRank(callname, filename, lineno)==curCheckpoints->levelCheckpoints->ckpRank);
}

void adProfileAdj_SNPRead(char* callname, char* filename, unsigned int lineno) {
  assert(getCheckpointLocationRank(callname, filename, lineno)==curCheckpoints->levelCheckpoints->ckpRank);
  curCheckpoints->levelCheckpoints->basisCkpT = mytime_();
}

void adProfileAdj_beginReverse(char* callname, char* filename, unsigned int lineno) {
  curCheckpoints->levelCheckpoints->ckpT += mytime_() - curCheckpoints->levelCheckpoints->basisCkpT;
  assert(getCheckpointLocationRank(callname, filename, lineno)==curCheckpoints->levelCheckpoints->ckpRank);
  ++(allCkpLocations[curCheckpoints->levelCheckpoints->ckpRank]->occurrences);
  curCheckpoints = openNewCkp(curCheckpoints);
}

void adProfileAdj_endReverse(char* callname, char* filename, unsigned int lineno) {
  assert(curCheckpoints->levelCheckpoints==NULL);
  CostBenefit *additionalCostBenefits = curCheckpoints->costBenefits;
  stack_size_t peakSize2 = curCheckpoints->stackPeak;
  stack_size_t peakSizeTurn2 = curCheckpoints->stackPeakTurn;
  // Special case if "callname" is an external, "callname_B" did not call adProfileAdj_turn,
  //  and curCheckpoints->stackPeak is still 0 ! Set it to the current stack size:
  if (peakSize2==0) peakSize2 = adStack_getCurrentStackSize();
  curCheckpoints = closeCkp(curCheckpoints);
  // Now merge the additional CostBenefit's from the checkpointed fragment into the current CostBenefit's:
  unsigned int localCkpRank = curCheckpoints->levelCheckpoints->ckpRank;
  cycle_time_t localCkpTimeCost = curCheckpoints->levelCheckpoints->ckpT;
  stack_size_t peakSize1 = curCheckpoints->stackPeak;
  stack_size_t peakSizeTurn1 = curCheckpoints->stackPeakTurn;
  int mergedLocal = 0;
  unsigned int additionalCkpRank;
  cycle_time_t additionalT;
  delta_size_t additionalPM;
  delta_size_t additionalPMTurn;
#ifdef _ADPROFILETRACE
  if (_ADPROFILETRACE==-1 || _ADPROFILETRACE==localCkpRank) {
    printf("ENDREVERSE OF %s::%i CALL %s, DEPTH:%i, SIZEBEFORESNP:%"PRId64", SIZEAFTERSNP:%"PRId64"\n",
           allCkpLocations[localCkpRank]->fname,
           allCkpLocations[localCkpRank]->line_nb,
           (allCkpLocations[localCkpRank]->callee_name ? allCkpLocations[localCkpRank]->callee_name : "MANUAL"),
           depthCkp(curCheckpoints),
           curCheckpoints->levelCheckpoints->sizeBeforeSnp,
           curCheckpoints->levelCheckpoints->sizeAfterSnp);
    printf("  BEFORE (peak:%"PRId64" bytes, peak at turn:%"PRId64" bytes) :", peakSize1, peakSizeTurn1);
    dumpCostBenefits(curCheckpoints->costBenefits);
    printf("\n  ADDING (peak:%"PRId64" bytes, peak at turn:%"PRId64" bytes) :", peakSize2, peakSizeTurn2);
    dumpCostBenefits(additionalCostBenefits);
    printf("\n");
  }
  int nbProfilesBefore = costBenefitsLength(curCheckpoints->costBenefits);
  int nbProfilesAdded = costBenefitsLength(additionalCostBenefits);
#endif
  CostBenefit **toCostBenefits = &(curCheckpoints->costBenefits);
  // compute merged peak size:
  if (peakSize2>peakSize1) curCheckpoints->stackPeak = peakSize2;
  // compute merged list of costs/benefits:
  while (additionalCostBenefits || !mergedLocal) {
#ifdef _ADPROFILETRACE
  if (_ADPROFILETRACE==-1 || _ADPROFILETRACE==localCkpRank) {
    printf("    mergedLocal:%i next additional:%u,%"PRIu64",%"PRId64",%"PRId64"\n",
           mergedLocal,(additionalCostBenefits?additionalCostBenefits->ckpRank:0),
           (additionalCostBenefits?additionalCostBenefits->deltaT:0),
           (additionalCostBenefits?additionalCostBenefits->deltaPM:0),
           (additionalCostBenefits?additionalCostBenefits->deltaPMTurn:0));
  }
#endif
    if (!mergedLocal
        && (additionalCostBenefits==NULL || additionalCostBenefits->ckpRank>localCkpRank)) {
      additionalCkpRank = localCkpRank;
      additionalT = localCkpTimeCost;
      additionalPM = 0;
      additionalPMTurn = 0;
      mergedLocal = 1;
    } else if (additionalCostBenefits->ckpRank==localCkpRank) {
      additionalCkpRank = localCkpRank;
      additionalT = additionalCostBenefits->deltaT + localCkpTimeCost;
      additionalPM = additionalCostBenefits->deltaPM ;
      additionalPMTurn = additionalCostBenefits->deltaPMTurn;
      CostBenefit* tmp = additionalCostBenefits->next;
      free(additionalCostBenefits);
      additionalCostBenefits = tmp;
      mergedLocal = 1;
    } else {
      additionalCkpRank = additionalCostBenefits->ckpRank;
      additionalT = additionalCostBenefits->deltaT;
      additionalPM = additionalCostBenefits->deltaPM;
      additionalPMTurn = additionalCostBenefits->deltaPMTurn;
      CostBenefit* tmp = additionalCostBenefits->next;
      free(additionalCostBenefits);
      additionalCostBenefits = tmp;
    }
    toCostBenefits = getSetCostBenefit(toCostBenefits, additionalCkpRank);
#ifdef _ADPROFILETRACE
  if (_ADPROFILETRACE==-1 || _ADPROFILETRACE==localCkpRank) {
    printf("    additionalCkpRank:%u deltaPM:%"PRId64" additionalPM:%"PRId64"\n", additionalCkpRank, (*toCostBenefits)->deltaPM, additionalPM);
  }
#endif
    if (additionalCkpRank==localCkpRank) {
      delta_size_t intermediate = peakSizeTurn2 + additionalPMTurn - curCheckpoints->levelCheckpoints->sizeAfterSnp;
#ifdef _ADPROFILETRACE
  if (_ADPROFILETRACE==-1 || _ADPROFILETRACE==localCkpRank) {
    printf("      %"PRIu64"  %"PRId64" %"PRIu64"\n", peakSizeTurn2, additionalPMTurn, curCheckpoints->levelCheckpoints->sizeAfterSnp) ;
    printf("      %"PRId64" += %"PRId64"\n", (*toCostBenefits)->deltaPMTurn, intermediate);
  }
#endif
      (*toCostBenefits)->deltaPMTurn += intermediate;
      stack_size_t costPeak1 = peakSize1 + (*toCostBenefits)->deltaPM + intermediate;
      stack_size_t costPeak2 = peakSize2 + additionalPM;
      (*toCostBenefits)->deltaPM = (costPeak1>costPeak2 ? costPeak1 : costPeak2) - curCheckpoints->stackPeak;
    } else {
      // Do a max on memory costs:
      stack_size_t costPeak1 = peakSize1 + (*toCostBenefits)->deltaPM;
      stack_size_t costPeak2 = peakSize2 + additionalPM;
      (*toCostBenefits)->deltaPM = (costPeak1>costPeak2 ? costPeak1 : costPeak2) - curCheckpoints->stackPeak;
    }
    (*toCostBenefits)->deltaT += additionalT;
    toCostBenefits = &((*toCostBenefits)->next);
  }
#ifdef _ADPROFILETRACE
  int nbProfilesAfter = costBenefitsLength(curCheckpoints->costBenefits);
  if (_ADPROFILETRACE==-1 || _ADPROFILETRACE==localCkpRank) {
    printf("  ==> MERGE %i COSTBENEFITS INTO %i --> %i\n",
           nbProfilesAdded, nbProfilesBefore, nbProfilesAfter);
    printf("  >AFTER (peak:%"PRId64" bytes, peak at turn:%"PRId64" bytes) :",
           curCheckpoints->stackPeak, curCheckpoints->stackPeakTurn);
    dumpCostBenefits(curCheckpoints->costBenefits);
    printf("\n");
  }
#endif
  // merging finished.
  assert(getCheckpointLocationRank(callname, filename, lineno)==curCheckpoints->levelCheckpoints->ckpRank);
  releaseLevelCkp(curCheckpoints);
}

void adProfileAdj_turn(char* callname, char* filename) {
  curCheckpoints->stackPeak = adStack_getCurrentStackSize();
  curCheckpoints->stackPeakTurn = curCheckpoints->stackPeak;
}

void adProfileAdj_startRepeat() {
  curCheckpoints->repeatLevel = curCheckpoints->levelCheckpoints;
}

void adProfileAdj_resetRepeat() {
  curCheckpoints->levelCheckpoints = curCheckpoints->repeatLevel;
}

void adProfileAdj_endRepeat() {
  while (curCheckpoints->repeatLevel != curCheckpoints->levelCheckpoints) {
    // Now that we are no longer in Repeat mode, free preserved LevelCkpStack's:
    LevelCkpStack *previousElem = curCheckpoints->repeatLevel->previous;
    free(curCheckpoints->repeatLevel);
    curCheckpoints->repeatLevel = previousElem;
  }
  curCheckpoints->repeatLevel = NULL;
}

typedef struct _SortedCostBenefit {
  unsigned int ckpRank; 
  cycle_time_t deltaT;
  delta_size_t deltaPM;
  float ratio;
  struct _SortedCostBenefit* next;
} SortedCostBenefit;

int callNameComesAfter(unsigned int ckpRank1, unsigned int ckpRank2) {
  int comparison = strcmp(allCkpLocations[ckpRank1]->callee_name, allCkpLocations[ckpRank2]->callee_name);
  return (comparison>0 || (comparison==0 && allCkpLocations[ckpRank1]->line_nb > allCkpLocations[ckpRank2]->line_nb));
}

void showOneCostBenefit(SortedCostBenefit* sortedCostBenefits) {
  CkpLocation* ckpLocation = allCkpLocations[sortedCostBenefits->ckpRank];
  cycle_time_t deltaT = (sortedCostBenefits->deltaT / 1000); // now milliSeconds!
  int deltaTsec = deltaT/1000;
  int deltaTmillisec = deltaT%1000;
  printf("  - Time gain -%2d.%03d s.", deltaTsec, deltaTmillisec);
  delta_size_t deltaPM = sortedCostBenefits->deltaPM;
  if (deltaPM<0) {
    printf(" and peak memory gain %"PRId64"b", deltaPM);
  } else if (deltaPM==0) {
    printf(" at peak memory cost zero");
  } else {
    int deltaPMMb = deltaPM/1000000;
    int deltaPMb = deltaPM%1000000;
    printf(" at peak memory cost %4d.%06d Mb", deltaPMMb, deltaPMb);
  }
  if (ckpLocation->callee_name) {
    printf(" for call %s (%u times), at", ckpLocation->callee_name, ckpLocation->occurrences);
  } else {
    printf(" for checkpoint (%u times) starting at", ckpLocation->occurrences);
  }
  printf(" location#%u: line %u of file %s\n", sortedCostBenefits->ckpRank, ckpLocation->line_nb, ckpLocation->fname);
  sortedCostBenefits = sortedCostBenefits->next;
}

void adProfileAdj_showProfiles() {
  SortedCostBenefit* sortedCostBenefitsN = NULL;
  SortedCostBenefit* sortedCostBenefitsZ = NULL;
  SortedCostBenefit* sortedCostBenefitsP = NULL;
  CostBenefit *costBenefits = initialCkpStack.costBenefits;
  SortedCostBenefit** toSortedCostBenefits;
  float ratio;
  while (costBenefits) {
    unsigned int newCkpRank = costBenefits->ckpRank;
    if (costBenefits->deltaPM <= 0) {
      // Sort negative memory costs and zero memory costs in two lists ordered by routine name then location rank:
      toSortedCostBenefits = (costBenefits->deltaPM==0 ? &sortedCostBenefitsZ : &sortedCostBenefitsN);
      ratio = 0.0;
      while (*toSortedCostBenefits && callNameComesAfter(newCkpRank, (*toSortedCostBenefits)->ckpRank)) {
        toSortedCostBenefits = &((*toSortedCostBenefits)->next);
      }
    } else {
      // Sort positive memory costs by decreasing time/memory benefit:
      toSortedCostBenefits = &sortedCostBenefitsP;
      ratio = ((float)costBenefits->deltaT) / ((float)costBenefits->deltaPM);
      while (*toSortedCostBenefits && ratio<(*toSortedCostBenefits)->ratio) {
        toSortedCostBenefits = &((*toSortedCostBenefits)->next);
      }
    }
    SortedCostBenefit *newSorted = (SortedCostBenefit*)malloc(sizeof(SortedCostBenefit));
    newSorted->ckpRank = newCkpRank;
    newSorted->deltaT = costBenefits->deltaT;
    newSorted->deltaPM = costBenefits->deltaPM;
    newSorted->ratio = ratio;
    newSorted->next = *toSortedCostBenefits;
    *toSortedCostBenefits = newSorted;
    costBenefits = costBenefits->next;
  }
  printf("CLOCKS_PER_SEC:%d\n",CLOCKS_PER_SEC);
  printf("PEAK STACK:%"PRId64" bytes\n", initialCkpStack.stackPeak);
  printf("SUGGESTED NOCHECKPOINTs:\n");
  printf(" * Peak memory gain:\n");
  while (sortedCostBenefitsN) {
    showOneCostBenefit(sortedCostBenefitsN);
    sortedCostBenefitsN = sortedCostBenefitsN->next;
  }
  printf(" * Peak memory neutral:\n");
  while (sortedCostBenefitsZ) {
    showOneCostBenefit(sortedCostBenefitsZ);
    sortedCostBenefitsZ = sortedCostBenefitsZ->next;
  }
  printf(" * Peak memory cost:\n");
  while (sortedCostBenefitsP) {
    showOneCostBenefit(sortedCostBenefitsP);
    sortedCostBenefitsP = sortedCostBenefitsP->next;
  }
}

/****************** INTERFACE CALLED FROM FORTRAN *******************/

void adprofileadj_snpwrite_(char* callname, char* filename, unsigned int *lineno) {
  adProfileAdj_SNPWrite(callname, filename, *lineno);
}

void adprofileadj_beginadvance_(char* callname, char* filename, unsigned int *lineno) {
  adProfileAdj_beginAdvance(callname, filename, *lineno);
}

void adprofileadj_endadvance_(char* callname, char* filename, unsigned int *lineno) {
  adProfileAdj_endAdvance(callname, filename, *lineno);
}

void adprofileadj_snpread_(char* callname, char* filename, unsigned int *lineno) {
  adProfileAdj_SNPRead(callname, filename, *lineno);
}

void adprofileadj_beginreverse_(char* callname, char* filename, unsigned int *lineno){
  adProfileAdj_beginReverse(callname, filename, *lineno);
}

void adprofileadj_endreverse_(char* callname, char* filename, unsigned int *lineno) {
  adProfileAdj_endReverse(callname, filename, *lineno);
}

void adprofileadj_turn_(char* callname, char* filename) {
  adProfileAdj_turn(callname, filename);
}

void adprofileadj_startrepeat_() {
  adProfileAdj_startRepeat();
}

void adprofileadj_resetrepeat_() {
  adProfileAdj_resetRepeat();
}

void adprofileadj_endrepeat_() {
  adProfileAdj_endRepeat();
}

void adprofileadj_showprofiles_() {
  adProfileAdj_showProfiles();
}
