/*
 * TAPENADE Automatic Differentiation Engine
 * Copyright (C) 1999-2024 Inria
 * See the LICENSE.md file in the project root for more information.
 *
 * 
 * This package gathers functions used to profile a piece of (adjoint differentiated) code in terms of memory and time. 
 */

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

// #define DEBUG_PROFILE 1

typedef uint64_t stack_size_t;
typedef int64_t  delta_size_t;

typedef clock_t cycle_time_t;
cycle_time_t mytime_() {
        return clock();
}

/** The maximum number of code checkpoint locations that we can handle */
#define MAX_NB_CHECKPOINTS 2000

/** Elementary cost/benefit data of a (yes/no) checkpointing decision on a static checkpoint location. */
typedef struct _CostBenefit{
  unsigned int ckpRank; // Rank of the checkpoint location this cost/benefit applies to. 
  cycle_time_t deltaT; // The (always >=0) extra runtime cost of checkpointing this location.
  delta_size_t deltaPM; // The extra peak memory cost of NOT checkpointing this location.
  struct _CostBenefit *next; // Following CostBenefit infos, about other checkpoint locations.
} CostBenefit;

/** One element in a list of runtime checkpoints, each with some time/peak memory cost data.
 * The checkpoint locations in the list are all directly contained in the same parent runtime checkpoint.
 * The list is ordered from the most recent (downtream) checkpoint till the most upsream. */
typedef struct _LevelCkpStack {
  unsigned int ckpRank;
  cycle_time_t ckpT;
  cycle_time_t basisCkpT;
  stack_size_t sizeAfterSnp;
  struct _LevelCkpStack *previous;
} LevelCkpStack;

/** One element in the list of nested runtime checkpoints at some present time during the profiling run.
 * The list is ordered from the deepest runtime checkpoint to the topmost enclosing runtime checkpoint. */
//TODO meilleur nom que CkpStack (stack deja pris!)
typedef struct _CkpStack {
  LevelCkpStack *levelCheckpoints;
  CostBenefit *costBenefits;
  stack_size_t stackPeak;
  LevelCkpStack *repeatLevel; // for management of start/reset/end Repeat
  struct _CkpStack *above;
  struct _CkpStack *below;
} CkpStack;

/** The main global stack of active checkpoints at some present time during the profiling run. */
CkpStack initialCkpStack = {.levelCheckpoints=NULL, .costBenefits=NULL,
                            .stackPeak=0, .repeatLevel=NULL,
                            .above=NULL, .below=NULL};
CkpStack *curCheckpoints = &initialCkpStack;

/************************* UTILITIES ************************/ 

/** Describes one (static) checkpoint location in the code source, with a unique rank attached to it.
 * A (yes/no) checkpointing decision applies to one static checkpoint location globally. */
typedef struct {
  char* fname; // Name of the file in which the checkpoint appears
  unsigned int line_nb; // Line number in the original code where the checkpoint appears
  char* callee_name; // Name of the procedure being called (when checkpoint is on a call)
  unsigned int occurrences;
  unsigned int ckpRank; // Unique integer rank (from 0 up) attached to this checkpoint location.
} CkpLocation;

/** The array of all static checkpoint locations encountered so far */
CkpLocation* allCkpLocations[MAX_NB_CHECKPOINTS] = {};

unsigned int getCheckpointLocationRank(char* callname, char* filename, unsigned int lineno) {
  unsigned int i = 0;
  int found = 0;
  while (i<MAX_NB_CHECKPOINTS && !found) {
    if (allCkpLocations[i]==NULL) {
      allCkpLocations[i] = (CkpLocation *)malloc(sizeof(CkpLocation));
      allCkpLocations[i]->ckpRank = i;
      allCkpLocations[i]->fname = (char*)malloc((1+strlen(filename))*sizeof(char));
      strcpy(allCkpLocations[i]->fname, filename);
      allCkpLocations[i]->line_nb = lineno ;
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
    if (!found) ++i ;
  }
  return i;
}

LevelCkpStack* stackNewLevelCkp(unsigned int ckpRank, LevelCkpStack* levelCheckpoints) {
  LevelCkpStack* additionalLevelCkp = (LevelCkpStack*)malloc(sizeof(LevelCkpStack));
  additionalLevelCkp->ckpRank = ckpRank;
  additionalLevelCkp->ckpT = 0;
  additionalLevelCkp->basisCkpT = 0;
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
    newCkp->above = NULL ;
  }
  newCkp->levelCheckpoints = NULL;
  newCkp->costBenefits = NULL;
  newCkp->stackPeak = 0;
  return newCkp;
}

CkpStack* closeCkp(CkpStack* checkpoints) {
  return checkpoints->below;
}

CostBenefit** getSetTCB(CostBenefit **toCostBenefits, unsigned int rank) {
  while (*toCostBenefits!=NULL && (*toCostBenefits)->ckpRank < rank) {
    toCostBenefits = &((*toCostBenefits)->next) ;
  }
  if (*toCostBenefits==NULL || (*toCostBenefits)->ckpRank > rank) {
    CostBenefit* additionalCostBenefits = (CostBenefit*)malloc(sizeof(CostBenefit)) ;
    additionalCostBenefits->ckpRank = rank ;
    additionalCostBenefits->deltaT = 0;
    additionalCostBenefits->deltaPM = 0;
    additionalCostBenefits->next = *toCostBenefits;
    *toCostBenefits = additionalCostBenefits ;
  }
  return toCostBenefits;
}

#ifdef DEBUG_PROFILE

int depthCkp(CkpStack* checkpoints) {
  int i = 0;
  while (checkpoints) {
    ++i;
    checkpoints = checkpoints->below ;
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
    printf(" %u:(%"PRIu64";%"PRId64")",
           costBenefits->ckpRank, costBenefits->deltaT, costBenefits->deltaPM) ;
    costBenefits = costBenefits->next;
  }
}

#endif

/****************** MAIN PUBLISHED PRIMITIVES ****************/

void adProfileAdj_SNPWrite(char* callname, char* filename, unsigned int lineno) {
  unsigned int ckpRank = getCheckpointLocationRank(callname, filename, lineno);
  curCheckpoints->levelCheckpoints = stackNewLevelCkp(ckpRank, curCheckpoints->levelCheckpoints);
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
  curCheckpoints = openNewCkp(curCheckpoints) ;
}

void adProfileAdj_endReverse(char* callname, char* filename, unsigned int lineno) {
  assert(curCheckpoints->levelCheckpoints==NULL);
  CostBenefit *additionalCostBenefits = curCheckpoints->costBenefits ;
  stack_size_t peakSize2 = curCheckpoints->stackPeak ;
  curCheckpoints = closeCkp(curCheckpoints);
  // Now merge the additional costs from the checkpointed fragment into the current costs:
  unsigned int localCkpRank = curCheckpoints->levelCheckpoints->ckpRank;
  cycle_time_t localCkpTimeCost = curCheckpoints->levelCheckpoints->ckpT;
  stack_size_t peakSize1 = curCheckpoints->stackPeak ;
  delta_size_t localCkpPeakSizeCost =
    (peakSize1>peakSize2 ? peakSize2 : peakSize1) - curCheckpoints->levelCheckpoints->sizeAfterSnp;
  int mergedLocal = 0 ;
  unsigned int additionalCkpRank;
  cycle_time_t additionalT;
  delta_size_t additionalPM;
#ifdef DEBUG_PROFILE
  printf("ENDREVERSE OF %s::%i CALL %s, DEPTH:%i, SIZEAFTERSNP:%"PRId64"\n",
         allCkpLocations[localCkpRank]->fname,
         allCkpLocations[localCkpRank]->line_nb,
         (allCkpLocations[localCkpRank]->callee_name ? allCkpLocations[localCkpRank]->callee_name : "MANUAL"),
         depthCkp(curCheckpoints),
         curCheckpoints->levelCheckpoints->sizeAfterSnp);
  printf("  BEFORE (peak:%"PRId64" bytes) :", peakSize1);
  dumpCostBenefits(curCheckpoints->costBenefits);
  printf("\n  ADDING (peak:%"PRId64" bytes) :", peakSize2);
  dumpCostBenefits(additionalCostBenefits);
  printf("\n") ;
  int nbProfilesBefore = costBenefitsLength(curCheckpoints->costBenefits);
  int nbProfilesAdded = costBenefitsLength(additionalCostBenefits);
#endif
  CostBenefit **toCostBenefits = &(curCheckpoints->costBenefits) ;
  // compute merged peak size:
  if (peakSize2>peakSize1) curCheckpoints->stackPeak = peakSize2 ;
  // compute merged list of costs/benefits:
  while (additionalCostBenefits || !mergedLocal) {
    if (!mergedLocal
        && (additionalCostBenefits==NULL || additionalCostBenefits->ckpRank>localCkpRank)) {
      additionalCkpRank = localCkpRank;
      additionalT = localCkpTimeCost;
      additionalPM = localCkpPeakSizeCost;
      mergedLocal = 1;
    } else if (additionalCostBenefits->ckpRank==localCkpRank) {
      additionalCkpRank = localCkpRank ;
      additionalT = additionalCostBenefits->deltaT + localCkpTimeCost ;
      additionalPM = additionalCostBenefits->deltaPM + localCkpPeakSizeCost ;
      CostBenefit* tmp = additionalCostBenefits->next ;
      free(additionalCostBenefits);
      additionalCostBenefits = tmp ;
      mergedLocal = 1;
    } else {
      additionalCkpRank = additionalCostBenefits->ckpRank;
      additionalT = additionalCostBenefits->deltaT;
      additionalPM = additionalCostBenefits->deltaPM;
      CostBenefit* tmp = additionalCostBenefits->next ;
      free(additionalCostBenefits);
      additionalCostBenefits = tmp ;
    }
    toCostBenefits = getSetTCB(toCostBenefits, additionalCkpRank) ;
    if (additionalCkpRank==localCkpRank) {
      // Do a sum on memory costs:
      (*toCostBenefits)->deltaPM += additionalPM;
    } else {
      // Do a max on memory costs:
      stack_size_t costPeak1 = peakSize1+(*toCostBenefits)->deltaPM ;
      stack_size_t costPeak2 = peakSize2+additionalPM ;
      (*toCostBenefits)->deltaPM = (costPeak1>costPeak2 ? costPeak1: costPeak2) - curCheckpoints->stackPeak ;
    }
    (*toCostBenefits)->deltaT += additionalT;
    toCostBenefits = &((*toCostBenefits)->next) ;
  }
#ifdef DEBUG_PROFILE
  int nbProfilesAfter = costBenefitsLength(curCheckpoints->costBenefits);
  printf("  ==> MERGE %i COSTBENEFITS INTO %i --> %i\n",
         nbProfilesAdded,
         nbProfilesBefore,
         nbProfilesAfter);
  printf("  >AFTER (peak:%"PRId64" bytes) :", curCheckpoints->stackPeak);
  dumpCostBenefits(curCheckpoints->costBenefits);
  printf("\n") ;
#endif
  // merging finished.
  assert(getCheckpointLocationRank(callname, filename, lineno)==curCheckpoints->levelCheckpoints->ckpRank);
  releaseLevelCkp(curCheckpoints) ;
}

void adProfileAdj_turn(char* callname, char* filename) {
  curCheckpoints->stackPeak = adStack_getCurrentStackSize();
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
    LevelCkpStack *previousElem = curCheckpoints->repeatLevel->previous ;
    free(curCheckpoints->repeatLevel) ;
    curCheckpoints->repeatLevel = previousElem ;
  }
  curCheckpoints->repeatLevel = NULL;
}

void adProfileAdj_showProfiles() {
  CostBenefit *costBenefits = initialCkpStack.costBenefits;
  printf("CLOCKS_PER_SEC:%d\n",CLOCKS_PER_SEC) ;
  printf("PEAK STACK:%"PRId64" bytes\n", initialCkpStack.stackPeak) ;
  printf("CHECKPOINTING COST/BENEFITS:\n");
  while (costBenefits) {
    CkpLocation* ckpLocation = allCkpLocations[costBenefits->ckpRank];
    cycle_time_t deltaT = costBenefits->deltaT;
    cycle_time_t deltaPM = costBenefits->deltaPM;
    if (ckpLocation->callee_name) {
      printf("  - for call %s (executed %u times) at line %u of file %s:\n", ckpLocation->callee_name, ckpLocation->occurrences, ckpLocation->line_nb, ckpLocation->fname);
    } else {
      printf("  - for checkpoint (occurs %u times) starting at line %u of file %s:\n", ckpLocation->occurrences, ckpLocation->line_nb, ckpLocation->fname);
    }
    printf("    time cost of checkpointing: %"PRIu64"\n", costBenefits->deltaT);
    printf("    peak memory cost of NOT checkpointing: %"PRId64" bytes\n", costBenefits->deltaPM);
    costBenefits = costBenefits->next;
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
